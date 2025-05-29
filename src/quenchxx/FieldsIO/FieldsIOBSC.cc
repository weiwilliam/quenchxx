/*
 * (C) Copyright 2025 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "quenchxx/FieldsIO/FieldsIOBSC.h"

#include <netcdf.h>

#include <vector>

#include "oops/util/Logger.h"

#include "quenchxx/Geometry.h"

#define ERR(e, msg) {std::string s(nc_strerror(e)); \
  throw eckit::Exception(s + " : " + msg, Here());}

namespace quenchxx {

// -----------------------------------------------------------------------------

static FieldsIOMaker<FieldsIOBSC> makerBSC_("bsc");

// -----------------------------------------------------------------------------

static std::vector<std::string> existingFiles_;

// -----------------------------------------------------------------------------

void FieldsIOBSC::read(const Geometry & geom,
                       const varns::Variables & vars,
                       const eckit::Configuration & conf,
                       atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::read starting" << std::endl;

  // Get function space
  const atlas::functionspace::StructuredColumns fs(geom.functionSpace());

  // Clear local fieldset
  fset.clear();

  // Create local fieldset
  for (const auto & var : vars) {
    atlas::Field field = fs.createField<double>(
      atlas::option::name(var.name()) | atlas::option::levels(var.getLevels()));
    fset.add(field);
  }

  // Initialize local fieldset
  for (auto & field : fset) {
    auto view = atlas::array::make_view<double, 2>(field);
    view.assign(0.0);
  }

  // Global data
  atlas::FieldSet globalData;
  for (const auto & var : vars) {
    atlas::Field field = fs.createField<double>(
      atlas::option::name(var.name()) | atlas::option::levels(var.getLevels())
      | atlas::option::global());
    globalData.add(field);
  }

  // Get filepath
  const std::string ncFilePath = conf.getString("filepath");

  // Get levels selection (from geometry section)
  const std::vector<size_t> levels = geom.io().getUnsignedVector("levels selection");

  // Get file initial time
  const util::DateTime initialTime(geom.io().getString("initial date"));

  // Get file final time
  const util::DateTime finalTime(geom.io().getString("final date"));

  // Get total number of hours
  ASSERT(finalTime >= initialTime);
  const size_t timeMax = (finalTime-initialTime).toSeconds()/3600+1;

  // Get read time
  const util::DateTime validTime(conf.getString("date"));

  // Difference in hours
  ASSERT(validTime >= initialTime);
  ASSERT(validTime <= finalTime);
  const size_t time = (validTime-initialTime).toSeconds()/3600;

  // NetCDF IDs
  int ncid, retval, time_id, var_id[vars.size()];

  if (geom.getComm().rank() == 0) {
    // Get grid
    const atlas::StructuredGrid grid = fs.grid();

    // Get sizes
    const size_t nx = grid.nxmax();
    const size_t ny = grid.ny();

    oops::Log::info() << "Info     : Reading file: " << ncFilePath << std::endl;

    // Open NetCDF file
    if ((retval = nc_open(ncFilePath.c_str(), NC_NOWRITE, &ncid))) ERR(retval, ncFilePath);

    // Check that initialTime is consistent with "time" metadata in file
    if ((retval = nc_inq_varid(ncid, "time", &time_id))) ERR(retval, "time");
    const std::string time_units_key = "units";
    size_t attlen;
    if (retval = nc_inq_attlen(ncid, time_id, time_units_key.c_str(), &attlen))
      ERR(retval, "time");
    char *time_units_char = reinterpret_cast<char*>(malloc(attlen+1));
    if (retval = nc_get_att_text(ncid, time_id, time_units_key.c_str(), time_units_char))
      ERR(retval, "time");
    const std::string time_units_value(time_units_char);
    const std::string initialTimeFromFileStr = time_units_value.substr(12, 10) + "T"
      + time_units_value.substr(23, 5) + ":00Z";
    const util::DateTime initialTimeFromFile(initialTimeFromFileStr);
    ASSERT(initialTime == initialTimeFromFile);

    // Check number of times
    size_t timeMaxFromFile;
    if ((retval = nc_inq_dimlen(ncid, time_id, &timeMaxFromFile))) ERR(retval, "time");
    ASSERT(timeMax == timeMaxFromFile);

    for (size_t jvar = 0; jvar < vars.size(); ++jvar) {
      // Get variables ID
      if ((retval = nc_inq_varid(ncid, vars[jvar].name().c_str(), &var_id[jvar])))
        ERR(retval, vars[jvar].name());

      // Get variable view
      auto varView = atlas::array::make_view<double, 2>(globalData[vars[jvar].name()]);

      for (size_t k = 0; k < levels.size(); ++k) {
        // Read level
        std::vector<double> zvar(nx*ny);
        const std::vector<size_t> startp({time, levels[k]-1, 0, 0});
        const std::vector<size_t> countp({1, 1, ny, nx});
        if ((retval = nc_get_vars_double(ncid, var_id[jvar], startp.data(), countp.data(), NULL,
          zvar.data()))) ERR(retval, vars[jvar].name());

        // Deserialize data to view
        for (atlas::idx_t j = 0; j < ny; ++j) {
          for (atlas::idx_t i = 0; i < grid.nx(ny-1-j); ++i) {
            atlas::gidx_t gidx = grid.index(i, ny-1-j);
            varView(gidx, k) = zvar[j*nx+i];
          }
        }
      }
    }

    // Close file
    if ((retval = nc_close(ncid))) ERR(retval, ncFilePath);
  }

  // Scatter data from main processor
  fs.scatter(globalData, fset);

  // Mark dirty to be safe
  fset.set_dirty();

  oops::Log::trace() << classname() << "::read done" << std::endl;
}

// -----------------------------------------------------------------------------

void FieldsIOBSC::write(const Geometry & geom,
                        const eckit::Configuration & conf,
                        const atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::write starting" << std::endl;

  // Get function space
  const atlas::functionspace::StructuredColumns fs(geom.functionSpace());

  // Define variables vector from fset
  const std::vector<std::string> vars = fset.field_names();

  // Get filepath
  const std::string ncFilePath = conf.getString("filepath");

  // Check if this file already exists
  const bool existingFile = std::find(existingFiles_.begin(), existingFiles_.end(), ncFilePath)
    != existingFiles_.end();
  if (!existingFile) {
    existingFiles_.push_back(ncFilePath);
  }

  // Get total number of levels (from geometry section)
  const size_t lmMax = geom.io().getUnsigned("total number of levels");

  // Get levels selection (from geometry section)
  const std::vector<size_t> levels = geom.io().getUnsignedVector("levels selection");

  // Get file initial time
  const util::DateTime initialTime(geom.io().getString("initial date"));

  // Get file final time
  const util::DateTime finalTime(geom.io().getString("final date"));

  // Get total number of hours
  ASSERT(finalTime >= initialTime);
  const size_t timeMax = (finalTime-initialTime).toSeconds()/3600+1;

  // Get write time
  const util::DateTime validTime(conf.getString("date"));

  // Difference in hours
  ASSERT(validTime >= initialTime);
  ASSERT(validTime <= finalTime);
  const size_t time = (validTime-initialTime).toSeconds()/3600;

  // NetCDF IDs
  int retval, ncid, rlon_id, rlat_id, lm_id, time_id,
    dRlon_id[1], dRlat_id[1], dLm_id[1], dTime_id[1], d2D_id[2], d4D_id[4],
    vRlon_id, vRlat_id, vLm_id, vTime_id,
    lon_id, lat_id, var_id[vars.size()];

  // Prepare local coordinates and data
  atlas::FieldSet localData;
  if (!existingFile) {
    atlas::Field lonLocal = fs.createField<double>(atlas::option::name("lon"));
    localData.add(lonLocal);
    atlas::Field latLocal = fs.createField<double>(atlas::option::name("lat"));
    localData.add(latLocal);
    const auto lonlatView = atlas::array::make_view<double, 2>(fs.lonlat());
    auto lonViewLocal = atlas::array::make_view<double, 1>(localData["lon"]);
    auto latViewLocal = atlas::array::make_view<double, 1>(localData["lat"]);
    for (atlas::idx_t jnode = 0; jnode < fs.lonlat().shape(0); ++jnode) {
       lonViewLocal(jnode) = lonlatView(jnode, 0);
       latViewLocal(jnode) = lonlatView(jnode, 1);
    }
  }
  for (size_t jvar = 0; jvar < vars.size(); ++jvar) {
    localData.add(fset[vars[jvar]]);
  }

  // Prepare global coordinates and data
  atlas::FieldSet globalData;
  if (!existingFile) {
    atlas::Field lonGlobal = fs.createField<double>(
      atlas::option::name("lon") | atlas::option::global());
    globalData.add(lonGlobal);
    atlas::Field latGlobal = fs.createField<double>(
      atlas::option::name("lat") | atlas::option::global());
    globalData.add(latGlobal);
  }
  for (size_t jvar = 0; jvar < vars.size(); ++jvar) {
    atlas::Field globalField = fs.createField<double>(atlas::option::name(vars[jvar]) |
      atlas::option::levels(fset[vars[jvar]].shape(1)) | atlas::option::global());
    globalData.add(globalField);
  }

  // Gather coordinates and data on main processor
  fs.gather(localData, globalData);

  if (geom.getComm().rank() == 0) {
    if (existingFile) {
      oops::Log::info() << "Info     : Updating file: " << ncFilePath << std::endl;
    } else {
      oops::Log::info() << "Info     : Writing file: " << ncFilePath << std::endl;
    }

    // Get grid
    const atlas::StructuredGrid grid = fs.grid();

    // Get sizes
    const size_t nx = grid.nxmax();
    const size_t ny = grid.ny();

    // Definition mode

    if (existingFile) {
      // Open NetCDF file
      if ((retval = nc_open(ncFilePath.c_str(),
        NC_64BIT_OFFSET | NC_WRITE, &ncid))) ERR(retval, ncFilePath);

      // Switch to definition mode
      if ((retval = nc_redef(ncid))) ERR(retval, ncFilePath);
    } else {
      // Create NetCDF file
      if ((retval = nc_create(ncFilePath.c_str(),
        NC_64BIT_OFFSET | NC_CLOBBER, &ncid))) ERR(retval, ncFilePath);
    }

    if (existingFile) {
      // Get dimension
      if ((retval = nc_inq_dimid(ncid, "rlon", &rlon_id))) ERR(retval, "rlon");
      if ((retval = nc_inq_dimid(ncid, "rlat", &rlat_id))) ERR(retval, "rlat");
      if ((retval = nc_inq_dimid(ncid, "lm", &lm_id))) ERR(retval, "lm");
      if ((retval = nc_inq_dimid(ncid, "time", &time_id))) ERR(retval, "time");
    } else {
      // Create dimensions
      if ((retval = nc_def_dim(ncid, "rlon", nx, &rlon_id))) ERR(retval, "rlon");
      if ((retval = nc_def_dim(ncid, "rlat", ny, &rlat_id))) ERR(retval, "rlat");
      if ((retval = nc_def_dim(ncid, "lm", lmMax, &lm_id))) ERR(retval, "lm");
      if ((retval = nc_def_dim(ncid, "time", NC_UNLIMITED, &time_id))) ERR(retval, "time");
    }

    // Dimensions arrays
    dRlon_id[0] = rlon_id;
    dRlat_id[0] = rlat_id;
    dLm_id[0] = lm_id;
    dTime_id[0] = time_id;
    d2D_id[0] = rlat_id;
    d2D_id[1] = rlon_id;
    d4D_id[0] = time_id;
    d4D_id[1] = lm_id;
    d4D_id[2] = rlat_id;
    d4D_id[3] = rlon_id;

    if (!existingFile) {
      // Define coordinates
      if ((retval = nc_def_var(ncid, "rlon", NC_FLOAT, 1, dRlon_id, &vRlon_id)))
        ERR(retval, "rlon");
      // TODO(Andrea): add metadata
      if ((retval = nc_def_var(ncid, "rlat", NC_FLOAT, 1, dRlat_id, &vRlat_id)))
        ERR(retval, "rlat");
      // TODO(Andrea): add metadata
      if ((retval = nc_def_var(ncid, "lm", NC_INT, 1, dLm_id, &vLm_id))) ERR(retval, "lm");
      // TODO(Andrea): add metadata
      if ((retval = nc_def_var(ncid, "time", NC_INT, 1, dTime_id, &vTime_id))) ERR(retval, "time");
      // TODO(Andrea): add metadata
      if ((retval = nc_def_var(ncid, "lon", NC_FLOAT, 2, d2D_id, &lon_id))) ERR(retval, "lon");
      // TODO(Andrea): add metadata
      if ((retval = nc_def_var(ncid, "lat", NC_FLOAT, 2, d2D_id, &lat_id))) ERR(retval, "lat");
      // TODO(Andrea): add metadata
    }
    for (size_t jvar = 0; jvar < vars.size(); ++jvar) {
      // Check whether this variable exists
      if (nc_inq_varid(ncid, vars[jvar].c_str(), &var_id[jvar]) != NC_NOERR) {
        // Define variable
        if ((retval = nc_def_var(ncid, vars[jvar].c_str(), NC_FLOAT, 4, d4D_id, &var_id[jvar])))
          ERR(retval, vars[jvar]);
        // TODO(Andrea): add metadata
      }
    }

    // End definition mode
    if ((retval = nc_enddef(ncid))) ERR(retval, ncFilePath);

    // Data mode

    if (!existingFile) {
      // Copy coordinates
      const auto lonViewGlobal = atlas::array::make_view<double, 1>(globalData["lon"]);
      const auto latViewGlobal = atlas::array::make_view<double, 1>(globalData["lat"]);
      std::vector<float> zlon(ny*nx);
      std::vector<float> zlat(ny*nx);
      for (atlas::idx_t j = 0; j < ny; ++j) {
        for (atlas::idx_t i = 0; i < grid.nx(ny-1-j); ++i) {
          atlas::gidx_t gidx = grid.index(i, ny-1-j);
          zlon[j*nx + i] = lonViewGlobal(gidx);
          zlat[j*nx + i] = latViewGlobal(gidx);
        }
      }

      // Create rlon
      std::vector<float> zRlon(nx);
      const float rlonStart = grid.spec().getFloat("xspace.start");
      const float rlonEnd = grid.spec().getFloat("xspace.end");
      for (size_t jLon = 0; jLon < nx; ++jLon) {
        zRlon[jLon] = rlonStart + static_cast<float>(jLon)*(rlonEnd-rlonStart)
          /static_cast<float>(nx-1);
      }

      // Create rlat
      std::vector<float> zRlat(ny);
      const float rlatStart = grid.spec().getFloat("yspace.start");
      const float rlatEnd = grid.spec().getFloat("yspace.end");
      for (size_t jLat = 0; jLat < ny; ++jLat) {
        zRlat[jLat] = rlatStart + static_cast<float>(jLat)*(rlatEnd-rlatStart)
          /static_cast<float>(ny-1);
      }

      // Create lm
      std::vector<int> zLm(lmMax);
      for (size_t jLm = 0; jLm < lmMax; ++jLm) {
        zLm[jLm] = jLm;
      }

      // Create time
      std::vector<int> zTime(timeMax);
      for (size_t jTime = 0; jTime < timeMax; ++jTime) {
        zTime[jTime] = jTime;
      }

      // Write coordinates
      if ((retval = nc_put_var_float(ncid, vRlon_id, zRlon.data()))) ERR(retval, "rlon");
      if ((retval = nc_put_var_float(ncid, vRlat_id, zRlat.data()))) ERR(retval, "rlat");
      if ((retval = nc_put_var_int(ncid, vLm_id, zLm.data()))) ERR(retval, "lm");
      const std::vector<size_t> startp({0});
      const std::vector<size_t> countp({timeMax});
      if ((retval = nc_put_vars_int(ncid, vTime_id, startp.data(), countp.data(), NULL,
        zTime.data()))) ERR(retval, "time");
      if ((retval = nc_put_var_float(ncid, lon_id, zlon.data()))) ERR(retval, "lon");
      if ((retval = nc_put_var_float(ncid, lat_id, zlat.data()))) ERR(retval, "lat");
    }

    for (size_t jvar = 0; jvar < vars.size(); ++jvar) {
      // Get variable view
      const auto varView = atlas::array::make_view<double, 2>(globalData[vars[jvar]]);

      for (atlas::idx_t k = 0; k < fset[vars[jvar]].shape(1); ++k) {
        // Copy data
        std::vector<float> zvar(ny*nx);
        for (atlas::idx_t j = 0; j < ny; ++j) {
          for (atlas::idx_t i = 0; i < grid.nx(ny-1-j); ++i) {
            atlas::gidx_t gidx = grid.index(i, ny-1-j);
            zvar[j*nx + i] = varView(gidx, k);
          }
        }

        // Write data
        const std::vector<size_t> startp({time, levels[k]-1, 0, 0});
        const std::vector<size_t> countp({1, 1, ny, nx});
        if ((retval = nc_put_vars_float(ncid, var_id[jvar], startp.data(), countp.data(), NULL,
          zvar.data()))) ERR(retval, vars[jvar]);
      }
    }

    // Close file
    if ((retval = nc_close(ncid))) ERR(retval, ncFilePath);
  }

  oops::Log::trace() << classname() << "::write done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace quenchxx
