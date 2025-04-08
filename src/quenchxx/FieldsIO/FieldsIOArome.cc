/*
 * (C) Copyright 2025 Meteorologisk Institutt
 *
 */

#include "quenchxx/FieldsIO/FieldsIOArome.h"

#include <netcdf.h>

#include <string>
#include <vector>

#include "oops/util/Logger.h"

#define ERR(e, msg) {std::string s(nc_strerror(e)); \
  throw eckit::Exception(s + " : " + msg, Here());}

namespace quenchxx {

// -----------------------------------------------------------------------------

void readArome(const Geometry & geom,
               const varns::Variables & vars,
               const eckit::Configuration & config,
               atlas::FieldSet & fset) {
  oops::Log::trace() << "quenchxx::readArome starting" << std::endl;

  // Build filepath
  std::string filepath = config.getString("filepath");
  if (config.has("member")) {
    std::ostringstream out;
    out << std::setfill('0') << std::setw(6) << config.getInt("member");
    filepath.append("_");
    filepath.append(out.str());
  }

  // NetCDF file path
  std::string ncFilePath = filepath + ".nc";

  // Clear local fieldset
  fset.clear();

  // Variables to copy / to read
  varns::Variables varsToRead;
  for (const auto & var : vars) {
    if (var.name() == "air_pressure" || var.name() == "air_pressure_half") {
      // Get surface pressure and retrieve air_pressure or air_pressure_half from ak/bk
      varsToRead.push_back("SURFPRESSION");
      varsToRead["SURFPRESSION"].setLevels(1);
    } else if (var.name() == "height_above_mean_sea_level_at_surface") {
      // Get surface geopotential and retrieve surface height
      varsToRead.push_back("SPECSURFGEOPOTEN");
      varsToRead["SPECSURFGEOPOTEN"].setLevels(1);
    } else {
      // Variable to read
      varsToRead.push_back(var);
    }
  }

  // Check that both wind components are required or none
  if (varsToRead.has("WIND.U.PHYS") || varsToRead.has("WIND.V.PHYS")) {
    ASSERT(varsToRead.has("WIND.U.PHYS") && varsToRead.has("WIND.V.PHYS"));
  }

  // Create local fieldset
  atlas::FieldSet fsetToRead;
  for (const auto & var : varsToRead) {
    atlas::Field field = geom.functionSpace().createField<double>(
      atlas::option::name(var.name()) | atlas::option::levels(var.getLevels()));
    fsetToRead.add(field);
  }

  // Initialize local fieldset
  for (auto & varField : fsetToRead) {
    auto view = atlas::array::make_view<double, 2>(varField);
    view.assign(0.0);
  }

  // File variables names
  size_t nVarLev = 0;
  std::vector<std::string> varLevName;
  for (const auto & var : varsToRead) {
    for (int jlevel = 0; jlevel < var.getLevels(); ++jlevel) {
      if (var.name() == "SURFPRESSION" || var.name() == "SPECSURFGEOPOTEN") {
        varLevName.push_back(var.name());
      } else {
        const std::string level = std::to_string(jlevel+1);
        varLevName.push_back("S" + std::string(3-level.length(), '0') + level + var.name());
      }
      ++nVarLev;
    }
  }

  // NetCDF IDs
  int ncid, retval, var_id[nVarLev];

  // Global data
  atlas::FieldSet globalData;
  for (const auto & var : varsToRead) {
    atlas::Field varField = geom.functionSpace().createField<double>(
      atlas::option::name(var.name())
      | atlas::option::levels(var.getLevels()) | atlas::option::global());
    globalData.add(varField);
  }

  // StructuredColumns
  atlas::functionspace::StructuredColumns fs(geom.functionSpace());

  if (geom.getComm().rank() == 0) {
    // Get grid
    atlas::StructuredGrid grid = fs.grid();

    // Get sizes
    int nx = grid.nxmax();
    int ny = grid.ny();

    oops::Log::info() << "Info     : Reading file: " << ncFilePath << std::endl;

    // Open NetCDF file
    if ((retval = nc_open(ncFilePath.c_str(), NC_NOWRITE, &ncid))) ERR(retval, ncFilePath);

    // Get variables
    for (size_t jVarLev = 0; jVarLev < nVarLev; ++jVarLev) {
      if ((retval = nc_inq_varid(ncid, varLevName[jVarLev].c_str(), &var_id[jVarLev]))) {
        ERR(retval, varLevName[jVarLev]);
      }
    }

    size_t iVarLev = 0;
    for (const auto & var : varsToRead) {
      auto varField = globalData[var.name()];
      auto varView = atlas::array::make_view<double, 2>(varField);
      for (int jlevel = 0; jlevel < var.getLevels(); ++jlevel) {
        // Read data
        std::vector<double> zvar(ny*nx);
        if ((retval = nc_get_var_double(ncid, var_id[iVarLev], zvar.data()))) {
          ERR(retval, varLevName[iVarLev]);
        }
        ++iVarLev;

        // Copy data
        for (int j = 0; j < ny; ++j) {
          for (int i = 0; i < grid.nx(j); ++i) {
            atlas::gidx_t gidx = grid.index(i, j);
            varView(gidx, jlevel) = zvar[j*nx+i];
          }
        }
      }
    }
  }

  // Scatter data from main processor
  fs.scatter(globalData, fsetToRead);

  // Get pressure from its logarithm
  if (fsetToRead.has("SURFPRESSION")) {
    // Get field
    auto psField = fsetToRead["SURFPRESSION"];

    // Get view
    auto psView = atlas::array::make_view<double, 2>(psField);

    // Apply exponential
    for (int jnode = 0; jnode < psField.shape(0); ++jnode) {
      psView(jnode, 0) = std::exp(psView(jnode, 0));
    }
  }

  // Get eastward and northward winds from local grid winds
  if (fsetToRead.has("WIND.U.PHYS") || fsetToRead.has("WIND.V.PHYS")) {
    // Get fields
    auto uField = fsetToRead["WIND.U.PHYS"];
    auto vField = fsetToRead["WIND.V.PHYS"];

    // Get views
    auto uView = atlas::array::make_view<double, 2>(uField);
    auto vView = atlas::array::make_view<double, 2>(vField);

    // Get lon/lat view
    const auto lonlatView = atlas::array::make_view<double, 2>(geom.functionSpace().lonlat());

    for (int jnode = 0; jnode < uField.shape(0); ++jnode) {
      // Get local point
      atlas::PointLonLat p({lonlatView(jnode, 0), lonlatView(jnode, 1)});

      // Get local Jacobian
      double dx_dlon = geom.grid().projection().jacobian(p).dx_dlon();
      double dx_dlat = geom.grid().projection().jacobian(p).dx_dlat();
      double dy_dlon = geom.grid().projection().jacobian(p).dy_dlon();
      double dy_dlat = geom.grid().projection().jacobian(p).dy_dlat();

      // Normalize Jacobian
      const double dlonNorm = 1.0/std::sqrt(dx_dlon*dx_dlon+dy_dlon*dy_dlon);
      const double dlatNorm = 1.0/std::sqrt(dx_dlat*dx_dlat+dy_dlat*dy_dlat);
      dx_dlon *= dlonNorm;
      dy_dlon *= dlonNorm;
      dx_dlat *= dlatNorm;
      dy_dlat *= dlatNorm;

      // Apply transform
      for (int jlevel = 0; jlevel < uField.shape(1); ++jlevel) {
        const double uPhys = uView(jnode, jlevel);
        const double vPhys = vView(jnode, jlevel);
        uView(jnode, jlevel) = uPhys*dx_dlon + vPhys*dy_dlon;
        vView(jnode, jlevel) = uPhys*dx_dlat + vPhys*dy_dlat;
      }
    }
  }

  // Processing
  for (const auto & var : vars) {
    if (var.name() == "air_pressure" || var.name() == "air_pressure_half") {
      // Retrieve air_pressure or air_pressure_half from ak/bk

      // Hybrid coordinates
      std::vector<double> ak(var.getLevels());
      std::vector<double> bk(var.getLevels());

      if (geom.getComm().rank() == 0) {
        // NetCDF IDs
        int ak_id, bk_id, dim_id;
        size_t nab;

        // Get hybrid coordinates IDs
        const std::string akName = config.getString("ak", "hybrid_coef_A");
        const std::string bkName = config.getString("bk", "hybrid_coef_B");
        if ((retval = nc_inq_varid(ncid, akName.c_str(), &ak_id))) ERR(retval, akName);
        if ((retval = nc_inq_varid(ncid, bkName.c_str(), &bk_id))) ERR(retval, bkName);

        // Get hybrid coordinates dimension
        if ((retval = nc_inq_vardimid(ncid, ak_id, &dim_id))) ERR(retval, akName);
        if ((retval = nc_inq_dimlen(ncid, dim_id, &nab))) ERR(retval, "nab");

        // Read data
        std::vector<double> akFromFile(nab);
        std::vector<double> bkFromFile(nab);
        if ((retval = nc_get_var_double(ncid, ak_id, akFromFile.data()))) ERR(retval, akName);
        if ((retval = nc_get_var_double(ncid, bk_id, bkFromFile.data()))) ERR(retval, bkName);

        if (var.name() == "air_pressure") {
          // Pressure at full levels
          ASSERT(static_cast<int>(nab) == var.getLevels()+1);
          for (int jlevel = 0; jlevel < var.getLevels(); ++jlevel) {
            ak[jlevel] = 0.5*(akFromFile[jlevel]+akFromFile[jlevel+1]);
            bk[jlevel] = 0.5*(bkFromFile[jlevel]+bkFromFile[jlevel+1]);
          }
        } else if (var.name() == "air_pressure_half") {
          // Pressure at half levels
          ASSERT(static_cast<int>(nab) == var.getLevels());
          for (int jlevel = 0; jlevel < var.getLevels(); ++jlevel) {
            ak[jlevel] = akFromFile[jlevel];
            bk[jlevel] = bkFromFile[jlevel];
          }
        }
      }

      // Broadcast hybrid coordinates
      geom.getComm().broadcast(ak.begin(), ak.end(), 0);
      geom.getComm().broadcast(bk.begin(), bk.end(), 0);

      // Create field
      atlas::Field varField = geom.functionSpace().createField<double>(
        atlas::option::name(var.name()) | atlas::option::levels(var.getLevels()));
      fset.add(varField);

      // Get views
      const auto psView = atlas::array::make_view<double, 2>(fsetToRead["SURFPRESSION"]);
      auto varView = atlas::array::make_view<double, 2>(varField);

      // Compute pressure
      for (int jnode = 0; jnode < varField.shape(0); ++jnode) {
        for (int jlevel = 0; jlevel < var.getLevels(); ++jlevel) {
          varView(jnode, jlevel) = ak[jlevel] + bk[jlevel]*psView(jnode, 0);
        }
      }
    } else if (var.name() == "height_above_mean_sea_level_at_surface") {
      // Retrieve surface height

      // Create field
      atlas::Field varField = geom.functionSpace().createField<double>(
        atlas::option::name(var.name()) | atlas::option::levels(1));
      fset.add(varField);

      // Get views
      const auto zsView = atlas::array::make_view<double, 2>(fsetToRead["SPECSURFGEOPOTEN"]);
      auto varView = atlas::array::make_view<double, 2>(varField);

      // Compute surface height
      const double gInv = 1.0/9.81;
      for (int jnode = 0; jnode < varField.shape(0); ++jnode) {
        varView(jnode, 0) = zsView(jnode, 0)*gInv;
      }
    } else {
      // Add fields
      fset.add(fsetToRead[var.name()]);
    }
  }

  if (geom.getComm().rank() == 0) {
    // Close file
    if ((retval = nc_close(ncid))) ERR(retval, ncFilePath);
  }

  // Code is too complicated, mark dirty to be safe
  fset.set_dirty();

  oops::Log::trace() << "quenchxx::readArome done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace quenchxx
