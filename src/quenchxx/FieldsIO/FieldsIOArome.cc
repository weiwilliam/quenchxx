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
    if (var.name() == "log_of_air_pressure_at_surface" || var.name() == "air_pressure_at_surface"
      || var.name() == "air_pressure" || var.name() == "air_pressure_at_half_levels") {
      // Get surface pressure
      varsToRead.push_back("SURFPRESSION");
      varsToRead["SURFPRESSION"].setLevels(1);
    } else if (var.name() == "height_above_mean_sea_level_at_surface") {
      // Get surface geopotential and retrieve surface height
      varsToRead.push_back("SPECSURFGEOPOTEN");
      varsToRead["SPECSURFGEOPOTEN"].setLevels(1);
    } else if (var.name() == "geographical_x_wind" || var.name() == "eastward_wind") {
      // Get u wind
      varsToRead.push_back("WIND.U.PHYS");
      varsToRead["WIND.U.PHYS"].setLevels(var.getLevels());
    } else if (var.name() == "geographical_y_wind" || var.name() == "northward_wind") {
      // Get v wind
      varsToRead.push_back("WIND.V.PHYS");
      varsToRead["WIND.V.PHYS"].setLevels(var.getLevels());
    } else if (var.name() == "air_temperature") {
      // Get temperature
      varsToRead.push_back("TEMPERATURE");
      varsToRead["TEMPERATURE"].setLevels(var.getLevels());
    } else if (var.name() == "water_vapor_mixing_ratio_wrt_moist_air") {
      // Get specific humidity
      varsToRead.push_back("HUMI.SPECIFI");
      varsToRead["HUMI.SPECIFI"].setLevels(var.getLevels());
    } else {
      // Unknown variable
      throw eckit::Exception("unknown variable", Here());
    }
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

  // Processing
  for (const auto & var : vars) {
    if (var.name() == "log_of_air_pressure_at_surface") {
      // Share field
      fset.add(fsetToRead["SURFPRESSION"]);
      fset["SURFPRESSION"].rename("log_of_air_pressure_at_surface");
    }

    if (var.name() == "air_pressure_at_surface"
      || var.name() == "air_pressure" || var.name() == "air_pressure_at_half_levels") {
      // Create field
      atlas::Field varField = geom.functionSpace().createField<double>(
        atlas::option::name(var.name()) | atlas::option::levels(var.getLevels()));
      fset.add(varField);

      // Get view
      auto varView = atlas::array::make_view<double, 2>(varField);

      // Get read view
      const auto logOfPsView = atlas::array::make_view<double, 2>(fsetToRead["SURFPRESSION"]);

      if (var.name() == "air_pressure_at_surface") {
        // Apply exp
        for (int jnode = 0; jnode < varField.shape(0); ++jnode) {
          varView(jnode, 0) = std::exp(logOfPsView(jnode, 0));
        }
      } else if (var.name() == "air_pressure" || var.name() == "air_pressure_at_half_levels") {
        // Retrieve air_pressure or air_pressure_at_half_levels from ak/bk

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
          } else if (var.name() == "air_pressure_at_half_levels") {
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

        // Compute pressure
        for (int jnode = 0; jnode < varField.shape(0); ++jnode) {
          for (int jlevel = 0; jlevel < var.getLevels(); ++jlevel) {
            varView(jnode, jlevel) = ak[jlevel] + bk[jlevel]*std::exp(logOfPsView(jnode, 0));
          }
        }
      }
    } 

    if (var.name() == "height_above_mean_sea_level_at_surface") {
      // Create field
      atlas::Field varField = geom.functionSpace().createField<double>(
        atlas::option::name(var.name()) | atlas::option::levels(var.getLevels()));
      fset.add(varField);

      // Get views
      const auto zsView = atlas::array::make_view<double, 2>(fsetToRead["SPECSURFGEOPOTEN"]);
      auto varView = atlas::array::make_view<double, 2>(varField);

      // Compute surface height
      const double gInv = 1.0/9.81;
      for (int jnode = 0; jnode < varField.shape(0); ++jnode) {
        varView(jnode, 0) = zsView(jnode, 0)*gInv;
      }
    } 

    if (var.name() == "geographical_x_wind") {
      // Share field
      fset.add(fsetToRead["WIND.U.PHYS"]);
      fset["WIND.U.PHYS"].rename("geographical_x_wind");
    }

    if (var.name() == "geographical_y_wind") {
      // Share field
      fset.add(fsetToRead["WIND.V.PHYS"]);
      fset["WIND.V.PHYS"].rename("geographical_y_wind");
    }

    if (var.name() == "eastward_wind" || var.name() == "northward_wind") {
      // Compute spherical winds

      // Create field
      atlas::Field varField = geom.functionSpace().createField<double>(
        atlas::option::name(var.name()) | atlas::option::levels(var.getLevels()));
      fset.add(varField);

      // Get views
      const auto uView = atlas::array::make_view<double, 2>(fsetToRead["WIND.U.PHYS"]);
      const auto vView = atlas::array::make_view<double, 2>(fsetToRead["WIND.V.PHYS"]);
      auto varView = atlas::array::make_view<double, 2>(varField);

      // Get lon/lat view
      const auto lonlatView = atlas::array::make_view<double, 2>(geom.functionSpace().lonlat());

      for (int jnode = 0; jnode < varField.shape(0); ++jnode) {
        // Get local point
        atlas::PointLonLat p({lonlatView(jnode, 0), lonlatView(jnode, 1)});

        if (var.name() == "eastward_wind") {
          // Get local Jacobian
          double dx_dlon = geom.grid().projection().jacobian(p).dx_dlon();
          double dy_dlon = geom.grid().projection().jacobian(p).dy_dlon();

          // Normalize Jacobian
          const double dlonNorm = 1.0/std::sqrt(dx_dlon*dx_dlon+dy_dlon*dy_dlon);
          dx_dlon *= dlonNorm;
          dy_dlon *= dlonNorm;

          // Apply transform
          for (int jlevel = 0; jlevel < varField.shape(1); ++jlevel) {
            varView(jnode, jlevel) = uView(jnode, jlevel)*dx_dlon + vView(jnode, jlevel)*dy_dlon;
          }
        } else if (var.name() == "northward_wind") {
          // Get local Jacobian
          double dx_dlat = geom.grid().projection().jacobian(p).dx_dlat();
          double dy_dlat = geom.grid().projection().jacobian(p).dy_dlat();

          // Normalize Jacobian
          const double dlatNorm = 1.0/std::sqrt(dx_dlat*dx_dlat+dy_dlat*dy_dlat);
          dx_dlat *= dlatNorm;
          dy_dlat *= dlatNorm;

          // Apply transform
          for (int jlevel = 0; jlevel < varField.shape(1); ++jlevel) {
            varView(jnode, jlevel) = uView(jnode, jlevel)*dx_dlat + vView(jnode, jlevel)*dy_dlat;
          }
        }
      }
    } 

    if (var.name() == "air_temperature") {
      // Share field
      fset.add(fsetToRead["TEMPERATURE"]);
      fset["TEMPERATURE"].rename("air_temperature");
    }

    if (var.name() == "water_vapor_mixing_ratio_wrt_moist_air") {
      // Share field
      fset.add(fsetToRead["HUMI.SPECIFI"]);
      fset["HUMI.SPECIFI"].rename("water_vapor_mixing_ratio_wrt_moist_air");
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
