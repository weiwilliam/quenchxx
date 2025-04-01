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
    if (var.name() == "air_pressure") {
      // Copy air_pressure from vertical coordinates
      ASSERT(geom.fields().has("air_pressure"));
      fset.add(geom.fields()["air_pressure"]);
    } else {
      // Variable to read
      varsToRead.push_back(var);
    }
  }

  // Create local fieldset
  atlas::FieldSet fsetToRead;
  for (size_t jvar = 0; jvar < varsToRead.size(); ++jvar) {
    atlas::Field field = geom.functionSpace().createField<double>(
      atlas::option::name(varsToRead[jvar].name()) |
      atlas::option::levels(varsToRead[jvar].getLevels()));
    fsetToRead.add(field);
  }

  // Initialize local fieldset
  for (auto & field : fsetToRead) {
    auto view = atlas::array::make_view<double, 2>(field);
    view.assign(0.0);
  }

  // File variables names
  size_t nVarLev = 0;
  std::vector<std::string> varLevName;
  for (size_t jvar = 0; jvar < varsToRead.size(); ++jvar) {
    for (int k = 0; k < varsToRead[jvar].getLevels(); ++k) {
      if (varsToRead[jvar].name() == "SURFPRESSION") {
        varLevName.push_back(varsToRead[jvar].name());
      } else {
        const std::string level = std::to_string(k+1);
        varLevName.push_back("S" + std::string(3-level.length(), '0') + level
          + varsToRead[jvar].name());
      }
      ++nVarLev;
    }
  }

  // NetCDF IDs
  int ncid, retval, var_id[nVarLev];

  // Global data
  atlas::FieldSet globalData;
  for (size_t jvar = 0; jvar < varsToRead.size(); ++jvar) {
    atlas::Field field = geom.functionSpace().createField<double>(
      atlas::option::name(varsToRead[jvar].name())
      | atlas::option::levels(varsToRead[jvar].getLevels()) | atlas::option::global());
    globalData.add(field);
  }

  // StructuredColumns
  atlas::functionspace::StructuredColumns fs(geom.functionSpace());

  if (geom.getComm().rank() == 0) {
    // Get grid
    atlas::StructuredGrid grid = fs.grid();

    // Get sizes
    atlas::idx_t nx = grid.nxmax();
    atlas::idx_t ny = grid.ny();

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
    for (size_t jvar = 0; jvar < varsToRead.size(); ++jvar) {
      auto varField = globalData[varsToRead[jvar].name()];
      auto varView = atlas::array::make_view<double, 2>(varField);
      for (int k = 0; k < varsToRead[jvar].getLevels(); ++k) {
        // Read data
        std::vector<double> zvar(ny*nx);
        if ((retval = nc_get_var_double(ncid, var_id[iVarLev], zvar.data()))) {
          ERR(retval, varLevName[iVarLev]);
        }
        ++iVarLev;

        // Copy data
        for (atlas::idx_t j = 0; j < ny; ++j) {
          for (atlas::idx_t i = 0; i < grid.nx(j); ++i) {
            atlas::gidx_t gidx = grid.index(i, j);
            varView(gidx, k) = zvar[j*nx+i];
          }
        }

        // Get pressure from its logarithm
        if (varsToRead[jvar].name() == "SURFPRESSION") {
          for (int jnodeGlb = 0; jnodeGlb < varField.shape(0); ++jnodeGlb) {
            varView(jnodeGlb, 0) = std::exp(varView(jnodeGlb, 0));
          }
        }
      }
    }

    // Close file
    if ((retval = nc_close(ncid))) ERR(retval, ncFilePath);
  }

  // Scatter data from main processor
  fs.scatter(globalData, fsetToRead);

  // Add fields
  for (const auto & field : fsetToRead) {
    fset.add(field);
  }

  // Code is too complicated, mark dirty to be safe
  fset.set_dirty();

  oops::Log::trace() << "quenchxx::readArome done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace quenchxx
