/*
 * (C) Copyright 2025 Meteorologisk Institutt
 *
 */

#include "quenchxx/FieldsIO/FieldsIOGmsh.h"

#include <string>

#include "atlas/output/Gmsh.h"
#include "atlas/util/Config.h"

#include "eckit/exception/Exceptions.h"

#include "oops/util/Logger.h"

namespace quenchxx {

// -----------------------------------------------------------------------------

void writeGmsh(const Geometry & geom,
               const eckit::Configuration & config,
               const atlas::FieldSet & fset) {
  oops::Log::trace() << "quenchxx::writeGmsh starting" << std::endl;

  if (geom.mesh().generated()) {
    // GMSH file path
    std::string gmshFilePath = config.getString("filepath");;
    gmshFilePath.append(".msh");
    oops::Log::info() << "Info     : Writing file: " << gmshFilePath << std::endl;

    // GMSH configuration
    const auto gmshConfig =
    atlas::util::Config("coordinates", "xyz") | atlas::util::Config("ghost", true) |
    atlas::util::Config("info", true);
    atlas::output::Gmsh gmsh(gmshFilePath, gmshConfig);

     // Write GMSH
    gmsh.write(geom.mesh());
    gmsh.write(fset, fset[0].functionspace());
  } else {
    throw eckit::Exception("mesh should be generated for GMSH output", Here());
  }

  oops::Log::trace() << "quenchxx::writeGmsh done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace quenchxx
