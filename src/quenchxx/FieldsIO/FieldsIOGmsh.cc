/*
 * (C) Copyright 2025 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "quenchxx/FieldsIO/FieldsIOGmsh.h"

#include <vector>

#include "atlas/meshgenerator/MeshGenerator.h"

#include "eckit/exception/Exceptions.h"

#include "oops/util/FieldSetHelpers.h"
#include "oops/util/Logger.h"

#include "quenchxx/Geometry.h"

namespace quenchxx {

// -----------------------------------------------------------------------------

static FieldsIOMaker<FieldsIOGmsh> makerGmsh_("gmsh");

// -----------------------------------------------------------------------------

void FieldsIOGmsh::read(const Geometry & geom,
                        const varns::Variables & vars,
                        const eckit::Configuration & conf,
                        atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::read starting" << std::endl;

  throw eckit::NotImplemented("GMSH input not implemented yet", Here());

  oops::Log::trace() << classname() << "::read done" << std::endl;
}

// -----------------------------------------------------------------------------

void FieldsIOGmsh::write(const Geometry & geom,
                         const eckit::Configuration & conf,
                         const atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::write starting" << std::endl;

  if (!geom.mesh().generated()) {
    const atlas::MeshGenerator gen("delaunay");
    geom.mesh() = gen(geom.grid(), geom.partitioner());
  }

  // GMSH file path
  std::string gmshfilepath = conf.getString("filepath");;
  gmshfilepath.append(".msh");
  oops::Log::info() << "Info     : Writing file: " << gmshfilepath << std::endl;

  // GMSH configuration
  const auto gmshConfig =
  atlas::util::Config("coordinates", "xyz") | atlas::util::Config("ghost", true) |
  atlas::util::Config("info", true);
  atlas::output::Gmsh gmsh(gmshfilepath, gmshConfig);

  // Write GMSH
  gmsh.write(geom.mesh());
  gmsh.write(fset, geom.functionSpace());

  oops::Log::trace() << classname() << "::write done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace quenchxx
