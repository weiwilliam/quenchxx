/*
 * (C) Copyright 2025 Meteorologisk Institutt
 *
 */

#include "quenchxx/FieldsIO/FieldsIODefault.h"

#include <vector>

#include "oops/util/FieldSetHelpers.h"
#include "oops/util/Logger.h"

namespace quenchxx {

// -----------------------------------------------------------------------------

void readDefault(const Geometry & geom,
                 const varns::Variables & vars_in_file,
                 const eckit::Configuration & config,
                 atlas::FieldSet & fset) {
  oops::Log::trace() << "quenchxx::readDefault starting" << std::endl;

  // Create variableSizes
  std::vector<size_t> variableSizes;
  for (const auto & var : vars_in_file) {
    variableSizes.push_back(var.getLevels());
  }

  // Update configuration
  eckit::LocalConfiguration conf(config);
  if (!conf.has("latitude south to north")) {
    conf.set("latitude south to north", geom.latSouthToNorth());
  }

  // Read fieldset
  util::readFieldSet(geom.getComm(),
                     geom.functionSpace(),
                     variableSizes,
                     vars_in_file.variables(),
                     conf,
                     fset);

  oops::Log::trace() << "quenchxx::readDefault done" << std::endl;
}

// -----------------------------------------------------------------------------

void writeDefault(const Geometry & geom,
                  const eckit::Configuration & config,
                  const atlas::FieldSet & fset) {
  oops::Log::trace() << "quenchxx::writeDefault starting" << std::endl;

  // Update configuration
  eckit::LocalConfiguration conf(config);
  if (!conf.has("latitude south to north")) {
    conf.set("latitude south to north", geom.latSouthToNorth());
  }

  // Write fieldset
  util::writeFieldSet(geom.getComm(), conf, fset);

  oops::Log::trace() << "quenchxx::writeDefault done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace quenchxx
