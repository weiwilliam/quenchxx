/*
 * (C) Copyright 2025 Meteorologisk Institutt
 *
 */

#pragma once

#include "atlas/field.h"

#include "eckit/config/Configuration.h"

#include "quenchxx/Geometry.h"
#include "quenchxx/VariablesSwitch.h"

namespace quenchxx {

// -----------------------------------------------------------------------------

void readGrib(const Geometry &,
              const varns::Variables &,
              const eckit::Configuration &,
              atlas::FieldSet &);

// -----------------------------------------------------------------------------

}  // namespace quenchxx
