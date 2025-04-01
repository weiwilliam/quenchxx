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

void writeGmsh(const Geometry &,
               const eckit::Configuration &,
               const atlas::FieldSet &);

// -----------------------------------------------------------------------------

}  // namespace quenchxx
