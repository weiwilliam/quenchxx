/*
 * (C) Copyright 2025 Meteorologisk Institutt
 *
 */

#pragma once

#include "atlas/field.h"
#include "atlas/functionspace.h"

#include "eckit/config/Configuration.h"
#include "eckit/mpi/Comm.h"

namespace quenchxx {

// -----------------------------------------------------------------------------

extern "C" {
  void fieldsio_arome_fa_f90(const eckit::Configuration &,
                             const eckit::mpi::Comm *,
                             const atlas::functionspace::FunctionSpaceImpl *,
                             const atlas::field::FieldSetImpl *,
                             const atlas::field::FieldSetImpl *);
}

// -----------------------------------------------------------------------------

}  // namespace quenchxx
