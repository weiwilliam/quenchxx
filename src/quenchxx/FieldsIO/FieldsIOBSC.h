/*
 * (C) Copyright 2025 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>

#include "atlas/field.h"

#include "eckit/config/Configuration.h"

#include "quenchxx/FieldsIO/FieldsIOBase.h"
#include "quenchxx/VariablesSwitch.h"

namespace quenchxx {
  class Geometry;

// -----------------------------------------------------------------------------

class FieldsIOBSC : public FieldsIOBase {
 public:
  static const std::string classname()
    {return "quenchxx::FieldsIOBSC";}

  // Constructor/destructor
  explicit FieldsIOBSC(const std::string & ioFormat)
    : FieldsIOBase(ioFormat) {}
  ~FieldsIOBSC() = default;

  // Read
  void read(const Geometry &,
            const varns::Variables &,
            const eckit::Configuration &,
            atlas::FieldSet &) const override;

  // Write
  void write(const Geometry &,
             const eckit::Configuration &,
             const atlas::FieldSet &) const override;
};

// -----------------------------------------------------------------------------

}  // namespace quenchxx
