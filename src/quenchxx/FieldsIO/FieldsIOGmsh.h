/*
 * (C) Copyright 2025 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>
#include <vector>

#include "atlas/field.h"
#include "atlas/output/Gmsh.h"

#include "quenchxx/FieldsIO/FieldsIOBase.h"
#include "quenchxx/VariablesSwitch.h"

namespace quenchxx {
  class Geometry;

// -----------------------------------------------------------------------------

class FieldsIOGmsh : public FieldsIOBase {
 public:
  static const std::string classname()
    {return "quenchxx::FieldsIOGmsh";}

  // Constructor/destructor
  explicit FieldsIOGmsh(const std::string & ioFormat)
    : FieldsIOBase(ioFormat) {}
  ~FieldsIOGmsh() = default;

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

