/*
 * (C) Copyright 2025 Meteorologisk Institutt
 *
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

class FieldsIOGrib : public FieldsIOBase {
 public:
  static const std::string classname()
    {return "quenchxx::FieldsIOGrib";}

  // Constructor/destructor
  explicit FieldsIOGrib(const std::string & ioFormat)
    : FieldsIOBase(ioFormat) {}
  ~FieldsIOGrib() = default;

  // Read
  void read(const Geometry &,
            const varns::Variables &,
            const eckit::Configuration &,
            atlas::FieldSet &) const override;
};

// -----------------------------------------------------------------------------

}  // namespace quenchxx

