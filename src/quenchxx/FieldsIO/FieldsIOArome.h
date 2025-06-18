/*
 * (C) Copyright 2025 Meteorologisk Institutt
 *
 */

#pragma once

#include <string>

#include "atlas/field.h"

#include "eckit/config/Configuration.h"

#ifdef READFA
#include "ectrans/transi.h"
#endif

#include "quenchxx/FieldsIO/FieldsIOBase.h"
#include "quenchxx/VariablesSwitch.h"

namespace quenchxx {
  class Geometry;

// -----------------------------------------------------------------------------

class FieldsIOArome : public FieldsIOBase {
 public:
  static const std::string classname()
    {return "quenchxx::FieldsIOArome";}

  // Constructor/destructor
  explicit FieldsIOArome(const std::string & ioFormat)
    : FieldsIOBase(ioFormat) {}
  ~FieldsIOArome() = default;

  // Read
  void read(const Geometry &,
            const varns::Variables &,
            const eckit::Configuration &,
            atlas::FieldSet &) const override;

 private:
#ifdef READFA
  std::unique_ptr<struct Trans_t> trans_;
#endif
};

// -----------------------------------------------------------------------------

}  // namespace quenchxx

