/*
 * (C) Copyright 2022 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>

#include "src/Fields.h"

namespace quenchxx {

// -----------------------------------------------------------------------------
/// Class to represent a Fields for the  model
class Fields : public quench::Fields {
  using quenchFields = quench::Fields;
  using quenchFields::quenchFields;

 public:
  static const std::string classname() {return "quenchxx::Fields";}
};
// -----------------------------------------------------------------------------

}  // namespace quenchxx
