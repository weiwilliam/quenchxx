/*
 * (C) Copyright 2024 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>

#include "src/LinearVariableChange.h"

namespace quenchxx {

// -----------------------------------------------------------------------------

class LinearVariableChange: public quench::LinearVariableChange {
  using quenchLinearVariableChange = quench::LinearVariableChange;
  using quenchLinearVariableChange::quenchLinearVariableChange;

 public:
  static const std::string classname() {return "quenchxx::LinearVariableChange";}
};

// -----------------------------------------------------------------------------

}  // namespace quenchxx
