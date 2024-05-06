/*
 * (C) Copyright 2022 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>

#include "src/Increment.h"

namespace quenchxx {

// -----------------------------------------------------------------------------

class Increment : public quench::Increment {
  using quenchIncrement = quench::Increment;
  using quenchIncrement::quenchIncrement;

 public:
  static const std::string classname() {return "quenchxx::Increment";}
};

// -----------------------------------------------------------------------------

}  // namespace quenchxx
