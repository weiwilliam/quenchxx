/*
 * (C) Copyright 2022 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>

#include "src/State.h"

namespace quenchxx {

// -----------------------------------------------------------------------------

class State : public quench::State {
  using quenchState = quench::State;
  using quenchState::quenchState;

 public:
  static const std::string classname() {return "quenchxx::State";}
};

// -----------------------------------------------------------------------------

}  // namespace quenchxx
