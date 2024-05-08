/*
 * (C) Copyright 2024 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>

#include "src/Covariance.h"

namespace quenchxx {

// -----------------------------------------------------------------------------

class Covariance : public quench::Covariance {
  using quenchCovariance = quench::Covariance;
  using quenchCovariance::quenchCovariance;

 public:
  static const std::string classname() {return "quenchxx::Covariance";}
};

// -----------------------------------------------------------------------------

}  // namespace quenchxx
