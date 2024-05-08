/*
 * (C) Copyright 2024 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>

#include "src/Geometry.h"

namespace quenchxx {

// -----------------------------------------------------------------------------

class Geometry : public quench::Geometry {
  using quenchGeometry = quench::Geometry;
  using quenchGeometry::quenchGeometry;

 public:
  static const std::string classname() {return "quenchxx::Geometry";}
};

// -----------------------------------------------------------------------------

}  // namespace quenchxx
