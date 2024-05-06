/*
 * (C) Copyright 2022 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>

#include "src/Geometry.h"

namespace quenchxx {

// -----------------------------------------------------------------------------
/// Geometry handles geometry for quenchxx model.

class Geometry : public quench::Geometry {
  using quenchGeometry = quench::Geometry;
  using quenchGeometry::quenchGeometry;

 public:
  static const std::string classname() {return "quenchxx::Geometry";}
};

// -----------------------------------------------------------------------------

}  // namespace quenchxx
