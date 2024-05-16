/*
 * (C) Copyright 2024 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>

#include "oops/base/LocalIncrement.h"

#include "quenchxx/GeometryIterator.h"

#include "src/Increment.h"

namespace quenchxx {

// -----------------------------------------------------------------------------

class Increment : public quench::Increment {
  using quenchIncrement = quench::Increment;
  using quenchIncrement::quenchIncrement;

 public:
  static const std::string classname() {return "quenchxx::Increment";}

  void ones();
  oops::LocalIncrement getLocal(const GeometryIterator & geometryIterator) const;
  void setLocal(const oops::LocalIncrement & localIncrement,
                const GeometryIterator & geometryIterator);
};

// -----------------------------------------------------------------------------

}  // namespace quenchxx
