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

  // Extended constructors
  Increment(const Geometry & geom, const oops::Variables & vars, const util::DateTime & vt) : quench::Increment(geom, vars, vt), geom_(new Geometry(geom)) {}
  Increment(const Geometry & geom, const Increment & other) : quench::Increment(geom, other), geom_(new Geometry(geom)) {}
  Increment(const Increment & other, const bool & copy) : quench::Increment(other, copy), geom_(other.geom_) {}

  void ones();
  oops::LocalIncrement getLocal(const GeometryIterator & geometryIterator) const;
  void setLocal(const oops::LocalIncrement & localIncrement,
                const GeometryIterator & geometryIterator);

  // Extended accessor
  std::shared_ptr<const Geometry> geometry() const {return geom_;} 

 private:
  // Extended members
  std::shared_ptr<const Geometry> geom_;
};

// -----------------------------------------------------------------------------

}  // namespace quenchxx
