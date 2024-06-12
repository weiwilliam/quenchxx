/*
 * (C) Copyright 2024 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>
#include <vector>

#include "src/Geometry.h"

namespace quenchxx {
  class GeometryIterator;

// -----------------------------------------------------------------------------

class Geometry : public quench::Geometry {
  using quenchGeometry = quench::Geometry;
  using quenchGeometry::quenchGeometry;

 public:
  static const std::string classname() {return "quenchxx::Geometry";}

  Geometry(const eckit::Configuration &,
           const eckit::mpi::Comm & comm = oops::mpi::world());
  Geometry(const Geometry &);

  GeometryIterator begin() const;
  GeometryIterator end() const;
  std::vector<double> verticalCoord(std::string &) const;

  const size_t & iteratorDimension() const {return iteratorDimension_;}
  const size_t & nnodes() const {return nnodes_;}
  const size_t & nlevs() const {return nlevs_;}

 private:
  size_t iteratorDimension_;
  size_t nnodes_;
  size_t nlevs_;
  std::vector<double> vert_coord_avg_;
};

// -----------------------------------------------------------------------------

}  // namespace quenchxx
