/*
 * (C) Copyright 2024 Meteorologisk Institutt
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "quenchxx/Geometry.h"

#include <string>
#include <vector>

#include "quenchxx/GeometryIterator.h"

// -----------------------------------------------------------------------------

namespace quenchxx {

// -----------------------------------------------------------------------------

Geometry::Geometry(const eckit::Configuration & config,
                   const eckit::mpi::Comm & comm)
  : quench::Geometry(config, comm)
{
  // Iterator dimension
  iteratorDimension_ = config.getInt("iterator dimension", 2);
  ASSERT((iteratorDimension_ == 2) || (iteratorDimension_ == 3));

  // Domain size
  nnodes_ = fields().field("vert_coord_0").shape(0);
  nlevs_ = fields().field("vert_coord_0").shape(1);

  // Averaged vertical coordinate
  const auto ghostView = atlas::array::make_view<int, 1>(functionSpace().ghost());
  const auto vert_coordView = atlas::array::make_view<double, 2>(fields().field("vert_coord_0"));
  for (atlas::idx_t jlevel = 0; jlevel < nlevs_; ++jlevel) {
    double vert_coord_avg = 0.0;
    double counter = 0.0;
    for (atlas::idx_t jnode = 0; jnode < nnodes_; ++jnode) {
      if (ghostView(jnode) == 0) {
        vert_coord_avg += vert_coordView(jnode, jlevel);
        counter += 1.0;
      }
    }
    comm.allReduceInPlace(vert_coord_avg, eckit::mpi::sum());
    comm.allReduceInPlace(counter, eckit::mpi::sum());
    if (counter > 0.0) {
      vert_coord_avg /= counter;
    }
    vert_coord_avg_.push_back(vert_coord_avg);
  }  
}

// -----------------------------------------------------------------------------

Geometry::Geometry(const Geometry & other)
  : quench::Geometry(other), iteratorDimension_(other.iteratorDimension_), nnodes_(other.nnodes_),
  nlevs_(other.nlevs_), vert_coord_avg_(other.vert_coord_avg_) 
{}

// -----------------------------------------------------------------------------

GeometryIterator Geometry::begin() const {
  return GeometryIterator(*this, 0, 0);
}

// -----------------------------------------------------------------------------

GeometryIterator Geometry::end() const {
  return GeometryIterator(*this, nnodes_, nlevs_);
}

// -----------------------------------------------------------------------------

std::vector<double> Geometry::verticalCoord(std::string & vcUnits) const {
  return vert_coord_avg_;
}

// -----------------------------------------------------------------------------

}  // namespace quenchxx

