/*
 * (C) Copyright 2024 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "quenchxx/ObsLocalizationEC.h"

#include <algorithm>
#include <cmath>
#include <vector>

#include "eckit/config/Configuration.h"

#include "oops/generic/gc99.h"
#include "oops/util/missingValues.h"

#include "quenchxx/GeometryIterator.h"

// -----------------------------------------------------------------------------

namespace quenchxx {

// -----------------------------------------------------------------------------

static oops::ObsLocalizationMaker<Traits, ObsLocalization>
  makerObsLocalizationEC_("default");

// -----------------------------------------------------------------------------

ObsLocalization::ObsLocalization(const eckit::Configuration & config,
                                 const ObsSpace & obsSpace)
  : locs_(obsSpace.locations()),
  locFunc_(config.getString("localization function", "gc99")),
  horScale_(config.getDouble("horizontal length-scale", 0.0)),
  verScale_(config.getDouble("vertical length-scale", 0.0)) {}

// -----------------------------------------------------------------------------

void ObsLocalization::computeLocalization(const GeometryIterator & geometryIterator,
                                          ObsVector & obsVector) const {
  oops::Log::trace() << "ObsLocalization::computeLocalization starting" << std::endl;

  // Get grid point coordinates
  eckit::geometry::Point3 gridPoint = *geometryIterator;
  const atlas::PointLonLat horGridPoint({gridPoint[0], gridPoint[1]});

  // Loop over observations
  const double missing = util::missingValue<double>();
  for (size_t jo = 0; jo < locs_.size(); ++jo) {
    // Compute normalized horizontal distance
    const atlas::PointLonLat horObsPoint({locs_[jo][0], locs_[jo][1]});
    double horDist = atlas::util::Earth().distance(horGridPoint, horObsPoint);
    if (horDist > 0.0) {
      if (horScale_ > 0.0) {
        horDist /= horScale_;
      } else {
        horDist = 1.0;
      }
    }

    // Compute normalized vertical distance
    double verDist = 0.0;
    if (geometryIterator.iteratorDimension() == 3) {
      verDist = std::abs(gridPoint[2] - locs_[jo][2]);
    }
    if (verDist > 0.0) {
      if (verScale_ > 0.0) {
        verDist /= verScale_;
      } else {
        verDist = 1.0;
      }
    }

    if ((horDist < 1.0) && (verDist < 1.0)) {
      // Compute localization as a product of horizontal and vertical components
      const double loc = locFunc(horDist)*locFunc(verDist);
      for (size_t jvar = 0; jvar < obsVector.nvars(); ++jvar) {
        if (obsVector(jvar, jo) != missing) {
          obsVector.set(jvar, jo, obsVector(jvar, jo)*loc);
        }
      }
    } else {
      for (size_t jvar = 0; jvar < obsVector.nvars(); ++jvar) {
        // Set at missing value
        obsVector.set(jvar, jo, missing);
      }
    }
  }

  oops::Log::trace() << "ObsLocalization::computeLocalization done" << std::endl;
}

// -----------------------------------------------------------------------------

double ObsLocalization::locFunc(const double & normDist) const {
  double res = 0.0;
  if (locFunc_ == "gc99") {
    res = oops::gc99(normDist);
  } else {
    throw eckit::UserError("Wrong localization function", Here());
  }

  return res;
}

// -----------------------------------------------------------------------------

void ObsLocalization::print(std::ostream & os) const {
  os << "ObsLocalization with length-scale: " << horScale_ << " / " << verScale_ << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace quenchxx
