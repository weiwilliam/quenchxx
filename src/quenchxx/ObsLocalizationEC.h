/*
 * (C) Copyright 2024 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <ostream>
#include <string>
#include <vector>

#include "eckit/config/Configuration.h"
#include "eckit/geometry/Point3.h"

#include "oops/base/ObsLocalizationBase.h"

#include "quenchxx/ObsSpace.h"
#include "quenchxx/ObsVector.h"
#include "quenchxx/Traits.h"

namespace quenchxx {
  class GeometryIterator;

// -----------------------------------------------------------------------------

class ObsLocalization: public oops::ObsLocalizationBase<Traits> {
 public:
  ObsLocalization(const eckit::Configuration &,
                  const ObsSpace &);

 protected:
  void computeLocalization(const GeometryIterator &,
                           ObsVector &) const override;

 private:
  void print(std::ostream &) const override;

  // Localization function
  double locFunc(const double &) const;

  // Observations coordinates
  std::vector<atlas::Point3> locs_;

  // Localization function and scales
  const std::string locFunc_;
  const double horScale_;
  const double verScale_;
};

// -----------------------------------------------------------------------------

}  // namespace quenchxx
