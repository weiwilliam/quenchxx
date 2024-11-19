/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#include "quenchxx/HybridWeight.h"

#include <cmath>
#include <vector>

#include "quenchxx/Increment.h"

#include "oops/util/Logger.h"

namespace quenchxx {

// -----------------------------------------------------------------------------

HybridWeight::HybridWeight(const eckit::Configuration & config)
  : wgt_(std::sqrt(config.getDouble("weight"))) {
  oops::Log::trace() << classname() << "::HybridWeight" << std::endl;
}


// -----------------------------------------------------------------------------

void HybridWeight::multiply(Increment & dx) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;

  dx *= wgt_;

  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace quenchxx
