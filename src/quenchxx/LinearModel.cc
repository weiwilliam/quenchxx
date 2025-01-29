/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#include "oops/util/Logger.h"

#include "quenchxx/Increment.h"
#include "quenchxx/LinearModel.h"
#include "quenchxx/TraitsFwd.h"

namespace quenchxx {

// -----------------------------------------------------------------------------

static oops::LinearModelMaker<Traits, LinearModel> makerLinearModelDefault_("default");

// -----------------------------------------------------------------------------

LinearModel::LinearModel(const Geometry &,
                         const eckit::Configuration & config)
  : timeResolution_(config.getString("tstep")) {
  oops::Log::trace() << classname() << "::LinearModel starting" << std::endl;

  oops::Log::info() << "Persistance linear model" << std::endl;

  oops::Log::trace() << classname() << "::LinearModel done" << std::endl;
}

// -----------------------------------------------------------------------------

void LinearModel::stepTL(Increment & dx,
                         const ModelAuxIncrement &) const {
  oops::Log::trace() << classname() << "::stepTL starting" << std::endl;

  dx.updateTime(timeResolution_);

  oops::Log::trace() << classname() << "::stepTL done" << std::endl;
}

// -----------------------------------------------------------------------------

void LinearModel::stepAD(Increment & dx,
                         ModelAuxIncrement &) const {
  oops::Log::trace() << classname() << "::stepAD starting" << std::endl;

  dx.updateTime(-timeResolution_);

  oops::Log::trace() << classname() << "::stepAD done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace quenchxx
