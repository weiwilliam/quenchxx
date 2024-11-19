/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#include "quenchxx/ObsAuxCtlVec.h"

#include "oops/interface/ObsAuxCovarianceBase.h"

#include "oops/util/Logger.h"

namespace quenchxx {

// -----------------------------------------------------------------------------

static oops::ObsAuxCtlVecMaker<Traits, ObsAuxCtlVec> makeObsAuxCtlVecDefault_("default");

// -----------------------------------------------------------------------------

}  // namespace quenchxx
