/*
 * (C) Copyright 2024 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "quenchxx/Model.h"

#include <vector>

#include "oops/mpi/mpi.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "quenchxx/State.h"

namespace quenchxx {

// -------------------------------------------------------------------------------------------------

static oops::interface::ModelMaker<Traits, ModelPseudo> makermodel_("pseudo");

// -------------------------------------------------------------------------------------------------

ModelPseudo::ModelPseudo(const Geometry & resol,
                         const eckit::Configuration & config)
  : tstep_(0) {
  oops::Log::trace() << "ModelPseudo::ModelPseudo starting" << std::endl;

  ModelPseudoParameters params;
  params.deserialize(config);

  tstep_ = util::Duration(config.getString("tstep"));

  oops::Log::trace() << "ModelPseudo::ModelPseudo done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void ModelPseudo::step(State & xx,
                       const ModelAuxControl &) const {
  // Update validTime
  xx.updateTime(tstep_);

  // Do nothing and print message
  if (oops::mpi::world().rank() == 0) {
    oops::Log::warning() << "Pseudo model: just ticking the clock." << std::endl;
  }

  oops::Log::trace() << "ModelPseudo::step" << xx.validTime() << std::endl;
}


// -------------------------------------------------------------------------------------------------

}  // namespace quenchxx
