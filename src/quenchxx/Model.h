/*
 * (C) Copyright 2024 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <ostream>
#include <string>

#include "oops/base/Variables.h"
#include "oops/interface/ModelBase.h"
#include "oops/util/Duration.h"
#include "oops/util/ObjectCounter.h"

#include "quenchxx/Traits.h"

namespace quenchxx {
  class Geometry;
  class ModelAuxControl;
  class Increment;
  class State;

// -------------------------------------------------------------------------------------------------

class ModelPseudoParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ModelPseudoParameters, Parameters)
 public:
  oops::RequiredParameter<util::Duration> tstep{ "tstep", this};
};

// -------------------------------------------------------------------------------------------------

class ModelPseudo: public oops::interface::ModelBase<Traits>,
                   private util::ObjectCounter<ModelPseudo> {
 public:
  static const std::string classname() {return "quenchxx::ModelPseudo";}

  ModelPseudo(const Geometry &,
              const eckit::Configuration &);
  ~ModelPseudo() {}

/// Prepare model integration
  void initialize(State &) const {}

/// Model integration
  void step(State &,
            const ModelAuxControl &) const;

/// Finish model integration
  void finalize(State &) const {}

/// Utilities
  const util::Duration & timeResolution() const {return tstep_;}

 private:
  void print(std::ostream & os) const {os << "ModelPseudo";}
  util::Duration tstep_;
};

// -------------------------------------------------------------------------------------------------

}  // namespace quenchxx
