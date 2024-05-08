/*
 * (C) Copyright 2021-2023 UCAR
 * (C) Copyright 2024 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "quenchxx/VariableChange.h"

#include <map>
#include <ostream>
#include <string>
#include <vector>

#include <typeinfo>

#include "oops/base/Variables.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/Logger.h"

#include "quenchxx/Geometry.h"
#include "quenchxx/ModelData.h"
#include "quenchxx/State.h"
#include "quenchxx/VaderCookbook.h"
#include "quenchxx/VariableChangeParameters.h"

namespace quenchxx {

// -------------------------------------------------------------------------------------------------

VariableChange::VariableChange(const eckit::Configuration & config,
                               const Geometry & geometry)
  : quench::VariableChange(config, geometry), alias_(geometry.alias()), vader_() {
  oops::Log::trace() << "VariableChange::VariableChange starting" << std::endl;

  // Deserialize configuration
  VariableChangeParameters params;
  params.deserialize(config);

  // Pass model data parameters to vader configuration
  ModelData modelData(geometry);
  eckit::LocalConfiguration vaderConfig;
  vaderConfig.set(vader::configCookbookKey,
    params.toConfiguration().getSubConfiguration("vader custom cookbook"));
  vaderConfig.set(vader::configModelVarsKey, modelData.modelData());

  // Create vader with quenchxx custom cookbook
  vader_.reset(new vader::Vader(params.vaderParam, vaderConfig));

  oops::Log::trace() << "VariableChange::VariableChange done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void VariableChange::changeVar(State & x,
                               const oops::Variables & vars_out) const {
  oops::Log::trace() << "VariableChange::changeVar starting" << std::endl;

  // State to FieldSet
  atlas::FieldSet fset;
  x.toFieldSet(fset);

  // Call vader
  oops::Variables varsCha(vars_out);
  vader_->changeVar(fset, varsCha);

  // FieldSet to State
  x.fromFieldSet(fset);

  oops::Log::trace() << "VariableChange::changeVar done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void VariableChange::changeVarInverse(State & x,
                                      const oops::Variables & vars_out) const {
  oops::Log::trace() << "VariableChange::changeVarInverse starting" << std::endl;

  // State to FieldSet
  atlas::FieldSet fset;
  x.toFieldSet(fset);

  // Call vader
  oops::Variables varsCha(vars_out);
  vader_->changeVar(fset, varsCha);

  // FieldSet to State
  x.fromFieldSet(fset);

  oops::Log::trace() << "VariableChange::changeVarInverse done" << std::endl;
}

// -------------------------------------------------------------------------------------------------


}  // namespace quenchxx
