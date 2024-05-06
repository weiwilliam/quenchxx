/*
 * (C) Copyright 2021-2023 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <map>
#include <ostream>
#include <string>
#include <vector>

#include <typeinfo>

#include "oops/base/VariableChangeParametersBase.h"
#include "oops/mpi/mpi.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "oops/util/FieldSetHelpers.h"

#include "src/Geometry.h"
#include "src/ModelData.h"
#include "src/State.h"
#include "src/Constants.h"
#include "src/VariableChange.h"

namespace genint {

// -------------------------------------------------------------------------------------------------

VariableChange::VariableChange(const eckit::Configuration & config, const Geometry & geometry)
//VariableChange::VariableChange(const Parameters_ & params, const Geometry & geometry)
    : mapVariables_(geometry.mapVariables()), vader_() {
//  : mapVariables_(geometry.mapVariables()), inputParam_(), vader_() {

  VariableChangeParameters params;
  params.deserialize(config);
  eckit::LocalConfiguration variableChangeConfig = params.toConfiguration();
  ModelData modelData{geometry};
  eckit::LocalConfiguration vaderConfig;
  vaderConfig.set(vader::configCookbookKey,
                                variableChangeConfig.getSubConfiguration("vader custom cookbook"));
  vaderConfig.set(vader::configModelVarsKey, modelData.modelData());

  // Create vader with genint custom cookbook
  vader_.reset(new vader::Vader(params.vaderParam,
                                vaderConfig));
  // Create the variable change
  //variableChange_.reset(VariableChangeFactory::create(geometry,
  //                                                    params.variableChangeParameters.value()));
}

// -------------------------------------------------------------------------------------------------

VariableChange::~VariableChange() {}

// -------------------------------------------------------------------------------------------------

void VariableChange::changeVar(State & x, const oops::Variables & vars_out) const {
  // Trace
  oops::Log::trace() << "VariableChange::changeVar starting" << std::endl;

  // Needed Variables and fieldsets copies
  oops::Variables varsCha = vars_out;
  oops::Variables varsState =  x.variables();
  oops::Variables varsAdd = x.variables();
  atlas::FieldSet xfs;
  x.toFieldSet(xfs);

  // Convert to jedi names using geometry map for variables.
  // *FIXME* due to current atlas bug the fiedset rename method
  // is not updating the index map, for this reason we need to create
  // a new fieldset xfsVarder to get the index map updated
  // https://github.com/ecmwf/atlas/issues/147
  std::map<std::string,std::string> mapVars = mapVariables_;
  const std::vector<std::string>& varsVec = xfs.field_names();
  for (auto &var : varsVec) {
    xfs.field(var).rename(mapVars[var]);
  }

  // Call vader and get the out variables names
  varsAdd += vader_->changeVar(xfs, varsCha);
  varsAdd -= varsState;

  // Create and update the output fieldset
  // *FIXME* this step is also necessary because of rename bug
  // in atlas. Once fixed state should use rename and jedi var names
  // should be used in the entire application
  atlas::FieldSet xfsOut;
  x.toFieldSet(xfsOut);
  util::removeFieldsFromFieldSet(xfsOut, varsAdd.variables());
  for (auto &var : varsAdd.variables()) {
    xfsOut.add(xfs.field(var));
  }
  x.fromFieldSet(xfsOut);

  oops::Log::trace() << "VariableChange::changeVar done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void VariableChange::changeVarInverse(State & x, const oops::Variables & vars_out) const {
  // Trace
  oops::Log::trace() << "VariableChange::changeVarInverse starting" << std::endl;

  // Needed Variables and fieldsets copies
  oops::Variables varsCha = vars_out;
  oops::Variables varsState =  x.variables();
  oops::Variables varsAdd = x.variables();
  atlas::FieldSet xfs;
  x.toFieldSet(xfs);

  // Convert to jedi names using geometry map for variables.
  // *FIXME* due to current atlas bug the fiedset rename method
  // is not updating the index map, for this reason we need to create
  // a new fieldset xfsVarder to get the index map updated
  // https://github.com/ecmwf/atlas/issues/147
  std::map<std::string,std::string> mapVars = mapVariables_;
  const std::vector<std::string>& varsVec = xfs.field_names();
  for (auto &var : varsVec) {
    xfs.field(var).rename(mapVars[var]);
  }

  // Call vader and get the out variables names
  varsAdd += vader_->changeVar(xfs, varsCha);
  varsAdd -= varsState;

  // Create and update the output fieldset
  // *FIXME* this step is also necessary because of rename bug
  // in atlas. Once fixed state should use rename and jedi var names
  // should be used in the entire application
  atlas::FieldSet xfsOut;
  x.toFieldSet(xfsOut);
  util::removeFieldsFromFieldSet(xfsOut, varsAdd.variables());
  for (auto &var : varsAdd.variables()) {
    xfsOut.add(xfs.field(var));
  }
  x.fromFieldSet(xfsOut);

  oops::Log::trace() << "VariableChange::changeVarInverse done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void VariableChange::print(std::ostream & os) const {
  os << *vader_;
}

// -------------------------------------------------------------------------------------------------


}  // namespace genint
