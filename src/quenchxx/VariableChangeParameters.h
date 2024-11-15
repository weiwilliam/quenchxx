/*
 * (C) Copyright 2022 UCAR.
 * (C) Copyright 2023-2024 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <map>
#include <string>
#include <vector>

#include "oops/base/VariableChangeParametersBase.h"

namespace quenchxx {

// -------------------------------------------------------------------------------------------------
/// VariableChange parameters class

class VariableChangeParameters : public oops::VariableChangeParametersBase {
  OOPS_CONCRETE_PARAMETERS(VariableChangeParameters, oops::VariableChangeParametersBase)
 public:
  oops::Parameter<std::map<std::string, std::vector<std::string>>> vaderCustomCookbook{
    "vader custom cookbook", vaderQuenchxxCustomCookbook(), this};
  oops::Parameter<vader::VaderParameters> vaderParam{"vader", {}, this};
};

// -------------------------------------------------------------------------------------------------

}  // namespace quenchxx
