/*
 * (C) Copyright 2022 UCAR
 * (C) Copyright 2024 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <map>
#include <string>
#include <vector>

#include "src/VariableChangeParameters.h"

namespace quenchxx {

// -------------------------------------------------------------------------------------------------

class VariableChangeParameters : public quench::VariableChangeParameters {
  OOPS_CONCRETE_PARAMETERS(VariableChangeParameters, quench::VariableChangeParameters)
 public:
  oops::Parameter<std::map<std::string, std::vector<std::string>>> vaderCustomCookbook{
    "vader custom cookbook", vaderQuenchxxCustomCookbook(), this};
  oops::Parameter<vader::VaderParameters> vaderParam{"vader", {}, this};
};

// -------------------------------------------------------------------------------------------------

}  // namespace quenchxx
