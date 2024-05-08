/*
 * (C) Copyright 2024 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>

#include "src/ModelData.h"

namespace quenchxx {

// -------------------------------------------------------------------------------------------------

class ModelData : public quench::ModelData {
  using quenchModelData = quench::ModelData;
  using quenchModelData::quenchModelData;

 public:
  static const std::string classname() {return "quenchxx::ModelData";}
};

// -------------------------------------------------------------------------------------------------

}  // namespace quenchxx
