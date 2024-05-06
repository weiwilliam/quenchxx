/*
 * (C) Copyright 2022 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>

namespace quenchxx {

class Covariance;
class Geometry;
class Increment;
class LinearVariableChange;
class ModelData;
class State;
class VariableChange;

struct Traits {
  static std::string name() {return "quenchxx";}
  static std::string nameCovar() {return "quenchxxCovariance";}

  typedef quenchxx::Covariance           Covariance;
  typedef quenchxx::Geometry             Geometry;
  typedef quenchxx::Increment            Increment;
  typedef quenchxx::LinearVariableChange LinearVariableChange;
  typedef quenchxx::ModelData            ModelData;
  typedef quenchxx::State                State;
  typedef quenchxx::VariableChange       VariableChange;
};

}  // namespace quenchxx
