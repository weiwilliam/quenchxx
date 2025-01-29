/*
 * (C) Copyright 2022 UCAR.
 * (C) Copyright 2023-2024 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>

namespace quenchxx {

class Covariance;
class Geometry;
class GeometryIterator;
class GeoVaLs;
class HorizScaleDecomposition;
class Increment;
class IncrEnsCtlVec;
class IncrModCtlVec;
class Interpolator;
class LinearVariableChange;
class LocalizationMatrix;
class Locations;
class Model;
class ModelAuxControl;
class ModelAuxControlEstimator;
class ModelAuxCovariance;
class ModelAuxCtlVec;
class ModelAuxIncrement;
class ModelData;
class ObsSpace;
class ObsVector;
class State;
class VariableChange;
#ifdef ECSABER
class Variables;
#endif

struct Traits {
  static std::string name()
    {return "quenchxx";}
  static std::string nameCovar()
    {return "quenchxxCovariance";}

  using Covariance = quenchxx::Covariance;
  using Geometry = quenchxx::Geometry;
  using GeometryIterator = quenchxx::GeometryIterator;
  using GeoVaLs = quenchxx::GeoVaLs;
  using HorizScaleDecomposition = quenchxx::HorizScaleDecomposition;
  using Increment = quenchxx::Increment;
  using IncrEnsCtlVec = quenchxx::IncrEnsCtlVec;
  using IncrModCtlVec = quenchxx::IncrModCtlVec;
  using Interpolator = quenchxx::Interpolator;
  using LinearVariableChange = quenchxx::LinearVariableChange;
  using LocalizationMatrix = quenchxx::LocalizationMatrix;
  using Locations = quenchxx::Locations;
  using Model = quenchxx::Model;
  using ModelAuxControl = quenchxx::ModelAuxControl;
  using ModelAuxControlEstimator = quenchxx::ModelAuxControlEstimator;
  using ModelAuxCovariance = quenchxx::ModelAuxCovariance;
  using ModelAuxCtlVec = quenchxx::ModelAuxCtlVec;
  using ModelAuxIncrement = quenchxx::ModelAuxIncrement;
  using ModelData = quenchxx::ModelData;
  using ObsSpace = quenchxx::ObsSpace;
  using ObsVector = quenchxx::ObsVector;
  using State = quenchxx::State;
  using VariableChange = quenchxx::VariableChange;
#ifdef ECSABER
  using Variables = quenchxx::Variables;
#endif
};

}  // namespace quenchxx
