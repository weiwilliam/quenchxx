/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#pragma once

#include <memory>
#include <ostream>
#include <string>

#include "oops/interface/LinearModelBase.h"

#include "oops/util/Duration.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "quenchxx/TraitsFwd.h"

namespace eckit {
  class Configuration;
}

namespace quenchxx {
  class State;

// -----------------------------------------------------------------------------
/// LinearModel class

class LinearModel: public oops::LinearModelBase<Traits>,
           private util::ObjectCounter<LinearModel> {
 public:
  static const std::string classname()
    {return "quenchxx::LinearModel";}

/// OOPS interface

// Constructors/destructor
  LinearModel()
    {}
  LinearModel(const Geometry &,
              const eckit::Configuration &);
  ~LinearModel()
    {}

// Model trajectory computation
  void setTrajectory(State &,
                     const Model &,
                     const ModelAuxControl &)
    {}

// Run TLM and its adjoint
  void initializeTL(Increment &,
                    const ModelAuxIncrement &) const
    {}
  void stepTL(Increment &,
              const ModelAuxIncrement &) const;
  void finalizeTL(Increment &) const
    {}

  void initializeAD(Increment &,
                    const ModelAuxIncrement &) const
    {}
  void stepAD(Increment &,
              ModelAuxIncrement &) const;
  void finalizeAD(Increment &) const
    {}

// Other utilities
  const util::Duration & timeResolution() const
    {return timeResolution_;}

 private:
  void print(std::ostream &) const
    {throw eckit::NotImplemented(Here());}

  const util::Duration timeResolution_;
};
// -----------------------------------------------------------------------------

}  // namespace quenchxx
