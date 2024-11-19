/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#pragma once

#include <ostream>
#include <string>
#include <vector>

#include "eckit/exception/Exceptions.h"

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"

namespace eckit {
  class Configuration;
}

namespace quenchxx {

// -----------------------------------------------------------------------------
///  Variables class

class Variables : public oops::JediVariables,
                  private util::ObjectCounter<Variables> {
  using quenchxxVariables = oops::JediVariables;
  using quenchxxVariables::quenchxxVariables;

 public:
  static const std::string classname()
    {return "quenchxx::Variables";}

// Extra constructor
  explicit Variables(const eckit::Configuration &);

// Extra accessor
  std::vector<std::string> variablesList() const
    {return this->variables();}
};

// -----------------------------------------------------------------------------

}  // namespace quenchxx
