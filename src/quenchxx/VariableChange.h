/*
 * (C) Copyright 2017-2021 UCAR.
 * (C) Copyright 2023-2024 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "eckit/config/Configuration.h"

#include "oops/util/Printable.h"

#include "vader/vader.h"

#ifdef ECSABER
#include "quenchxx/Variables.h"
namespace varns = quenchxx;
#else
#include "oops/base/Variables.h"
namespace varns = oops;
#endif

namespace quenchxx {
  class Geometry;
  class State;

// -----------------------------------------------------------------------------

class VariableChange : public util::Printable {
 public:
  static const std::string classname()
    {return "quenchxx::VariableChange";}

  // Constructor/destructor
  VariableChange(const eckit::Configuration &,
                 const Geometry &);

  void changeVar(State &,
                 const varns::Variables &) const;
  void changeVarInverse(State &,
                        const varns::Variables &) const;

 private:
  void print(std::ostream & os) const override
    {os << *vader_;};
  std::vector<eckit::LocalConfiguration> alias_;
  std::unique_ptr<vader::Vader> vader_;
};

// -----------------------------------------------------------------------------

}  // namespace quenchxx
