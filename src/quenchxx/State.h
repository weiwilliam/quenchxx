/*
 * (C) Copyright 2022 UCAR.
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

#include "oops/util/DateTime.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "quenchxx/Fields.h"
#include "quenchxx/GeoVaLs.h"
#include "quenchxx/LinearModel.h"
#include "quenchxx/Locations.h"
#include "quenchxx/Model.h"
#include "quenchxx/VariablesSwitch.h"

namespace eckit {
  class Configuration;
}

namespace quenchxx {
  class Geometry;
  class Increment;

// -----------------------------------------------------------------------------
/// State class

class State : public util::Printable,
              private util::ObjectCounter<State> {
 public:
  static const std::string classname()
    {return "quenchxx::State";}

  // Constructors
  State(const Geometry &,
        const varns::Variables &,
        const util::DateTime &);
  State(const Geometry &,
        const eckit::Configuration &);
  State(const Geometry & resol,
        const State & other)
    : fields_(new Fields(*other.fields_, resol)) {}
  State(const varns::Variables & vars,
        const State & other)
    : fields_(new Fields(*other.fields_)) {}
  State(const State & other)
    : fields_(new Fields(*other.fields_)) {}
  State(const Geometry & resol,
        const Model &,
        const eckit::Configuration & conf)
    : State(resol, conf) {}
  State(const Geometry & resol,
        const LinearModel &,
        const eckit::Configuration & conf)
    : State(resol, conf) {}
  State(const Geometry & resol,
        const Model &,
        const State & other)
    : State(resol, other) {}
  State(const Geometry & resol,
        const Model &,
        const State & other,
        const eckit::Configuration &)
    : State(resol, other) {}

  // Assignment
  State & operator=(const State &);

  // Interactions with Increment
  State & operator+=(const Increment &);

  // I/O and diagnostics
  void read(const eckit::Configuration & config)
    {fields_->read(config);}
  void write(const eckit::Configuration & config) const
    {fields_->write(config);}
  double norm() const
    {return fields_->norm();}
  const util::DateTime & validTime() const
    {return fields_->time();}
  util::DateTime & validTime()
    {return fields_->time();}
  void updateTime(const util::Duration & dt)
    {fields_->time() += dt;}

  // Access to fields
  Fields & fields()
    {return *fields_;}
  const Fields & fields() const
    {return *fields_;}
  std::shared_ptr<const Geometry> geometry() const
    {return fields_->geometry();}

  // ATLAS FieldSet accessor
  void toFieldSet(atlas::FieldSet & fset) const
    {fields_->toFieldSet(fset);}
  void fromFieldSet(const atlas::FieldSet & fset)
    {fields_->fromFieldSet(fset);}
  const atlas::FieldSet & fieldSet() const
    {return fields_->fieldSet();}
  atlas::FieldSet & fieldSet()
    {return fields_->fieldSet();}
  void synchronizeFields()
    {fields_->synchronizeFields();}

  // Other
  void zero()
    {fields_->zero();}
  void accumul(const double & zz,
               const State & xx)
    {fields_->axpy(zz, xx.fields());}
  const varns::Variables & variables() const
    {return fields_->variables();}
  void interpolate(const Locations & locs,
                   GeoVaLs & gv) const
    {fields_->interpolate(locs, gv);}
  void forceWith(const State & other,
                 const varns::Variables & vars)
    {fields_->forceWith(*(other.fields_), vars);}

  // Serialization
  size_t serialSize() const
    {return fields_->serialSize();}
  void serialize(std::vector<double> & vect) const
    {fields_->serialize(vect);}
  void deserialize(const std::vector<double> & vect,
                   size_t & index)
    {fields_->deserialize(vect, index);}
  void transpose(const State & /*FCState*/, const eckit::mpi::Comm & /*global*/,
     const int /*ensNum*/, const int /*transNum*/) {
     throw eckit::NotImplemented("quenchxx::State::transpose not implemented", Here());
  }

 private:
  // Print
  void print(std::ostream &) const;

  // Fields
  std::unique_ptr<Fields> fields_;
};

// -----------------------------------------------------------------------------

}  // namespace quenchxx
