/*
 * (C) Copyright 2024 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/generic/instantiateModelFactory.h"
#include "oops/runs/Run.h"
#include "oops/runs/Variational.h"
#include "quenchxx/Traits.h"
#include "saber/oops/instantiateCovarFactory.h"
#include "ufo/instantiateObsErrorFactory.h"
#include "ufo/instantiateObsFilterFactory.h"
#include "ufo/ObsTraits.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  saber::instantiateCovarFactory<quenchxx::Traits>();
  ufo::instantiateObsErrorFactory();
  ufo::instantiateObsFilterFactory();
  oops::instantiateModelFactory<quenchxx::Traits>();
  oops::Variational<quenchxx::Traits, ufo::ObsTraits> var;
  return run.execute(var);
}
