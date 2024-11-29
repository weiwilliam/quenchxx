/*
 * (C) Copyright 2024 Meteorologisk Institutt.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/assimilation/instantiateEvilMinFactory.h"
#include "oops/runs/Run.h"
#include "oops/runs/Variational.h"
#include "quenchxx/instantiateQuenchMatrices.h"
#include "quenchxx/Logbook.h"
#include "quenchxx/Traits.h"
#include "saber/oops/instantiateCovarFactory.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  oops::instantiateEvilMinFactory<quenchxx::Traits>();
  saber::instantiateCovarFactory<quenchxx::Traits>();
  quenchxx::instantiateQuenchMatrices();
  oops::Variational<quenchxx::Traits> var;
  quenchxx::Logbook::start();
  run.execute(var);
  quenchxx::Logbook::stop();
  return 0;
}
