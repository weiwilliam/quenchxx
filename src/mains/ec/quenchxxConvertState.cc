/*
 * (C) Copyright 2024 Meteorologisk Institutt.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/runs/ConvertState.h"
#include "oops/runs/Run.h"
#include "quenchxx/instantiateQuenchMatrices.h"
#include "quenchxx/Logbook.h"
#include "quenchxx/Traits.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  quenchxx::instantiateQuenchMatrices();
  oops::ConvertState<quenchxx::Traits> cs;
  quenchxx::Logbook::start();
  run.execute(cs);
  quenchxx::Logbook::stop();
  return 0;
}
