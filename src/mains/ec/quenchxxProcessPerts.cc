/*
 * (C) Copyright 2024 Meteorologisk Institutt.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <limits>

#include "oops/runs/Run.h"
#include "oops/util/Logger.h"
#include "quenchxx/instantiateQuenchMatrices.h"
#include "quenchxx/Logbook.h"
#include "quenchxx/Traits.h"
#include "saber/oops/ProcessPerts.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  quenchxx::instantiateQuenchMatrices();
  saber::ProcessPerts<quenchxx::Traits> pp;
  quenchxx::Logbook::start();
  run.execute(pp);
  quenchxx::Logbook::stop();
  return 0;
}
