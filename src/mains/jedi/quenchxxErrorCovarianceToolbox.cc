/*
 * (C) Copyright 2024 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/runs/Run.h"
#include "quenchxx/Traits.h"
#include "saber/oops/ErrorCovarianceToolbox.h"
#include "saber/oops/instantiateCovarFactory.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  saber::instantiateCovarFactory<quenchxx::Traits>();
  saber::ErrorCovarianceToolbox<quenchxx::Traits> ect;
  return run.execute(ect);
}
