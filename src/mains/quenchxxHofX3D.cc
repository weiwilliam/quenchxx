

#include "oops/runs/HofX3D.h"
#include "oops/runs/Run.h"
#include "quenchxx/Traits.h"
#include "ufo/instantiateObsErrorFactory.h"
#include "ufo/instantiateObsFilterFactory.h"
#include "ufo/ObsTraits.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  ufo::instantiateObsErrorFactory();
  ufo::instantiateObsFilterFactory();
  oops::HofX3D<quenchxx::Traits, ufo::ObsTraits> hofx;
  return run.execute(hofx);
}
