#include "oops/runs/Run.h"
#include "quenchxx/Traits.h"
#include "test/interface/Geometry.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  test::Geometry<quenchxx::Traits> tests;
  return run.execute(tests);
}
