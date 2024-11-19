/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#include "quenchxx/LinearObsOperator.h"

#include <string>

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"

#include "quenchxx/GeoVaLs.h"
#include "quenchxx/ObsAuxIncrement.h"
#include "quenchxx/ObsSpace.h"
#include "quenchxx/ObsVector.h"
#include "quenchxx/Variables.h"

namespace quenchxx {

// -----------------------------------------------------------------------------

static oops::LinearObsOpMaker<Traits, LinearObsOperator> makerLinearObsOpDefault_("default");

// -----------------------------------------------------------------------------

LinearObsOperator::LinearObsOperator(const ObsSpace &,
                                     const eckit::Configuration & config)
  : inputs_(new Variables(config)) {
  oops::Log::trace() << classname() << "::LinearObsOperator" << std::endl;
}

// -----------------------------------------------------------------------------

void LinearObsOperator::obsEquivTL(const GeoVaLs & gv,
                                   ObsVector & ovec,
                                   const ObsAuxIncrementPtrMap_ & bias) const {
  oops::Log::trace() << classname() << "::obsEquivTL starting" << std::endl;

  // Check number of GeoVaLs variables
  ASSERT(gv.fieldSet().size() == 1);

  // Get GeoVaLs view
  const auto gvField = gv.fieldSet()[0];
  const auto gvView = atlas::array::make_view<double, 1>(gvField);

  // Get bias
  double bias_ = 0.0;  // TODO(Benjamin): bias should also be an atlas fieldset
  using icst_ = typename ObsAuxIncrementPtrMap_::const_iterator;
  icst_ it = bias.find("ObsAuxControl");
  if (it != bias.end()) {
    const ObsAuxIncrement * pbias = dynamic_cast<const ObsAuxIncrement*>(it->second.get());
    ASSERT(pbias != nullptr);
    bias_ = (*pbias).value();
  }

  // Compute observation equivalent
  for (int jo = 0; jo < gvField.shape(0); ++jo) {
    const int ii = gv.obsIndex(jo);
    ovec(ii) = gvView(jo)+bias_;
  }

  oops::Log::trace() << classname() << "::obsEquivTL done" << std::endl;
}

// -----------------------------------------------------------------------------

void LinearObsOperator::obsEquivAD(GeoVaLs & gv,
                                   const ObsVector & ovec,
                                   ObsAuxIncrementPtrMap_ & bias) const {
  oops::Log::trace() << classname() << "::obsEquivAD starting" << std::endl;

  // Check number of GeoVaLs variables
  ASSERT(gv.fieldSet().size() == 1);

  // Get GeoVaLs view
  auto gvField = gv.fieldSet()[0];
  auto gvView = atlas::array::make_view<double, 1>(gvField);

  // Get bias
  using iter_ = typename ObsAuxIncrementPtrMap_::iterator;
  iter_ it = bias.find("ObsAuxControl");
  if (it != bias.end()) {
    ObsAuxIncrement * pbias = dynamic_cast<ObsAuxIncrement*>(it->second.get());
    ASSERT(pbias != nullptr);
    for (int jo = 0; jo < gvField.shape(0); ++jo) {
      const int ii = gv.obsIndex(jo);
      (*pbias).value() += ovec(ii);
    }
  }

  // Compute observation equivalent
  for (int jo = 0; jo < gvField.shape(0); ++jo) {
    const int ii = gv.obsIndex(jo);
    gvView(jo) = ovec(ii);
  }

  oops::Log::trace() << classname() << "::obsEquivAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void LinearObsOperator::print(std::ostream & os) const {
  oops::Log::trace() << classname() << "::print starting" << std::endl;

  os << "LinearObsOperator: quenchxx LinearObsOperator";

  oops::Log::trace() << classname() << "::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace quenchxx
