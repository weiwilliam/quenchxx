/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#pragma once

#include "oops/interface/GenericMatrix.h"

#include "quenchxx/HybridWeight.h"
#include "quenchxx/TraitsFwd.h"

namespace quenchxx {

// -----------------------------------------------------------------------------

void instantiateQuenchMatrices() {
  static oops::GenericMatrixMaker<quenchxx::Traits, HybridWeight> makerEnsWgt_("hybrid_weight");
}

// -----------------------------------------------------------------------------

}  // namespace quenchxx
