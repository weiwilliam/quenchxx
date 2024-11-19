/*
 * (C) Copyright 2023 Meteorologisk Institutt
 *
 */

#include "quenchxx/Locations.h"

#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"

#include "quenchxx/ObsSpace.h"

namespace quenchxx {

// -----------------------------------------------------------------------------

Locations::Locations(const ObsSpace & obsSpace,
                     const util::DateTime & t1,
                     const util::DateTime & t2)
  : obsSpace_(obsSpace) {
  oops::Log::trace() << classname() << "::Locations starting" << std::endl;

  // Get observations coordinates
  locs_ = obsSpace_.locations(t1, t2);

  // Get vector of observation lon/lat
  std::vector<double> obsLonLatLoc;
  for (size_t jo = 0; jo < locs_.size(); ++jo) {
    obsLonLatLoc.push_back(locs_[jo][0]);
    obsLonLatLoc.push_back(locs_[jo][1]);
  }

  // Gather number of observations
  const int nobsLoc = locs_.size();
  nobsLocVec_.resize(obsSpace_.getComm().size());
  obsSpace_.getComm().allGather(nobsLoc, nobsLocVec_.begin(), nobsLocVec_.end());
  size_t nobsGlb = 0;
  for (size_t jt = 0; jt < obsSpace_.getComm().size(); ++jt) {
    nobsGlb += nobsLocVec_[jt];
  }

  // Allocation
  std::vector<double> obsLonLatGlb(2*nobsGlb);

  // Define counts and displacements
  std::vector<int> counts(obsSpace_.getComm().size());
  for (size_t jt = 0; jt < obsSpace_.getComm().size(); ++jt) {
    counts[jt] = 2*nobsLocVec_[jt];
  }
  std::vector<int> displs;
  displs.push_back(0);
  for (size_t jt = 0; jt < obsSpace_.getComm().size()-1; ++jt) {
    displs.push_back(displs[jt]+counts[jt]);
  }

  // AllGatherv lon/lat
  obsSpace_.getComm().allGatherv(obsLonLatLoc.begin(), obsLonLatLoc.end(), obsLonLatGlb.begin(),
    counts.data(), displs.data());

  // Create observation grid
  eckit::LocalConfiguration gridConf;
  gridConf.set("type", "unstructured");
  gridConf.set("xy", obsLonLatGlb);
  grid_ = atlas::Grid(gridConf);

  oops::Log::trace() << classname() << "::Locations done" << std::endl;
}

// -----------------------------------------------------------------------------

void Locations::print(std::ostream & os) const {
  oops::Log::trace() << classname() << "::print starting" << std::endl;

  os << "Locations::print not implemented";

  oops::Log::trace() << classname() << "::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace quenchxx
