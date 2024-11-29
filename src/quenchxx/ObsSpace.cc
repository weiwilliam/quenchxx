/*
 * (C) Copyright 2023 Meteorologisk Institutt
 *
 */

#include "quenchxx/ObsSpace.h"

#include <netcdf.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include "atlas/functionspace.h"
#include "atlas/grid/Grid.h"
#include "atlas/util/KDTree.h"

#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/mpi/mpi.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"

#include "quenchxx/Geometry.h"
#include "quenchxx/GeoVaLs.h"
#include "quenchxx/ObsVector.h"

#define ERR(e) {throw eckit::Exception(nc_strerror(e), Here());}

namespace quenchxx {

// -----------------------------------------------------------------------------

ObsSpace::ObsSpace(const eckit::Configuration & config,
                   const Geometry & geom,
                   const util::DateTime & bgn,
                   const util::DateTime & end,
                   const bool lscreened)
  : winbgn_(bgn), winend_(end), lscreened_(lscreened), comm_(geom.getComm()),
    geom_(new Geometry(geom)), nobsOwn_(0), nobsLoc_(0), nobsGlb_(0),
    vars_(config.getStringVector("variables")) {
  oops::Log::trace() << classname() << "::ObsSpace starting" << std::endl;

  nameIn_.clear();
  nameOut_.clear();
  if (config.has("ObsData")) {
    const eckit::LocalConfiguration dataConfig(config, "ObsData");
    if (dataConfig.has("distribution")) {
      distribution_ = dataConfig.getSubConfiguration("distribution");
      if (config.has("obs localizations")) {
        const std::vector<eckit::LocalConfiguration> configs =
          config.getSubConfigurations("obs localizations");
        ASSERT(configs.size() == 1);  // No use for more than one localization so far
        const double horLengthScale = configs[0].getDouble("horizontal length-scale");
        distribution_.set("horizontal length-scale", horLengthScale);
      }
    }
    if (lscreened) {
      if (dataConfig.has("ObsDataInScreened")) {
        nameIn_ = dataConfig.getString("ObsDataInScreened.filepath");
        oops::Log::trace() << classname() << "::ObsSpace reading screened observations from "
          << nameIn_ << std::endl;
        read(nameIn_);
      }
      if (dataConfig.has("ObsDataOutScreened")) {
        nameOut_ = dataConfig.getString("ObsDataOutScreened.filepath");
      }
    } else {
      if (dataConfig.has("ObsDataIn")) {
        nameIn_ = dataConfig.getString("ObsDataIn.filepath");
        oops::Log::trace() << classname() << "::ObsSpace reading observations from " << nameIn_
          << std::endl;
        read(nameIn_);
      }
      if (dataConfig.has("ObsDataOut")) {
        nameOut_ = dataConfig.getString("ObsDataOut.filepath");
      }
    }
  }

  oops::Log::trace() << classname() << "::ObsSpace done" << std::endl;
}

// -----------------------------------------------------------------------------

ObsSpace::~ObsSpace() {
  oops::Log::trace() << classname() << "::~ObsSpace starting" << std::endl;

  if (!nameOut_.empty()) {
    oops::Log::trace() << classname() << "::~ObsSpace saving nameOut = " << nameOut_ << std::endl;
    write(nameOut_, Write::Original);
  }

  oops::Log::trace() << classname() << "::~ObsSpace done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsSpace::putdb(const atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::putdb starting" << std::endl;

  // Check Fields size
  for (const auto & field : fset) {
    ASSERT(field.shape(0) == nobsLoc_);
    ASSERT(field.shape(1) == 1);
  }

  bool existingFieldSet = false;
  for (auto & dataFset : data_) {
    if (dataFset.name() == fset.name()) {
      // Existing FieldSet, copy fields
      for (auto & dataField : dataFset) {
        const atlas::Field field = fset[dataField.name()];
        auto dataView = atlas::array::make_view<double, 2>(dataField);
        const auto view = atlas::array::make_view<double, 2>(field);
        dataView.assign(view);
      }
      existingFieldSet = true;
    }
  }

  if (!existingFieldSet) {
    // FieldSet name not found, inserting new FieldSet
    atlas::FieldSet dataFset;
    util::copyFieldSet(fset, dataFset);
    data_.push_back(dataFset);
  }

  oops::Log::trace() << classname() << "::putdb done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsSpace::getdb(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::getdb starting" << std::endl;

  bool existingFieldSet = false;
  for (const auto & dataFset : data_) {
    if (dataFset.name() == fset.name()) {
      // FieldSet found, copy fields
      fset.clear();
      util::copyFieldSet(dataFset, fset);
      existingFieldSet = true;
    }
  }

  if (!existingFieldSet) {
    // FieldSet name not found
    std::string message = "Cannot find group " + fset.name() + " in observation database,"
      + "existing groups are: ";
    for (const auto & dataFset : data_) {
      message += dataFset.name() + " ";
    }
    throw eckit::Exception(message, Here());
  }

  oops::Log::trace() << classname() << "::getdb done" << std::endl;
}

// -----------------------------------------------------------------------------

std::vector<atlas::Point3> ObsSpace::locations(const util::DateTime & t1,
                                               const util::DateTime & t2) const {
  oops::Log::trace() << classname() << "::locations starting" << std::endl;

  std::vector<size_t> olist = timeSelect(t1, t2);
  const size_t nobs = olist.size();
  std::vector<atlas::Point3> locs(nobs);
  for (size_t jo = 0; jo < nobs; ++jo) {
    locs[jo] = locs_[olist[jo]];
  }

  oops::Log::trace() << classname() << "::locations done" << std::endl;
  return locs;
}

// -----------------------------------------------------------------------------

std::vector<size_t> ObsSpace::timeSelect(const util::DateTime & t1,
                                         const util::DateTime & t2) const {
  oops::Log::trace() << classname() << "::timeSelect starting" << std::endl;

  std::vector<size_t> mask;
  for (size_t jo = 0; jo < nobsLoc_; ++jo) {
    if (t1 == t2) {
      if (times_[jo] == t1) {
        mask.push_back(jo);
      }
    } else {
      if (times_[jo] > t1 && times_[jo] <= t2) {
        mask.push_back(jo);
      }
    }
  }

  oops::Log::trace() << classname() << "::timeSelect done" << std::endl;
  return mask;
}

// -----------------------------------------------------------------------------

void ObsSpace::generateDistribution(const eckit::Configuration & config) {
  oops::Log::trace() << classname() << "::generateDistribution starting" << std::endl;

  // TODO(Benjamin): same as JEDI

  // No forecast: assimilation window start time must be equal to final time
  ASSERT(winend_ == winbgn_);

  // Parameters
  nobsGlb_ = config.getInt("density");
  const std::vector<double> vert_coord_avg = geom_->vert_coord_avg(config.getString("variable"));
  const double bottom = *std::min_element(vert_coord_avg.begin(), vert_coord_avg.end());
  const double top = *std::max_element(vert_coord_avg.begin(), vert_coord_avg.end());

  // Random number generator
  static std::mt19937 generator(winbgn_.timestamp());
  static std::uniform_real_distribution<double> randomDistrib(0.0, 1.0);

  // Global vectors
  std::vector<double> locsGlb;

  // Define source grid horizontal distribution
  const atlas::grid::Distribution distribution = geom_->partitioner().partition(geom_->grid());

  // Local number of observations
  nobsOwnVec_.resize(comm_.size());

  if (comm_.rank() == 0) {
    // Resize vectors
    std::vector<double> locsGlbTmp(3*nobsGlb_);
    std::vector<size_t> orderGlbTmp(nobsGlb_);

    // Generate locations
    size_t iobs = 0;
    while (iobs < nobsGlb_) {
      // Generate location
      const double lon = 360.0*randomDistrib(generator);
      const double lat = 90.0-std::acos(1.0-2.0*randomDistrib(generator))*180.0/M_PI;

      // Check point validity
      bool validPoint = true;
      if (!geom_->grid().domain().global()) {
        const atlas::functionspace::StructuredColumns fs(geom_->functionSpace());
        const atlas::RegularGrid grid(fs.grid());
        atlas::PointLonLat p({lon, lat});
        grid.projection().lonlat2xy(p);
        validPoint = grid.domain().contains(p);
      }

      // Process valid point
      if (validPoint) {
        locsGlbTmp[3*iobs] = lon;
        locsGlbTmp[3*iobs+1] = lat;
        locsGlbTmp[3*iobs+2] = bottom+(top-bottom)*randomDistrib(generator);
        orderGlbTmp[iobs] = iobs;
        ++iobs;
      }
    }

    // Split observations between tasks
    std::vector<int> partition(nobsGlb_);
    if (true) {  // TODO(Benjamin): other distributions
      // Using nearest neighbor
      atlas::util::IndexKDTree search;
      search.reserve(geom_->grid().size());
      size_t jnode = 0;
      for (const auto & lonLat : geom_->grid().lonlat()) {
        atlas::PointLonLat pointLonLat(lonLat);
        pointLonLat.normalise();
        atlas::PointXY point(pointLonLat);
        search.insert(point, jnode);
        ++jnode;
      }
      search.build();

      for (size_t jo = 0; jo < nobsGlb_; ++jo) {
        // Find MPI task
        atlas::PointLonLat pointLonLat(locsGlbTmp[3*jo], locsGlbTmp[3*jo+1]);
        pointLonLat.normalise();

        // Search nearest neighbor
        atlas::util::IndexKDTree::ValueList neighbor = search.closestPoints(pointLonLat, 1);

        // Define partition
        partition[jo] = distribution.partition(neighbor[0].payload());
        ++nobsOwnVec_[partition[jo]];
      }
    }

    // Reorder vectors
    locsGlb.resize(3*nobsGlb_);
    order_.resize(nobsGlb_);
    std::vector<int> iobsOwnVec(comm_.size(), 0);
    for (size_t jo = 0; jo < nobsGlb_; ++jo) {
      size_t offset = 0;
      for (int jt = 0; jt < partition[jo]; ++jt) {
        offset += nobsOwnVec_[jt];
      }
      offset += iobsOwnVec[partition[jo]];
      locsGlb[3*offset+0] = locsGlbTmp[3*jo+0];
      locsGlb[3*offset+1] = locsGlbTmp[3*jo+1];
      locsGlb[3*offset+2] = locsGlbTmp[3*jo+2];
      order_[offset] = orderGlbTmp[jo];
      ++iobsOwnVec[partition[jo]];
    }
  }

  // Broadcast number of observations
  comm_.broadcast(nobsOwnVec_, 0);
  nobsOwn_ = nobsOwnVec_[comm_.rank()];

  // Allocation of local vector
  std::vector<double> locsOwn(3*nobsOwn_);

  // Define counts and displacements
  std::vector<int> locsCounts;
  std::vector<int> locsDispls;
  for (size_t jt = 0; jt < comm_.size(); ++jt) {
    locsCounts.push_back(3*nobsOwnVec_[jt]);
  }
  locsDispls.push_back(0);
  for (size_t jt = 0; jt < comm_.size()-1; ++jt) {
    locsDispls.push_back(locsDispls[jt]+locsCounts[jt]);
  }

  // Scatter locations
  comm_.scatterv(locsGlb.begin(), locsGlb.end(), locsCounts, locsDispls,
    locsOwn.begin(), locsOwn.end(), 0);

  // Format data
  for (size_t jo = 0; jo < nobsOwn_; ++jo) {
    locs_.push_back(atlas::Point3(locsOwn[3*jo], locsOwn[3*jo+1], locsOwn[3*jo+2]));
  }

  // Generate times
  for (size_t jo = 0; jo < nobsOwn_; ++jo) {
    times_.push_back(winbgn_);
  }

  // Generate observations error
  const std::vector<double> err = config.getDoubleVector("error");
  atlas::FieldSet fset;
  fset.name() = config.getString("obserror");
  for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
    atlas::Field field(vars_[jvar].name(), atlas::array::make_datatype<double>(),
      atlas::array::make_shape(nobsOwn_, 1));
    auto view = atlas::array::make_view<double, 2>(field);
    for (size_t jo = 0; jo < nobsOwn_; ++jo) {
      view(jo, 0) = err[jvar];
    }
    fset.add(field);
  }

  // Save FieldSet
  nobsLoc_ = nobsOwn_;
  this->putdb(fset);

  // Fill halo
  this->fillHalo();

  oops::Log::trace() << classname() << "::generateDistribution done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsSpace::printJo(const ObsVector &,
                       const ObsVector &) {
  oops::Log::trace() << classname() << "::printJo starting" << std::endl;

  oops::Log::info() << "ObsSpace::printJo not implemented" << std::endl;

  oops::Log::trace() << classname() << "::printJo starting" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsSpace::screenObservations(const ObsVector & dep,
                                  const GeoVaLs & gv) const {
  oops::Log::trace() << classname() << "::screenObservations starting" << std::endl;

  // Define whether observations are screened or not
  // TODO(Benjamin): more clever screening
  std::vector<bool> validObs(nobsLoc_);
  for (size_t jo = 0; jo < nobsLoc_; ++jo) {
    validObs[jo] = true;
  }
  const int nobsLocScr = std::count(validObs.cbegin(), validObs.cend(), true);

  // Apply screening on time, location and data
  for (size_t jo = 0; jo < nobsLoc_; ++jo) {
    if (validObs[jo]) {
      screenedTimes_.push_back(times_[jo]);
      screenedLocations_.push_back(locs_[jo]);
    }
  }
  for (const auto & dataFset : data_) {
    atlas::FieldSet fset;
    fset.name() = dataFset.name();
    for (const auto & dataField : dataFset) {
      atlas::Field field(dataField.name(), atlas::array::make_datatype<double>(),
        atlas::array::make_shape(nobsLocScr, 1));
      const auto dataView = atlas::array::make_view<double, 2>(dataField);
      auto view = atlas::array::make_view<double, 2>(field);
      size_t joScr = 0;
      for (size_t jo = 0; jo < nobsLoc_; ++jo) {
        if (validObs[jo]) {
          view(joScr, 0) = dataView(jo, 0);
          ++joScr;
        }
      }
      fset.add(field);
    }
    screenedData_.push_back(fset);
  }

  // Write screened observations
  write(nameOut_, Write::Screened);

  oops::Log::trace() << classname() << "::screenObservations done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsSpace::read(const std::string & filePath) {
  oops::Log::trace() << classname() << "::read starting" << std::endl;

  // Global size and vectors
  int ngrp;
  std::vector<int> partition;
  std::vector<int> timesGlb;
  std::vector<double> locsGlb;
  std::vector<double> dataGlb;

  // Counts and displacements
  nobsOwnVec_.resize(comm_.size());

  // Define source grid horizontal distribution
  const atlas::grid::Distribution distribution = geom_->partitioner().partition(geom_->grid());

  // NetCDF IDs
  int retval, ncid, nobs_id, dateTime_id, order_id, longitude_id, latitude_id, height_id, data_id;
  std::vector<int> group_ids;

  if (comm_.rank() == 0) {
    // Open NetCDF file
    std::string ncFilePath = filePath + ".nc";
    oops::Log::info() << "Reading file: " << ncFilePath << std::endl;
    if (retval = nc_open(ncFilePath.c_str(), NC_NOWRITE, &ncid)) ERR(retval);

    // Get dimension
    if (retval = nc_inq_dimid(ncid, "Location", &nobs_id)) ERR(retval);
    if (retval = nc_inq_dimlen(ncid, nobs_id, &nobsGlb_)) ERR(retval);

    // Get groups list
    if (retval = nc_inq_grps(ncid, &ngrp, NULL)) ERR(retval);
    group_ids.resize(ngrp);
    if (retval = nc_inq_grps(ncid, NULL, group_ids.data())) ERR(retval);

    // Get order if available
    order_.resize(nobsGlb_);
    retval = nc_inq_varid(ncid, "order", &order_id);
    if (retval == NC_NOERR) {
      if (retval = nc_get_var_int(ncid, order_id, order_.data())) ERR(retval);
    } else {
      for (size_t jo = 0; jo < nobsGlb_; ++jo) {
        order_[jo] = jo;
      }
    }

    // Get MetaData
    int meta_group_id;
    if (retval = nc_inq_grp_ncid(ncid, "MetaData", &meta_group_id)) ERR(retval);

    // Get dateTime from MetaData
    std::vector<int64_t> dateTime(nobsGlb_);
    if (retval = nc_inq_varid(meta_group_id, "dateTime", &dateTime_id)) ERR(retval);
    if (retval = nc_get_var_long(meta_group_id, dateTime_id, dateTime.data())) ERR(retval);
    const std::string dateTime_units_key = "units";
    size_t attlen;
    if (retval = nc_inq_attlen(meta_group_id, dateTime_id, dateTime_units_key.c_str(),
      &attlen)) ERR(retval);
    char **dateTime_units_char = reinterpret_cast<char**>(malloc(attlen*sizeof(char*)));
    memset(dateTime_units_char, 0, attlen*sizeof(char*));
    if (retval = nc_get_att_string(meta_group_id, dateTime_id, dateTime_units_key.c_str(),
      dateTime_units_char)) ERR(retval);
    std::string dateTime_units_value(*dateTime_units_char);
    util::DateTime start(dateTime_units_value.substr(14, 20));

    // Get longitude, latitude and height from MetaData
    std::vector<float> longitude(nobsGlb_);
    std::vector<float> latitude(nobsGlb_);
    std::vector<float> height(nobsGlb_);
    if (retval = nc_inq_varid(meta_group_id, "longitude", &longitude_id)) ERR(retval);
    if (retval = nc_get_var_float(meta_group_id, longitude_id, longitude.data()))
      ERR(retval);
    if (retval = nc_inq_varid(meta_group_id, "latitude", &latitude_id)) ERR(retval);
    if (retval = nc_get_var_float(meta_group_id, latitude_id, latitude.data())) ERR(retval);
    if (retval = nc_inq_varid(meta_group_id, "height", &height_id)) ERR(retval);
    if (retval = nc_get_var_float(meta_group_id, height_id, height.data())) ERR(retval);

    // Split observations between tasks
    // TODO(Benjamin): other distributions
    partition.resize(nobsGlb_);
    if (true) {
      // Using nearest neighbor
      atlas::util::IndexKDTree search;
      search.reserve(geom_->grid().size());
      size_t jnode = 0;
      for (const auto & lonLat : geom_->grid().lonlat()) {
        atlas::PointLonLat pointLonLat(lonLat);
        pointLonLat.normalise();
        atlas::PointXY point(pointLonLat);
        search.insert(point, jnode);
        ++jnode;
      }
      search.build();

      for (size_t jo = 0; jo < nobsGlb_; ++jo) {
        // Find MPI task
        atlas::PointLonLat pointLonLat(longitude[jo], latitude[jo]);
        pointLonLat.normalise();

        // Search nearest neighborass
        atlas::util::IndexKDTree::ValueList neighbor = search.closestPoints(pointLonLat, 1);

        // Define partition
        partition[jo] = distribution.partition(neighbor[0].payload());
        ++nobsOwnVec_[partition[jo]];
      }
    }

    // Get ordered time and location
    timesGlb.resize(6*nobsGlb_);
    locsGlb.resize(3*nobsGlb_);
    std::vector<int> iobsOwnVec(comm_.size(), 0);
    for (size_t jo = 0; jo < nobsGlb_; ++jo) {
      size_t offset = 0;
      for (int jt = 0; jt < partition[jo]; ++jt) {
        offset += nobsOwnVec_[jt];
      }
      offset += iobsOwnVec[partition[jo]];
      util::DateTime startTest = start + util::Duration(dateTime[jo]);
      startTest.toYYYYMMDDhhmmss(timesGlb[6*offset+0], timesGlb[6*offset+1], timesGlb[6*offset+2],
        timesGlb[6*offset+3], timesGlb[6*offset+4], timesGlb[6*offset+5]);
      locsGlb[3*offset+0] = static_cast<double>(longitude[jo]);
      locsGlb[3*offset+1] = static_cast<double>(latitude[jo]);
      locsGlb[3*offset+2] = static_cast<double>(height[jo]);
      ++iobsOwnVec[partition[jo]];
    }
  }

  // Broadcast number of observations for each task
  comm_.broadcast(nobsOwnVec_, 0);
  nobsOwn_ = nobsOwnVec_[comm_.rank()];
  nobsGlb_ = 0;
  for (size_t jt = 0; jt < comm_.size(); ++jt) {
    nobsGlb_ += nobsOwnVec_[jt];
  }

  // Allocation of local vectors
  std::vector<int> timesOwn(6*nobsOwn_);
  std::vector<double> locsOwn(3*nobsOwn_);

  // Define counts and displacements
  std::vector<int> timesCounts;
  std::vector<int> locsCounts;
  std::vector<int> dataCounts;
  std::vector<int> timesDispls;
  std::vector<int> locsDispls;
  std::vector<int> dataDispls;
  for (size_t jt = 0; jt < comm_.size(); ++jt) {
    timesCounts.push_back(6*nobsOwnVec_[jt]);
    locsCounts.push_back(3*nobsOwnVec_[jt]);
    dataCounts.push_back(vars_.size()*nobsOwnVec_[jt]);
  }
  timesDispls.push_back(0);
  locsDispls.push_back(0);
  dataDispls.push_back(0);
  for (size_t jt = 0; jt < comm_.size()-1; ++jt) {
    timesDispls.push_back(timesDispls[jt]+timesCounts[jt]);
    locsDispls.push_back(locsDispls[jt]+locsCounts[jt]);
    dataDispls.push_back(dataDispls[jt]+dataCounts[jt]);
  }

  // Scatter times and locations
  comm_.scatterv(timesGlb.begin(), timesGlb.end(), timesCounts, timesDispls,
    timesOwn.begin(), timesOwn.end(), 0);
  comm_.scatterv(locsGlb.begin(), locsGlb.end(), locsCounts, locsDispls,
    locsOwn.begin(), locsOwn.end(), 0);

  // Format local times and locations
  for (size_t jo = 0; jo < nobsOwn_; ++jo) {
    times_.push_back(util::DateTime(timesOwn[6*jo], timesOwn[6*jo+1], timesOwn[6*jo+2],
      timesOwn[6*jo+3], timesOwn[6*jo+4], timesOwn[6*jo+5]));
    locs_.push_back(atlas::Point3(locsOwn[3*jo], locsOwn[3*jo+1], locsOwn[3*jo+2]));
  }

  // Broadcast number of groups
  comm_.broadcast(ngrp, 0);

  for (int jgrp = 0; jgrp < ngrp; ++jgrp) {
    // Get group name
    std::string grpName;
    size_t grpNameLen = 0;
    if (comm_.rank() == 0) {
      if (retval = nc_inq_grpname_len(group_ids[jgrp], &grpNameLen)) ERR(retval);
      grpName.resize(grpNameLen);
      if (retval = nc_inq_grpname(group_ids[jgrp], &grpName[0])) ERR(retval);
    }
    oops::mpi::broadcastString(comm_, grpName, 0);

    if (grpName.substr(0, grpNameLen-1) != "MetaData") {
      // Non-MetaData group
      if (comm_.rank() == 0) {
        std::vector<float> dataVar(nobsGlb_);
        dataGlb.resize(vars_.size()*nobsGlb_);
        for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
          // Get data
          if (retval = nc_inq_varid(group_ids[jgrp], vars_[jvar].name().c_str(), &data_id))
            ERR(retval);
          if (retval = nc_get_var_float(group_ids[jgrp], data_id, dataVar.data())) ERR(retval);

          // Get ordered data
          std::vector<int> iobsOwnVec(comm_.size(), 0);
          for (size_t jo = 0; jo < nobsGlb_; ++jo) {
            size_t offset = 0;
            for (int jt = 0; jt < partition[jo]; ++jt) {
              offset += nobsOwnVec_[jt];
            }
            offset += iobsOwnVec[partition[jo]];
            dataGlb[vars_.size()*offset+jvar] = static_cast<double>(dataVar[jo]);
            ++iobsOwnVec[partition[jo]];
          }
        }
      }

      // Allocation of local vector
      std::vector<double> dataOwn(vars_.size()*nobsOwn_);

      // Scatter data
      comm_.scatterv(dataGlb.begin(), dataGlb.end(), dataCounts, dataDispls,
        dataOwn.begin(), dataOwn.end(), 0);

      // Format data
      atlas::FieldSet fset;
      fset.name() = grpName.substr(0, grpNameLen-1);
      for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
        atlas::Field field(vars_[jvar].name(), atlas::array::make_datatype<double>(),
          atlas::array::make_shape(nobsOwn_, 1));
        auto view = atlas::array::make_view<double, 2>(field);
        for (size_t jo = 0; jo < nobsOwn_; ++jo) {
          view(jo, 0) = dataOwn[jo*vars_.size()+jvar];
        }
        fset.add(field);
      }
      data_.push_back(fset);
    }
  }

  if (comm_.rank() == 0) {
    // Close file
    if (retval = nc_close(ncid)) ERR(retval);
  }

  // Add one second if time_ = winbgn_ and  winend_ > winbgn_
  if (winend_ > winbgn_) {
    for (auto & time : times_) {
      if (time == winbgn_) {
        time += util::Duration("PT1S");
      }
    }
  }

  // Halo expansion
  this->fillHalo();

  oops::Log::trace() << classname() << "::read done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsSpace::write(const std::string & filePath,
                     const bool & writeScreened) const {
  oops::Log::trace() << classname() << "::write starting" << std::endl;

  std::vector<util::DateTime> * times = 0;
  std::vector<atlas::Point3> * locs = 0;
  std::vector<atlas::FieldSet> * data = 0;

  // Set pointers
  if (writeScreened) {
    times = &screenedTimes_;
    locs = &screenedLocations_;
    data = &screenedData_;
  } else {
    times = &times_;
    locs = &locs_;
    data = &data_;
  }
  ASSERT(times);
  ASSERT(locs);
  ASSERT(data);


  // Global size check
  ASSERT(nobsOwnVec_.size() == comm_.size());

  // Allocation of global vectors
  std::vector<int> timesGlb;
  std::vector<double> locsGlb;
  std::vector<double> dataGlb;
  if (comm_.rank() == 0) {
    timesGlb.resize(6*nobsGlb_);
    locsGlb.resize(3*nobsGlb_);
    dataGlb.resize(vars_.size()*nobsGlb_);
  }

  // Allocation of local vectors
  std::vector<int> timesOwn(6*nobsOwn_);
  std::vector<double> locsOwn(3*nobsOwn_);
  std::vector<double> dataOwn(vars_.size()*nobsOwn_);

  // Define counts and displacements
  std::vector<int> timesCounts;
  std::vector<int> locsCounts;
  std::vector<int> dataCounts;
  std::vector<int> timesDispls;
  std::vector<int> locsDispls;
  std::vector<int> dataDispls;
  for (size_t jt = 0; jt < comm_.size(); ++jt) {
    timesCounts.push_back(6*nobsOwnVec_[jt]);
    locsCounts.push_back(3*nobsOwnVec_[jt]);
    dataCounts.push_back(vars_.size()*nobsOwnVec_[jt]);
  }
  timesDispls.push_back(0);
  locsDispls.push_back(0);
  dataDispls.push_back(0);
  for (size_t jt = 0; jt < comm_.size()-1; ++jt) {
    timesDispls.push_back(timesDispls[jt]+timesCounts[jt]);
    locsDispls.push_back(locsDispls[jt]+locsCounts[jt]);
    dataDispls.push_back(dataDispls[jt]+dataCounts[jt]);
  }

  // Format time and location
  for (size_t jo = 0; jo < nobsOwn_; ++jo) {
    (*times)[jo].toYYYYMMDDhhmmss(timesOwn[6*jo], timesOwn[6*jo+1], timesOwn[6*jo+2],
      timesOwn[6*jo+3], timesOwn[6*jo+4], timesOwn[6*jo+5]);
    locsOwn[3*jo+0] = (*locs)[jo][0];
    locsOwn[3*jo+1] = (*locs)[jo][1];
    locsOwn[3*jo+2] = (*locs)[jo][2];
  }

  // Gather times, locations and data
  comm_.gatherv(timesOwn, timesGlb, timesCounts, timesDispls, 0);
  comm_.gatherv(locsOwn, locsGlb, locsCounts, locsDispls, 0);

  // NetCDF IDs
  int retval, ncid, nobs_id, d_id[1], Locations_id, order_id, metaData_id, dateTime_id,
    longitude_id, latitude_id, height_id;
  std::vector<int> group_ids;
  std::vector<int> groupVar_ids;

  if (comm_.rank() == 0) {
    // Define locations vector
    std::vector<int> Locations(nobsGlb_);
    for (size_t jo = 0; jo < nobsGlb_; ++jo) {
      Locations[jo] = jo;
    }

    // Find start dateTime
    util::DateTime start(timesGlb[0], timesGlb[1], timesGlb[2], timesGlb[3], timesGlb[4],
      timesGlb[5]);
    for (size_t jo = 1; jo < nobsGlb_; ++jo) {
      util::DateTime startTest(timesGlb[6*jo], timesGlb[6*jo+1], timesGlb[6*jo+2], timesGlb[6*jo+3],
        timesGlb[6*jo+4], timesGlb[6*jo+5]);
      if (startTest < start) {
        start = startTest;
      }
    }

    // Define dateTime
    std::vector<int64_t> dateTime(nobsGlb_);
    for (size_t jo = 0; jo < nobsGlb_; ++jo) {
      util::DateTime startTest(timesGlb[6*jo], timesGlb[6*jo+1], timesGlb[6*jo+2], timesGlb[6*jo+3],
        timesGlb[6*jo+4], timesGlb[6*jo+5]);
      dateTime[jo] = (startTest-start).toSeconds();
    }

    // Define longitude
    std::vector<float> longitude(nobsGlb_);
    for (size_t jo = 0; jo < nobsGlb_; ++jo) {
      longitude[jo] = locsGlb[3*jo];
    }

    // Define latitude
    std::vector<float> latitude(nobsGlb_);
    for (size_t jo = 0; jo < nobsGlb_; ++jo) {
      latitude[jo] = locsGlb[3*jo+1];
    }

    // Define height
    std::vector<float> height(nobsGlb_);
    for (size_t jo = 0; jo < nobsGlb_; ++jo) {
      height[jo] = locsGlb[3*jo+2];
    }

    // Create NetCDF file
    const std::string ncFilePath = writeScreened ? filePath + "_screened.nc" : filePath + ".nc";
    if (retval = nc_create(ncFilePath.c_str(), NC_CLOBBER|NC_NETCDF4, &ncid)) ERR(retval);

    // Global attributes
    const std::string ioda_layout_key = "_ioda_layout";
    const std::string ioda_layout_value = "ObsGroup";
    const char *ioda_layout_char[1] = {ioda_layout_value.c_str()};
    if (retval = nc_put_att_string(ncid, NC_GLOBAL, ioda_layout_key.c_str(), 1, ioda_layout_char))
      ERR(retval);
    const std::string ioda_layout_version_key = "_ioda_layout_version";
    const int64_t ioda_layout_version_value = 0;
    if (retval = nc_put_att_long(ncid, NC_GLOBAL, ioda_layout_version_key.c_str(), NC_INT64, 1,
      &ioda_layout_version_value)) ERR(retval);

    // Create dimension
    if (retval = nc_def_dim(ncid, "Location", NC_UNLIMITED, &nobs_id)) ERR(retval);
    d_id[0] = nobs_id;

    // Missing value
    const std::string fillValue_key = "_FillValue";

    // Define global variables
    if (retval = nc_def_var(ncid, "Location", NC_INT, 1, d_id, &Locations_id)) ERR(retval);
    if (retval = nc_put_att_int(ncid, Locations_id, fillValue_key.c_str(), NC_INT, 1,
      &util::missingValue<int>())) ERR(retval);
    if (retval = nc_def_var(ncid, "order", NC_INT, 1, d_id, &order_id)) ERR(retval);
    if (retval = nc_put_att_int(ncid, order_id, fillValue_key.c_str(), NC_INT, 1,
      &util::missingValue<int>())) ERR(retval);

    // Define metadata group
    if (retval = nc_def_grp(ncid, "MetaData", &metaData_id)) ERR(retval);
    if (retval = nc_def_var(metaData_id, "dateTime", NC_INT64, 1, d_id, &dateTime_id))
      ERR(retval);
    if (retval = nc_put_att_long(metaData_id, dateTime_id, fillValue_key.c_str(), NC_INT64, 1,
      &util::missingValue<int64_t>())) ERR(retval);
    const std::string dateTime_units_key = "units";
    const std::string dateTime_units_value = "seconds since " + start.toString();
    const char *dateTime_units_char[1] = {dateTime_units_value.c_str()};
    if (retval = nc_put_att_string(metaData_id, dateTime_id, dateTime_units_key.c_str(),
      1, dateTime_units_char)) ERR(retval);
    if (retval = nc_def_var(metaData_id, "longitude", NC_FLOAT, 1, d_id, &longitude_id))
      ERR(retval);
    if (retval = nc_put_att_float(metaData_id, longitude_id, fillValue_key.c_str(), NC_FLOAT, 1,
      &util::missingValue<float>())) ERR(retval);
    const std::string longitude_units_key = "units";
    const std::string longitude_units_value = "degrees_east";
    const char *longitude_units_char[1] = {longitude_units_value.c_str()};
    if (retval = nc_put_att_string(metaData_id, longitude_id, longitude_units_key.c_str(), 1,
      longitude_units_char)) ERR(retval);
    if (retval = nc_def_var(metaData_id, "latitude", NC_FLOAT, 1, d_id, &latitude_id))
      ERR(retval);
    if (retval = nc_put_att_float(metaData_id, latitude_id, fillValue_key.c_str(), NC_FLOAT, 1,
      &util::missingValue<float>())) ERR(retval);
    const std::string latitude_units_key = "units";
    const std::string latitude_units_value = "degrees_north";
    const char *latitude_units_char[1] = {latitude_units_value.c_str()};
    if (retval = nc_put_att_string(metaData_id, latitude_id, latitude_units_key.c_str(), 1,
      latitude_units_char)) ERR(retval);
    if (retval = nc_def_var(metaData_id, "height", NC_FLOAT, 1, d_id, &height_id))
      ERR(retval);
    if (retval = nc_put_att_float(metaData_id, height_id, fillValue_key.c_str(), NC_FLOAT, 1,
      &util::missingValue<float>())) ERR(retval);
    const std::string height_units_key = "units";
    const std::string height_units_value = "m";
    const char *height_units_char[1] = {height_units_value.c_str()};
    if (retval = nc_put_att_string(metaData_id, height_id, height_units_key.c_str(), 1,
      height_units_char)) ERR(retval);

    // Define data groups
    for (auto const & fset : *data) {
      int group_id;
      if (retval = nc_def_grp(ncid, fset.name().c_str(), &group_id)) ERR(retval);
      group_ids.push_back(group_id);
      for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
        int groupVar_id;
        if (retval = nc_def_var(group_id, vars_[jvar].name().c_str(), NC_FLOAT, 1, d_id,
          &groupVar_id)) ERR(retval);
        if (retval = nc_put_att_float(group_id, groupVar_id, fillValue_key.c_str(), NC_FLOAT, 1,
          &util::missingValue<float>())) ERR(retval);
        groupVar_ids.push_back(groupVar_id);
      }
    }

    // End definition mode
    if (retval = nc_enddef(ncid)) ERR(retval);

    // Write metadata
    const size_t init = 0;
    const size_t nobs = Locations.size();
    if (retval = nc_put_vara_int(ncid, Locations_id, &init, &nobs, Locations.data())) ERR(retval);
    if (retval = nc_put_vara_int(ncid, order_id, &init, &nobs, order_.data())) ERR(retval);
    if (retval = nc_put_vara_long(metaData_id, dateTime_id, &init, &nobs, dateTime.data()))
      ERR(retval);
    if (retval = nc_put_vara_float(metaData_id, longitude_id, &init, &nobs, longitude.data()))
      ERR(retval);
    if (retval = nc_put_vara_float(metaData_id, latitude_id, &init, &nobs, latitude.data()))
      ERR(retval);
    if (retval = nc_put_vara_float(metaData_id, height_id, &init, &nobs, height.data()))
      ERR(retval);
  }

  size_t igrp = 0;
  for (auto const & fset : *data) {
    // Format data
    for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
      const atlas::Field field = fset[vars_[jvar].name()];
      const auto view = atlas::array::make_view<double, 2>(field);
      for (size_t jo = 0; jo < nobsOwn_; ++jo) {
        dataOwn[jo*vars_.size()+jvar] = view(jo, 0);
      }
    }

    // Gather data
    comm_.gatherv(dataOwn, dataGlb, dataCounts, dataDispls, 0);

    if (comm_.rank() == 0) {
      // Write data
      for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
        std::vector<float> dataVar(nobsGlb_);
        for (size_t jo = 0; jo < nobsGlb_; ++jo) {
          dataVar[jo] = static_cast<float>(dataGlb[jo*vars_.size()+jvar]);
        }
        if (retval = nc_put_var_float(group_ids[igrp], groupVar_ids[igrp*vars_.size()+jvar],
          dataVar.data())) ERR(retval);
      }
    }
    ++igrp;
  }

  if (comm_.rank() == 0) {
    // Close file
    if (retval = nc_close(ncid)) ERR(retval);
  }

  // Reset pointers
  times = 0;
  locs = 0;
  data = 0;

  oops::Log::trace() << classname() << "::write done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsSpace::fillHalo() {
  oops::Log::trace() << classname() << "::fillHalo starting" << std::endl;

  if (distribution_.empty()) {
    // No halo
    nobsLoc_ = nobsOwn_;
  } else {
    // Get center, radius and length-scale
    const std::vector<double> center = distribution_.getDoubleVector("center");
    const double radius = distribution_.getDouble("radius");
    const double horLengthScale = distribution_.getDouble("horizontal length-scale");

    // Allgather center and radius
    std::vector<double> centerLonVec(comm_.size());
    std::vector<double> centerLatVec(comm_.size());
    std::vector<double> radiusVec(comm_.size());
    comm_.allGather(center[0], centerLonVec.begin(), centerLonVec.end());
    comm_.allGather(center[1], centerLatVec.begin(), centerLatVec.end());
    comm_.allGather(radius, radiusVec.begin(), radiusVec.end());

    // Prepare send buffer index
    std::vector<int> sendCounts(comm_.size(), 0);
    for (size_t jt = 0; jt < comm_.size(); ++jt) {
      if (jt != comm_.rank()) {
        const atlas::PointLonLat centerPoint({centerLonVec[jt], centerLatVec[jt]});
        const double haloSize = radiusVec[jt]+horLengthScale;

        for (size_t jo = 0; jo < nobsOwn_; ++jo) {
          // Compute distance to task center
          const atlas::PointLonLat obsPoint({locs_[jo][0], locs_[jo][1]});
          const double dist = atlas::util::Earth().distance(centerPoint, obsPoint);
          if (dist <= haloSize) {
            sendBufIndex_.push_back(jo);
            ++sendCounts[jt];
          }
        }
      }
    }
    nSend_ = sendBufIndex_.size();

    // Communicate sendCounts to get recvCounts
    std::vector<int> recvCounts(comm_.size());
    std::vector<int> tmpCounts(comm_.size(), 1);
    std::vector<int> tmpDispls;
    tmpDispls.push_back(0);
    for (size_t jt = 0; jt < comm_.size()-1; ++jt) {
      tmpDispls.push_back(tmpDispls[jt]+tmpCounts[jt]);
    }
    comm_.allToAllv(sendCounts.data(), tmpCounts.data(), tmpDispls.data(), recvCounts.data(),
      tmpCounts.data(), tmpDispls.data());

    // Number of received observations, local number of observations
    nRecv_ = 0;
    for (const auto & item : recvCounts) {
      nRecv_ += item;
    }
    nobsLoc_ = nobsOwn_+nRecv_;

    // Define displs
    std::vector<int> sendDispls;
    std::vector<int> recvDispls;
    sendDispls.push_back(0);
    recvDispls.push_back(0);
    for (size_t jt = 0; jt < comm_.size()-1; ++jt) {
      sendDispls.push_back(sendDispls[jt]+sendCounts[jt]);
      recvDispls.push_back(recvDispls[jt]+recvCounts[jt]);
    }

    // Prepare counts and displs for time, locations and data
    std::vector<int> timeSendCounts(comm_.size());
    std::vector<int> timeRecvCounts(comm_.size());
    std::vector<int> timeSendDispls(comm_.size());
    std::vector<int> timeRecvDispls(comm_.size());
    std::vector<int> locsSendCounts(comm_.size());
    std::vector<int> locsRecvCounts(comm_.size());
    std::vector<int> locsSendDispls(comm_.size());
    std::vector<int> locsRecvDispls(comm_.size());
    dataSendCounts_.resize(comm_.size());
    dataRecvCounts_.resize(comm_.size());
    dataSendDispls_.resize(comm_.size());
    dataRecvDispls_.resize(comm_.size());
    for (size_t jt = 0; jt < comm_.size(); ++jt) {
      timeSendCounts[jt] = 6*sendCounts[jt];
      timeRecvCounts[jt] = 6*recvCounts[jt];
      timeSendDispls[jt] = 6*sendDispls[jt];
      timeRecvDispls[jt] = 6*recvDispls[jt];
      locsSendCounts[jt] = 3*sendCounts[jt];
      locsRecvCounts[jt] = 3*recvCounts[jt];
      locsSendDispls[jt] = 3*sendDispls[jt];
      locsRecvDispls[jt] = 3*recvDispls[jt];
      dataSendCounts_[jt] = vars_.size()*sendCounts[jt];
      dataRecvCounts_[jt] = vars_.size()*recvCounts[jt];
      dataSendDispls_[jt] = vars_.size()*sendDispls[jt];
      dataRecvDispls_[jt] = vars_.size()*recvDispls[jt];
    }

    // Format time and locations
    std::vector<int> timeSendBuf(6*nSend_);
    std::vector<double> locsSendBuf(3*nSend_);
    for (size_t jj = 0; jj < nSend_; ++jj) {
      size_t jo = sendBufIndex_[jj];
      times_[jo].toYYYYMMDDhhmmss(timeSendBuf[6*jj], timeSendBuf[6*jj+1], timeSendBuf[6*jj+2],
        timeSendBuf[6*jj+3], timeSendBuf[6*jj+4], timeSendBuf[6*jj+5]);
      locsSendBuf[3*jj+0] = locs_[jo][0];
      locsSendBuf[3*jj+1] = locs_[jo][1];
      locsSendBuf[3*jj+2] = locs_[jo][2];
    }

    // Communicate time and locations
    std::vector<int> timeRecvBuf(6*nRecv_);
    std::vector<double> locsRecvBuf(3*nRecv_);
    comm_.allToAllv(timeSendBuf.data(), timeSendCounts.data(), timeSendDispls.data(),
      timeRecvBuf.data(), timeRecvCounts.data(), timeRecvDispls.data());
    comm_.allToAllv(locsSendBuf.data(), locsSendCounts.data(), locsSendDispls.data(),
      locsRecvBuf.data(), locsRecvCounts.data(), locsRecvDispls.data());

    // Append halo to owned time and locations
    for (size_t jj = 0; jj < nRecv_; ++jj) {
      times_.push_back(util::DateTime(timeRecvBuf[6*jj], timeRecvBuf[6*jj+1], timeRecvBuf[6*jj+2],
        timeRecvBuf[6*jj+3], timeRecvBuf[6*jj+4], timeRecvBuf[6*jj+5]));
      locs_.push_back(atlas::Point3(locsRecvBuf[3*jj], locsRecvBuf[3*jj+1], locsRecvBuf[3*jj+2]));
    }

    for (auto & fset : data_) {
      // Fill halo of FieldSet
      this->fillHalo(fset);
    }
  }

  oops::Log::trace() << classname() << "::fillHalo done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsSpace::fillHalo(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::fillHalo starting" << std::endl;

  if (!distribution_.empty()) {
    // Format data
    std::vector<double> dataSendBuf(vars_.size()*nSend_);
    for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
      const atlas::Field field = fset[vars_[jvar].name()];
      ASSERT(field.shape(0) == nobsOwn_);
      const auto view = atlas::array::make_view<double, 2>(field);
      for (size_t jj = 0; jj < nSend_; ++jj) {
        size_t jo = sendBufIndex_[jj];
        dataSendBuf[vars_.size()*jj+jvar] = view(jo, 0);
      }
    }

    // Communicate time and locations
    std::vector<double> dataRecvBuf(vars_.size()*nRecv_);
    comm_.allToAllv(dataSendBuf.data(), dataSendCounts_.data(), dataSendDispls_.data(),
      dataRecvBuf.data(), dataRecvCounts_.data(), dataRecvDispls_.data());

    // Append halo to owned data
    for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
      atlas::Field field = fset[vars_[jvar].name()];
      field.resize(atlas::array::make_shape(nobsLoc_, 1));
      auto view = atlas::array::make_view<double, 2>(field);
      for (size_t jj = 0; jj < nRecv_; ++jj) {
        view(nobsOwn_+jj, 0) = dataRecvBuf[vars_.size()*jj+jvar];
      }
    }
  }

  oops::Log::trace() << classname() << "::fillHalo done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsSpace::print(std::ostream & os) const {
  oops::Log::trace() << classname() << "::print starting" << std::endl;

  os << "ObsSpace: assimilation window = " << winbgn_ << " to " << winend_ << std::endl;
  os << "ObsSpace: file in = " << nameIn_ << ", file out = " << nameOut_ << std::endl;

  oops::Log::trace() << classname() << "::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace quenchxx
