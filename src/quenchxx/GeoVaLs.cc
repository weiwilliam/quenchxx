/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#include "quenchxx/GeoVaLs.h"

#include <netcdf.h>
#include <limits>

#include "eckit/exception/Exceptions.h"

#include "oops/util/Logger.h"
#include "oops/util/RandomField.h"

#include "quenchxx/ObsSpace.h"

#define ERR(e) {throw eckit::Exception(nc_strerror(e), Here());}

namespace quenchxx {

// -----------------------------------------------------------------------------

GeoVaLs::GeoVaLs(const ObsSpace & obsSpace,
                 const varns::Variables & vars,
                 const Increment &,
                 const util::DateTime & t1,
                 const util::DateTime & t2)
  : comm_(obsSpace.getComm()), obsSpace_(obsSpace), iobs_(), gvFieldSet_() {
  oops::Log::trace() << classname() << "::GeoVaLs starting" << std::endl;

  // Active observations indices
  iobs_ = obsSpace_.timeSelect(t1, t2);

  // Create GeoVaLs fieldset
  for (const auto & var : vars.variables()) {
    atlas::Field gvField(var, atlas::array::make_datatype<double>(),
      atlas::array::make_shape(iobs_.size()));
    gvFieldSet_.add(gvField);
  }

  oops::Log::trace() << classname() << "::GeoVaLs done" << std::endl;
}

// -----------------------------------------------------------------------------

GeoVaLs::GeoVaLs(const ObsSpace & obsSpace,
                 const varns::Variables & vars,
                 const State &,
                 const util::DateTime & t1,
                 const util::DateTime & t2)
  : comm_(obsSpace.getComm()), obsSpace_(obsSpace), iobs_(), gvFieldSet_() {
  oops::Log::trace() << classname() << "::GeoVaLs starting" << std::endl;

  // Active observations indices
  iobs_ = obsSpace_.timeSelect(t1, t2);

  // Create GeoVaLs fieldset
  for (const auto & var : vars.variables()) {
    atlas::Field gvField(var, atlas::array::make_datatype<double>(),
      atlas::array::make_shape(iobs_.size()));
    gvFieldSet_.add(gvField);
  }

  oops::Log::trace() << classname() << "::GeoVaLs done" << std::endl;
}

// -----------------------------------------------------------------------------

void GeoVaLs::zero() {
  oops::Log::trace() << classname() << "::zero starting" << std::endl;

  for (auto gvField : gvFieldSet_) {
    auto gvView = atlas::array::make_view<double, 1>(gvField);
    for (size_t jj = 0; jj < gvField.size(); ++jj) {
      gvView(jj) = 0.0;
    }
  }

  oops::Log::trace() << classname() << "::zero done" << std::endl;
}

// -----------------------------------------------------------------------------

void GeoVaLs::random() {
  oops::Log::trace() << classname() << "::random starting" << std::endl;

  for (auto gvField : gvFieldSet_) {
    util::NormalDistributionField dist(gvField.size(), 0.0, 1.0);
    auto gvView = atlas::array::make_view<double, 1>(gvField);
    for (size_t jj = 0; jj < gvField.size(); ++jj) {
      gvView(jj) = dist[jj];
    }
  }

  oops::Log::trace() << classname() << "::random done" << std::endl;
}

// -----------------------------------------------------------------------------

double GeoVaLs::dot_product_with(const GeoVaLs & gv) const {
  oops::Log::trace() << classname() << "::dot_product_with starting" << std::endl;

  ASSERT(gvFieldSet_.field_names() == gv.gvFieldSet_.field_names());
  double zz = 0.0;
  for (const auto & var : gvFieldSet_.field_names()) {
    const atlas::Field gvField1 = gvFieldSet_[var];
    const atlas::Field gvField2 = gv.gvFieldSet_[var];
    ASSERT(gvField1.size() == gvField2.size());
    const auto gvView1 = atlas::array::make_view<double, 1>(gvField1);
    const auto gvView2 = atlas::array::make_view<double, 1>(gvField2);
    for (size_t jj = 0; jj < gvField1.size(); ++jj) {
      zz += gvView1(jj)*gvView2(jj);
    }
  }
  comm_.allReduceInPlace(zz, eckit::mpi::sum());

  oops::Log::trace() << classname() << "::dot_product_with done" << std::endl;
  return zz;
}

// -----------------------------------------------------------------------------

void GeoVaLs::read(const eckit::Configuration & conf) {
  oops::Log::trace() << classname() << "::read starting " << std::endl;

  // TODO(Benjamin): make it parallel

  // NetCDF IDs
  int retval, ncid, var_id;

  // Open NetCDF file
  const std::string filePath(conf.getString("filepath"));
  std::string ncFilePath = filePath + ".nc";
  if ((retval = nc_open(ncFilePath.c_str(), NC_NOWRITE, &ncid))) ERR(retval);

  for (size_t jvar = 0; jvar < iobs_.size(); ++jvar) {
    // Get field
    atlas::Field gvField = gvFieldSet_[static_cast<int>(jvar)];

    // Get variables ids
    if ((retval = nc_inq_varid(ncid, gvField.name().c_str(), &var_id))) ERR(retval);

    // Get data
    std::vector<double> var(iobs_.size());

    // Get data
    if ((retval = nc_get_var_double(ncid, var_id, var.data()))) ERR(retval);

    // Fill field
    auto gvView = atlas::array::make_view<double, 1>(gvField);
    for (size_t jo = 0; jo < iobs_.size(); ++jo) {
      gvView(jo) = var[jo];
    }
  }

  // Close file
  if ((retval = nc_close(ncid))) ERR(retval);

  oops::Log::trace() << classname() << "::read done" << std::endl;
}

// -----------------------------------------------------------------------------

void GeoVaLs::write(const eckit::Configuration & conf) const {
  oops::Log::trace() << classname() << "::write starting" << std::endl;

  // TODO(Benjamin): make it parallel

  // NetCDF IDs
  int retval, ncid, nobs_id, dvar_id[1], var_id[gvFieldSet_.size()];

  // Create NetCDF file
  const std::string filePath(conf.getString("filepath"));
  const std::string ncFilePath = filePath + ".nc";
  if ((retval = nc_create(ncFilePath.c_str(), NC_CLOBBER, &ncid))) ERR(retval);

  // Create dimensions
  if ((retval = nc_def_dim(ncid, "nobs", iobs_.size(), &nobs_id))) ERR(retval);

  // Dimension array
  dvar_id[0] = nobs_id;

  for (int jvar = 0; jvar < gvFieldSet_.size(); ++jvar) {
    // Get field
    const atlas::Field gvField = gvFieldSet_[jvar];

    // Define variables
    if ((retval = nc_def_var(ncid, gvField.name().c_str(), NC_DOUBLE, 1, dvar_id,
      &var_id[jvar]))) ERR(retval);
  }

  // End definition mode
  if ((retval = nc_enddef(ncid))) ERR(retval);

  for (int jvar = 0; jvar < gvFieldSet_.size(); ++jvar) {
    // Get field
    const atlas::Field gvField = gvFieldSet_[jvar];

    // Copy data
    const auto gvView = atlas::array::make_view<double, 1>(gvField);
    std::vector<double> var(iobs_.size());
    for (size_t jo = 0; jo < iobs_.size(); ++jo) {
      var[jo] = gvView[jo];
    }

    // Write data
    if ((retval = nc_put_var_double(ncid, var_id[jvar], var.data()))) ERR(retval);
  }

  if ((retval = nc_close(ncid))) ERR(retval);

  oops::Log::trace() << classname() << "::write done" << std::endl;
}

// -----------------------------------------------------------------------------

void GeoVaLs::print(std::ostream & os) const {
  oops::Log::trace() << classname() << "::print starting" << std::endl;

  for (const auto & gvField : gvFieldSet_) {
    const auto gvView = atlas::array::make_view<double, 1>(gvField);
    double zmin = std::numeric_limits<double>::max();
    double zmax = -std::numeric_limits<double>::max();
    double zavg = 0.0;
    for (size_t jj = 0; jj < gvField.size(); ++jj) {
      if (gvView(jj) < zmin) zmin = gvView(jj);
      if (gvView(jj) > zmax) zmax = gvView(jj);
      zavg += gvView(jj);
    }
    comm_.allReduceInPlace(zmin, eckit::mpi::min());
    comm_.allReduceInPlace(zmax, eckit::mpi::max());
    if (obsSpace_.size() > 0) {
      comm_.allReduceInPlace(zavg, eckit::mpi::sum());
      zavg /= obsSpace_.size();
    }
    os << "GeoVaLs for " << gvField.name() << "[" << obsSpace_.size() << "]: "
      << "Min=" << zmin << ", Max=" << zmax << ", Average=" << zavg << std::endl;
  }

  oops::Log::trace() << classname() << "::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace quenchxx
