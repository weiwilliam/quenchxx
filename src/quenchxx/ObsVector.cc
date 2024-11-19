/*
 * (C) Copyright 2023 Meteorologisk Institutt
 *
 */

#include "quenchxx/ObsVector.h"

#include <math.h>
#include <limits>
#include <random>

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

#include "quenchxx/GeoVaLs.h"

namespace quenchxx {

// -----------------------------------------------------------------------------

ObsVector::ObsVector(const ObsSpace & obsSpace)
  : comm_(obsSpace.getComm()), obsSpace_(obsSpace), data_(obsSpace.sizeLoc()) {
  oops::Log::trace() << classname() << "::ObsVector starting" << std::endl;

  zero();

  oops::Log::trace() << classname() << "::ObsVector done" << std::endl;
}

// -----------------------------------------------------------------------------

ObsVector::ObsVector(const ObsVector & other,
                     const bool copy)
  : comm_(other.comm_), obsSpace_(other.obsSpace_), data_(other.data_.size()) {
  oops::Log::trace() << classname() << "::ObsVector starting" << std::endl;

  if (copy) {
    data_ = other.data_;
  } else {
    zero();
  }

  oops::Log::trace() << classname() << "::ObsVector done" << std::endl;
}

// -----------------------------------------------------------------------------

ObsVector & ObsVector::operator= (const ObsVector & rhs) {
  oops::Log::trace() << classname() << "::operator= starting" << std::endl;

  ASSERT(data_.size() == rhs.data_.size());
  data_ = rhs.data_;

  oops::Log::trace() << classname() << "::operator= done" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------

ObsVector & ObsVector::operator*= (const double & zz) {
  oops::Log::trace() << classname() << "::operator*= starting" << std::endl;

  for (double & val : data_) {
    val *= zz;
  }

  oops::Log::trace() << classname() << "::operator*= done" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------

ObsVector & ObsVector::operator+= (const ObsVector & rhs) {
  oops::Log::trace() << classname() << "::operator+= starting" << std::endl;

  ASSERT(data_.size() == rhs.data_.size());
  for (size_t jj = 0; jj < data_.size(); ++jj) {
    data_[jj] += rhs.data_[jj];
  }

  oops::Log::trace() << classname() << "::operator+= done" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------

ObsVector & ObsVector::operator-= (const ObsVector & rhs) {
  oops::Log::trace() << classname() << "::operator-= starting" << std::endl;

  ASSERT(data_.size() == rhs.data_.size());
  for (size_t jj = 0; jj < data_.size(); ++jj) {
    data_[jj] -= rhs.data_[jj];
  }

  oops::Log::trace() << classname() << "::operator-= done" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------

ObsVector & ObsVector::operator*= (const ObsVector & rhs) {
  oops::Log::trace() << classname() << "::operator*= starting" << std::endl;

  ASSERT(data_.size() == rhs.data_.size());
  for (size_t jj = 0; jj < data_.size(); ++jj) {
    data_[jj] *= rhs.data_[jj];
  }

  oops::Log::trace() << classname() << "::operator*= done" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------

ObsVector & ObsVector::operator/= (const ObsVector & rhs) {
  oops::Log::trace() << classname() << "::operator/= starting" << std::endl;

  ASSERT(data_.size() == rhs.data_.size());
  for (size_t jj = 0; jj < data_.size(); ++jj) {
    data_[jj] /= rhs.data_[jj];
  }

  oops::Log::trace() << classname() << "::operator/= done" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------

void ObsVector::zero() {
  oops::Log::trace() << classname() << "::zero starting" << std::endl;

  for (double & val : data_) {
    val = 0.0;
  }

  oops::Log::trace() << classname() << "::zero done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsVector::invert() {
  oops::Log::trace() << classname() << "::invert starting" << std::endl;

  for (double & val : data_) {
    val = 1.0/val;
  }

  oops::Log::trace() << classname() << "::invert done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsVector::axpy(const double & zz,
                     const ObsVector & rhs) {
  oops::Log::trace() << classname() << "::axpy starting" << std::endl;

  ASSERT(data_.size() == rhs.data_.size());
  for (size_t jj = 0; jj < data_.size(); ++jj) {
    data_[jj] += zz*rhs.data_[jj];
  }

  oops::Log::trace() << classname() << "::axpy done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsVector::random() {
  oops::Log::trace() << classname() << "::random starting" << std::endl;

  // Global vector
  std::vector<double> randomPert;

  if (comm_.rank() == 0) {
    // Random numbers generator
    static std::mt19937 generator(2);
    static std::normal_distribution<double> randomDistrib(0.0, 1.0);

    // Generate perturbations
    std::vector<double> randomPertTmp(obsSpace_.size());
    for (size_t jj = 0; jj < obsSpace_.size(); ++jj) {
      randomPertTmp[jj] = randomDistrib(generator);
    }

    // Reorder perturbations
    randomPert.resize(obsSpace_.size());
    for (size_t jj = 0; jj < obsSpace_.size(); ++jj) {
      randomPert[jj] = randomPertTmp[obsSpace_.order()[jj]];
    }
  }

  // Define counts and displacements
  std::vector<int> dataCounts;
  std::vector<int> dataDispls;
  for (size_t jt = 0; jt < comm_.size(); ++jt) {
    dataCounts.push_back(obsSpace_.sizeVec()[jt]);
  }
  dataDispls.push_back(0);
  for (size_t jt = 0; jt < comm_.size()-1; ++jt) {
    dataDispls.push_back(dataDispls[jt]+dataCounts[jt]);
  }

  // Scatter perturbation
  comm_.scatterv(randomPert.begin(), randomPert.end(), dataCounts, dataDispls,
    data_.begin(), data_.end(), 0);

  oops::Log::trace() << classname() << "::random done" << std::endl;
}

// -----------------------------------------------------------------------------

double ObsVector::dot_product_with(const ObsVector & other) const {
  oops::Log::trace() << classname() << "::dot_product_with starting" << std::endl;

  ASSERT(data_.size() == other.data_.size());
  double zz = 0.0;
  for (size_t jj = 0; jj < data_.size(); ++jj) {
    zz += data_[jj]*other.data_[jj];
  }
  comm_.allReduceInPlace(zz, eckit::mpi::sum());

  oops::Log::trace() << classname() << "::dot_product_with done" << std::endl;
  return zz;
}

// -----------------------------------------------------------------------------

double ObsVector::rms() const {
  oops::Log::trace() << classname() << "::rms starting" << std::endl;

  double zz = 0.0;
  for (size_t jj = 0; jj < data_.size(); ++jj) {
    zz += data_[jj]*data_[jj];
  }
  if (obsSpace_.size() > 0) {
    comm_.allReduceInPlace(zz, eckit::mpi::sum());
    zz = sqrt(zz/obsSpace_.size());
  }

  oops::Log::trace() << classname() << "::rms done" << std::endl;
  return zz;
}

// -----------------------------------------------------------------------------

void ObsVector::print(std::ostream & os) const {
  oops::Log::trace() << classname() << "::print starting" << std::endl;

  double zmin = std::numeric_limits<double>::max();
  double zmax = -std::numeric_limits<double>::max();
  double zavg = 0.0;
  for (const double & val : data_) {
    if (val < zmin) zmin = val;
    if (val > zmax) zmax = val;
    zavg += val;
  }
  comm_.allReduceInPlace(zmin, eckit::mpi::min());
  comm_.allReduceInPlace(zmax, eckit::mpi::max());
  if (obsSpace_.size() > 0) {
    comm_.allReduceInPlace(zavg, eckit::mpi::sum());
    zavg /= obsSpace_.size();
  }

  os << "quenchxx[" << obsSpace_.size() << "]: Min=" << zmin << ", Max=" << zmax
    << ", Average=" << zavg;

  oops::Log::trace() << classname() << "::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace quenchxx
