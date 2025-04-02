/*
 * (C) Copyright 2022 UCAR.
 * (C) Copyright 2023-2024 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "quenchxx/Geometry.h"

#include <netcdf.h>

#include <cmath>
#include <sstream>

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"
#include "atlas/meshgenerator.h"
#include "atlas/util/Geometry.h"
#include "atlas/util/KDTree.h"
#include "atlas/util/Point.h"

#include "eckit/exception/Exceptions.h"
#include "eckit/mpi/Comm.h"

#include "oops/generic/gc99.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/FunctionSpaceHelpers.h"
#include "oops/util/Logger.h"

#include "quenchxx/Fields.h"
#include "quenchxx/GeometryIterator.h"

#define ERR(e, msg) {std::string s(nc_strerror(e)); throw eckit::Exception(s + ": " + msg, Here());}

namespace quenchxx {

// -----------------------------------------------------------------------------

Geometry::Geometry(const eckit::Configuration & config,
                   const eckit::mpi::Comm & comm)
  : comm_(comm), groups_() {
  oops::Log::trace() << classname() << "::Geometry starting" << std::endl;

  GeometryParameters params;
  params.deserialize(config);

  // Setup atlas geometric data structures
  atlas::FieldSet fieldsetOwnedMask;
  util::setupFunctionSpace(comm_, config, grid_, partitioner_, mesh_, functionSpace_,
    fieldsetOwnedMask);
  halo_ = params.halo.value();
  gridType_ = params.grid.value().getString("type", "no_type");

  // Setup geometry fields
  fields_ = atlas::FieldSet();

  // Add owned points mask -- this mask does not depend on the group so was precomputed
  fields_->add(fieldsetOwnedMask.field("owned"));

  // Levels direction
  levelsAreTopDown_ = params.levelsAreTopDown.value();

  // Model data
  modelData_ = params.modelData.value();

  // Latitudes from south to north in files
  latSouthToNorth_ = params.latSouthToNorth.value();

  // Variable name alias
  setupAlias(params);

  // Groups
  size_t groupIndex = 0;
  for (const auto & groupParams : params.groups.value()) {
    // Use this group index for all the group variables
    for (const auto & var : groupParams.variables.value()) {
      if (groupIndex_.find(var) != groupIndex_.end()) {
        throw eckit::UserError(
          "Same variable present in distinct groups " + var, Here());
      } else {
        groupIndex_[var] = groupIndex;
      }
    }

    // Define group
    groupData group;

    // Number of levels
    group.levels_ = groupParams.levels.value();

    // Corresponding level for 2D variables (first or last)
    group.lev2d_ = groupParams.lev2d.value();

    // Vertical coordinate
    setupVertCoord(groupParams, groupIndex, group);

    // Setup mask
    setupMask(groupParams, groupIndex, group);

    // Save group
    groups_.push_back(group);

    // Increment group index
    groupIndex++;
  }

  // Interpolation
  const auto &interpParams = params.interpolation.value();
  if (interpParams != boost::none) {
    interpolation_ = interpParams->toConfiguration();
  } else {
    interpolation_ = eckit::LocalConfiguration();
      interpolation_.set("interpolation type", "unstructured");
  }

  // GeometryData
  if (interpolation_.getString("interpolation type") == "unstructured") {
    geomData_.reset(new oops::GeometryData(functionSpace_, fields_, levelsAreTopDown_, comm_));
  }

  // Check for duplicate points
  const auto ghostView = atlas::array::make_view<int, 1>(functionSpace_.ghost());
  const auto ownedView = atlas::array::make_view<int, 2>(fields_.field("owned"));
  size_t duplicatedPointsCount = 0;
  for (atlas::idx_t jnode = 0; jnode < fields_.field("owned").shape(0); ++jnode) {
    // Duplicate point = owned==0 and ghost==0 (see util::setupFunctionSpace in oops)
    if (ghostView(jnode) == 0 && ownedView(jnode, 0) == 0) {
      ++duplicatedPointsCount;
    }
  }
  comm_.allReduceInPlace(duplicatedPointsCount, eckit::mpi::sum());
  duplicatePoints_ = (duplicatedPointsCount > 0);

  // Check lon/lat from files
  const auto &checkLonLatConf = params.checkLonLat.value();
  if (checkLonLatConf != boost::none) {
    checkLonLat(*checkLonLatConf);
  }

  // Iterator dimension
  iteratorDimension_ = config.getInt("iterator dimension", 2);
  ASSERT((iteratorDimension_ == 2) || (iteratorDimension_ == 3));

  // First group vertical coordinate field
  const auto vertCoord = groups_[0].vertCoord_;

  // Domain size
  nnodes_ = vertCoord.shape(0);
  nlevs_ = vertCoord.shape(1);

  // Averaged vertical coordinate
  const auto vertCoordView = atlas::array::make_view<double, 2>(vertCoord);
  for (atlas::idx_t jlevel = 0; jlevel < nlevs_; ++jlevel) {
    double vertCoordAvg = 0.0;
    double counter = 0.0;
    for (atlas::idx_t jnode = 0; jnode < nnodes_; ++jnode) {
      if (ghostView(jnode) == 0) {
        vertCoordAvg += vertCoordView(jnode, jlevel);
        counter += 1.0;
      }
    }
    comm.allReduceInPlace(vertCoordAvg, eckit::mpi::sum());
    comm.allReduceInPlace(counter, eckit::mpi::sum());
    if (counter > 0.0) {
      vertCoordAvg /= counter;
    }
    vertCoordAvg_.push_back(vertCoordAvg);
  }

  // Print summary
  if (params.printSummary.value()) {
    this->print(oops::Log::info());
  }

  oops::Log::trace() << classname() << "::Geometry done" << std::endl;
}

// -----------------------------------------------------------------------------

Geometry::Geometry(const Geometry & other)
  : comm_(other.comm_), halo_(other.halo_), grid_(other.grid_), gridType_(other.gridType_),
  partitioner_(other.partitioner_), mesh_(other.mesh_), groupIndex_(other.groupIndex_),
  levelsAreTopDown_(other.levelsAreTopDown_), modelData_(other.modelData_), alias_(other.alias_),
  latSouthToNorth_(other.latSouthToNorth_), interpolation_(other.interpolation_),
  duplicatePoints_(other.duplicatePoints_), iteratorDimension_(other.iteratorDimension_),
  nnodes_(other.nnodes_), nlevs_(other.nlevs_), vertCoordAvg_(other.vertCoordAvg_) {
  oops::Log::trace() << classname() << "::Geometry starting" << std::endl;

  // Copy function space
  if (other.functionSpace_.type() == "StructuredColumns") {
    // StructuredColumns
    functionSpace_ = atlas::functionspace::StructuredColumns(other.functionSpace_);
  } else if (other.functionSpace_.type() == "NodeColumns") {
    // NodeColumns
    if (grid_.name().compare(0, 2, std::string{"CS"}) == 0) {
      // CubedSphere
      functionSpace_ = atlas::functionspace::CubedSphereNodeColumns(other.functionSpace_);
    } else {
      // Other NodeColumns
      functionSpace_ = atlas::functionspace::NodeColumns(other.functionSpace_);
    }
  } else if (other.functionSpace_.type() == "PointCloud") {
    throw eckit::NotImplemented(other.functionSpace_.type() + " function space not supported",
      Here());
  } else {
    throw eckit::NotImplemented(other.functionSpace_.type() + " function space not supported yet",
      Here());
  }

  // Copy geometry fields
  fields_ = util::shareFields(other.fields_);

  // Copy groups
  for (size_t groupIndex = 0; groupIndex < other.groups_.size(); ++groupIndex) {
    // Define group
    groupData group;

    // Copy number of levels
    group.levels_ = other.groups_[groupIndex].levels_;

    // Copy corresponding level for 2D variables (first or last)
    group.lev2d_ = other.groups_[groupIndex].lev2d_;

    // Copy vertical coordinate
    group.vertCoord_ = other.groups_[groupIndex].vertCoord_;

    // Copy averaged vertical coordinate
    group.vertCoordAvg_ = other.groups_[groupIndex].vertCoordAvg_;

    // Copy mask size
    group.gmaskSize_ = other.groups_[groupIndex].gmaskSize_;

    // Save group
    groups_.push_back(group);
  }

  // Geometry data
  if (interpolation_.getString("interpolation type") == "unstructured") {
    geomData_.reset(new oops::GeometryData(functionSpace_, fields_, levelsAreTopDown_, comm_));
  }

  oops::Log::trace() << classname() << "::Geometry done" << std::endl;
}

// -----------------------------------------------------------------------------

size_t Geometry::groupIndex(const std::string & var) const {
  oops::Log::trace() << classname() << "::groupIndex starting" << std::endl;

  if (groupIndex_.find(var) == groupIndex_.end()) {
    throw eckit::Exception("cannot find group index for variable " + var, Here());
  }
  const size_t groupIndex = groupIndex_.at(var);

  oops::Log::trace() << classname() << "::groupIndex done" << std::endl;
  return groupIndex;
}

// -----------------------------------------------------------------------------

std::vector<size_t> Geometry::variableSizes(const varns::Variables & vars) const {
  oops::Log::trace() << classname() << "::variableSizes starting" << std::endl;

  std::vector<size_t> sizes;
  for (const auto & var : vars) {
    sizes.push_back(levels(var.name()));
  }

  oops::Log::trace() << classname() << "::variableSizes done" << std::endl;
  return sizes;
}

// -----------------------------------------------------------------------------

std::vector<size_t> Geometry::variableSizes(const std::vector<std::string> & varNames) const {
  oops::Log::trace() << classname() << "::variableSizes starting" << std::endl;

  // Create variables
  const varns::Variables vars(varNames);

  oops::Log::trace() << classname() << "::variableSizes done" << std::endl;
  return variableSizes(vars);
}

// -----------------------------------------------------------------------------

void Geometry::print(std::ostream & os) const {
  oops::Log::trace() << classname() << "::print starting" << std::endl;

  std::string prefix;
  if (os.rdbuf() == oops::Log::info().rdbuf()) {
    prefix = "Info     : ";
  }
  os << prefix <<  "Quench geometry grid:" << std::endl;
  os << prefix << "- name: " << grid_.name() << std::endl;
  os << prefix << "- size: " << grid_.size() << std::endl;
  if (!grid_.domain().global()) {
    os << prefix << "Regional grid detected" << std::endl;
  }
  if (partitioner_) {
    os << prefix << "Partitioner:" << std::endl;
    os << prefix << "- type: " << partitioner_.type() << std::endl;
  }
  os << prefix << "Function space:" << std::endl;
  os << prefix << "- type: " << functionSpace_.type() << std::endl;
  os << prefix << "- halo: " << halo_ << std::endl;
  os << prefix << "Groups: " << std::endl;
  for (size_t groupIndex = 0; groupIndex < groups_.size(); ++groupIndex) {
    os << prefix << "- Group " << groupIndex << ":" << std::endl;
    os << prefix << "  Vertical levels: " << std::endl;
    os << prefix << "  - number: " << levels(groupIndex) << std::endl;
    os << prefix << "  - vertCoord: " << groups_[groupIndex].vertCoordAvg_ << std::endl;
    os << prefix << "  Mask size: " << static_cast<int>(groups_[groupIndex].gmaskSize_*100.0)
       << "%" << std::endl;
  }

  if (!modelData_.empty()) {
    os << prefix << "Model data: " << modelData_ << std::endl;
  }

  oops::Log::trace() << classname() << "::print done" << std::endl;
}

// -----------------------------------------------------------------------------

void Geometry::setupAlias(const GeometryParameters & params) {
  oops::Log::trace() << classname() << "::setupAlias starting" << std::endl;

  for (const auto & item : params.alias.value()) {
    eckit::LocalConfiguration confItem;
    item.serialize(confItem);
    alias_.push_back(confItem);
  }

  if (params.checkAliasConsistency.value()) {
    // Check alias consistency
    std::vector<std::string> vars;
    for (const auto & groupParams : params.groups.value()) {
      const std::vector<std::string> grpVars = groupParams.variables.value();
      vars.insert(vars.end(), grpVars.begin(), grpVars.end());
    }
    for (const auto & item : alias_) {
      const std::string codeVar = item.getString("in code");
      if (std::find(vars.begin(), vars.end(), codeVar) == vars.end()) {
        // Code variable not available in the list of variables anymore
        throw eckit::UserError("Alias error: code variable not available anymore", Here());
     } else {
        // Remove code variable from the list of available variables
        vars.erase(std::remove(vars.begin(), vars.end(), codeVar), vars.end());
      }
    }
    for (const auto & item : alias_) {
      const std::string fileVar = item.getString("in file");
      if (std::find(vars.begin(), vars.end(), fileVar) == vars.end()) {
        // Add file variable to the list of variables
        vars.push_back(fileVar);
      } else {
        // File variable is already present in the list of variables
        throw eckit::UserError("Alias error: duplicated file variable", Here());
      }
    }
  }

  oops::Log::trace() << classname() << "::setupAlias done" << std::endl;
}

// -----------------------------------------------------------------------------

void Geometry::setupVertCoord(const GroupParameters & groupParams,
                              const size_t & groupIndex,
                              groupData & group) {
  oops::Log::trace() << classname() << "::setupVertCoord starting" << std::endl;

  // Get optional parameters
  const auto &vertCoordConf = groupParams.vertCoordConf.value();

  // Get vertical coordinate name
  std::string vertCoordName = "vert_coord_" + std::to_string(groupIndex);
  if (vertCoordConf != boost::none) {
    if (vertCoordConf->has("name")) {
      vertCoordName = vertCoordConf->getString("name");
    }
  }

  // Create vertical coordinate
  group.vertCoord_ = functionSpace_.createField<double>(
    atlas::option::name(vertCoordName) | atlas::option::levels(group.levels_));

  // Set interpolation metadata
  group.vertCoord_.metadata().set("interp_type", "default");

  // Get view
  auto vertCoordView = atlas::array::make_view<double, 2>(group.vertCoord_);

  if (vertCoordConf != boost::none) {
    // Vertical coordinate from a configuration
    if (vertCoordConf->has("profile")) {
      // From a vector of doubles (one for each level)
      const std::vector<double> profile = vertCoordConf->getDoubleVector("profile");
      if (profile.size() != group.levels_) {
        throw eckit::UserError("Wrong number of levels in the user-specified vertical coordinate",
          Here());
      }
      for (atlas::idx_t jnode = 0; jnode < group.vertCoord_.shape(0); ++jnode) {
        for (size_t jlevel = 0; jlevel < group.levels_; ++jlevel) {
          vertCoordView(jnode, jlevel) = profile[jlevel];
        }
      }
    } else {
      // From a file
      const bool hybridVertCoord = (vertCoordConf->has("ak") && vertCoordConf->has("bk"));

      // Get variable to read
      const std::string varName = vertCoordConf->getString("variable");
      const oops::Variables vertCoordVars(std::vector<std::string>({varName}));

      // Create field
      Fields field(*this, vertCoordVars, util::DateTime());

      // Read field
      field.read(*vertCoordConf);

      // Get view
      const auto view = atlas::array::make_view<double, 2>(field.fieldSet()[varName]);

      if (hybridVertCoord) {
        // Hybrid coordinates
        std::vector<double> ak(group.levels_);
        std::vector<double> bk(group.levels_);
        if (comm_.rank() == 0) {
          // NetCDF file path
          const std::string ncFilePath = vertCoordConf->getString("filepath") + ".nc";

          // NetCDF IDs
          int ncid, retval, dim_id, ak_id, bk_id;
          size_t nab;

          // Open NetCDF file
          if ((retval = nc_open(ncFilePath.c_str(), NC_NOWRITE, &ncid))) ERR(retval, ncFilePath);

          // Get hybrid coordinates IDs
          const std::string akName = vertCoordConf->getString("ak");
          const std::string bkName = vertCoordConf->getString("bk");
          if ((retval = nc_inq_varid(ncid, akName.c_str(), &ak_id))) ERR(retval, akName);
          if ((retval = nc_inq_varid(ncid, bkName.c_str(), &bk_id))) ERR(retval, bkName);

          // Get hybrid coordinates dimension
          if ((retval = nc_inq_vardimid(ncid, ak_id, &dim_id))) ERR(retval, akName);
          if ((retval = nc_inq_dimlen(ncid, dim_id, &nab))) ERR(retval, "nab");

          // Read data
          std::vector<double> akFromFile(nab);
          std::vector<double> bkFromFile(nab);
          if ((retval = nc_get_var_double(ncid, ak_id, akFromFile.data()))) ERR(retval, akName);
          if ((retval = nc_get_var_double(ncid, bk_id, bkFromFile.data()))) ERR(retval, bkName);

          // Close file
          if ((retval = nc_close(ncid))) ERR(retval, ncFilePath);

          if (nab == group.levels_) {
            // Copy hybrid coefficients
            for (size_t jlevel = 0; jlevel < group.levels_; ++jlevel) {
              ak[jlevel] = akFromFile[jlevel];
              bk[jlevel] = bkFromFile[jlevel+1];
            }
          } else if (nab == group.levels_+1) {
            // Assuming field levels at hybrid coefficients half levels
            for (size_t jlevel = 0; jlevel < group.levels_; ++jlevel) {
              ak[jlevel] = 0.5*(akFromFile[jlevel]+akFromFile[jlevel+1]);
              bk[jlevel] = 0.5*(bkFromFile[jlevel]+bkFromFile[jlevel+1]);
            }
          } else {
            throw eckit::Exception("wrong number of levels in hybrid vertical coordinates", Here());
          }
        }

        // Broadcast hybrid coordinates
        comm_.broadcast(ak.begin(), ak.end(), 0);
        comm_.broadcast(bk.begin(), bk.end(), 0);

        // Compute hybrid vertical coordinate
        for (atlas::idx_t jnode = 0; jnode < group.vertCoord_.shape(0); ++jnode) {
          for (size_t jlevel = 0; jlevel < group.levels_; ++jlevel) {
            vertCoordView(jnode, jlevel) = ak[jlevel] + bk[jlevel]*view(jnode, 0);
          }
        }
      } else {
        // Copy 3D field
        for (atlas::idx_t jnode = 0; jnode < group.vertCoord_.shape(0); ++jnode) {
          for (size_t jlevel = 0; jlevel < group.levels_; ++jlevel) {
            vertCoordView(jnode, jlevel) = view(jnode, jlevel);
          }
        }
      }
    }
  } else {
    // From level index (default)
    for (atlas::idx_t jnode = 0; jnode < group.vertCoord_.shape(0); ++jnode) {
      for (size_t jlevel = 0; jlevel < group.levels_; ++jlevel) {
        vertCoordView(jnode, jlevel) = static_cast<double>(jlevel+1);
      }
    }
  }

  // Get ghost and owned views
  const auto ghostView = atlas::array::make_view<int, 1>(functionSpace_.ghost());
  const auto ownedView = atlas::array::make_view<int, 2>(fields_.field("owned"));

  // Average vertical coordinate
  for (size_t jlevel = 0; jlevel < group.levels_; ++jlevel) {
    // Initialization
    double avg = 0.0;
    double counter = 0.0;

    // Loop over owned points
    for (atlas::idx_t jnode = 0; jnode < group.vertCoord_.shape(0); ++jnode) {
      if (ghostView(jnode) == 0 && ownedView(jnode, 0) == 1) {
        avg += vertCoordView(jnode, jlevel);
        counter += 1.0;
      }
    }

    // Communication
    comm_.allReduceInPlace(avg, eckit::mpi::sum());
    comm_.allReduceInPlace(counter, eckit::mpi::sum());

    // Normalization
    if (counter > 0.0) {
      avg /= counter;
    }

    // Update profile
    group.vertCoordAvg_.push_back(avg);
  }

  // Add orography (mountain) on bottom level
  const auto &orographyParams = groupParams.orography.value();
  if (orographyParams != boost::none) {
    // Get top latitude value
    const atlas::PointLonLat topPoint({orographyParams->topLon.value(),
      orographyParams->topLat.value()});

    // Get lon/lat view
    const auto lonlatView = atlas::array::make_view<double, 2>(functionSpace_.lonlat());

    for (atlas::idx_t jnode = 0; jnode < lonlatView.shape(0); ++jnode) {
      // Get delta
      const double delta = (group.levels_ == 1) ? 1.0 :
        vertCoordView(jnode, group.levels_-2)-vertCoordView(jnode, group.levels_-1);

      // Get x and y points
      const atlas::PointLonLat xPoint({lonlatView(jnode, 0), orographyParams->topLat.value()});
      const atlas::PointLonLat yPoint({orographyParams->topLon.value(), lonlatView(jnode, 1)});

      // Compute normalization
      const double dxNorm = atlas::util::Earth().distance(xPoint, topPoint)
        /orographyParams->zonalLength.value();
      const double dyNorm = atlas::util::Earth().distance(yPoint, topPoint)
        /orographyParams->meridionalLength.value();
      const double distNorm = std::sqrt(dxNorm*dxNorm+dyNorm*dyNorm);

      // Define orography
      const double orography = delta*orographyParams->height.value()*oops::gc99(distNorm);

      // Add orography to existing vertical coordinate
      vertCoordView(jnode, group.levels_-1) += orography;
    }
  }

  // Add vertical coordinate in Geometry fields
  fields_->add(group.vertCoord_);

  oops::Log::trace() << classname() << "::setupVertCoord starting" << std::endl;
}


// -----------------------------------------------------------------------------

void Geometry::setupMask(const GroupParameters & groupParams,
                         const size_t & groupIndex,
                         groupData & group) {
  oops::Log::trace() << classname() << "::setupMask starting" << std::endl;

  // Default mask, set to 1 (true)
  const std::string gmaskName = "gmask_" + std::to_string(groupIndex);
  atlas::Field gmask = functionSpace_.createField<int>(
    atlas::option::name(gmaskName) | atlas::option::levels(group.levels_));
  auto maskView = atlas::array::make_view<int, 2>(gmask);
  maskView.assign(1);

  // Ghost view
  auto ghostView = atlas::array::make_view<int, 1>(functionSpace_.ghost());

  // Specific mask
  if (groupParams.maskType.value() == "none") {
    // No mask
  } else if (groupParams.maskType.value() == "sea") {
    // Read sea mask

    // Lon/lat sizes
    size_t nlon = 0;
    size_t nlat = 0;

    // File path
    ASSERT(groupParams.maskPath.value());
    const std::string ncFilePath = *groupParams.maskPath.value();

    // NetCDF IDs
    int ncid, retval, nlon_id, nlat_id, lon_id, lat_id, lsm_id;

    if (comm_.rank() == 0) {
      // Open NetCDF file
      if ((retval = nc_open(ncFilePath.c_str(), NC_NOWRITE, &ncid))) ERR(retval, ncFilePath);

      // Get lon/lat sizes
      if ((retval = nc_inq_dimid(ncid, "lon", &nlon_id))) ERR(retval, "lon");
      if ((retval = nc_inq_dimid(ncid, "lat", &nlat_id))) ERR(retval, "lat");
      if ((retval = nc_inq_dimlen(ncid, nlon_id, &nlon))) ERR(retval, "lon");
      if ((retval = nc_inq_dimlen(ncid, nlat_id, &nlat))) ERR(retval, "lat");
    }

    // Broadcast lon/lat sizes
    comm_.broadcast(nlon, 0);
    comm_.broadcast(nlat, 0);

    // Coordinates and land-sea mask
    std::vector<double> lon(nlon);
    std::vector<double> lat(nlat);
    std::vector<int> lsm(nlat*nlon);

    if (comm_.rank() == 0) {
      // Get lon/lat
      if ((retval = nc_inq_varid(ncid, "lon", &lon_id))) ERR(retval, "lon");
      if ((retval = nc_inq_varid(ncid, "lat", &lat_id))) ERR(retval, "lat");
      if ((retval = nc_inq_varid(ncid, "LSMASK", &lsm_id))) ERR(retval, "LMASK");

      // Read data
      std::vector<float> zlon(nlon);
      std::vector<float> zlat(nlat);
      std::vector<uint8_t> zlsm(nlat*nlon);
      if ((retval = nc_get_var_float(ncid, lon_id, zlon.data()))) ERR(retval, "lon");
      if ((retval = nc_get_var_float(ncid, lat_id, zlat.data()))) ERR(retval, "lat");
      if ((retval = nc_get_var_ubyte(ncid, lsm_id, zlsm.data()))) ERR(retval, "LMASK");

      // Copy data
      for (size_t ilon = 0; ilon < nlon; ++ilon) {
        lon[ilon] = static_cast<double>(zlon[ilon]);
      }
      for (size_t ilat = 0; ilat < nlat; ++ilat) {
        lat[ilat] = static_cast<double>(zlat[ilat]);
      }
      for (size_t ilat = 0; ilat < nlat; ++ilat) {
       for (size_t ilon = 0; ilon < nlon; ++ilon) {
          lsm[ilat*nlon+ilon] = static_cast<int>(zlsm[ilat*nlon+ilon]);
        }
      }

      // Close file
      if ((retval = nc_close(ncid))) ERR(retval, ncFilePath);
    }

    // Broadcast coordinates and land-sea mask
    comm_.broadcast(lon.begin(), lon.end(), 0);
    comm_.broadcast(lat.begin(), lat.end(), 0);
    comm_.broadcast(lsm.begin(), lsm.end(), 0);

    // Build KD-tree
    atlas::Geometry geometry(atlas::util::Earth::radius());
    atlas::util::IndexKDTree2D search(geometry);
    search.reserve(nlat*nlon);
    std::vector<double> lon2d;
    std::vector<double> lat2d;
    std::vector<size_t> payload2d;
    int jnode = 0;
    for (size_t ilat = 0; ilat < nlat; ++ilat) {
      for (size_t ilon = 0; ilon < nlon; ++ilon) {
        lon2d.push_back(lon[ilon]);
        lat2d.push_back(lat[ilat]);
        payload2d.push_back(jnode);
        ++jnode;
      }
    }
    search.build(lon2d, lat2d, payload2d);

    if (functionSpace_.type() == "StructuredColumns") {
      // StructuredColumns
      atlas::functionspace::StructuredColumns fs(functionSpace_);
      auto lonlatView = atlas::array::make_view<double, 2>(fs.xy());
      auto maskView = atlas::array::make_view<int, 2>(gmask);
      for (atlas::idx_t jnode = 0; jnode < fs.xy().shape(0); ++jnode) {
        if (ghostView(jnode) == 0) {
          // Find nearest neighbor
          size_t nn = search.closestPoint(atlas::PointLonLat{lonlatView(jnode, 0),
            lonlatView(jnode, 1)}).payload();

          // Ocean points for all levels
          for (size_t jlevel = 0; jlevel < group.levels_; ++jlevel) {
            if (lsm[nn] == 0) {
               maskView(jnode, jlevel) = 1;
             } else {
               maskView(jnode, jlevel) = 0;
             }
           }

          // Ocean + small islands for:
          // - the first level of 3D fields,
          // - the 2D fields if lev2d = "first"
          if (lsm[nn] == 3) {
            if ((group.levels_ > 1) || (group.lev2d_ == "first")) {
              maskView(jnode, 0) = 1;
            }
          }
        }
      }
    } else {
      throw eckit::NotImplemented("Sea mask not supported for " + functionSpace_.type() + " yet",
        Here());
    }
  } else {
    throw eckit::UserError("Wrong mask type", Here());
  }
  fields_->add(gmask);

  // Mask size
  group.gmaskSize_ = 0.0;
  size_t domainSize = 0.0;
  for (atlas::idx_t jnode = 0; jnode < gmask.shape(0); ++jnode) {
    for (atlas::idx_t jlevel = 0; jlevel < gmask.shape(1); ++jlevel) {
      if (ghostView(jnode) == 0) {
        if (maskView(jnode, jlevel) == 1) {
          group.gmaskSize_ += 1.0;
        }
        domainSize++;
      }
    }
  }
  comm_.allReduceInPlace(group.gmaskSize_, eckit::mpi::sum());
  comm_.allReduceInPlace(domainSize, eckit::mpi::sum());
  if (domainSize > 0) {
    group.gmaskSize_ = group.gmaskSize_/static_cast<double>(domainSize);
  }

  oops::Log::trace() << classname() << "::setupMask done" << std::endl;
}

// -----------------------------------------------------------------------------

void Geometry::checkLonLat(const eckit::Configuration & checkLonLatConf) const {
  oops::Log::trace() << classname() << "::checkLonLat starting" << std::endl;

  // Return if configuration is empty
  if (checkLonLatConf.empty()) {
    return;
  }

  // Get variable to read
  const std::string lonName = checkLonLatConf.getString("longitude", "longitude");
  const std::string latName = checkLonLatConf.getString("latitude", "latitude");
  const oops::Variables lonLatVars(std::vector<std::string>({lonName, latName}));

  // Create field
  Fields field(*this, lonLatVars, util::DateTime());

  // Read field
  field.read(checkLonLatConf);

  // Get views
  const auto lonView = atlas::array::make_view<double, 2>(field.fieldSet()[lonName]);
  const auto latView = atlas::array::make_view<double, 2>(field.fieldSet()[latName]);

  // Get lon/lat view
  const auto lonlatView = atlas::array::make_view<double, 2>(functionSpace_.lonlat());

  // Get ghost view
  const auto ghostView = atlas::array::make_view<int, 1>(functionSpace_.ghost());

  // Check lon/lat
  for (atlas::idx_t jnode = 0; jnode < functionSpace_.lonlat().shape(0); ++jnode) {
    if (ghostView(jnode) == 0) {
      if (std::abs(lonView(jnode, 0)-lonlatView(jnode, 0)) > 1.0e-6) {
        std::cout << lonView(jnode, 0) << " = " << lonlatView(jnode, 0) << std::endl;
        throw eckit::Exception("inaccurate longitude", Here());
      }
      if (std::abs(latView(jnode, 0)-lonlatView(jnode, 1)) > 1.0e-6) {
        std::cout << latView(jnode, 0) << " = " << lonlatView(jnode, 1) << std::endl;
        throw eckit::Exception("inaccurate latitude", Here());
      }
    }
  }

  oops::Log::trace() << classname() << "::checkLonLat starting" << std::endl;
}

// -----------------------------------------------------------------------------

GeometryIterator Geometry::begin() const {
  return GeometryIterator(*this, 0, 0);
}

// -----------------------------------------------------------------------------

GeometryIterator Geometry::end() const {
  return GeometryIterator(*this, nnodes_, nlevs_);
}

// -----------------------------------------------------------------------------

std::vector<double> Geometry::verticalCoord(std::string & vcUnits) const {
  return vertCoordAvg_;
}

// -----------------------------------------------------------------------------

}  // namespace quenchxx
