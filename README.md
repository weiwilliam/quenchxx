# QUENCHXX

### Licence

(C) Copyright 2024 Meteorologisk Institutt

[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

### Description

Extension of the saber/quench, a file-based and ATLAS-based model interfaced with OOPS.

### Installation

ECBUILD bundle:

```cmake
# (C) Copyright 2022 UCAR
#

cmake_minimum_required( VERSION 3.23 FATAL_ERROR )

find_package( ecbuild 3.6 REQUIRED HINTS ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/../ecbuild)

project( quenchxx-bundle VERSION 8.0.0 LANGUAGES C CXX Fortran )

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

include( ecbuild_bundle )

# Default release mode
set( ECBUILD_DEFAULT_BUILD_TYPE Release )

# Enable OpenMP and MPI
set( ENABLE_MPI ON CACHE BOOL "Compile with MPI" )
set( ENABLE_OMP ON CACHE BOOL "Compile with OpenMP" )
set( ENABLE_CUDA OFF CACHE BOOL "Compile with CUDA" )

# Depend path for non-ecbuild packages
set(DEPEND_LIB_ROOT ${CMAKE_CURRENT_BINARY_DIR}/Depends)
list(APPEND CMAKE_PREFIX_PATH ${DEPEND_LIB_ROOT})

# Library path for non-ecbuild packages
link_directories(${CMAKE_CURRENT_BINARY_DIR}/lib)

include( GNUInstallDirs )
if(APPLE)
  list( APPEND CMAKE_INSTALL_RPATH $ENV{llvm_openmp_ROOT}/lib )
endif()
list( APPEND CMAKE_INSTALL_RPATH ${CMAKE_CURRENT_BINARY_DIR}/fv3 )

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)

# when building, already use the install RPATH
set(CMAKE_BUILD_WITH_INSTALL_RPATH ON)

# Define bundle
ecbuild_bundle_initialize()

# Use external jedi-cmake or build in bundle
if(DEFINED ENV{jedi_cmake_ROOT})
  include( $ENV{jedi_cmake_ROOT}/share/jedicmake/Functions/git_functions.cmake )
else()
  ecbuild_bundle( PROJECT jedicmake GIT "https://github.com/jcsda-internal/jedi-cmake.git" BRANCH develop UPDATE RECURSIVE )
  include( jedicmake/cmake/Functions/git_functions.cmake )
endif()

# ECMWF tools
ecbuild_bundle( PROJECT eckit    GIT "https://github.com/ecmwf/eckit.git"                 TAG 1.26.0 )
ecbuild_bundle( PROJECT fckit    GIT "https://github.com/ecmwf/fckit.git"                 TAG 0.11.0 )
ecbuild_bundle( PROJECT fiat     GIT "https://github.com/ecmwf-ifs/fiat.git"              BRANCH main )
ecbuild_bundle( PROJECT ectrans  GIT "https://github.com/ecmwf-ifs/ectrans.git"           BRANCH main )
ecbuild_bundle( PROJECT atlas    GIT "https://github.com/ecmwf/atlas.git"                 BRANCH develop UPDATE)
ecbuild_bundle( PROJECT gsw      GIT "https://github.com/jcsda-internal/GSW-Fortran.git"  BRANCH develop UPDATE )

# OOPS / VADER / SABER
ecbuild_bundle( PROJECT oops     GIT "https://github.com/jcsda-internal/oops.git"         BRANCH develop UPDATE )
#ecbuild_bundle( PROJECT vader    GIT "https://github.com/jcsda-internal/vader.git"       BRANCH develop UPDATE )
ecbuild_bundle( PROJECT vader    GIT "https://github.com/jcsda-internal/vader.git"        BRANCH feature/recipes_hmr_at_ap UPDATE )
ecbuild_bundle( PROJECT saber    GIT "https://github.com/jcsda-internal/saber.git"        BRANCH develop UPDATE )

# Observation operators
ecbuild_bundle( PROJECT crtm     GIT "https://github.com/jcsda/CRTMv3.git"                BRANCH develop UPDATE )

option(ENABLE_IODA_DATA "Obtain ioda test data from ioda-data repository (vs tarball)" ON)
ecbuild_bundle( PROJECT ioda-data GIT "https://github.com/jcsda-internal/ioda-data.git"   BRANCH develop UPDATE )
ecbuild_bundle( PROJECT ioda     GIT "https://github.com/jcsda-internal/ioda.git"         BRANCH develop UPDATE )

option(ENABLE_UFO_DATA "Obtain ufo test data from ufo-data repository (vs tarball)" ON)
ecbuild_bundle( PROJECT ufo-data GIT "https://github.com/jcsda-internal/ufo-data.git"     BRANCH develop UPDATE )
ecbuild_bundle( PROJECT ufo      GIT "https://github.com/jcsda-internal/ufo.git"          BRANCH develop UPDATE )

# QUENCHXX
ecbuild_bundle( PROJECT quenchxx GIT "https://github.com/benjaminmenetrier/quenchxx.git   BRANCH develop UPDATE)

ecbuild_bundle_finalize()
```
