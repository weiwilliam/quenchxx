# QUENCH model

### Licence

(C) Copyright 2024 Meteorologisk Institutt

[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

### Description

File-based and ATLAS-based model interfaced with OOPS.

### Compilation and test

To use SABER with OOPS, the bundle should be as follows:

```bash
# Environment --- Edit as needed
OOPS_BUNDLE_SRC=$(pwd)
OOPS_BUNDLE_BUILD=build
OOPS_BUNDLE_INSTALL=$HOME/local

# 1. Create a CMakeLitst.txt file for the OOPS bundle
cd $OOPS_BUNDLE_SRC
cat > ./CMakeLists.txt <<EOF

cmake_minimum_required( VERSION 3.3.2 FATAL_ERROR )
find_package( ecbuild 3.6 REQUIRED )
project( oops-bundle VERSION 1.0 LANGUAGES C CXX Fortran )
set( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} )

# Enable OpenMP and MPI
set( ENABLE_MPI ON CACHE BOOL "Compile with MPI" )
set( ENABLE_OMP ON CACHE BOOL "Compile with OpenMP" )

ecbuild_bundle_initialize()
ecbuild_requires_macro_version( 2.7 )

ecbuild_bundle( PROJECT eckit    GIT "ssh://git@git.ecmwf.int:7999/ecsdk/eckit.git"    BRANCH develop )
ecbuild_bundle( PROJECT fckit    GIT "ssh://git@git.ecmwf.int:7999/ecsdk/fckit.git"    BRANCH develop )
ecbuild_bundle( PROJECT fiat     GIT "ssh://git@git.ecmwf.int:7999/ecsdk/fiat.git"     BRANCH main )
ecbuild_bundle( PROJECT ectrans  GIT "ssh://git@git.ecmwf.int:7999/ecsdk/ectrans.git"  BRANCH main )
ecbuild_bundle( PROJECT atlas    GIT "ssh://git@git.ecmwf.int:7999/ecsdk/atlas.git"    BRANCH develop )
ecbuild_bundle( PROJECT oops     GIT "ssh://git@git.ecmwf.int:7999/oops/oops.git"      BRANCH feature/develop )
ecbuild_bundle( PROJECT ecsaber  GIT "ssh://git@git.ecmwf.int:7999/sab/ecsaber.git"    BRANCH develop RECURSIVE )

ecbuild_bundle_finalize()

EOF

# 2. Create the build directory:
mkdir -p $OOPS_BUNDLE_BUILD
cd $OOPS_BUNDLE_BUILD

# 3. Run CMake
ecbuild --prefix=$OOPS_BUNDLE_INSTALL -DECSABER=ON $OOPS_BUNDLE_SRC

# 4. Compile / Install
make -j4
make install

# 5. Run CI tests
cd $OOPS_BUNDLE_BUILD
ctest
```
