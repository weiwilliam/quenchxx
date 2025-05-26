# (C) Copyright 2024 Meteorologisk Institutt
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

################################################################################
# QUENCHXX
################################################################################

cmake_minimum_required( VERSION 3.23 )

project( quenchxx VERSION 0.0.1 LANGUAGES C CXX Fortran )

################################################################################
# Update source from saber/quench
################################################################################

# Update flag
set( QUENCHXX_UPDATE OFF )
if( DEFINED ENV{QUENCHXX_UPDATE} )
    set( QUENCHXX_UPDATE $ENV{QUENCHXX_UPDATE} )
endif()

if( QUENCHXX_UPDATE )
    message( STATUS "Update QUENCHXX source from saber/quench" )

    # Remove update script
    file( REMOVE ${CMAKE_CURRENT_SOURCE_DIR}/patch_update.sh )

    # Create update script
    file( WRITE ${CMAKE_CURRENT_SOURCE_DIR}/patch_update.sh "#!/usr/bin/env bash\n" )

    # Define source directory
    if( ECSABER )
        set( QUENCH_SRC ${CMAKE_SOURCE_DIR}/ecsaber/saber-jedi/quench/src )
    else()
        set( QUENCH_SRC ${CMAKE_SOURCE_DIR}/saber/quench/src )
    endif()

    # Run check_diff_recursively script
    execute_process( COMMAND bash
                             ${CMAKE_BINARY_DIR}/bin/quenchxx_check_diff_recursively.sh
                             ${CMAKE_BINARY_DIR}/bin/quenchxx_check_diff.sh
                             ${QUENCH_SRC} 
                             ${CMAKE_CURRENT_SOURCE_DIR}/src/quenchxx
                             "QUENCHXX"
                             ${CMAKE_CURRENT_SOURCE_DIR}/patch_update.sh )

    # Run update script
    execute_process( COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/patch_update.sh
                     COMMAND_ERROR_IS_FATAL ANY )

    # Remove update script
    file( REMOVE ${CMAKE_CURRENT_SOURCE_DIR}/patch_update.sh )
else()
    message( STATUS "Do not update QUENCHXX source from saber/quench" )
endif()

################################################################################

option(OPENMP "Build quenchxx with OpenMP support" ON)

set( CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake;${CMAKE_MODULE_PATH})

set( CMAKE_DIRECTORY_LABELS "quenchxx" )

find_package( ecbuild 3.6 REQUIRED )

set( ECBUILD_DEFAULT_BUILD_TYPE Release )
set( ENABLE_OS_TESTS           OFF CACHE BOOL "Disable OS tests" FORCE )
set( ENABLE_LARGE_FILE_SUPPORT OFF CACHE BOOL "Disable testing of large file support" FORCE )
set( ENABLE_MPI ON CACHE BOOL "Compile with MPI" )
set( ENABLE_OBS ON CACHE BOOL "Enable observations" )

include( ecbuild_system NO_POLICY_SCOPE )

ecbuild_requires_macro_version( 2.5 )

################################################################################
# Project
################################################################################

ecbuild_declare_project()

ecbuild_enable_fortran( REQUIRED )

set( QUENCHXX_LINKER_LANGUAGE CXX )

################################################################################
# Dependencies
################################################################################

# jedi-cmake
find_package( jedicmake QUIET )

# NetCDF
find_package( NetCDF REQUIRED COMPONENTS C Fortran )

# atlas
find_package( atlas 0.35.0 REQUIRED )

# eccodes
find_package( eccodes QUIET )

# ectrans and falfilfa
find_package( ECTRANS 1.6.0 COMPONENTS etrans double QUIET )
find_package( falfilfa QUIET )
if( ( ECTRANS_FOUND AND ectrans_HAVE_ETRANS ) AND falfilfa_FOUND )
    add_definitions( -DREADFA )
endif()

# oops
if( ECSABER )
    find_package( oops REQUIRED )
    find_package( oops_patch REQUIRED )
else()
    find_package( oops 1.10.0 REQUIRED )
endif()

# saber
find_package( saber 1.10.0 REQUIRED )

if( ( NOT ECSABER ) AND ENABLE_OBS )
    # crtm
    find_package( crtm 2.2.3 REQUIRED )

    # ufo
    find_package( ufo 1.10.0 REQUIRED )
endif()

# vader
find_package( vader 1.7.0 REQUIRED )

# OpenMP
if(OPENMP)
  find_package( OpenMP COMPONENTS Fortran )
endif()

# MPI
ecbuild_find_mpi( COMPONENTS CXX Fortran REQUIRED )

# QUENCHXX required public libraries
if( ECSABER )
    list( APPEND QUENCHXX_PUBLIC_LIBS "oops" "oops_patch" "saber" "vader" )
else()
    list( APPEND QUENCHXX_PUBLIC_LIBS "oops" "saber" "vader" )
    if( ENABLE_OBS )
        list( APPEND QUENCHXX_PUBLIC_LIBS "ufo" )
    endif()
endif()

################################################################################
# Export package info
################################################################################

list( APPEND QUENCHXX_TPLS atlas atlas_f oops saber )

set( QUENCHXX_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/src
                           ${CMAKE_CURRENT_BINARY_DIR}/src
                           ${DEPEND_LIB_ROOT}/include
                           ${DEPEND_LIB_ROOT}/mod
                           ${DEPEND_LIB_ROOT}/include/quenchxx
                           ${CMAKE_Fortran_MODULE_DIRECTORY})
set( QUENCHXX_LIBRARIES quenchxx )

get_directory_property( QUENCH_DEFINITIONS COMPILE_DEFINITIONS )

foreach( _tpl ${QUENCHXX_TPLS} )
  string( TOUPPER ${_tpl} TPL )
  list( APPEND QUENCHXX_EXTRA_DEFINITIONS   ${${TPL}_DEFINITIONS}  ${${TPL}_TPL_DEFINITIONS}  )
  list( APPEND QUENCHXX_EXTRA_INCLUDE_DIRS  ${${TPL}_INCLUDE_DIRS} ${${TPL}_TPL_INCLUDE_DIRS} )
  list( APPEND QUENCHXX_EXTRA_LIBRARIES     ${${TPL}_LIBRARIES}    ${${TPL}_TPL_LIBRARIES}    )
endforeach()

################################################################################
# Sources
################################################################################

if( ECSABER )
    add_definitions(-DECSABER)
endif()

add_subdirectory( tools )
add_subdirectory( src )
add_subdirectory( test )

################################################################################
# Finalise configuration
################################################################################

# prepares a tar.gz of the sources and/or binaries
ecbuild_install_project( NAME quenchxx )

# print the summary of the configuration
ecbuild_print_summary()
