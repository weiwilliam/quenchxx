# (C) Copyright 2022 UCAR.
# (C) Copyright 2023-2024 Meteorologisk Institutt
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

list( APPEND quenchxx_src_files

Constants.cc
Constants.h
Covariance.h
Fields.cc
Fields.h
Geometry.cc
Geometry.h
GeometryIterator.cc
GeometryIterator.h
GeoVaLs.cc
GeoVaLs.h
Increment.cc
Increment.h
Interpolation.cc
Interpolation.h
IncrModCtlVec.h
LinearVariableChange.cc
LinearVariableChange.h
LinearVariableChangeParameters.h
ModelData.h
ObsSpace.cc
ObsSpace.h
ObsVector.cc
ObsVector.h
State.cc
State.h
Traits.h
TraitsFwd.h
Utilities.cc
Utilities.h
VaderCookbook.h
VariableChange.cc
VariableChange.h
VariableChangeParameters.h
VariablesSwitch.h

)

if( ECSABER )
    list( APPEND quenchxx_src_files

    HorizScaleDecomposition.h
    HybridWeight.cc
    HybridWeight.h
    IncrEnsCtlVec.h
    instantiateQuenchMatrices.h
    Interpolator.h
    LinearObsOperator.cc
    LinearObsOperator.h
    LocalizationMatrix.h
    Locations.cc
    Locations.h
    Logbook.cc
    Logbook.h
    LinearModel.cc
    LinearModel.h
    Model.cc
    Model.h
    ModelAuxControl.h
    ModelAuxCovariance.h
    ModelAuxCtlVec.h
    ModelAuxControlEstimator.h
    ModelAuxIncrement.h
    ObsAuxControl.cc
    ObsAuxControl.h
    ObsAuxCovariance.cc
    ObsAuxCovariance.h
    ObsAuxCtlVec.cc
    ObsAuxCtlVec.h
    ObsAuxIncrement.cc
    ObsAuxIncrement.h
    ObsError.cc
    ObsError.h
    ObsLocalizationEC.cc
    ObsLocalizationEC.h
    ObsOperator.cc
    ObsOperator.h
    Variables.h
    Variables.cc

    )
else()
    list( APPEND quenchxx_src_files

    ObsLocalizationJEDI.cc
    ObsLocalizationJEDI.h

    )
endif()

ecbuild_add_library( TARGET   quenchxx
                     SOURCES  ${quenchxx_src_files}
                     PUBLIC_LIBS  ${QUENCHXX_PUBLIC_LIBS}
                     HEADER_DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/quenchxx
                     INSTALL_HEADERS LISTED
                     LINKER_LANGUAGE CXX )

target_include_directories( quenchxx PUBLIC "$<BUILD_INTERFACE:${QUENCHXX_INCLUDE_DIRS}>" "$<BUILD_INTERFACE:${QUENCHXX_EXTRA_INCLUDE_DIRS}>" )

target_link_libraries( quenchxx PUBLIC NetCDF::NetCDF_Fortran NetCDF::NetCDF_C )
if( eccodes_FOUND )
    target_link_libraries( quenchxx PUBLIC eccodes )
endif()

#Configure include directory layout for build-tree to match install-tree
set(QUENCHXX_BUILD_DIR_INCLUDE_PATH ${CMAKE_BINARY_DIR}/${PROJECT_NAME}/include)
add_custom_target(quenchxx_headers ALL COMMAND ${CMAKE_COMMAND} -E make_directory "${QUENCHXX_BUILD_DIR_INCLUDE_PATH}"
                                   COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR} ${QUENCHXX_BUILD_DIR_INCLUDE_PATH}/quenchxx )
add_dependencies(quenchxx quenchxx_headers)
target_include_directories(quenchxx PUBLIC $<BUILD_INTERFACE:${QUENCHXX_BUILD_DIR_INCLUDE_PATH}>
                                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
                                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../../..>)
