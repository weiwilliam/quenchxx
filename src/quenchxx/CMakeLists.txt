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
Increment.cc
Increment.h
Interpolation.cc
Interpolation.h
LinearVariableChange.cc
LinearVariableChange.h
LinearVariableChangeParameters.h
Model.cc
Model.h
ModelAuxControl.h
ModelAuxCovariance.h
ModelAuxIncrement.h
ModelData.h
ObsLocalization.cc
ObsLocalization.h
State.cc
State.h
Traits.h
TraitsFwd.h
VaderCookbook.h
VariableChange.cc
VariableChange.h
VariableChangeParameters.h

)

ecbuild_add_library( TARGET   quenchxx
                     SOURCES  ${quenchxx_src_files}
                     PUBLIC_LIBS  oops saber ufo vader
                     HEADER_DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/quenchxx
                     INSTALL_HEADERS LISTED
                     LINKER_LANGUAGE CXX )

target_include_directories( quenchxx PUBLIC "$<BUILD_INTERFACE:${QUENCHXX_INCLUDE_DIRS}>" "$<BUILD_INTERFACE:${QUENCHXX_EXTRA_INCLUDE_DIRS}>" )

target_link_libraries( quenchxx PUBLIC NetCDF::NetCDF_Fortran NetCDF::NetCDF_C )
