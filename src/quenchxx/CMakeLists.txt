# (C) Copyright 2024 Meteorologisk Institutt
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

list( APPEND quenchxx_src_files

Constants.h
Constants.cc
Covariance.h
Fields.h
Geometry.h
Increment.h
LinearVariableChange.h
Model.h
Model.cc
ModelAuxControl.h
ModelAuxCovariance.h
ModelAuxIncrement.h
ModelData.h
State.h
VariableChange.cc
VariableChange.h

)

ecbuild_add_library( TARGET   quenchxx
                     SOURCES  ${quenchxx_src_files}
                     PUBLIC_LIBS  oops quench saber ufo vader
                     HEADER_DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/quenchxx
                     INSTALL_HEADERS LISTED
                     LINKER_LANGUAGE CXX )

target_include_directories( quenchxx PUBLIC "$<BUILD_INTERFACE:${QUENCHXX_INCLUDE_DIRS}>" "$<BUILD_INTERFACE:${QUENCHXX_EXTRA_INCLUDE_DIRS}>" )

target_link_libraries( quenchxx PUBLIC NetCDF::NetCDF_Fortran NetCDF::NetCDF_C )
