# (C) Copyright 2017-2019 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 
# VADER libraries
if( vader_FOUND )
    set( vader_LIBRARIES "vader" )
else()
    set( vader_LIBRARIES "" )
endif()

# Executables
ecbuild_add_executable( TARGET  quenchxx_geometry.x
                        SOURCES quenchxxGeometry.cc
                        LIBS    quenchxx )

ecbuild_add_executable( TARGET  quenchxx_error_covariance_toolbox.x
                        SOURCES quenchxxErrorCovarianceToolbox.cc
                        LIBS    quenchxx
                                ${vader_LIBRARIES}
                                ${gsibec_LIBRARIES}
                                saber )

ecbuild_add_executable( TARGET  quenchxx_process_perts.x
                        SOURCES quenchxxProcessPerts.cc
                        LIBS    quenchxx
                                ${vader_LIBRARIES}
                                ${gsibec_LIBRARIES}
                                saber )

ecbuild_add_executable( TARGET  quenchxx_convertstate.x
                        SOURCES quenchxxConvertState.cc
                        LIBS    quenchxx
                                ${vader_LIBRARIES}
                                saber )

ecbuild_add_executable( TARGET  quenchxx_hofx3d.x
                        SOURCES quenchxxHofX3D.cc
                        LIBS    quenchxx
                                ${vader_LIBRARIES}
                                ufo )

ecbuild_add_executable( TARGET  quenchxx_variational.x
                        SOURCES quenchxxVariational.cc
                        LIBS    quenchxx
                                ${vader_LIBRARIES}
                                ${gsibec_LIBRARIES}
                                saber
                                ufo )

ecbuild_add_executable( TARGET  quenchxx_letkf.x
                        SOURCES quenchxxLETKF.cc
                        LIBS    quenchxx
                                ${vader_LIBRARIES}
                                ${gsibec_LIBRARIES}
                                ufo )
