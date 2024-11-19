# (C) Copyright 2017-2019 UCAR
# (C) Copyright 2024 Meteorologisk Institutt
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Select version
if( ECSABER )
    set( version "ec" )
    list( APPEND libs "quenchxx" "saber" )
else()
    set( version "jedi" )
    list( APPEND libs "quenchxx" "saber" "ufo")
endif()

# Executables
if( NOT ECSABER )
    ecbuild_add_executable( TARGET  quenchxx_geometry.x
                            SOURCES ${version}/quenchxxGeometry.cc
                            LIBS    ${libs} )
endif()

ecbuild_add_executable( TARGET  quenchxx_error_covariance_toolbox.x
                        SOURCES ${version}/quenchxxErrorCovarianceToolbox.cc
                        LIBS    ${libs} )

ecbuild_add_executable( TARGET  quenchxx_process_perts.x
                        SOURCES ${version}/quenchxxProcessPerts.cc
                        LIBS    ${libs} )

ecbuild_add_executable( TARGET  quenchxx_convertstate.x
                        SOURCES ${version}/quenchxxConvertState.cc
                        LIBS    ${libs} )

if( ECSABER )
    ecbuild_add_executable( TARGET  quenchxx_hofx.x
                            SOURCES ${version}/quenchxxHofX.cc
                            LIBS    ${libs} )
else()
    ecbuild_add_executable( TARGET  quenchxx_hofx3d.x
                            SOURCES ${version}/quenchxxHofX3D.cc
                            LIBS    ${libs} )
endif()

ecbuild_add_executable( TARGET  quenchxx_variational.x
                        SOURCES ${version}/quenchxxVariational.cc
                        LIBS    ${libs} )

if( NOT ECSABER )
    ecbuild_add_executable( TARGET  quenchxx_letkf.x
                            SOURCES ${version}/quenchxxLETKF.cc
                            LIBS    ${libs} )
endif()

if( ECSABER )
    ecbuild_add_executable( TARGET  quenchxx_makeobs.x
                            SOURCES ${version}/quenchxxMakeObs.cc
                            LIBS    ${libs} )

    ecbuild_add_executable( TARGET  quenchxx_evil_update.x
                            SOURCES ${version}/quenchxxEvilUpdate.cc
                            LIBS    ${libs} )

    ecbuild_add_executable( TARGET  quenchxx_assim_ensemble.x
                            SOURCES ${version}/quenchxxAssimEnsemble.cc
                            LIBS    ${libs} )


    ecbuild_add_executable( TARGET  quenchxx_ensemble_variance.x
                            SOURCES ${version}/quenchxxEnsembleVariance.cc
                            LIBS    ${libs} )

    ecbuild_add_executable( TARGET  quenchxx_bgos.x
                            SOURCES ${version}/quenchxxBGOS.cc
                            LIBS    ${libs} )

    ecbuild_add_executable( TARGET  quenchxx_synthetic_data.x
                            SOURCES ${version}/quenchxxSyntheticData.cc
                            LIBS    ${libs} )
endif()
