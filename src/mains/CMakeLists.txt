# (C) Copyright 2017-2019 UCAR
# (C) Copyright 2024 Meteorologisk Institutt
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Select version
if( ECSABER )
    set( version "ec" )
else()
    set( version "jedi" )
endif()

# Libraries
set( LIBS ${QUENCHXX_PUBLIC_LIBS} )
list( APPEND LIBS "quenchxx" )

# OOPS-only executables
ecbuild_add_executable( TARGET  quenchxx_convertstate.x
                        SOURCES ${version}/quenchxxConvertState.cc
                        LIBS    ${LIBS} )

if( ECSABER )
    ecbuild_add_executable( TARGET  quenchxx_ensemble_variance.x
                            SOURCES ${version}/quenchxxEnsembleVariance.cc
                            LIBS    ${LIBS} )
else()
    ecbuild_add_executable( TARGET  quenchxx_geometry.x
                            SOURCES ${version}/quenchxxGeometry.cc
                            LIBS    ${LIBS} )

    ecbuild_add_executable( TARGET  quenchxx_state.x
                            SOURCES ${version}/quenchxxTestState.cc
                            LIBS    ${LIBS} )

    ecbuild_add_executable( TARGET  quenchxx_increment.x
                            SOURCES ${version}/quenchxxTestIncrement.cc
                            LIBS    ${LIBS} )

    ecbuild_add_executable( TARGET  quenchxx_ens_mean_variance.x
                            SOURCES ${version}/quenchxxEnsMeanVariance.cc
                            LIBS    ${LIBS} )
endif()

# SABER-only executables
if( ENABLE_SABER )
    ecbuild_add_executable( TARGET  quenchxx_error_covariance_toolbox.x
                            SOURCES ${version}/quenchxxErrorCovarianceToolbox.cc
                            LIBS    ${LIBS} )

    ecbuild_add_executable( TARGET  quenchxx_process_perts.x
                            SOURCES ${version}/quenchxxProcessPerts.cc
                            LIBS    ${LIBS} )
endif()

# OBS executables
if( ENABLE_OBS )
    if( ECSABER )
        ecbuild_add_executable( TARGET  quenchxx_hofx.x
                                SOURCES ${version}/quenchxxHofX.cc
                                LIBS    ${LIBS} )

        ecbuild_add_executable( TARGET  quenchxx_makeobs_patched.x
                                SOURCES ${version}/quenchxxMakeObsPatched.cc
                                LIBS    ${LIBS} )
    else()
        ecbuild_add_executable( TARGET  quenchxx_hofx3d.x
                                SOURCES ${version}/quenchxxHofX3D.cc
                                LIBS    ${LIBS} )

        ecbuild_add_executable( TARGET  quenchxx_letkf.x
                                SOURCES ${version}/quenchxxLETKF.cc
                                LIBS    ${LIBS} )
    endif()
endif()

# SABER+OBS executable
if( ENABLE_SABER AND ENABLE_OBS )
    ecbuild_add_executable( TARGET  quenchxx_variational.x
                            SOURCES ${version}/quenchxxVariational.cc
                            LIBS    ${LIBS} )

    if( ECSABER )
        ecbuild_add_executable( TARGET  quenchxx_evil_update.x
                                SOURCES ${version}/quenchxxEvilUpdate.cc
                                LIBS    ${LIBS} )

        ecbuild_add_executable( TARGET  quenchxx_assim_ensemble_patched.x
                                SOURCES ${version}/quenchxxAssimEnsemblePatched.cc
                                LIBS    ${LIBS} )

        ecbuild_add_executable( TARGET  quenchxx_bgos.x
                                SOURCES ${version}/quenchxxBGOS.cc
                                LIBS    ${LIBS} )

        ecbuild_add_executable( TARGET  quenchxx_synthetic_data.x
                                SOURCES ${version}/quenchxxSyntheticData.cc
                                LIBS    ${LIBS} )
    endif()
endif()
