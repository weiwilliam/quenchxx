# (C) Copyright 2024 Meteorlogisk Institutt
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Macro to link files
macro( link_files src_dir dest_dir )
    foreach( _f IN ITEMS ${ARGN} )
        execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink ${src_dir}/${_f} ${dest_dir}/${_f} )
    endforeach()
endmacro()

# Use external jedi-cmake or bundle build
if( DEFINED ENV{jedi_cmake_ROOT} )
  include( $ENV{jedi_cmake_ROOT}/share/jedicmake/Functions/git_functions.cmake )
else()
  include( ${CMAKE_SOURCE_DIR}/jedicmake/cmake/Functions/git_functions.cmake )
endif()

# Test data
# ---------
list( APPEND quenchxx_data
Data/genint/background/geos_cf_N10_subset.20210801_0000z.nc
Data/genint/background/camchem_N30x20_subset.nc
Data/genint/background/wrfout_d01_subset_20230519_000000.nc
Data/genint/observation/tropomi_no2_tropo_2020090318_m.nc4
Data/genint/observation/tropomi_no2_tropo_2020090318_m_singleobs.nc4
)

# Create the Data/ directory and subdirectories
# ---------------------------------------------
# GENINT
file( MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/genint/background )
file( MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/genint/observation )
# SYNTHETIC
file( MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/synthetic )

# Link test data
# --------------
link_files( ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${quenchxx_data} )

# Create the input/output directories
# -----------------------------------
# GENINT
file( MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testinput/genint )
file( MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testoutput/genint )
# SYNTHETIC
file( MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testinput/synthetic )
file( MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testoutput/synthetic )

# Test input scripts
# ------------------
list( APPEND quenchxx_testinput
# GENINT
testinput/genint/geometry_regular_lonlat.yaml
testinput/genint/geometry_lambertCC.yaml
testinput/genint/convertstate_regular_lonlat.yaml
testinput/genint/convertstate_nonuni_lonlat.yaml
testinput/genint/convertstate_varcha_regular_lonlat.yaml
testinput/genint/convertstate_varcha_nonuni_lonlat.yaml
testinput/genint/hofx3d_regular_lonlat.yaml
testinput/genint/hofx3d_nonuni_lonlat.yaml
testinput/genint/hofx3d_lambertCC.yaml

# SYNTHETIC
testinput/synthetic/glb_stddev.yaml
testinput/synthetic/glb_truth_06.yaml
testinput/synthetic/glb_truth_12.yaml
testinput/synthetic/glb_truth_18.yaml
testinput/synthetic/glb_background_06.yaml
testinput/synthetic/glb_background_12.yaml
testinput/synthetic/glb_background_18.yaml
testinput/synthetic/glb_ensemble_06.yaml
testinput/synthetic/glb_ensemble_12.yaml
testinput/synthetic/glb_ensemble_18.yaml
testinput/synthetic/glb_makeobs.yaml
testinput/synthetic/glb_3dvar.yaml
testinput/synthetic/glb_3densvar.yaml
testinput/synthetic/glb_4densvar.yaml
testinput/synthetic/glb_letkf.yaml
testinput/synthetic/reg_stddev.yaml
testinput/synthetic/reg_truth_06.yaml
testinput/synthetic/reg_truth_12.yaml
testinput/synthetic/reg_truth_18.yaml
testinput/synthetic/reg_background_06.yaml
testinput/synthetic/reg_background_12.yaml
testinput/synthetic/reg_background_18.yaml
testinput/synthetic/reg_ensemble_06.yaml
testinput/synthetic/reg_ensemble_12.yaml
testinput/synthetic/reg_ensemble_18.yaml
testinput/synthetic/reg_makeobs.yaml
testinput/synthetic/reg_3dvar.yaml
testinput/synthetic/reg_3densvar.yaml
testinput/synthetic/reg_4densvar.yaml
testinput/synthetic/reg_letkf.yaml
)

# Test output references
# ----------------------
list( APPEND quenchxx_testoutput
# GENINT
testoutput/genint/geometry_regular_lonlat.ref
testoutput/genint/geometry_lambertCC.ref
testoutput/genint/convertstate_regular_lonlat.ref
testoutput/genint/convertstate_nonuni_lonlat.ref
testoutput/genint/convertstate_varcha_regular_lonlat.ref
testoutput/genint/convertstate_varcha_nonuni_lonlat.ref
testoutput/genint/hofx3d_regular_lonlat.ref
testoutput/genint/hofx3d_nonuni_lonlat.ref
testoutput/genint/hofx3d_lambertCC.ref

# SYNTHETIC
testoutput/synthetic/glb_stddev.ref
testoutput/synthetic/glb_truth_06.ref
testoutput/synthetic/glb_truth_12.ref
testoutput/synthetic/glb_truth_18.ref
testoutput/synthetic/glb_background_06.ref
testoutput/synthetic/glb_background_12.ref
testoutput/synthetic/glb_background_18.ref
testoutput/synthetic/glb_ensemble_06.ref
testoutput/synthetic/glb_ensemble_12.ref
testoutput/synthetic/glb_ensemble_18.ref
testoutput/synthetic/glb_makeobs.ref
testoutput/synthetic/glb_3dvar.ref
testoutput/synthetic/glb_3densvar.ref
testoutput/synthetic/glb_4densvar.ref
testoutput/synthetic/glb_letkf.ref
testoutput/synthetic/reg_stddev.ref
testoutput/synthetic/reg_truth_06.ref
testoutput/synthetic/reg_truth_12.ref
testoutput/synthetic/reg_truth_18.ref
testoutput/synthetic/reg_background_06.ref
testoutput/synthetic/reg_background_12.ref
testoutput/synthetic/reg_background_18.ref
testoutput/synthetic/reg_ensemble_06.ref
testoutput/synthetic/reg_ensemble_12.ref
testoutput/synthetic/reg_ensemble_18.ref
testoutput/synthetic/reg_makeobs.ref
testoutput/synthetic/reg_3dvar.ref
testoutput/synthetic/reg_3densvar.ref
testoutput/synthetic/reg_4densvar.ref
testoutput/synthetic/reg_letkf.ref
)

# Link input configurations and output references
# -----------------------------------------------
link_files( ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${quenchxx_testinput} )
link_files( ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${quenchxx_testoutput} )

# Interface tests
# ---------------

if( crtm_FOUND )
    # GENINT
    ecbuild_add_test( TARGET   quenchxx_test_genint_geometry_regular_lonlat
                      MPI      6
                      ARGS     testinput/genint/geometry_regular_lonlat.yaml
                      COMMAND  quenchxx_geometry.x
                      DEPENDS  quenchxx_geometry.x )

    ecbuild_add_test( TARGET   quenchxx_test_genint_geometry_lambertCC
                      MPI      6
                      ARGS     testinput/genint/geometry_lambertCC.yaml
                      COMMAND  quenchxx_geometry.x
                      DEPENDS  quenchxx_geometry.x )
endif()

# Integration tests
# -----------------

if( crtm_FOUND )
    # GENINT
    ecbuild_add_test( TARGET   quenchxx_test_genint_convertstate_regular_lonlat
                      MPI      6
                      ARGS     testinput/genint/convertstate_regular_lonlat.yaml
                      COMMAND  quenchxx_convertstate.x
                      DEPENDS  quenchxx_convertstate.x )

    ecbuild_add_test( TARGET   quenchxx_test_genint_convertstate_nonuni_lonlat
                      MPI      6
                      ARGS     testinput/genint/convertstate_nonuni_lonlat.yaml
                      COMMAND  quenchxx_convertstate.x
                      DEPENDS  quenchxx_convertstate.x )

    ecbuild_add_test( TARGET   quenchxx_test_genint_convertstate_varcha_regular_lonlat
                      MPI      6
                      ARGS     testinput/genint/convertstate_varcha_regular_lonlat.yaml
                      COMMAND  quenchxx_convertstate.x
                      DEPENDS  quenchxx_convertstate.x )

    ecbuild_add_test( TARGET   quenchxx_test_genint_convertstate_varcha_nonuni_lonlat
                      MPI      6
                      ARGS     testinput/genint/convertstate_varcha_nonuni_lonlat.yaml
                      COMMAND  quenchxx_convertstate.x
                      DEPENDS  quenchxx_convertstate.x )

    ecbuild_add_test( TARGET   quenchxx_test_genint_hofx3d_regular_lonlat
                      MPI      6
                      ARGS     testinput/genint/hofx3d_regular_lonlat.yaml
                      COMMAND  quenchxx_hofx3d.x
                      DEPENDS  quenchxx_hofx3d.x )

    ecbuild_add_test( TARGET   quenchxx_test_genint_hofx3d_nonuni_lonlat
                      MPI      6
                      ARGS     testinput/genint/hofx3d_nonuni_lonlat.yaml
                      COMMAND  quenchxx_hofx3d.x
                      DEPENDS  quenchxx_hofx3d.x )

    ecbuild_add_test( TARGET   quenchxx_test_genint_hofx3d_lambertCC
                      MPI      6
                      ARGS     testinput/genint/hofx3d_lambertCC.yaml
                      COMMAND  quenchxx_hofx3d.x
                      DEPENDS  quenchxx_hofx3d.x )
endif()

# SYNTHETIC
foreach( DOMAIN "glb" "reg" )
  ecbuild_add_test( TARGET   quenchxx_test_synthetic_${DOMAIN}_stddev
                    MPI      1
                    ARGS     testinput/synthetic/${DOMAIN}_stddev.yaml
                    COMMAND  quenchxx_convertstate.x
                    DEPENDS  quenchxx_convertstate.x )
  
  foreach( HH "06" "12" "18" )
    ecbuild_add_test( TARGET   quenchxx_test_synthetic_${DOMAIN}_truth_${HH}
                      MPI      1
                      ARGS     testinput/synthetic/${DOMAIN}_truth_${HH}.yaml
                      COMMAND  quenchxx_convertstate.x
                      DEPENDS  quenchxx_convertstate.x )
  
    ecbuild_add_test( TARGET   quenchxx_test_synthetic_${DOMAIN}_background_${HH}
                      MPI      1
                      ARGS     testinput/synthetic/${DOMAIN}_background_${HH}.yaml
                      COMMAND  quenchxx_error_covariance_toolbox.x
                      DEPENDS  quenchxx_error_covariance_toolbox.x
                      TEST_DEPENDS synthetic_${DOMAIN}_stddev
                                   synthetic_${DOMAIN}_truth_${HH} )
  
    ecbuild_add_test( TARGET   quenchxx_test_synthetic_${DOMAIN}_ensemble_${HH}
                      MPI      1
                      ARGS     testinput/synthetic/${DOMAIN}_ensemble_${HH}.yaml
                      COMMAND  quenchxx_error_covariance_toolbox.x
                      DEPENDS  quenchxx_error_covariance_toolbox.x
                      TEST_DEPENDS synthetic_${DOMAIN}_stddev
                                   synthetic_${DOMAIN}_background_${HH} )
  endforeach()
  
  ecbuild_add_test( TARGET   quenchxx_test_synthetic_${DOMAIN}_makeobs
                    MPI      1
                    ARGS     testinput/synthetic/${DOMAIN}_makeobs.yaml
                    COMMAND  quenchxx_hofx3d.x
                    DEPENDS  quenchxx_hofx3d.x
                    TEST_DEPENDS synthetic_${DOMAIN}_truth_12 )
  
  ecbuild_add_test( TARGET   quenchxx_test_synthetic_${DOMAIN}_3dvar
                    MPI      1
                    ARGS     testinput/synthetic/${DOMAIN}_3dvar.yaml
                    COMMAND  quenchxx_variational.x
                    DEPENDS  quenchxx_variational.x
                    TEST_DEPENDS synthetic_${DOMAIN}_background_12
                                 synthetic_${DOMAIN}_makeobs )
  
  ecbuild_add_test( TARGET   quenchxx_test_synthetic_${DOMAIN}_3densvar
                    MPI      1
                    ARGS     testinput/synthetic/${DOMAIN}_3densvar.yaml
                    COMMAND  quenchxx_variational.x
                    DEPENDS  quenchxx_variational.x
                    TEST_DEPENDS synthetic_${DOMAIN}_background_12
                                 synthetic_${DOMAIN}_ensemble_12
                                 synthetic_${DOMAIN}_makeobs )
  
  ecbuild_add_test( TARGET   quenchxx_test_synthetic_${DOMAIN}_4densvar
                    MPI      1
                    ARGS     testinput/synthetic/${DOMAIN}_4densvar.yaml
                    COMMAND  quenchxx_variational.x
                    DEPENDS  quenchxx_variational.x
                    TEST_DEPENDS synthetic_${DOMAIN}_background_06
                                 synthetic_${DOMAIN}_background_12
                                 synthetic_${DOMAIN}_background_18
                                 synthetic_${DOMAIN}_ensemble_06
                                 synthetic_${DOMAIN}_ensemble_12
                                 synthetic_${DOMAIN}_ensemble_18
                                 synthetic_${DOMAIN}_makeobs )

  ecbuild_add_test( TARGET   quenchxx_test_synthetic_${DOMAIN}_letkf
                    MPI      1
                    ARGS     testinput/synthetic/${DOMAIN}_letkf.yaml
                    COMMAND  quenchxx_letkf.x
                    DEPENDS  quenchxx_letkf.x
                    TEST_DEPENDS synthetic_${DOMAIN}_background_12
                                 synthetic_${DOMAIN}_ensemble_12
                                 synthetic_${DOMAIN}_makeobs )
endforeach()
