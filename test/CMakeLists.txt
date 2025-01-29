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

# Macro to create test
macro( create_test name mpi exec )
    if( ECSABER )
        ecbuild_add_test( TARGET quenchxx_test_${name}_${mpi}
                          COMMAND ${CMAKE_BINARY_DIR}/bin/quenchxx_compare.py
                          ARGS ${mpi}
                               1
                               ${CMAKE_BINARY_DIR}/bin/quenchxx_${exec}.x
                               testinput/ec/${name}.json
                          DEPENDS quenchxx_${exec}.x )
    else()
        ecbuild_add_test( TARGET   quenchxx_test_${name}_${mpi}
                          MPI      ${mpi}
                          COMMAND  quenchxx_${exec}.x
                          ARGS     testinput/jedi/${name}.yaml
                          DEPENDS  quenchxx_${exec}.x )
    endif()
endmacro()

# Test data
# ---------
list( APPEND quenchxx_data
testdata/geos_cf_N10_subset.20210801_0000z.nc
testdata/camchem_N30x20_subset.nc
testdata/wrfout_d01_subset_20230519_000000.nc
testdata/tropomi_no2_tropo_2020090318_m.nc4
testdata/tropomi_no2_tropo_2020090318_m_singleobs.nc4
)

# Create the testdata directory
# -----------------------------
file( MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testdata )

# Link test data
# --------------
link_files( ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${quenchxx_data} )

# Create the input/output directories
# -----------------------------------
file( MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testinput )
file( MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testinput/ec )
file( MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testinput/jedi )
file( MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testref )
file( MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testref/ec )
file( MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testref/jedi )

# Test input scripts
# ------------------
list( APPEND quenchxx_testinput

# EC

# BGOS
testinput/ec/bgos_background.json
testinput/ec/bgos_ensemble.json
testinput/ec/bgos_linear_ensemble.json
testinput/ec/bgos_linear_randomized.json
testinput/ec/bgos_makeobs.json
testinput/ec/bgos_nonlinear_ensemble.json
testinput/ec/bgos_nonlinear_randomized.json
testinput/ec/bgos_stddev.json
testinput/ec/bgos_truth.json

# EVIL
testinput/ec/evil_assim_ensemble.json
testinput/ec/evil_background.json
testinput/ec/evil_control_3dvar.json
testinput/ec/evil_control_update_D.json
testinput/ec/evil_control_update_R.json
testinput/ec/evil_control_update_S.json
testinput/ec/evil_dual_3dvar.json
testinput/ec/evil_dual_update_D.json
testinput/ec/evil_dual_update_R.json
testinput/ec/evil_dual_update_S.json
testinput/ec/evil_ensemble.json
testinput/ec/evil_ensemble_variance.json
testinput/ec/evil_makeobs.json
testinput/ec/evil_primal_3dvar.json
testinput/ec/evil_primal_update_D.json
testinput/ec/evil_primal_update_R.json
testinput/ec/evil_primal_update_S.json
testinput/ec/evil_stddev.json
testinput/ec/evil_truth.json

# GLOBAL
testinput/ec/glb_3dvar.json
testinput/ec/glb_3densvar.json
testinput/ec/glb_4densvar.json
testinput/ec/glb_background_06.json
testinput/ec/glb_background_12.json
testinput/ec/glb_background_18.json
testinput/ec/glb_ensemble_06.json
testinput/ec/glb_ensemble_12.json
testinput/ec/glb_ensemble_18.json
testinput/ec/glb_hofx.json
testinput/ec/glb_interfaces.json
testinput/ec/glb_letkf_linear.json
testinput/ec/glb_letkf_linear_4d.json
testinput/ec/glb_letkf_nonlinear.json
testinput/ec/glb_letkf_nonlinear_4d.json
testinput/ec/glb_letkf_read_members.json
testinput/ec/glb_makeobs_06.json
testinput/ec/glb_makeobs_12.json
testinput/ec/glb_makeobs_18.json
testinput/ec/glb_stddev.json
testinput/ec/glb_truth_06.json
testinput/ec/glb_truth_12.json
testinput/ec/glb_truth_18.json

# REGIONAL
testinput/ec/reg_3dvar.json
testinput/ec/reg_3densvar.json
testinput/ec/reg_4densvar.json
testinput/ec/reg_background_06.json
testinput/ec/reg_background_12.json
testinput/ec/reg_background_18.json
testinput/ec/reg_ensemble_06.json
testinput/ec/reg_ensemble_12.json
testinput/ec/reg_ensemble_18.json
testinput/ec/reg_hofx.json
testinput/ec/reg_interfaces.json
testinput/ec/reg_letkf_linear.json
testinput/ec/reg_letkf_linear_4d.json
testinput/ec/reg_letkf_nonlinear.json
testinput/ec/reg_letkf_nonlinear_4d.json
testinput/ec/reg_letkf_read_members.json
testinput/ec/reg_makeobs_06.json
testinput/ec/reg_makeobs_12.json
testinput/ec/reg_makeobs_18.json
testinput/ec/reg_stddev.json
testinput/ec/reg_truth_06.json
testinput/ec/reg_truth_12.json
testinput/ec/reg_truth_18.json

# SYNTHETIC DATA
testinput/ec/synthetic_data_generate.json
testinput/ec/synthetic_data_makeobs.json
testinput/ec/synthetic_data_stddev.json
testinput/ec/synthetic_data_truth.json

# JEDI

# GENINT
testinput/jedi/genint_geometry_regular_lonlat.yaml
testinput/jedi/genint_geometry_lambertCC.yaml
testinput/jedi/genint_convertstate_regular_lonlat.yaml
testinput/jedi/genint_convertstate_nonuni_lonlat.yaml
testinput/jedi/genint_convertstate_varcha_regular_lonlat.yaml
testinput/jedi/genint_convertstate_varcha_nonuni_lonlat.yaml
testinput/jedi/genint_hofx3d_regular_lonlat.yaml
testinput/jedi/genint_hofx3d_nonuni_lonlat.yaml
testinput/jedi/genint_hofx3d_lambertCC.yaml

# GLOBAL
testinput/jedi/glb_stddev.yaml
testinput/jedi/glb_truth_06.yaml
testinput/jedi/glb_truth_12.yaml
testinput/jedi/glb_truth_18.yaml
testinput/jedi/glb_background_06.yaml
testinput/jedi/glb_background_12.yaml
testinput/jedi/glb_background_18.yaml
testinput/jedi/glb_ensemble_06.yaml
testinput/jedi/glb_ensemble_12.yaml
testinput/jedi/glb_ensemble_18.yaml
testinput/jedi/glb_makeobs_06.yaml
testinput/jedi/glb_makeobs_12.yaml
testinput/jedi/glb_makeobs_18.yaml
testinput/jedi/glb_3dvar.yaml
testinput/jedi/glb_3densvar.yaml
testinput/jedi/glb_4densvar.yaml
testinput/jedi/glb_letkf_linear.yaml
testinput/jedi/glb_letkf_linear_4d.yaml
testinput/jedi/glb_letkf_nonlinear.yaml
testinput/jedi/glb_letkf_nonlinear_4d.yaml
testinput/jedi/glb_letkf_read_members.yaml

# REGIONAL
testinput/jedi/reg_stddev.yaml
testinput/jedi/reg_truth_06.yaml
testinput/jedi/reg_truth_12.yaml
testinput/jedi/reg_truth_18.yaml
testinput/jedi/reg_background_06.yaml
testinput/jedi/reg_background_12.yaml
testinput/jedi/reg_background_18.yaml
testinput/jedi/reg_ensemble_06.yaml
testinput/jedi/reg_ensemble_12.yaml
testinput/jedi/reg_ensemble_18.yaml
testinput/jedi/reg_makeobs_06.yaml
testinput/jedi/reg_makeobs_12.yaml
testinput/jedi/reg_makeobs_18.yaml
testinput/jedi/reg_3dvar.yaml
testinput/jedi/reg_3densvar.yaml
testinput/jedi/reg_4densvar.yaml
testinput/jedi/reg_letkf_linear.yaml
testinput/jedi/reg_letkf_linear_4d.yaml
testinput/jedi/reg_letkf_nonlinear.yaml
testinput/jedi/reg_letkf_nonlinear_4d.yaml
testinput/jedi/reg_letkf_read_members.yaml
)

# Test output references
# ----------------------
list( APPEND quenchxx_testref

# COMMON

# GLOBAL
testref/glb_background_06.ref
testref/glb_background_12.ref
testref/glb_background_18.ref
testref/glb_ensemble_06.ref
testref/glb_ensemble_12.ref
testref/glb_ensemble_18.ref
testref/glb_makeobs_06.ref
testref/glb_makeobs_12.ref
testref/glb_makeobs_18.ref
testref/glb_stddev.ref
testref/glb_truth_06.ref
testref/glb_truth_12.ref
testref/glb_truth_18.ref

# REGIONAL
testref/reg_background_06.ref
testref/reg_background_12.ref
testref/reg_background_18.ref
testref/reg_ensemble_06.ref
testref/reg_ensemble_12.ref
testref/reg_ensemble_18.ref
testref/reg_makeobs_06.ref
testref/reg_makeobs_12.ref
testref/reg_makeobs_18.ref
testref/reg_stddev.ref
testref/reg_truth_06.ref
testref/reg_truth_12.ref
testref/reg_truth_18.ref

# EC

# BGOS
testref/ec/bgos_background.ref
testref/ec/bgos_ensemble.ref
testref/ec/bgos_linear_ensemble.ref
testref/ec/bgos_linear_randomized.ref
testref/ec/bgos_makeobs.ref
testref/ec/bgos_nonlinear_ensemble.ref
testref/ec/bgos_nonlinear_randomized.ref
testref/ec/bgos_stddev.ref
testref/ec/bgos_truth.ref

# EVIL
testref/ec/evil_assim_ensemble.ref
testref/ec/evil_background.ref
testref/ec/evil_control_3dvar.ref
testref/ec/evil_control_update_D.ref
testref/ec/evil_control_update_R.ref
testref/ec/evil_control_update_S.ref
testref/ec/evil_dual_3dvar.ref
testref/ec/evil_dual_update_D.ref
testref/ec/evil_dual_update_R.ref
testref/ec/evil_dual_update_S.ref
testref/ec/evil_ensemble.ref
testref/ec/evil_ensemble_variance.ref
testref/ec/evil_makeobs.ref
testref/ec/evil_primal_3dvar.ref
testref/ec/evil_primal_update_D.ref
testref/ec/evil_primal_update_R.ref
testref/ec/evil_primal_update_S.ref
testref/ec/evil_stddev.ref
testref/ec/evil_truth.ref

# GLOBAL
testref/ec/glb_3dvar.ref
testref/ec/glb_3densvar.ref
testref/ec/glb_4densvar.ref
testref/ec/glb_hofx.ref
testref/ec/glb_letkf_linear.ref
testref/ec/glb_letkf_linear_4d.ref
testref/ec/glb_letkf_nonlinear.ref
testref/ec/glb_letkf_nonlinear_4d.ref
testref/ec/glb_letkf_read_members.ref

# REGIONAL
testref/ec/reg_3dvar.ref
testref/ec/reg_3densvar.ref
testref/ec/reg_4densvar.ref
testref/ec/reg_hofx.ref
testref/ec/reg_letkf_linear.ref
testref/ec/reg_letkf_linear_4d.ref
testref/ec/reg_letkf_nonlinear.ref
testref/ec/reg_letkf_nonlinear_4d.ref
testref/ec/reg_letkf_read_members.ref

# SYNTHETIC DATA
testref/ec/synthetic_data_generate.ref
testref/ec/synthetic_data_makeobs.ref
testref/ec/synthetic_data_stddev.ref
testref/ec/synthetic_data_truth.ref

# JEDI

# GENINT
testref/jedi/genint_geometry_regular_lonlat.ref
testref/jedi/genint_geometry_lambertCC.ref
testref/jedi/genint_convertstate_regular_lonlat.ref
testref/jedi/genint_convertstate_nonuni_lonlat.ref
testref/jedi/genint_convertstate_varcha_regular_lonlat.ref
testref/jedi/genint_convertstate_varcha_nonuni_lonlat.ref
testref/jedi/genint_hofx3d_regular_lonlat.ref
testref/jedi/genint_hofx3d_nonuni_lonlat.ref
testref/jedi/genint_hofx3d_lambertCC.ref

# GLOBAL
testref/jedi/glb_3dvar.ref
testref/jedi/glb_3densvar.ref
testref/jedi/glb_4densvar.ref
testref/jedi/glb_letkf_linear.ref
testref/jedi/glb_letkf_linear_4d.ref
testref/jedi/glb_letkf_nonlinear.ref
testref/jedi/glb_letkf_nonlinear_4d.ref
testref/jedi/glb_letkf_read_members.ref

# REGIONAL
testref/jedi/reg_3dvar.ref
testref/jedi/reg_3densvar.ref
testref/jedi/reg_4densvar.ref
testref/jedi/reg_letkf_linear.ref
testref/jedi/reg_letkf_linear_4d.ref
testref/jedi/reg_letkf_nonlinear.ref
testref/jedi/reg_letkf_nonlinear_4d.ref
testref/jedi/reg_letkf_read_members.ref
)

# Link input configurations and output references
# -----------------------------------------------
link_files( ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${quenchxx_testinput} )
link_files( ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${quenchxx_testref} )

# Interface tests
# ---------------

if( crtm_FOUND AND NOT ECSABER )
    # GENINT
    create_test( genint_geometry_regular_lonlat 6 geometry )
    create_test( genint_geometry_lambertCC 6 geometry )
endif()

# Integration tests
# -----------------

if( crtm_FOUND AND NOT ECSABER )
    # GENINT tests
    create_test( genint_convertstate_regular_lonlat 6 convertstate )
    create_test( genint_convertstate_nonuni_lonlat 6 convertstate )
    create_test( genint_convertstate_varcha_regular_lonlat 6 convertstate )
    create_test( genint_convertstate_varcha_nonuni_lonlat 6 convertstate )
    create_test( genint_hofx3d_regular_lonlat 6 hofx3d )
    create_test( genint_hofx3d_nonuni_lonlat 6 hofx3d )
    create_test( genint_hofx3d_lambertCC 6 hofx3d )
endif()

foreach( mpi "1" "4" )
    # GLOBAL and REGIONAL tests
    foreach( domain "glb" "reg" )
        create_test( ${domain}_stddev ${mpi} convertstate )
        foreach( hh "06" "12" "18" )
            create_test( ${domain}_truth_${hh} ${mpi} convertstate )
            create_test( ${domain}_background_${hh} ${mpi} error_covariance_toolbox )
            create_test( ${domain}_ensemble_${hh} ${mpi} error_covariance_toolbox )
            if( ECSABER )
                create_test( ${domain}_makeobs_${hh} ${mpi} makeobs_patched )
            else()
                create_test( ${domain}_makeobs_${hh} ${mpi} hofx3d )
            endif()
        endforeach()

        create_test( ${domain}_3dvar ${mpi} variational )
        create_test( ${domain}_3densvar ${mpi} variational )
        create_test( ${domain}_4densvar ${mpi} variational )
        create_test( ${domain}_letkf_linear ${mpi} letkf )
        create_test( ${domain}_letkf_linear_4d ${mpi} letkf )
        create_test( ${domain}_letkf_nonlinear ${mpi} letkf )
        create_test( ${domain}_letkf_nonlinear_4d ${mpi} letkf )
        create_test( ${domain}_letkf_read_members ${mpi} letkf )
    endforeach()

    if( ECSABER )
        # BGOS tests
        create_test( bgos_truth ${mpi} convertstate )
        create_test( bgos_stddev ${mpi} convertstate )
        create_test( bgos_background ${mpi} error_covariance_toolbox )
        create_test( bgos_ensemble ${mpi} error_covariance_toolbox )
        create_test( bgos_makeobs ${mpi} makeobs_patched )
        create_test( bgos_linear_ensemble ${mpi} bgos )
        create_test( bgos_linear_randomized ${mpi} bgos )
        create_test( bgos_nonlinear_ensemble ${mpi} bgos )
        create_test( bgos_nonlinear_randomized ${mpi} bgos )

        # EVIL tests
        create_test( evil_truth ${mpi} convertstate )
        create_test( evil_stddev ${mpi} convertstate )
        create_test( evil_background ${mpi} error_covariance_toolbox )
        create_test( evil_ensemble ${mpi} error_covariance_toolbox )
        create_test( evil_makeobs ${mpi} makeobs_patched )
        create_test( evil_control_3dvar ${mpi} variational )
        create_test( evil_control_update_D ${mpi} evil_update )
        create_test( evil_control_update_R ${mpi} evil_update )
        create_test( evil_control_update_S ${mpi} evil_update )
        create_test( evil_dual_3dvar ${mpi} variational )
        create_test( evil_dual_update_D ${mpi} evil_update )
        create_test( evil_dual_update_R ${mpi} evil_update )
        create_test( evil_dual_update_S ${mpi} evil_update )
        create_test( evil_primal_3dvar ${mpi} variational )
        create_test( evil_primal_update_D ${mpi} evil_update )
        create_test( evil_primal_update_R ${mpi} evil_update )
        create_test( evil_primal_update_S ${mpi} evil_update )
        create_test( evil_assim_ensemble ${mpi} assim_ensemble_patched )
        create_test( evil_ensemble_variance ${mpi} ensemble_variance )

        # SYNTHETIC DATA tests
        create_test( synthetic_data_truth ${mpi} convertstate )
        create_test( synthetic_data_stddev ${mpi} convertstate )
        create_test( synthetic_data_makeobs ${mpi} makeobs_patched )
        create_test( synthetic_data_generate ${mpi} synthetic_data )
    endif()
endforeach()
