add_library(fc2d_hps_f OBJECT
    fortran/fc2d_hps_tag4refinement.f90
    fortran/fc2d_hps_tag4coarsening.f90
    fortran/fc2d_hps_apply_bc_default.f90
    fortran/fc2d_hps_eval_bc_default.f90
)

add_library(fc2d_hps STATIC
    fc2d_hps.cpp
    fc2d_hps_diagnostics.c
    fc2d_hps_options.c
    fc2d_hps_output_ascii.c
    fc2d_hps_physical_bc.c
    $<TARGET_OBJECTS:fc2d_hps_f>
)
target_include_directories(fc2d_hps PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${EXTERNAL_PKGS_INCLUDE})
target_link_directories(fc2d_hps PUBLIC ${EXTERNAL_PKGS_LIB})
target_link_libraries(fc2d_hps PUBLIC z sc p4est petsc forestclaw clawpatch EllipticForest blas lapack fishpack90)