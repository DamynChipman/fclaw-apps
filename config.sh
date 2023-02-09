#!/bin/sh

# HOME : Path to home directory
HOME=/Users/damynchipman

# P4EST_PATH : Path to p4est install (i.e., ${P4EST_PATH}/include, ${P4EST_PATH}/lib, ...)
P4EST_PATH=${HOME}/packages/p4est/p4est_source_git/build/local

# PETSC_PATH : Path to PETSc install (i.e., ${PETSC_PATH}/include, ${PETSC_PATH}/lib, ...) 
PETSC_PATH=${HOME}/packages/petsc/petsc-build

# FCLAW_APP : Absolute path to source code for fclaw-apps
FCLAW_APPS=${HOME}/packages/fclaw-apps

# FORESTCLAW_PATH : Path to ForestClaw install (i.e., ${FORESTCLAW_PATH}/include, ${FORESTCLAW_PATH}/lib, ...)
FORESTCLAW_PATH=${HOME}/packages/forestclaw/forestclaw-build/local

# ELLIPTICFOREST_PATH : Path to EllipticForest install (i.e., ${ELLIPTICFOREST_PATH}/include ${ELLIPTICFOREST_PATH}/lib, ...)
# ELLIPTICFOREST_PATH=${HOME}/packages/EllipticForest/build-develop/local
ELLIPTICFOREST_PATH=${HOME}/packages/EllipticForest/build-feature-hps-restructure/local

# --=== Create build directory ===--
BUILD_DIR=build-$(git branch --show-current)
mkdir -p ${BUILD_DIR}
cd ${BUILD_DIR}

# --=== CMake configure ===--
cmake ${FCLAW_APPS} \
    -DCMAKE_BUILD_TYPE=Debug \
    -DCMAKE_CXX_COMPILER=mpic++ \
    -DCMAKE_C_COMPILER=mpicc \
    -DP4EST_PATH=${P4EST_PATH} \
    -DPETSC_PATH=${PETSC_PATH} \
    -DFORESTCLAW_PATH=${FORESTCLAW_PATH} \
    -DELLIPTICFOREST_PATH=${ELLIPTICFOREST_PATH}

echo "Now cd to $(pwd)/${BUILD_DIR} and run make to compile"