#!/bin/bash -l

# Temporary example cmake configure script

set -e

MODULES=modules
COMPILER=gcc-7.4.0

module purge
module unuse ${MODULEPATH}
module use /opt/binaries/${MODULES}
module use /opt/compilers/${MODULES}
module use /opt/utilities/${MODULES}
module use /opt/software/${MODULES}/${COMPILER}
module load unzip
module load patch
module load bzip2
module load cmake
module load git
module load flex
module load bison
module load wget
module load bc
module load binutils
module load python/3.7.3
module load gcc
module load openmpi

cmake -DCMAKE_INSTALL_PREFIX:PATH=./install \
      -DCMAKE_BUILD_TYPE:STRING=Release \
      -DCMAKE_CXX_COMPILER:STRING=mpicxx \
      -DCMAKE_C_COMPILER:STRING=mpicc \
      -DCMAKE_Fortran_COMPILER:STRING=mpifort \
      -DAMR_WIND_ENABLE_EB:BOOL=OFF \
      -DAMR_WIND_ENABLE_MPI:BOOL=ON \
      -DAMR_WIND_ENABLE_FCOMPARE:BOOL=OFF \
      -DAMR_WIND_ENABLE_FEXTREMA:BOOL=OFF \
      -DAMR_WIND_ENABLE_TESTS:BOOL=ON \
      -DAMR_WIND_TEST_WITH_FCOMPARE:BOOL=OFF \
      -DAMR_WIND_TEST_WITH_FEXTREMA:BOOL=OFF \
      -DAMR_WIND_DIM:STRING=3 \
      -DPYTHON_EXECUTABLE=$(which python3) \
      -DCMAKE_EXPORT_COMPILE_COMMANDS:BOOL=TRUE \
      .. && make -j32
