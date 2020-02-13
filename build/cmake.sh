#!/bin/bash -l

# Temporary example cmake configure script

cmake -DCMAKE_INSTALL_PREFIX:PATH=./install \
      -DCMAKE_BUILD_TYPE:STRING=Release \
      -DCMAKE_CXX_COMPILER:STRING=mpicxx \
      -DCMAKE_C_COMPILER:STRING=mpicc \
      -DCMAKE_Fortran_COMPILER:STRING=mpifort \
      -DAMR_WIND_ENABLE_MPI:BOOL=ON \
      -DAMR_WIND_ENABLE_FCOMPARE:BOOL=OFF \
      -DAMR_WIND_ENABLE_TESTS:BOOL=OFF \
      -DAMR_WIND_TEST_WITH_FCOMPARE:BOOL=OFF \
      -DAMR_WIND_ENABLE_ALL_WARNINGS:BOOL=OFF \
      -DPYTHON_EXECUTABLE=$(which python3) \
      -DCMAKE_EXPORT_COMPILE_COMMANDS:BOOL=TRUE \
      .. && make -j8
