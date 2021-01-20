FROM exawind/exw-dev-deps:latest as base

WORKDIR /amr-wind
COPY . /amr-wind

ARG ENABLE_MPI=ON
ARG ENABLE_CUDA=OFF
ARG ENABLE_HIP=OFF
ARG ENABLE_DPCPP=OFF

RUN (\
    cmake \
        -Bbuild \
        -DCMAKE_INSTALL_PREFIX=/opt/exawind \
        -DAMR_WIND_ENABLE_MPI=${ENABLE_MPI} \
        -DAMR_WIND_ENABLE_CUDA=${ENABLE_CUDA} \
        -DAMR_WIND_ENABLE_HIP=${ENABLE_HIP} \
        -DAMR_WIND_ENABLE_DPCPP=${ENABLE_DPCPP} \
        -DAMR_WIND_ENABLE_OPENMP=OFF \
        -DAMR_WIND_ENABLE_NETCDF=ON \
        -DAMR_WIND_ENABLE_HYPRE=ON \
        -DAMR_WIND_ENABLE_MASA=OFF \
        -DAMR_WIND_ENABLE_TESTS=OFF . \
    && cd build \
    && make -j$(nproc) \
    && ./amr_wind_unit_tests \
    )
