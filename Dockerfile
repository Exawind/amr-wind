FROM exawind/exw-cuda-dev:latest as base

WORKDIR /amr-wind
COPY . /amr-wind

ARG ENABLE_CUDA=ON
ARG ENABLE_MPI=OFF

RUN (\
    cmake \
        -Bbuild \
        -DCMAKE_INSTALL_PREFIX=/opt/exawind \
        -DAMR_WIND_ENABLE_MPI=${ENABLE_MPI} \
        -DAMR_WIND_ENABLE_CUDA=${ENABLE_CUDA} \
        -DAMR_WIND_ENABLE_OPENMP=OFF \
        -DAMR_WIND_ENABLE_MASA=OFF \
        -DAMR_WIND_ENABLE_TESTS=OFF . \
    && cd build \
    && make -j$(nproc) \
    )
