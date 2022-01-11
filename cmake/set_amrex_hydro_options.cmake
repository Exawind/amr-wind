# Set amrex hydro options

set(HYDRO_SPACEDIM 3)
set(HYDRO_EB ${AMReX_EB})
set(HYDRO_MPI ${AMR_WIND_ENABLE_MPI})
set(HYDRO_OMP ${AMR_WIND_ENABLE_OPENMP})

if (AMR_WIND_ENABLE_CUDA)
  set(HYDRO_GPU_BACKEND CUDA CACHE STRING "AMReX GPU type" FORCE)
elseif(AMR_WIND_ENABLE_ROCM)
  set(HYDRO_GPU_BACKEND HIP CACHE STRING "AMReX GPU type" FORCE)
elseif(AMR_WIND_ENABLE_DPCPP)
  set(HYDRO_GPU_BACKEND SYCL CACHE STRING "AMReX GPU type" FORCE)
else()
  set(HYDRO_GPU_BACKEND NONE CACHE STRING "AMReX GPU type" FORCE)
endif()

