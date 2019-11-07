#User-specific source files
#list(APPEND AMR_WIND_EXTRA_SOURCES ${CMAKE_SOURCE_DIR}/Build/probdata.f90)
#list(APPEND AMR_WIND_EXTRA_SOURCES ${CMAKE_SOURCE_DIR}/Build/bc_fill_nd.F90)

#Compile-time options for executable
set(AMR_WIND_DIM 3)
set(AMR_WIND_ENABLE_EB ON)
set(AMR_WIND_ENABLE_MASA OFF)
