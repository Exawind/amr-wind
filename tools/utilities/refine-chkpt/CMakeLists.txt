set(tool_exe_name amr_wind_refine_chkpt)

add_executable(${tool_exe_name})
target_sources(${tool_exe_name}
  PRIVATE
  RefineCheckpt.cpp
  refine_chkpt.cpp)

target_include_directories(${tool_exe_name} PRIVATE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
target_link_libraries(${tool_exe_name} PUBLIC ${amr_wind_lib_name})
set_cuda_build_properties(${tool_exe_name})

install(TARGETS ${tool_exe_name}
  RUNTIME DESTINATION bin
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)
