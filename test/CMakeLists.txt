#=============================================================================
# AMR-Wind Testing
#=============================================================================

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/CTestCustom.cmake ${CMAKE_BINARY_DIR}/CTestCustom.cmake)

#Build GoogleTest (GoogleTest submodule must be checked out in this directory)
add_subdirectory(googletest)

add_subdirectory(test_files)

include(${CMAKE_CURRENT_SOURCE_DIR}/CTestList.cmake)

