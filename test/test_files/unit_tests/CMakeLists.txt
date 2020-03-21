add_executable(${amr_wind_unit_test_exe_name} "")

target_sources(${amr_wind_unit_test_exe_name}
   PRIVATE
   utest_main.cpp
   )


#AMR-Wind include directories
target_include_directories(${amr_wind_unit_test_exe_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

#Needed for AMReX_BuildInfo.H
target_include_directories(${amr_wind_unit_test_exe_name} SYSTEM PRIVATE ${CMAKE_SOURCE_DIR}/submods/amrex/Tools/C_scripts)
  
#Link our unit test executable with GoogleTest
#target_link_libraries(${utest_exe_name} PRIVATE GTest::GTest GTest::Main)
target_link_libraries(${amr_wind_unit_test_exe_name} PRIVATE gtest_main)
target_include_directories(${amr_wind_unit_test_exe_name} SYSTEM PRIVATE ${CMAKE_SOURCE_DIR}/submods/googletest/googletest/include)

target_link_libraries(${amr_wind_unit_test_exe_name} PRIVATE amrex)

#Define what we want to be installed during a make install 
install(TARGETS ${amr_wind_unit_test_exe_name}
        RUNTIME DESTINATION bin
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib)
