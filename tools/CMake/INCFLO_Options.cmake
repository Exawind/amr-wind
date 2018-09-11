#
# Here we define the default config options   
# that can be overwritten by the user         
# This file provides the followign variables  
#
# DEBUG
# INCFLO_BUILD_TYPE
# INCFLO_FFLAGS_OVERRIDES
# INCFLO_CXXFLAGS_OVERRIDES
# ENABLE_FPE
#

if (DEFINED __INCFLO_OPTIONS__)
   return ()
endif ()

# Define the following variable
# so that other included file can check if this file has been
# run already
set (__INCFLO_OPTIONS__ "")

# Creates `compile_commands.json` in the build build directory
#  `compile_commands.json` contains compiler flags used by plugins like YouCompleteMe
set ( CMAKE_EXPORT_COMPILE_COMMANDS ON )

#
# Populate the cache and check the value of the user-definable options 
#
option ( DEBUG "Build in debug mode" OFF )

if ( DEBUG )
   set ( CMAKE_BUILD_TYPE "Debug" )
else ()
   set ( CMAKE_BUILD_TYPE "Release" )
endif ()

string ( TOUPPER ${CMAKE_BUILD_TYPE} INCFLO_BUILD_TYPE )  

option ( ENABLE_FPE "Enable Floating Point Exceptions checks" OFF )

option ( ENABLE_PROJCC "Enable Approximate Projection" OFF)

option ( ENABLE_PTESTS "Include tests for projection method in Ctest suite " OFF )

option ( ENABLE_STESTS "Include tests for SIMPLE method in Ctest suite " ON )

