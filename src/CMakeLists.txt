#
# Define a macro to add sources.
# This will also accumulate include directories if
# a header file is spotted
# WARNING: this automagically prepends the correct
# path (relative to the file invoking the macro)
# to the source file.
#
set ( INCFLO_Fortran_MOD_PATH ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/mod )
set ( INCFLO_INCLUDE_PATHS  ${PROJECT_SOURCE_DIR}/src/include ${INCFLO_Fortran_MOD_PATH})

macro ( add_sources )
  foreach ( item IN ITEMS ${ARGN} )
    set ( source_name  ${CMAKE_CURRENT_LIST_DIR}/${item} )
    string ( REPLACE "//" "/" source_name ${source_name})
    target_sources ( ${INCFLO_LIBNAME} PRIVATE ${source_name} )

    # Add folder to list of includes if necessary
    get_filename_component ( filetype ${source_name} EXT )
    if ( (${filetype} MATCHES ".h") OR (${filetype} MATCHES ".H" ) )
      list (APPEND INCFLO_INCLUDE_PATHS ${CMAKE_CURRENT_LIST_DIR})
    endif()
  endforeach()
endmacro()

#
# Fist, define the library we want to add
#
add_library ( ${INCFLO_LIBNAME} "" )

set_target_properties ( ${INCFLO_LIBNAME}
   PROPERTIES
   Fortran_MODULE_DIRECTORY
   ${INCFLO_Fortran_MOD_PATH}
   )

#
# Now, one by one, we add all the sources from all the subdirs
# to the library defined above
# 

#
# Uncategorized
# 
add_sources ( incflo_advance.cpp )
add_sources ( incflo_arrays.cpp )
add_sources ( incflo_level.cpp )
add_sources ( incflo_regrid.cpp )

add_sources ( incflo_level.H )

add_sources ( main.cpp )

include ( convection/CMakeLists.txt )
include ( derive/CMakeLists.txt )
include ( diffusion/CMakeLists.txt )
include ( divop/CMakeLists.txt )
include ( embedded_boundaries/CMakeLists.txt )
include ( init_and_boundary_conds/CMakeLists.txt )
include ( io/CMakeLists.txt )
include ( modules/CMakeLists.txt )
include ( projection/CMakeLists.txt )


# Remove duplicates from include lists
list ( REMOVE_DUPLICATES INCFLO_INCLUDE_PATHS )
target_include_directories ( ${INCFLO_LIBNAME} PUBLIC ${INCFLO_INCLUDE_PATHS} )

#
# Link to imported AMReX
#
target_link_libraries ( ${INCFLO_LIBNAME} AMReX::amrex )


#
# Set the default flags 
# 
if ( NOT CMAKE_Fortran_FLAGS )
   target_compile_options ( ${INCFLO_LIBNAME}
      PUBLIC
      # GNU Debug
      $<BUILD_INTERFACE:$<$<CONFIG:Debug>:$<$<C_COMPILER_ID:GNU>:$<$<COMPILE_LANGUAGE:Fortran>:
      -O0 -ggdb -fbounds-check -fbacktrace -Wuninitialized -Wunused -finit-real=snan -finit-integer=2147483647>>>>
      # GNU Release
      $<BUILD_INTERFACE:$<$<CONFIG:Release>:$<$<C_COMPILER_ID:GNU>:$<$<COMPILE_LANGUAGE:Fortran>:
      >>>>
      # Intel Debug
      $<BUILD_INTERFACE:$<$<CONFIG:Debug>:$<$<C_COMPILER_ID:Intel>:$<$<COMPILE_LANGUAGE:Fortran>:
      -O0 -traceback -check bounds,uninit,pointers>>>>
      # Intel Release
      $<BUILD_INTERFACE:$<$<CONFIG:Release>:$<$<C_COMPILER_ID:Intel>:$<$<COMPILE_LANGUAGE:Fortran>:
      -ip -qopt-report=5 -qopt-report-phase=vec>>>>
      # Cray Debug
      $<BUILD_INTERFACE:$<$<CONFIG:Debug>:$<$<C_COMPILER_ID:Cray>:$<$<COMPILE_LANGUAGE:Fortran>:
      -O0 -e i>>>>
      # Cray Release 
      $<BUILD_INTERFACE:$<$<CONFIG:Release>:$<$<C_COMPILER_ID:Cray>:$<$<COMPILE_LANGUAGE:Fortran>:
      >>>>
      # PGI Debug
      $<BUILD_INTERFACE:$<$<CONFIG:Debug>:$<$<C_COMPILER_ID:PGI>:$<$<COMPILE_LANGUAGE:Fortran>:
      -O0 -Mbounds -Mchkptr>>>>
      # PGI Release
      $<BUILD_INTERFACE:$<$<CONFIG:Release>:$<$<C_COMPILER_ID:PGI>:$<$<COMPILE_LANGUAGE:Fortran>:
      -gopt -fast>>>>
      )	  
endif ()

if ( NOT CMAKE_CXX_FLAGS )
   target_compile_options ( ${INCFLO_LIBNAME}
      PUBLIC
      # GNU Debug
      $<BUILD_INTERFACE:$<$<CONFIG:Debug>:$<$<C_COMPILER_ID:GNU>:$<$<COMPILE_LANGUAGE:CXX>:
      -O0 -fno-inline -ggdb -Wall -Wno-sign-compare>>>>
      # GNU Release
      $<BUILD_INTERFACE:$<$<CONFIG:Release>:$<$<C_COMPILER_ID:GNU>:$<$<COMPILE_LANGUAGE:CXX>:
      >>>>
      # Intel Debug
      $<BUILD_INTERFACE:$<$<CONFIG:Debug>:$<$<C_COMPILER_ID:Intel>:$<$<COMPILE_LANGUAGE:CXX>:
      -O0 -traceback -Wcheck>>>>
      # Intel Release
      $<BUILD_INTERFACE:$<$<CONFIG:Release>:$<$<C_COMPILER_ID:Intel>:$<$<COMPILE_LANGUAGE:CXX>:
      -ip -qopt-report=5 -qopt-report-phase=vec>>>>
      # Cray Debug
      $<BUILD_INTERFACE:$<$<CONFIG:Debug>:$<$<C_COMPILER_ID:Cray>:$<$<COMPILE_LANGUAGE:CXX>:
      -O0>>>>
      # Cray Release 
      $<BUILD_INTERFACE:$<$<CONFIG:Release>:$<$<C_COMPILER_ID:Cray>:$<$<COMPILE_LANGUAGE:CXX>:
      >>>>
      # PGI Debug
      $<BUILD_INTERFACE:$<$<CONFIG:Debug>:$<$<C_COMPILER_ID:PGI>:$<$<COMPILE_LANGUAGE:CXX>:
      -O0 -Mbounds>>>>
      # PGI Release
      $<BUILD_INTERFACE:$<$<CONFIG:Release>:$<$<C_COMPILER_ID:PGI>:$<$<COMPILE_LANGUAGE:CXX>:
      -gopt -fast>>>>
      )	  
endif ()

if (DEFINED ENABLE_FPE)
   if (ENABLE_FPE)
      target_compile_options ( ${INCFLO_LIBNAME}
   	 PUBLIC
	 # GNU
	 $<$<C_COMPILER_ID:GNU>:$<$<COMPILE_LANGUAGE:Fortran>:
	 -ffpe-trap=invalid,zero -ftrapv>>
	 $<$<C_COMPILER_ID:GNU>:$<$<COMPILE_LANGUAGE:CXX>:
	 -ftrapv>>
	 # Intel
	 $<$<C_COMPILER_ID:GNU>:$<$<COMPILE_LANGUAGE:Fortran>:
	 -fpe3>>
	 $<$<C_COMPILER_ID:GNU>:$<$<COMPILE_LANGUAGE:CXX>:
	 -fpe3>>
	 # Cray
	 $<$<C_COMPILER_ID:GNU>:$<$<COMPILE_LANGUAGE:Fortran>:
	 -K trap=fp>>
	 $<$<C_COMPILER_ID:GNU>:$<$<COMPILE_LANGUAGE:CXX>:
	 -K trap=fp>>
	 #  PGI
	 $<$<C_COMPILER_ID:GNU>:$<$<COMPILE_LANGUAGE:Fortran>:
	 -Ktrap=divz,inv>>
	 $<$<C_COMPILER_ID:GNU>:$<$<COMPILE_LANGUAGE:CXX>:
	 >> )
   endif ()
else ()
   message (AUTHOR_WARNING "Variable ENABLE_FPE is not defined")
endif ()
