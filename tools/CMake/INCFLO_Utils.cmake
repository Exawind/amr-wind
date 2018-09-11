#
# Check if dir or file given by path exists and issue a warning or
# error if not
#
function ( check_path  path  message_type )
  if ( EXISTS ${path} )
  else ()
    message(${message_type} ${path} " does not exist!")
  endif ( EXISTS ${path} )
endfunction ()

#
# This function turns a list into a string
#
function ( list_to_string list )
  string (REPLACE ";" " " tmp "${${list}}")
  set ( ${list} "${tmp}" PARENT_SCOPE)
endfunction ()


#
# Append new_var to all_var
#
function ( append new_var all_var )
  if ( ${new_var} )
    set ( tmp  "${${all_var}} ${${new_var}}" )

    # Since I am OCD, remove the double spaces.
    string ( REPLACE "  " " " tmp ${tmp} )
    set ( ${all_var}  ${tmp} PARENT_SCOPE )
  endif ()
endfunction ()

#
# Function to append to link line
#
function ( append_to_link_line libs link_line )

  string ( STRIP "${${libs}}" libs )

  if ( ${ARGC} EQUAL 3 )  # Third one is optional flags
    set ( flags  ${${ARGV2}} )
    string ( STRIP "${flags}" flags )
    set (tmp "${flags} ${libs}")
  else ()
    set ( flags "")
    set (tmp "${libs}")
  endif ()

  if (tmp)
    list (APPEND ${link_line} ${tmp})
    set ( ${link_line} ${${link_line}} PARENT_SCOPE )
  endif ()

endfunction ()

#
# Function to accumulate preprocessor directives
#
function ( add_define new_define all_defines )

  set ( condition  1 )

  if ( ${ARGC} EQUAL 3 ) #
    set ( condition ${${ARGV2}} )
  elseif ( ${ARGC} GREATER 3 )
    message ( AUTHOR_WARNING "Function add_define accept AT MOST 3 args" )
  endif ()

  if ( ${condition} )
    set ( ${all_defines} "${${all_defines}} -D${new_define}" PARENT_SCOPE )
    #set ( ${all_defines} ${${all_defines}} -D${new_define} PARENT_SCOPE )
  endif ()

endfunction ()

#
# Stop if in-source build
#
macro ( check_build_tree_path )
  string ( FIND "${CMAKE_BINARY_DIR}" "${PROJECT_BINARY_DIR}/src" RESULT )
  if ( NOT "${RESULT}" STREQUAL "-1")
    message ( FATAL_ERROR  "ERROR: in-source builds are not allowed!")
  endif ()
endmacro( check_build_tree_path )


#
# Print variable (useful for debug)
#
function (print var)
  message (STATUS "   ${var} = ${${var}}")
endfunction ()


#
# Print option
#
function (print_option name value)
  message (STATUS "   ${name} = ${value}")
endfunction ()

#
# Find git info
#
macro ( get_git_info ) # EXTRA ARGS: branch commit

  # Find branch
  execute_process (
    COMMAND git branch
    COMMAND grep \*
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    OUTPUT_VARIABLE out
    ERROR_VARIABLE  err
    )

  if (err)
    message (WARNING "Failing to retrieve incflo Git branch")
  else ()
    string ( REPLACE "*" "" out ${out} )
    string ( STRIP ${out} out)
    message (STATUS "incflo branch: ${out}" )
    if (${ARGC} GREATER 0 ) # branch
      set ( ${ARGV0} ${out} )
    endif ()
  endif ()

  # Find commit
  execute_process (
    COMMAND git rev-parse HEAD
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    OUTPUT_VARIABLE out
    ERROR_VARIABLE  err
    )

  if (err)
    message (WARNING "Failing to retrieve incflo Git commit")
  else ()
    string (STRIP ${out} out)
    message (STATUS "incflo commit: ${out}" )
    if (${ARGC} EQUAL 2 ) # commit
      set ( ${ARGV1} ${out} )
    endif ()
  endif ()

  unset (out)

endmacro ()



#
# Function to prepend path to list items
#
function (prepend list prefix)

  set ( tmp "" )
  foreach (item ${${list}})
    set ( name   ${prefix}/${item} )
    string ( REPLACE "//" "/" name ${name})
    list ( APPEND tmp ${name} )
  endforeach ()

  set ( ${list} ${tmp}  PARENT_SCOPE )

endfunction ()

#
# This sets CMake_<LANG>_FLAGS_<CONFIG> to default values
# if the variable is empty
#
macro ( set_default_config_flags )

  if ( NOT CMAKE_Fortran_FLAGS_DEBUG )
    set (CMAKE_Fortran_FLAGS_DEBUG "-g")
  endif ()

  if ( NOT CMAKE_Fortran_FLAGS_RELEASE )
    set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
  endif ()

  if ( NOT CMAKE_CXX_FLAGS_DEBUG )
    set (CMAKE_CXX_FLAGS_DEBUG "-g")
  endif ()

  if ( NOT CMAKE_CXX_FLAGS_RELEASE )
    set (CMAKE_CXX_FLAGS_RELEASE "-O2 -DNDEBUG")
  endif ()

endmacro ()


#
# Print Configuration Summary
#
function (print_incflo_configuration_summary incflo_libname )

  if (NOT TARGET ${incflo_libname})
    message (AUTHOR_WARNING "Target ${incflo_libname} is not defined.")
    return ()
  endif ()

  if (NOT TARGET AMReX::amrex)
    message (AUTHOR_WARNING "Target ${incflo_libname} is not defined.")
    return ()
  endif ()

  string (TOUPPER "${CMAKE_BUILD_TYPE}" INCFLO_BUILD_TYPE)

  #
  # Get preprocessor flags
  #
  get_target_property ( INCFLO_DEFINES AMReX::amrex INTERFACE_COMPILE_DEFINITIONS )
  replace_genex ( INCFLO_DEFINES INCFLO_Fortran_DEFINES LANGUAGE Fortran )
  replace_genex ( INCFLO_DEFINES INCFLO_CXX_DEFINES LANGUAGE CXX )
  string (REPLACE " " ";-D" INCFLO_Fortran_DEFINES "-D${INCFLO_Fortran_DEFINES}")
  string (REPLACE " " ";-D" INCFLO_CXX_DEFINES "-D${INCFLO_CXX_DEFINES}")

  #
  # Get compiler flags flags
  #
  get_target_property ( INCFLO_FLAGS ${incflo_libname} INTERFACE_COMPILE_OPTIONS )
  replace_genex ( INCFLO_FLAGS INCFLO_Fortran_FLAGS LANGUAGE Fortran )
  replace_genex ( INCFLO_FLAGS INCFLO_CXX_FLAGS LANGUAGE CXX )

  if (NOT "${INCFLO_Fortran_FLAGS}")
    set (INCFLO_Fortran_FLAGS ${CMAKE_Fortran_FLAGS})
  endif ()

  if (NOT "${INCFLO_CXX_FLAGS}")
    set (INCFLO_CXX_FLAGS ${CMAKE_CXX_FLAGS})
  endif ()

  #
  # Get extra includes
  #
  get_target_property( INCFLO_INCLUDES AMReX::amrex INTERFACE_INCLUDE_DIRECTORIES )
  list ( REMOVE_DUPLICATES INCFLO_INCLUDES )

  #
  # Get extra libraries
  #
  get_target_property ( TMP AMReX::amrex INTERFACE_LINK_LIBRARIES )
  replace_genex ( TMP INCFLO_LINK_LINE )
  string (REPLACE ";" " " INCFLO_LINK_LINE "${INCFLO_LINK_LINE}")

  #
  # Config summary
  #
  message( STATUS "INCFLO configuration summary: ")
  message( STATUS "   C++ defines           = ${INCFLO_CXX_DEFINES}")
  message( STATUS "   Fortran defines       = ${INCFLO_Fortran_DEFINES}")
  message( STATUS "   C++ compiler          = ${CMAKE_CXX_COMPILER}")
  message( STATUS "   Fortran compiler      = ${CMAKE_Fortran_COMPILER}")
  message( STATUS "   C++ flags             = ${CMAKE_CXX_FLAGS_${INCFLO_BUILD_TYPE}} ${INCFLO_CXX_FLAGS}")
  message( STATUS "   Fortran flags         = ${CMAKE_Fortran_FLAGS_${INCFLO_BUILD_TYPE}} ${INCFLO_Fortran_FLAGS}")
  message( STATUS "   incflo includes         = ${INCFLO_INCLUDES}")
  message( STATUS "   incflo extra link line  = ${INCFLO_LINK_LINE}")

endfunction()

#
# Replace regex
#
macro (replace_genex input_list output_list )

  cmake_parse_arguments ( ARG "" "LANGUAGE" "" ${ARGN} )

  set (${output_list} "")

  # If input variables is NOT FOUND or empty, just return
  if (${input_list})

   set (tmp_list ${${input_list}})

    # Replace all ; with a place holder (*)
    string ( REPLACE ";" "*" tmp_list "${tmp_list}" )

    # Add tmp_list delimiter only where it suits us
    string ( REPLACE ">*" ">;" tmp_list "${tmp_list}" )
    string ( REPLACE "*$" ";$" tmp_list "${tmp_list}" )
    string ( REPLACE "*/" ";/" tmp_list "${tmp_list}" )
    string ( REPLACE "*" " "   tmp_list "${tmp_list}" )

    #
    # First remove entries related to:
    # 1) a compiler other than the one currently in use
    # 2) a build type other than the current one
    #
    foreach ( item IN ITEMS ${tmp_list} )
      string (REPLACE "$<" "" item ${item} )
      string (REPLACE ">" "" item ${item} )
      string (REPLACE ":" "" item ${item} )

      # Accept build interface generator expressions
      string (REPLACE "BUILD_INTERFACE" "" item ${item})

      # Skip genex for compilers other than the one in use
      string ( FIND ${item} "C_COMPILER_ID" idx1 )
      if ( ${idx1} GREATER -1 )
   	string ( FIND ${item} "${CMAKE_C_COMPILER_ID}" idx2 )
   	if ( ${idx2} GREATER -1 )
   	  string (REPLACE "C_COMPILER_ID${CMAKE_C_COMPILER_ID}" "" item ${item} )
   	else ()
   	  continue ()
   	endif ()
      endif ()

      string (FIND ${item} "CONFIG" idx3 )
      if ( ${idx3} GREATER -1 )
   	string (FIND ${item} "${CMAKE_BUILD_TYPE}" idx4)
   	if ( ${idx4} GREATER -1 )
   	  string (REPLACE "CONFIG${CMAKE_BUILD_TYPE}" "" item ${item} )
   	else ()
   	  continue ()
   	endif ()
      endif ()

      # Extract by Language part
      if ( ARG_LANGUAGE )
	string ( FIND ${item} "COMPILE_LANGUAGE" idx1 )
	if (${idx1} GREATER -1)
	  if (${ARG_LANGUAGE} STREQUAL Fortran )
	    string ( FIND ${item} "Fortran" idx2 )
	    if ( ${idx2} GREATER -1)
	      string (REPLACE "COMPILE_LANGUAGEFortran" "" item ${item} )
	    else()
	      continue ()
	    endif ()
	  elseif (${ARG_LANGUAGE} STREQUAL CXX)
	    string ( FIND ${item} "CXX" idx2 )
	    if ( ${idx2} GREATER -1)
	      string (REPLACE "COMPILE_LANGUAGECXX" "" item ${item} )
	    else()
	      continue ()
	    endif ()
	  endif ()
	endif ()
      endif ()

      # Now item should be ok to be added to final list
      list ( APPEND ${output_list} ${item})

    endforeach ()
  endif ()

  if (${output_list})
    list (REMOVE_DUPLICATES ${output_list} )
  endif ()

endmacro ()


#
# Print list
#
function ( print_list list )

  list ( LENGTH ${list} len )

  if ( ${len} GREATER 0 )
    message ("")
    message ( STATUS " LIST NAME:  ${list}")
    foreach ( item ${${list}})
      message ( STATUS "  ${item}")
    endforeach ()
    message ("")
  endif ()

endfunction ()


