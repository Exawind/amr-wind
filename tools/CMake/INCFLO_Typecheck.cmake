# 
# FUNCTION: add_typecheck_target
# 
# Adds a target to typecheck C++ calls to Fortran routines.
# Works only with GNU compiler so it returns if the compiler id is not GNU.
## Checks for Fortran/C++ headers in the dir of inclusion, and in all dirs below.
#
function ( add_typecheck_target )

   if (NOT (CMAKE_Fortran_COMPILER_ID MATCHES GNU))
      return ()
   endif ()

   if ( (NOT (TARGET AMReX::amrex)) OR (NOT (TARGET ${INCFLO_LIBNAME} )) )
      message (AUTHOR_WARNING "add_typecheck_target() can be called only after targets
AMReX::amrex and INCFLO_LIBNAME are defined ")
      return ()
   endif ()

   #
   # Directory for typecheck 
   #
   set ( TYPECHECK_DIR  ${CMAKE_BINARY_DIR}/TypeCheckTemp )


   # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   # STEP 1: create fortran modules
   # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

   #
   # Define the typecheckobjs library.
   # Basically, we only want to generate the module files
   # associated with fortran otherwise the other two steps
   # of the type check (see below) will fail because
   # it will not be able to include symbols
   # if modules are not there.
   # We use a "static library" rather than an "object library"
   # because this way we can use AMReX transitive dependencies
   # 
   add_library ( typecheckobjs STATIC EXCLUDE_FROM_ALL "" ) 
   target_link_libraries ( typecheckobjs AMReX::amrex )
   set_target_properties ( typecheckobjs
      PROPERTIES
      Fortran_MODULE_DIRECTORY   ${TYPECHECK_DIR} )


   #
   # Get the fortran sources and the fortran headers 
   #
   get_target_property ( INCFLO_ALLSRC ${INCFLO_LIBNAME} SOURCES )

   set ( F90SRC )
   set ( F90HEADERS )
   set ( HEXT ".H" )
   set ( FEXT ".f90;.F90;.f;.F")

   foreach ( item ${INCFLO_ALLSRC} )
      
      # Get the file extension
      get_filename_component ( FILETYPE ${item} EXT )

      # Add to F90 sources
      if ( FILETYPE IN_LIST FEXT )
	 list ( APPEND F90SRC ${item} )
      endif()

      # Add to F90 Headers
      if ( FILETYPE IN_LIST HEXT)
	 string ( REGEX MATCH "_f.H" COND1 ${item})
	 string ( REGEX MATCH "_F.H" COND2 ${item})

	 if ( COND1 OR COND2 )
	    list ( APPEND F90HEADERS ${item})	
	 endif ()

      endif ()
      
   endforeach ()

   # Set sources
   target_sources (typecheckobjs PRIVATE ${F90SRC} ${F90HEADERS} )


   #
   # Find includes needed for typecheck.
   # Must be done manually since we will use them in a custom command
   #
   set (TYPECHECK_INCLUDES)

   get_target_property ( INCFLO_INCLUDE_PATHS ${INCFLO_LIBNAME} INCLUDE_DIRECTORIES )
   get_target_property ( AMREX_INCLUDE_PATHS AMReX::amrex   INTERFACE_INCLUDE_DIRECTORIES )

   string (REPLACE ";" ";-I" TYPECHECK_INCLUDES
      "-I${INCFLO_INCLUDE_PATHS};${AMREX_INCLUDE_PATHS};${TYPECHECK_DIR}")
   list ( REMOVE_DUPLICATES TYPECHECK_INCLUDES )

   # incflo includes are needed by typecheckobjs
   target_include_directories ( typecheckobjs PUBLIC ${INCFLO_INCLUDE_PATHS} )


   # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   # STEP 2: create CPPD files from C++ headers
   # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   if (AMREX_ENABLE_DP)
      set (AMREX_REAL double)
   else ()
      set (AMREX_REAL  float)
   endif ()

   if (AMREX_ENABLE_DP_PARTICLES)
      set (AMREX_PARTICLE_REAL double)
   else ()
      set (AMREX_PARTICLE_REAL  float)
   endif ()


   #
   # Find AMReX defines 
   #
   get_target_property ( AMREX_DEFINES AMReX::amrex INTERFACE_COMPILE_DEFINITIONS )


   # Get rid of genex in defines list. Since they should be
   # only related to fortran stuff, we should be fine
   string ( GENEX_STRIP "${AMREX_DEFINES}" AMREX_DEFINES )
   string ( REPLACE ";" ";-D" TYPECHECK_DEFINES "-D${AMREX_DEFINES}" )

   set (CPPDHEADERS)
   foreach ( file ${F90HEADERS} )
      get_filename_component ( fname ${file} NAME ) # This strips away the path
      set ( CPPD_FILE ${fname}-cppd.h )
      get_filename_component ( fullname ${file} ABSOLUTE ) # This add the absolute path to fname
      add_custom_command ( OUTPUT  ${CPPD_FILE} COMMAND ${CMAKE_C_COMPILER}
	 ARGS ${TYPECHECK_DEFINES} ${TYPECHECK_INCLUDES} -E -P -x c -std=c99 ${fullname} > ${CPPD_FILE}
	 COMMAND sed
	 ARGS -i -e 's/amrex::Real/${AMREX_REAL}/g' ${CPPD_FILE} 
	 COMMAND sed
	 ARGS -i -e 's/amrex_real/${AMREX_REAL}/g' ${CPPD_FILE} 
	 COMMAND sed
	 ARGS -i -e 's/amrex_particle_real/${AMREX_PARTICLE_REAL}/g' ${CPPD_FILE}
	 COMMAND sed
	 ARGS -i -e '/typedef\\s*${AMREX_REAL}/d' ${CPPD_FILE} 
	 COMMAND sed
	 ARGS -i -e 's/\\&/*/g' ${CPPD_FILE} 
	 DEPENDS ${file} typecheckobjs # Leave dependency to typecheck so typecheckdir is created
	 WORKING_DIRECTORY ${TYPECHECK_DIR}  
	 COMMENT "Generating ${CPPD_FILE} " )
      list (APPEND CPPDHEADERS ${CPPD_FILE})
   endforeach ()



   # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   # STEP 3: generate origin files from fortran sources
   # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   set (F90ORIG)
   foreach ( file ${F90SRC} )
      get_filename_component ( fname ${file} NAME ) # This strips away the path
      set ( ORIG_FILE ${fname}.orig )
      get_filename_component ( fullname ${file} ABSOLUTE ) # This add the absolute path to fname
      add_custom_command (
	 OUTPUT   ${ORIG_FILE}
	 COMMAND  ${CMAKE_Fortran_COMPILER}
	 ARGS
	 ${TYPECHECK_DEFINES} ${TYPECHECK_INCLUDES} -fsyntax-only -fdump-fortran-original ${fullname} > ${ORIG_FILE}
	 DEPENDS   ${file} typecheckobjs
	 WORKING_DIRECTORY    ${TYPECHECK_DIR} 
	 COMMENT  "Generating ${ORIG_FILE} " )
      list (APPEND F90ORIG ${ORIG_FILE}) 
   endforeach ()


   # 
   # Add typecheck target
   #
   add_custom_target ( typecheck
      COMMAND python3  ${AMREX_TYPECHECKER}
      --workdir ${TYPECHECK_DIR} --output ${TYPECHECK_DIR}/amrex_typecheck.ou
      DEPENDS ${F90ORIG} ${CPPDHEADERS}
      WORKING_DIRECTORY ${TYPECHECK_DIR}
      )
   
endfunction ()
