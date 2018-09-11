#
# Macro to create tags
#
macro (add_tags_targets )

   if (  ( NOT DEFINED AMREX_INSTALL_DIR ) OR
	 ( NOT PROJECT_SOURCE_DIR ) )
      message (AUTHOR_WARNING "Some paths are not defined")
      return ()
   endif ()


   set (AMREX_SRC_DIR ${AMREX_INSTALL_DIR}/../sourcedir/Src/)
   set (INCFLO_SRC_DIR  ${PROJECT_SOURCE_DIR}/src/)

   # Add rule to generate TAGS
   # on macOS ctags-exuberant is just ctags

   find_program(CTAGS_EXU "ctags-exuberant")
   find_program(CTAGS_CMD "ctags")

   if(CTAGS_EXU)
      # Avoid "collisions" with AMReX => focus entirely on incflo sources:
      add_custom_target ( incflo_ctags
         COMMAND ctags-exuberant -R    --fortran-kinds=+i ${INCFLO_SRC_DIR}
         COMMENT "Generating only ctags file for incflo sources exclusively"
         WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/ )
   elseif(CTAGS_CMD)
      # Avoid "collisions" with AMReX => focus entirely on incflo sources:
      add_custom_target ( incflo_ctags
         COMMAND ctags -R    --fortran-kinds=+i ${INCFLO_SRC_DIR}
         COMMENT "Generating only ctags file for incflo sources exclusively"
         WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/ )
   endif()

   unset (INCFLO_SRC_DIR)
   unset (AMREX_SRC_DIR)

endmacro ()
