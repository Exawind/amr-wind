find_program(DOXYSPHINX_EXECUTABLE NAMES doxysphinx
    DOC "Doxysphinx Documentation Builder"
)

if(DOXYSPHINX_EXECUTABLE)
    execute_process(COMMAND ${DOXYSPHINX_EXECUTABLE} --version OUTPUT_VARIABLE DOXYSPHINX_VERSION_OUTPUT)
    if("${DOXYSPHINX_VERSION_OUTPUT}" MATCHES "version ([0-9]+.[0-9]+.[0-9]+)")
      set(DOXYSPHINX_VERSION "${CMAKE_MATCH_1}")
    endif()
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Doxysphinx REQUIRED_VARS DOXYSPHINX_EXECUTABLE
    VERSION_VAR DOXYSPHINX_VERSION
)

mark_as_advanced(DOXYSPHINX_EXECUTABLE)
