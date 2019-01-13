cmake_minimum_required(VERSION 3.7)
set(CMAKE_CXX_STANDARD 14)

project(POWannier)
set(VERSION_MAJOR   0   CACHE STRING "Project major version number.")
set(VERSION_MINOR   1   CACHE STRING "Project minor version number.")
set(VERSION_PATCH   0   CACHE STRING "Project patch version number.")

ENABLE_TESTING()

set(SCRIPT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/src/scripts)

# Set and create output directories
set(OUT_BIN_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)
set(OUT_LIB_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)
set(OUT_INC_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include)
file(MAKE_DIRECTORY ${OUT_BIN_DIRECTORY})
file(MAKE_DIRECTORY ${OUT_LIB_DIRECTORY})
file(MAKE_DIRECTORY ${OUT_INC_DIRECTORY})

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(src/util)
add_subdirectory(src/test)
add_subdirectory(src/powannier)

# Building documentation

set(DOC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/docs)

option(BUILD_DOCUMENTATION "Build documentation (requires Doxygen)" OFF)

if (BUILD_DOCUMENTATION)
  set(doxyfile_in ${DOC_DIR}/Doxyfile.in)
  set(doxyfile ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

  configure_file(${doxyfile_in} ${doxyfile} @ONLY)

  add_custom_target(doc
    COMMAND doxygen ${doxyfile}
    WORKING_DIRECTORY ${DOC_DIR}
    COMMENT "Generating API documentation with Doxygen"
    VERBATIM)

  install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html DESTINATION share/doc)
endif()
