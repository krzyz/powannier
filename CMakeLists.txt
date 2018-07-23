cmake_minimum_required(VERSION 3.7)
set(CMAKE_CXX_STANDARD 14)

project(POWannier)
set(VERSION_MAJOR   0   CACHE STRING "Project major version number.")
set(VERSION_MINOR   0   CACHE STRING "Project minor version number.")
set(VERSION_PATCH   1   CACHE STRING "Project patch version number.")

ENABLE_TESTING()
set(BASE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
set(THIS_LIBDIR ${BASE_DIRECTORY}/lib)
set(THIS_INCLUDEDIR ${BASE_DIRECTORY}/include)
set(THIRD_PARTY_LIBDIR ${BASE_DIRECTORY}/third-party/lib)
set(THIRD_PARTY_INCLUDEDIR ${BASE_DIRECTORY}/third-party/include)

add_subdirectory(src/util)
add_subdirectory(src/test)

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
