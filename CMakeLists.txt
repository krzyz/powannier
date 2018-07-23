# General options

cmake_minimum_required(VERSION 3.7)
set(CMAKE_CXX_STANDARD 14)

#set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)
#set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)
#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)
#
project(POWannier)
set(VERSION_MAJOR   0   CACHE STRING "Project major version number.")
set(VERSION_MINOR   0   CACHE STRING "Project minor version number.")
set(VERSION_PATCH   1   CACHE STRING "Project patch version number.")
#mark_as_advanced(VERSION_MAJOR VERSION_MINOR VERSION_PATCH)
#
#set(CMAKE_CXX_STANDARD 14)
#
#find_package(Armadillo REQUIRED)
#
#include_directories(${ARMADILLO_INCLUDE_DIRS})
#
#add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/third-party/cubature)
#add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/third-party/Catch2)
#
## Building the library
#
#file(GLOB UTIL_SOURCES 
#    src/util/*.h
#    src/util/*.cpp
#)
#
## Building tests
#
#file(GLOB TEST_SOURCES
#    src/test/*.cpp
#)
#
#add_executable(tests ${TEST_SOURCES} ${UTIL_SOURCES})
#
#target_include_directories(tests PUBLIC 
#  ${CMAKE_CURRENT_SOURCE_DIR}/third-party/cubature
#)
#
#target_link_libraries(tests Catch2::Catch2 cubature)
#
#list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/third-party/Catch2/contrib)
#
#include(CTest)
#include(Catch)
#catch_discover_tests(tests)

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
