cmake_minimum_required(VERSION 3.7)
set(CMAKE_CXX_STANDARD 14)

project(POWannier)
set(VERSION_MAJOR   0   CACHE STRING "Project major version number.")
set(VERSION_MINOR   2   CACHE STRING "Project minor version number.")
set(VERSION_PATCH   1   CACHE STRING "Project patch version number.")
mark_as_advanced(VERSION_MAJOR VERSION_MINOR VERSION_PATCH)

if (CMAKE_COMPILER_IS_GNUCC)
    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall")
endif()
if (MSVC)
    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} /W4")
endif()

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

ENABLE_TESTING()

set(SCRIPT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/src/scripts)

# Set and create output directories
set(OUT_BIN_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)
set(OUT_LIB_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)
set(OUT_INC_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(THIRD_PARTY_LIB_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/third_party/lib)
set(THIRD_PARTY_INC_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/third_party/include)

file(MAKE_DIRECTORY ${OUT_BIN_DIRECTORY})
file(MAKE_DIRECTORY ${OUT_LIB_DIRECTORY})
file(MAKE_DIRECTORY ${OUT_INC_DIRECTORY})
file(MAKE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/third_party)
file(MAKE_DIRECTORY ${THIRD_PARTY_LIB_DIRECTORY})
file(MAKE_DIRECTORY ${THIRD_PARTY_INC_DIRECTORY})

file(COPY src/powannier.h DESTINATION ${OUT_INC_DIRECTORY})

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
