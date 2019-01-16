cmake_minimum_required(VERSION 3.7)
set(CMAKE_CXX_STANDARD 14)

project(honeycomb)

set(POWANNIER_DIRECTORY ${CMAKE_SOURCE_DIR}/../)
set(THIRD_PARTY_DIRECTORY ${CMAKE_SOURCE_DIR}/../third-party/)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

find_package(Armadillo REQUIRED)

add_library(powannier STATIC IMPORTED)
set_target_properties(powannier PROPERTIES
    IMPORTED_LOCATION "${POWANNIER_DIRECTORY}/lib/libpowannier.a"
    INTERFACE_INCLUDE_DIRECTORIES "${POWANNIER_DIRECTORY}/include"
)

add_library(cubature SHARED IMPORTED)
set_target_properties(cubature PROPERTIES
  IMPORTED_LOCATION "${THIRD_PARTY_DIRECTORY}/lib/libcubature.so"
  INTERFACE_INCLUDE_DIRECTORIES "${THIRD_PARTY_DIRECTORY}/include"
)


add_executable(honeycomb honeycomb.cpp)
target_link_libraries(honeycomb cubature powannier armadillo)

#target_include_directories(honeycomb PUBLIC  ${ARMADILLO_INCLUDE_DIRS})
#target_link_libraries(honeycomb cubature powannier ${ARMADILLO_LIBRARIES})

set_target_properties(honeycomb
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}"
)
