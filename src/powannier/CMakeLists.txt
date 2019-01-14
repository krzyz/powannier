set(CMAKE_CXX_STANDARD 14)

include_directories(../utility)

find_package(Armadillo REQUIRED)

file(GLOB HEADERS *.h)
file(GLOB SOURCES *.cpp)

add_library(powannier SHARED ${SOURCES})
target_include_directories(powannier PUBLIC ${ARMADILLO_INCLUDE_DIRS}) 
target_include_directories(utility PUBLIC ${ARMADILLO_INCLUDE_DIRS})
target_link_libraries(powannier utility ${ARMADILLO_LIBRARIES})

set_target_properties(powannier
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${OUT_LIB_DIRECTORY}"
  ARCHIVE_OUTPUT_DIRECTORY "${OUT_LIB_DIRECTORY}"
  LIBRARY_OUTPUT_DIRECTORY "${OUT_LIB_DIRECTORY}"
)

file(COPY ${HEADERS} DESTINATION ${OUT_INC_DIRECTORY})
