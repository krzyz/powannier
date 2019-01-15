set(CMAKE_CXX_STANDARD 14)

include_directories(../utility)

find_package(Armadillo REQUIRED)

file(GLOB HEADERS *.h)
file(GLOB BUILDER_HEADERS builders/*.h)
file(GLOB SOURCES *.cpp)
file(GLOB BUILDER_SOURCES builders/*.cpp)

add_library(powannier STATIC ${SOURCES} ${BUILDER_SOURCES})
target_include_directories(powannier PUBLIC ${ARMADILLO_INCLUDE_DIRS}) 
target_include_directories(utility PUBLIC ${ARMADILLO_INCLUDE_DIRS})
target_link_libraries(powannier utility ${ARMADILLO_LIBRARIES})

set_target_properties(powannier
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${OUT_LIB_DIRECTORY}"
  ARCHIVE_OUTPUT_DIRECTORY "${OUT_LIB_DIRECTORY}"
  LIBRARY_OUTPUT_DIRECTORY "${OUT_LIB_DIRECTORY}"
)

if (WIN32)
add_custom_command(TARGET powannier POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        "${OUT_LIB_DIRECTORY}"
        "${OUT_BIN_DIRECTORY}")
endif()

file(MAKE_DIRECTORY ${OUT_INC_DIRECTORY}/powannier)
file(MAKE_DIRECTORY ${OUT_INC_DIRECTORY}/powannier/builders)
file(COPY ${HEADERS} DESTINATION ${OUT_INC_DIRECTORY}/powannier)
file(COPY ${BUILDER_HEADERS} DESTINATION ${OUT_INC_DIRECTORY}/powannier/builders)
