set(CMAKE_CXX_STANDARD 14)

include_directories(../utility)

find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

file(GLOB HEADERS *.h)
file(GLOB SOURCES *.cpp)

add_library(powannier ${SOURCES})
target_link_libraries(powannier utility ${ARMADILLO_LIBRARIES})

install(TARGETS powannier DESTINATION ${BASE_DIRECTORY}/lib)
install(FILES HEADERS DESTINATION ${BASE_DIRECTORY}/include)
