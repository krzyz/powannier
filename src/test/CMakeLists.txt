set(CMAKE_CXX_STANDARD 14)

include_directories(../utility)

file(GLOB TEST_SOURCES *.cpp)

add_subdirectory(${BASE_DIRECTORY}/third-party/Catch2 ${BASE_DIRECTORY}/build/third-party/Catch2 EXCLUDE_FROM_ALL)
add_executable(tests ${TEST_SOURCES})
target_link_libraries(tests utility Catch2::Catch2)

list(APPEND CMAKE_MODULE_PATH ${BASE_DIRECTORY}/third-party/Catch2/contrib)

include(CTest)
include(Catch)
catch_discover_tests(tests)

install(TARGETS tests DESTINATION ${BASE_DIRECTORY}/bin)
