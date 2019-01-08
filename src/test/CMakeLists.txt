set(CMAKE_CXX_STANDARD 14)

include_directories(../util)
include_directories(../powannier)

file(GLOB TEST_SOURCES *.cpp)

include(FetchContent)

FetchContent_Declare(
  catch
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
)

FetchContent_GetProperties(catch)
if(NOT catch_POPULATED)
  FetchContent_Populate(catch)
  add_subdirectory(${catch_SOURCE_DIR} ${catch_BINARY_DIR})
endif()


add_executable(tests ${TEST_SOURCES})
target_link_libraries(tests utility powannier Catch2::Catch2)

list(APPEND CMAKE_MODULE_PATH ${catch_SOURCE_DIR}/contrib)

include(CTest)
include(Catch)
catch_discover_tests(tests)

# copy test data to build directory
set(TEST_DATA_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testdata/)
file(MAKE_DIRECTORY ${TEST_DATA_DIRECTORY})
file(GLOB TEST_DATA ${CMAKE_CURRENT_SOURCE_DIR}/testdata/*)
file(COPY ${TEST_DATA} DESTINATION ${TEST_DATA_DIRECTORY})

install(TARGETS tests DESTINATION ${BASE_DIRECTORY}/bin)
