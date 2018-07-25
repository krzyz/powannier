set(CMAKE_CXX_STANDARD 14)

include_directories(../utility)

file(GLOB TEST_SOURCES *.cpp)

#add_subdirectory(${BASE_DIRECTORY}/third-party/Catch2 ${BASE_DIRECTORY}/build/third-party/Catch2 EXCLUDE_FROM_ALL)

#include(ExternalProject)
#
#ExternalProject_Add(
#  Catch
#  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
#  PREFIX ${EXTERNAL_INSTALL_LOCATION}
#)

include(FetchContent)

FetchContent_Declare(
  Catch
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
)

FetchContent_GetProperties(Catch)
if(NOT Catch_POPULATED)
  FetchContent_Populate(cubature)
  add_subdirectory(${Catch_SOURCE_DIR} ${Catch_BINARY_DIR})
endif()


add_executable(tests ${TEST_SOURCES})
target_link_libraries(tests utility Catch2::Catch2)

list(APPEND CMAKE_MODULE_PATH ${BASE_DIRECTORY}/third-party/Catch2/contrib)

include(CTest)
include(Catch)
catch_discover_tests(tests)

install(TARGETS tests DESTINATION ${BASE_DIRECTORY}/bin)