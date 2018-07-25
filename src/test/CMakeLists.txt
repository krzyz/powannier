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
  catch
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
)

FetchContent_GetProperties(catch)
if(NOT catch_POPULATED)
  FetchContent_Populate(catch)
  add_subdirectory(${catch_SOURCE_DIR} ${catch_BINARY_DIR})
endif()


add_executable(tests ${TEST_SOURCES})
target_link_libraries(tests utility Catch2::Catch2)

list(APPEND CMAKE_MODULE_PATH ${catch_SOURCE_DIR}/contrib)

include(CTest)
include(Catch)
catch_discover_tests(tests)

install(TARGETS tests DESTINATION ${BASE_DIRECTORY}/bin)
