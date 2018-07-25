set(CMAKE_INSTALL_LIBDIR ${THIRD_PARTY_LIBDIR})
set(CMAKE_INSTALL_INCLUDEDIR ${THIRD_PARTY_INCLUDEDIR})

file(GLOB HEADERS *.h)


#include(ExternalProject)

#ExternalProject_Add(
#  cubature
#  GIT_REPOSITORY https://github.com/stevengj/cubature.git
#  SOURCE_DIR     ${EXTERNAL_INSTALL_LOCATION}
#)

#include(ExternalProject)

#ExternalProject_Add(
    #cubature
    #DOWNLOAD_DIR "${EXTERNAL_INSTALL_LOCATION}/cubature" # The only dir option which is required
    #GIT_REPOSITORY https://github.com/stevengj/cubature.git
    #PREFIX ${EXTERNAL_INSTALL_LOCATION}
    #TIMEOUT 10
    #LOG_DOWNLOAD ON
    #CONFIGURE_COMMAND ""
    #BUILD_COMMAND ""
    #INSTALL_COMMAND ""
#)

#add_subdirectory(${EXTERNAL_INSTALL_LOCATION}/cubature ${BASE_DIRECTORY}/build/third-party/cubature)

#list(APPEND CMAKE_MODULE_PATH ${BASE_DIRECTORY})
include(FetchContent)

FetchContent_Declare(
  cubature
  GIT_REPOSITORY https://github.com/stevengj/cubature.git
)

FetchContent_GetProperties(cubature)
if(NOT cubature_POPULATED)
  FetchContent_Populate(cubature)
  add_subdirectory(${cubature_SOURCE_DIR} ${cubature_BINARY_DIR})
endif()


add_library(utility STATIC ${HEADERS})
add_dependencies(utility cubature)
set_target_properties(utility PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(utility cubature)

target_include_directories(utility PUBLIC 
  ${cubature_SOURCE_DIR}
)

message(STATUS ${THIS_INCLUDEDIR})
install(TARGETS utility DESTINATION ${THIS_LIBDIR})
install(FILES ${HEADERS} DESTINATION ${THIS_INCLUDEDIR})
