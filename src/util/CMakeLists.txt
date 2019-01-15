file(GLOB HEADERS *.h)

include(FetchContent)

FetchContent_Declare(
  cubature
  GIT_REPOSITORY https://github.com/stevengj/cubature.git
)

FetchContent_GetProperties(cubature)
if(NOT cubature_POPULATED)
  FetchContent_Populate(cubature)

  if (WIN32)
    execute_process(
      COMMAND powershell ${SCRIPT_DIRECTORY}/modifyCubatureCMakeLists.ps1
      WORKING_DIRECTORY ${cubature_SOURCE_DIR}
    )
  endif (WIN32)

  add_subdirectory(${cubature_SOURCE_DIR} ${cubature_BINARY_DIR})
endif()

add_library(utility STATIC ${HEADERS})
add_dependencies(utility cubature)
set_target_properties(utility PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(utility cubature)

target_include_directories(utility PUBLIC ${ARMADILLO_INCLUDE_DIRS})
target_include_directories(utility PUBLIC ${cubature_SOURCE_DIR})

file(GLOB CUBATURE_HEADERS ${cubature_SOURCE_DIR}/*.h)

file(MAKE_DIRECTORY ${OUT_INC_DIRECTORY}/util)
file(COPY ${HEADERS} DESTINATION ${OUT_INC_DIRECTORY}/util)
file(COPY ${CUBATURE_HEADERS} DESTINATION ${THIRD_PARTY_INC_DIRECTORY})

set_target_properties(cubature
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${THIRD_PARTY_LIB_DIRECTORY}"
  ARCHIVE_OUTPUT_DIRECTORY "${THIRD_PARTY_LIB_DIRECTORY}"
  LIBRARY_OUTPUT_DIRECTORY "${THIRD_PARTY_LIB_DIRECTORY}"
)