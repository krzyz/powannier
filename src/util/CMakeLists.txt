set(CMAKE_INSTALL_LIBDIR ${THIRD_PARTY_LIBDIR})
set(CMAKE_INSTALL_INCLUDEDIR ${THIRD_PARTY_INCLUDEDIR})

file(GLOB HEADERS *.h)

add_subdirectory(${BASE_DIRECTORY}/third-party/cubature ${BASE_DIRECTORY}/build/third-party/cubature)

add_library(utility STATIC ${HEADERS})
set_target_properties(utility PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(utility cubature)

target_include_directories(utility PUBLIC 
  ${BASE_DIRECTORY}/third-party/cubature
)

message(STATUS ${THIS_INCLUDEDIR})
install(TARGETS utility DESTINATION ${THIS_LIBDIR})
install(FILES ${HEADERS} DESTINATION ${THIS_INCLUDEDIR})
