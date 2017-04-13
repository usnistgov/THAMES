cmake_minimum_required(VERSION 2.6)

project (ThamesLibrary_project)

add_subdirectory(GEM3K)

set (EXTRA_LIBS ${EXTRA_LIBS} gem3k)

file (GLOB SOURCES "${CMAKE_SOURCE_DIR}/src/thameslib/*.cc")

#find_package (LibXml2)
#if (LIBXML2_FOUND)
#    include_directories(${LIBXML2_INCLUDE_DIR})
#    message ("Found libxml2 include files at " ${LIBXML2_INCLUDE_DIR})
#    message ("libxml2 libraries are " ${LIBXML2_LIBRARIES})
#    set (EXTRA_LIBS ${EXTRA_LIBS} ${LIBXML2_LIBRARIES})
#    set (XML_INCLUDE_PATH "-I${LIBXML2_INCLUDE_DIR}")
#else (LIBXML2_FOUND)
#    message (FATAL_ERROR "Library libxml not installed on system.  Please install it first.")
#endif (LIBXML2_FOUND)

#set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${XML_INCLUDE_PATH}")

message ("THAMES lib directory, CMAKE_CXX_FLAGS: " ${CMAKE_CXX_FLAGS})

add_library (thameslib STATIC ${SOURCES})

target_link_libraries(thameslib ${EXTRA_LIBS})

install (TARGETS thameslib DESTINATION ${CMAKE_SOURCE_DIR}/lib)
