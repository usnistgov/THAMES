cmake_minimum_required(VERSION 2.6)

project (thames)

set (THAMES_VERSION 2.0)
set (THAMES_VERSION_MAJOR 2)
set (THAMES_VERSION_MINOR 0)

# configure a header file to pass some CMake settings to source code
configure_file (
    "${CMAKE_SOURCE_DIR}/src/version.h.in"
    "${CMAKE_SOURCE_DIR}/src/version.h"
    )

# look for doxygen capability
find_package (Doxygen)
if (DOXYGEN_FOUND)
    configure_file ("${CMAKE_SOURCE_DIR}/Doxyfile.in" "${CMAKE_BINARY_DIR}/Doxyfile" @ONLY)
    add_custom_target (doc
        ${DOXYGEN_EXECUTABLE} ${CMAKE_BINARY_DIR}/Doxyfile
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        COMMENT "Generating API documentation with Doxygen" VERBATIM
        )
endif (DOXYGEN_FOUND)

find_package (LibXml2)
if (LIBXML2_FOUND)
    include_directories(${LIBXML2_INCLUDE_DIR})
    message ("Found libxml2 include files at " ${LIBXML2_INCLUDE_DIR})
    message ("libxml2 libraries are " ${LIBXML2_LIBRARIES})
    set (EXTRA_LIBS ${EXTRA_LIBS} ${LIBXML2_LIBRARIES})
    set (XML_INCLUDE_PATH "-I${LIBXML2_INCLUDE_DIR}")
else (LIBXML2_FOUND)
    message (FATAL_ERROR "Library libxml not installed on system.  Please install it first.")
endif (LIBXML2_FOUND)

set (THAMES_INCLUDE_PATH "-I${CMAKE_SOURCE_DIR}/src")
set (CMAKE_CXX_FLAGS "-O2 -DIPMGEMPLUGIN ${THAMES_INCLUDE_PATH} ${XML_INCLUDE_PATH}")

message ("Top level CMAKE_CXX_FLAGS: " ${CMAKE_CXX_FLAGS})

add_subdirectory (${CMAKE_SOURCE_DIR}/src/thameslib)

file (GLOB SOURCES "${CMAKE_SOURCE_DIR}/src/*.cc")

set (EXTRA_LIBS ${EXTRA_LIBS} z)
set (EXTRA_LIBS ${EXTRA_LIBS} m)
#set (EXTRA_LIBS ${EXTRA_LIBS} xml)

add_executable (thames ${SOURCES})

target_link_libraries (thames thameslib ${EXTRA_LIBS})

install (TARGETS thames DESTINATION ${CMAKE_SOURCE_DIR}/bin)
