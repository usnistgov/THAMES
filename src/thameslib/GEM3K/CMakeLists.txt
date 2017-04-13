cmake_minimum_required(VERSION 2.6)

project (gem3k_project)

file (GLOB SOURCES "${CMAKE_SOURCE_DIR}/src/thameslib/GEM3K/*.cpp")

add_library (gem3k STATIC ${SOURCES})

install (TARGETS gem3k DESTINATION ${CMAKE_SOURCE_DIR}/lib)
