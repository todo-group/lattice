cmake_minimum_required(VERSION 3.14)
include(cmake/prefix.cmake)
project(lattice CXX)
include(cmake/postfix.cmake)

add_library(lattice INTERFACE)
target_include_directories(lattice INTERFACE ${PROJECT_SOURCE_DIR})

include(eigen3)
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  include(googletest)
endif()
FetchContent_MakeAvailable(${FetchContents})

find_package(Boost)

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  add_subdirectory(test)
  add_subdirectory(example)
endif()
