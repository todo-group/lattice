cmake_minimum_required(VERSION 3.14)
include(cmake/prefix.cmake)
project(lattice CXX)
include(cmake/postfix.cmake)

option(LATTICE_INSTALL "Install lattice." ON)
option(LATTICE_BUILD_TESTS "Build all lattice's tests." ON)
option(LATTICE_BUILD_EXAMPLES "Build all lattice's examples." ON)

include_directories(${lattice_SOURCE_DIR})

include(googletest)
FetchContent_MakeAvailable(googletest)

# Eigen3  
find_package(Eigen3 REQUIRED)
include(${EIGEN3_USE_FILE})

# Boost
find_package(Boost)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
endif(Boost_FOUND)

if(LATTICE_BUILD_TESTS)
  add_subdirectory(test)
endif(LATTICE_BUILD_TESTS)
if(LATTICE_BUILD_EXAMPLES)
  add_subdirectory(example)
endif(LATTICE_BUILD_EXAMPLES)
