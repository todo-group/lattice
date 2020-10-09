cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
project(lattice NONE)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/config)
include(core)

include_directories(${PROJECT_SOURCE_DIR})

# Eigen3
find_package(Eigen3 REQUIRED)
include(${EIGEN3_USE_FILE})

# Boost
find_package(Boost)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
endif(Boost_FOUND)

# Google Test
set(GOOGLETEST_VERSION 1)
add_subdirectory(googletest)
include_directories(googletest/include)
enable_testing()

add_subdirectory(test)
add_subdirectory(example)
