# Copyright (C) 2019 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>

cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
project(standards NONE)

include(config/core.cmake)

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
add_subdirectory(tool/googletest/googletest)
include_directories(tools/googletest/googletest/include)
enable_testing()

add_subdirectory(test)
add_subdirectory(example)
