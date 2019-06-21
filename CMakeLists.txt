cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
project(standards NONE)

include(config/core.cmake)

include_directories(${PROJECT_SOURCE_DIR})

# Eigen3
find_package(Eigen3)
include(${EIGEN3_USE_FILE})

# Google Test
set(GOOGLETEST_VERSION 1)
add_subdirectory(tool/googletest)
include_directories(tools/googletest/include)
enable_testing()

add_subdirectory(test)
