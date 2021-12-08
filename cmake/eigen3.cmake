find_package(Eigen3)
if(Eigen3_FOUND)
  include(${EIGEN3_USE_FILE})
else(Eigen3_FOUND)
  include(FetchContent)
  FetchContent_Declare(
    eigen3
    URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.bz2
    SOURCE_SUBDIR none
  )
  include_directories(${eigen3_SOURCE_DIR})
  list(APPEND FetchContents eigen3)
endif(Eigen3_FOUND)
