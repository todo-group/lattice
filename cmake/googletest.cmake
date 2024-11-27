option(BUILD_GMOCK OFF)
option(INSTALL_GTEST OFF)

include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        v1.15.2
)
list(APPEND FetchContents googletest)
