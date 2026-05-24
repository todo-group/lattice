cmake_minimum_required(VERSION 3.14)
include(cmake/prefix.cmake)
project(lattice CXX)
include(cmake/postfix.cmake)
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
include(FetchContent)

set(LATTICE_PKG_VERSION 0.1.0)
set(LATTICE_INSTALL_CMAKEDIR ${CMAKE_INSTALL_LIBDIR}/cmake/lattice)
option(LATTICE_BUILD_TESTS "Build C++ tests" OFF)

add_library(lattice INTERFACE)
target_include_directories(lattice INTERFACE ${PROJECT_SOURCE_DIR})

set(LATTICE_RUST_FFI_RELEASE_LIB
  ${CMAKE_BINARY_DIR}/rust-target/release/${CMAKE_SHARED_LIBRARY_PREFIX}lattice_ffi${CMAKE_SHARED_LIBRARY_SUFFIX}
)
set(LATTICE_RUST_FFI_DEBUG_LIB
  ${CMAKE_BINARY_DIR}/rust-target/debug/${CMAKE_SHARED_LIBRARY_PREFIX}lattice_ffi${CMAKE_SHARED_LIBRARY_SUFFIX}
)
set(LATTICE_RUST_CARGO_TARGET_DIR ${CMAKE_BINARY_DIR}/rust-target)
set(LATTICE_RUST_WORKSPACE_STAGING_DIR ${CMAKE_BINARY_DIR}/rust-workspace)

find_program(LATTICE_CARGO_EXECUTABLE cargo)

if(LATTICE_CARGO_EXECUTABLE)
  add_custom_target(lattice_rust_ffi_build
    COMMAND ${CMAKE_COMMAND}
      -DLATTICE_CARGO_EXECUTABLE=${LATTICE_CARGO_EXECUTABLE}
      -DLATTICE_RUST_FFI_DEBUG_LIB=${LATTICE_RUST_FFI_DEBUG_LIB}
      -DLATTICE_RUST_FFI_RELEASE_LIB=${LATTICE_RUST_FFI_RELEASE_LIB}
      -DLATTICE_RUST_WORK_DIR=${PROJECT_SOURCE_DIR}/rust
      -DLATTICE_RUST_MANIFEST=${PROJECT_SOURCE_DIR}/rust/Cargo.toml
      -DLATTICE_RUST_CARGO_TARGET_DIR=${LATTICE_RUST_CARGO_TARGET_DIR}
      -DLATTICE_RUST_WORKSPACE_STAGING_DIR=${LATTICE_RUST_WORKSPACE_STAGING_DIR}
      -P ${PROJECT_SOURCE_DIR}/cmake/build_rust_ffi_if_needed.cmake
    COMMENT "Ensuring Rust FFI library exists"
  )
elseif(NOT EXISTS ${LATTICE_RUST_FFI_RELEASE_LIB} AND NOT EXISTS ${LATTICE_RUST_FFI_DEBUG_LIB})
  message(FATAL_ERROR
    "Rust XML backend requires rust/lattice-ffi and cargo is not available in PATH. Install cargo or prebuild rust/lattice-ffi."
  )
endif()

add_library(lattice_rust_ffi SHARED IMPORTED GLOBAL)
set_target_properties(lattice_rust_ffi PROPERTIES
  IMPORTED_LOCATION ${LATTICE_RUST_FFI_DEBUG_LIB}
  IMPORTED_LOCATION_DEBUG ${LATTICE_RUST_FFI_DEBUG_LIB}
  IMPORTED_LOCATION_RELEASE ${LATTICE_RUST_FFI_RELEASE_LIB}
  IMPORTED_LOCATION_RELWITHDEBINFO ${LATTICE_RUST_FFI_RELEASE_LIB}
  IMPORTED_LOCATION_MINSIZEREL ${LATTICE_RUST_FFI_RELEASE_LIB}
)
if(TARGET lattice_rust_ffi_build)
  add_dependencies(lattice_rust_ffi lattice_rust_ffi_build)
endif()
target_link_libraries(lattice INTERFACE lattice_rust_ffi)

include(eigen3)
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME AND LATTICE_BUILD_TESTS)
  include(googletest)
endif()
FetchContent_MakeAvailable(${FetchContents})

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  if(LATTICE_BUILD_TESTS)
    add_subdirectory(test)
  endif()
  add_subdirectory(example)
endif()

if(LATTICE_CARGO_EXECUTABLE)
  configure_file(
    ${PROJECT_SOURCE_DIR}/cmake/install_ensure_rust_ffi.cmake.in
    ${CMAKE_BINARY_DIR}/install_ensure_rust_ffi.cmake
    @ONLY
  )
  install(SCRIPT ${CMAKE_BINARY_DIR}/install_ensure_rust_ffi.cmake)
endif()

configure_package_config_file(
  ${PROJECT_SOURCE_DIR}/cmake/latticeConfig.cmake.in
  ${CMAKE_BINARY_DIR}/latticeConfig.cmake
  INSTALL_DESTINATION ${LATTICE_INSTALL_CMAKEDIR}
  PATH_VARS CMAKE_INSTALL_INCLUDEDIR CMAKE_INSTALL_LIBDIR
)

write_basic_package_version_file(
  ${CMAKE_BINARY_DIR}/latticeConfigVersion.cmake
  VERSION ${LATTICE_PKG_VERSION}
  COMPATIBILITY SameMajorVersion
)

configure_file(
  ${PROJECT_SOURCE_DIR}/cmake/lattice.pc.in
  ${CMAKE_BINARY_DIR}/lattice.pc
  @ONLY
)

install(DIRECTORY ${PROJECT_SOURCE_DIR}/lattice/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/lattice
  FILES_MATCHING
  PATTERN "*.hpp"
  PATTERN "*.h"
)

install(FILES
  $<$<CONFIG:Debug>:${LATTICE_RUST_FFI_DEBUG_LIB}>
  $<$<NOT:$<CONFIG:Debug>>:${LATTICE_RUST_FFI_RELEASE_LIB}>
  DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

install(FILES ${PROJECT_SOURCE_DIR}/LICENSE-2.0.txt
  DESTINATION ${CMAKE_INSTALL_DATADIR}/lattice
)

install(FILES
  ${CMAKE_BINARY_DIR}/latticeConfig.cmake
  ${CMAKE_BINARY_DIR}/latticeConfigVersion.cmake
  DESTINATION ${LATTICE_INSTALL_CMAKEDIR}
)

install(FILES
  ${CMAKE_BINARY_DIR}/lattice.pc
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig
)
