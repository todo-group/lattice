cmake_minimum_required(VERSION 3.14)
include(cmake/prefix.cmake)
project(lattice CXX)
include(cmake/postfix.cmake)

option(LATTICE_BUILD_RUST_FFI "Build Rust FFI shared library target via cargo" OFF)

add_library(lattice INTERFACE)
target_include_directories(lattice INTERFACE ${PROJECT_SOURCE_DIR})

set(LATTICE_RUST_FFI_RELEASE_LIB
  ${PROJECT_SOURCE_DIR}/rust/target/release/${CMAKE_SHARED_LIBRARY_PREFIX}lattice_ffi${CMAKE_SHARED_LIBRARY_SUFFIX}
)
set(LATTICE_RUST_FFI_DEBUG_LIB
  ${PROJECT_SOURCE_DIR}/rust/target/debug/${CMAKE_SHARED_LIBRARY_PREFIX}lattice_ffi${CMAKE_SHARED_LIBRARY_SUFFIX}
)

if(LATTICE_BUILD_RUST_FFI)
  find_program(LATTICE_CARGO_EXECUTABLE cargo)
  if(NOT LATTICE_CARGO_EXECUTABLE)
    message(FATAL_ERROR "LATTICE_BUILD_RUST_FFI=ON requires cargo in PATH")
  endif()
  add_custom_target(lattice_rust_ffi_build
    COMMAND ${LATTICE_CARGO_EXECUTABLE} build --manifest-path ${PROJECT_SOURCE_DIR}/rust/Cargo.toml -p lattice-ffi
    COMMAND ${LATTICE_CARGO_EXECUTABLE} build --release --manifest-path ${PROJECT_SOURCE_DIR}/rust/Cargo.toml -p lattice-ffi
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/rust
    COMMENT "Building Rust FFI library for lattice XML compatibility"
  )
endif()

if(EXISTS ${LATTICE_RUST_FFI_RELEASE_LIB})
  set(LATTICE_RUST_FFI_LOCATION ${LATTICE_RUST_FFI_RELEASE_LIB})
elseif(EXISTS ${LATTICE_RUST_FFI_DEBUG_LIB})
  set(LATTICE_RUST_FFI_LOCATION ${LATTICE_RUST_FFI_DEBUG_LIB})
elseif(LATTICE_BUILD_RUST_FFI)
  # Link to the expected location; build target above creates this artifact.
  set(LATTICE_RUST_FFI_LOCATION ${LATTICE_RUST_FFI_DEBUG_LIB})
else()
  message(FATAL_ERROR
    "Rust XML backend requires lattice-ffi. Either set LATTICE_BUILD_RUST_FFI=ON or build rust/lattice-ffi first."
  )
endif()

add_library(lattice_rust_ffi SHARED IMPORTED GLOBAL)
set_target_properties(lattice_rust_ffi PROPERTIES IMPORTED_LOCATION ${LATTICE_RUST_FFI_LOCATION})
if(TARGET lattice_rust_ffi_build)
  add_dependencies(lattice_rust_ffi lattice_rust_ffi_build)
endif()
target_link_libraries(lattice INTERFACE lattice_rust_ffi)

include(eigen3)
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  include(googletest)
endif()
FetchContent_MakeAvailable(${FetchContents})

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  add_subdirectory(test)
  add_subdirectory(example)
endif()
