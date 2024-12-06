set(PROGS ising construct1 construct2 construct3 construct4)
foreach(name ${PROGS})
  set(target_name ${name})
  add_executable(${target_name} ${name}.cpp)
  target_compile_options(${target_name} PRIVATE $<$<CXX_COMPILER_ID:GNU>: -Wall -Wextra -Wno-class-memaccess>)
  target_link_libraries(${target_name} lattice Eigen3::Eigen)
endforeach(name)

if(Boost_FOUND)
  set(PROGS construct_xml)
  foreach(name ${PROGS})
    set(target_name ${name})
    add_executable(${target_name} ${name}.cpp)
    target_compile_options(${target_name} PRIVATE $<$<CXX_COMPILER_ID:GNU>: -Wall -Wextra -Wno-class-memaccess>)
    target_link_libraries(${target_name} lattice Boost::boost Eigen3::Eigen)
  endforeach(name)
  file(COPY lattices.xml DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endif()
