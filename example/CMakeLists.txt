# Copyright (C) 2019 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>

set(PROGS ising construct1 construct2 construct3 construct4)
if(Boost_FOUND)
  list(APPEND PROGS construct_xml)
  file(COPY lattices.xml DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endif(Boost_FOUND)
foreach(name ${PROGS})
  add_executable(${name} ${name}.cpp)
endforeach(name)
