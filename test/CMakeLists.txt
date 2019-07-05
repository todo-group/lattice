# Copyright (C) 2019 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>

set(PROGS basis extent unitcell supercell graph)
if(Boost_FOUND)
  list(APPEND PROGS basis_xml unitcell_xml graph_xml)
endif(Boost_FOUND)
foreach(name ${PROGS})
  set(target_name test_${name})
  add_executable(${target_name} ${name}.cpp)
  set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
  target_link_libraries(${target_name} gtest_main)
  add_test(${name} ${name})
endforeach(name)
