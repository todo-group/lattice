# Copyright (C) 2019 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>

set(PROGS ising construct1 construct2 construct3 construct4)
foreach(name ${PROGS})
  add_executable(${name} ${name}.cpp)
endforeach(name)
