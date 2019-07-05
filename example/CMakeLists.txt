# Copyright (C) 2019 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>

set(PROGS ising)
foreach(name ${PROGS})
  add_executable(${name} ${name}.cpp)
endforeach(name)
