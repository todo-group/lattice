/*
   Copyright (C) 2019 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/

#include <iostream>
#include <random>
#include <vector>
#include "lattice/graph.hpp"
#include "observable.hpp"

using namespace lattice;

int main() {
  std::cout << "Metropolis Algorithm for Classical Ferromagnetic Ising Model\n";
  
  // inverse temperature
  double beta = 1 / 2.269;

  // square lattice
  std::size_t dim = 2;
  std::size_t length = 32;
  lattice::graph lat = lattice::graph::simple(dim, length);

  // random number generators
  int seed = 12345;
  std::mt19937 eng(seed);
  std::uniform_real_distribution<> uniform_01;

  // spin configuration
  std::vector<double> spins(lat.num_sites(), 1); // all up

  // measurements
  observable energy, mag2;
  
  std::size_t sweeps = 65536;
  std::size_t therm = sweeps / 8;
  for (std::size_t mcs = 0; mcs < therm + sweeps; ++mcs) {
    // metropolis update
    for (std::size_t s = 0; s < lat.num_sites(); ++s) {
      double diff = 0.0;
      for (std::size_t k = 0; k < lat.num_neighbors(s); ++k)
        diff += 2 * spins[s] * spins[lat.neighbor(s, k)];
      if (uniform_01(eng) < std::exp(- beta * diff)) spins[s] = - spins[s];
    }
    if (mcs > therm) {
      double ene = 0;
      for (std::size_t b = 0; b < lat.num_bonds(); ++b)
        ene -= spins[lat.source(b)] * spins[lat.target(b)];
      energy << ene / lat.num_sites();
      double mag = std::accumulate(spins.begin(), spins.end(), 0.0) / lat.num_sites();
      mag2 << mag * mag;
    }
  }
  
  std::cout << "dimension  = " << dim << std::endl
            << "length = " << length << std::endl
            << "beta = " << beta << std::endl
            << "energy density = " << energy.mean() << " +/- " << energy.error() << std::endl
            << "magnetization density^2 = " << mag2.mean() << " +/- " << mag2.error() << std::endl;
}
