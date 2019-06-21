#ifndef LATTICE_TYPES_HPP
#define LATTICE_TYPES_HPP

#include <string>
#include <Eigen/Dense>

namespace lattice {

enum boundary_t { open, periodic };

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> basis_t;

// set of spanning vectors
typedef Eigen::Matrix<long, Eigen::Dynamic, Eigen::Dynamic> span_t;

// extent vector = diagonal spanning matrix
typedef Eigen::Matrix<long, Eigen::Dynamic, 1> extent_t;

// coordinate vector
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> coordinate_t;

// unitcell offset
typedef Eigen::Matrix<long, Eigen::Dynamic, 1> offset_t;

} // end namespace lattice

#endif
