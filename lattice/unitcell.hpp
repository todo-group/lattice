#ifndef LATTICE_UNITCELL_HPP
#define LATTICE_UNITCELL_HPP

#include "basis.hpp"
#include "types.hpp"

namespace lattice {

class unitcell {
public:
  struct vertex_t {
    vertex_t() {}
    vertex_t(const coordinate_t& pos, int tp) :
      coordinate(pos), type(tp), neighbors(0), neighbor_edges(0) {}
    coordinate_t coordinate;
    int type;
    std::vector<std::size_t> neighbors;
    std::vector<std::size_t> neighbor_edges;
  };
  
  struct edge_t {
    edge_t() {}
    edge_t(std::size_t s, std::size_t t, offset_t os, int tp) :
      source(s), target(t), target_offset(os), type(tp) {}
    std::size_t source, target;
    offset_t target_offset;
    int type;
  };

  unitcell() {}
  unitcell(const basis& bs) : basis_(bs) {}
  unitcell(const unitcell& cell, const extent_t& extent);
  unitcell(const unitcell& cell, const span_t& span);

  std::size_t add_vertex(const coordinate_t& pos, int tp) {
    if (pos.size() != dimension())
      throw std::invalid_argument("vertex coordinate dimension mismatch");
    for (std::size_t i = 0; i < dimension(); ++i) {
      if (pos[i] < 0 || pos[i] >= 1.0)
        throw std::invalid_argument("vertex coordinate out of range");
    }
    std::size_t s = vertices_.size();
    vertices_.push_back(vertex_t(pos, tp));
    return s;
  }
    
  std::size_t add_edge(std::size_t s, std::size_t t, const offset_t& os, int tp) {
    if (s >= num_vertices() || t >= num_vertices())
      throw std::invalid_argument("vertex index out of range");
    if (os.size() != dimension())
      throw std::invalid_argument("unitcell offset dimension mismatch");
    std::size_t b = edges_.size();
    edges_.push_back(edge_t(s, t, os, tp));
    vertices_[s].neighbor_edges.push_back(b);
    vertices_[t].neighbor_edges.push_back(b);
    vertices_[s].neighbors.push_back(t);
    vertices_[t].neighbors.push_back(s);
    return b;
  }

  std::size_t dimension() const { return basis_.dimension(); }
  std::size_t num_vertices() const { return vertices_.size(); }
  std::size_t num_edges() const { return edges_.size(); }
  const vertex_t& vertex(std::size_t s) const { return vertices_[s]; }
  const edge_t& edge(std::size_t b) const { return edges_[b]; }
  double volume() const { return basis_.volume(); }
  
private:
  basis basis_;
  std::vector<vertex_t> vertices_;
  std::vector<edge_t> edges_;
};

std::size_t dimension(const unitcell& cell) {
  return cell.dimension();
}
  
} // end namespace lattice

#endif
