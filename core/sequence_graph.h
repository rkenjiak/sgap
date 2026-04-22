#pragma once

#include "types.h"
#include "sequence.h"
#include <vector>
#include <span>
#include <ostream>

namespace seq {

class SequenceGraph {
public:
  explicit SequenceGraph(std::size_t vertex_count);

  void add_arc(VertexId u, VertexId v);  

  void   set_label(VertexId v, Symbol s);
  Symbol label(VertexId v) const;

  std::span<const VertexId> neighbors(VertexId v) const;

  std::size_t vertex_count() const noexcept;
  std::size_t arc_count() const noexcept;

  Sequence induced_sequence(std::span<const VertexId> walk) const;

  bool is_valid_walk(std::span<const VertexId> walk) const;
  bool is_valid_path(std::span<const VertexId> path) const;

  friend std::ostream& operator<<(std::ostream& os, const SequenceGraph& g);

private:
  struct Vertex
  {
    Symbol label = '\0';
    std::vector<VertexId> adj;
  };

  bool has_arc(VertexId u, VertexId v) const;

  std::vector<Vertex> vertices_;
  std::size_t arc_count_ = 0;
};

} // namespace seq
