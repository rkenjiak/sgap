#pragma once

#include "types.h"
#include "sequence.h"
#include "sequence_graph.h"
#include <vector>
#include <span>
#include <optional>
#include <ostream>

namespace seq {

struct AlignmentCosts {
  int substitution = 1;
  int insertion    = 1;
  int deletion     = 1;
};

class AlignmentGraph {
public:
  struct Arc {
    VertexId target;
    int weight;
  };

  struct LayerPosition {
    std::size_t layer;
    std::size_t local_id;
  };

  AlignmentGraph(const SequenceGraph& graph,
                 const Sequence& query,
                 AlignmentCosts costs = {});

  VertexId source_id() const noexcept;
  VertexId sink_id() const noexcept;

  VertexId to_vertex_id(std::size_t layer, std::size_t local_v) const;
  LayerPosition from_vertex_id(VertexId id) const;

  bool is_source(VertexId id) const noexcept;
  bool is_sink(VertexId id) const noexcept;
  bool is_dummy(VertexId id) const;

  std::size_t vertex_count() const noexcept;
  std::size_t arc_count() const noexcept;
  std::size_t layer_count() const noexcept;
  std::size_t graph_vertex_count() const noexcept;

  std::span<const Arc> neighbors(VertexId v) const;

  AlignmentCosts costs() const noexcept;

  std::vector<VertexId> to_graph_walk(std::span<const VertexId> ag_walk) const;

  friend std::ostream& operator<<(std::ostream& os, const AlignmentGraph& g);

private:
  struct Vertex {
    std::vector<Arc> adj;
  };

  std::vector<Vertex> vertices_;
  std::size_t arc_count_ = 0;
  std::size_t layers_;
  std::size_t stride_;
  AlignmentCosts costs_;

  void build(const SequenceGraph& graph, const Sequence& query);
  void add_arc(VertexId u, VertexId v, int weight);
  int match_cost(Symbol query_char, Symbol graph_label) const;
};

struct AlignmentCost { int cost; };
struct AlignmentWalk { int cost; std::vector<VertexId> vertices; };
struct AlignmentPath : AlignmentWalk {};

} // namespace seq
