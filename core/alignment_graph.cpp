#include "alignment_graph.h"
#include <algorithm>
#include <cassert>

namespace seq {

AlignmentGraph::AlignmentGraph(const SequenceGraph& graph,
                               const Sequence& query,
                               AlignmentCosts costs)
  : layers_(query.size()),
    stride_(graph.vertex_count() + 1),
    costs_(costs)
{
  vertices_.resize(layers_ * stride_ + 2);
  build(graph, query);
}

void
AlignmentGraph::build(const SequenceGraph& graph, const Sequence& query)
{
  const std::size_t n = graph.vertex_count();
  const std::size_t m = layers_;
  const std::size_t dummy = n; // local id of dummy within a layer
  const VertexId s = source_id();
  const VertexId t = sink_id();

  // A. Source edges (layer 0)
  for (std::size_t v = 0; v < n; ++v)
    add_arc(s, to_vertex_id(0, v), match_cost(query[0], graph.label(v)));

  add_arc(s, to_vertex_id(0, dummy), costs_.deletion);

  // B. Inter-layer edges (layers 1..m-1)
  for (std::size_t i = 1; i < m; ++i) {
    // Substitution/Match: for each arc (u,v) in E
    for (std::size_t u = 0; u < n; ++u) {
      for (VertexId v : graph.neighbors(u)) {
        add_arc(to_vertex_id(i - 1, u),
                to_vertex_id(i, v),
                match_cost(query[i], graph.label(v)));
      }
    }

    // Deletion: for each v in V U {d}
    for (std::size_t v = 0; v <= n; ++v)
      add_arc(to_vertex_id(i - 1, v), to_vertex_id(i, v), costs_.deletion);

    // Dummy-to-vertex: (i-1, d) -> (i, v) for each v in V
    for (std::size_t v = 0; v < n; ++v)
      add_arc(to_vertex_id(i - 1, dummy),
              to_vertex_id(i, v),
              match_cost(query[i], graph.label(v)));
  }

  // C. Intra-layer insertion edges (each layer)
  for (std::size_t i = 0; i < m; ++i) {
    for (std::size_t u = 0; u < n; ++u) {
      for (VertexId v : graph.neighbors(u)) {
        add_arc(to_vertex_id(i, u),
                to_vertex_id(i, v),
                costs_.insertion);
      }
    }
  }

  // D. Sink edges (from last layer)
  for (std::size_t v = 0; v <= n; ++v)
    add_arc(to_vertex_id(m - 1, v), t, 0);
}

void
AlignmentGraph::add_arc(VertexId u, VertexId v, int weight)
{
  vertices_[u].adj.push_back({v, weight});
  ++arc_count_;
}

int
AlignmentGraph::match_cost(Symbol query_char, Symbol graph_label) const
{
  return (query_char == graph_label) ? 0 : costs_.substitution;
}

VertexId
AlignmentGraph::source_id() const noexcept
{
  return layers_ * stride_;
}

VertexId
AlignmentGraph::sink_id() const noexcept
{
  return layers_ * stride_ + 1;
}

VertexId
AlignmentGraph::to_vertex_id(std::size_t layer, std::size_t local_v) const
{
  return layer * stride_ + local_v;
}

AlignmentGraph::LayerPosition
AlignmentGraph::from_vertex_id(VertexId id) const
{
  assert(id < layers_ * stride_);
  return {id / stride_, id % stride_};
}

bool
AlignmentGraph::is_source(VertexId id) const noexcept
{
  return id == source_id();
}

bool
AlignmentGraph::is_sink(VertexId id) const noexcept
{
  return id == sink_id();
}

bool
AlignmentGraph::is_dummy(VertexId id) const
{
  if (id >= layers_ * stride_) return false;
  return (id % stride_) == (stride_ - 1);
}

std::size_t
AlignmentGraph::vertex_count() const noexcept
{
  return vertices_.size();
}

std::size_t
AlignmentGraph::arc_count() const noexcept
{
  return arc_count_;
}

std::size_t
AlignmentGraph::layer_count() const noexcept
{
  return layers_;
}

std::size_t
AlignmentGraph::graph_vertex_count() const noexcept
{
  return stride_ - 1;
}

std::span<const AlignmentGraph::Arc>
AlignmentGraph::neighbors(VertexId v) const
{
  assert(v < vertices_.size());
  return vertices_[v].adj;
}

AlignmentCosts
AlignmentGraph::costs() const noexcept
{
  return costs_;
}

std::vector<VertexId>
AlignmentGraph::to_graph_walk(std::span<const VertexId> ag_walk) const
{
  std::vector<VertexId> result;
  for (VertexId v : ag_walk) {
    if (is_source(v) || is_sink(v) || is_dummy(v))
      continue;
    auto [layer, local_id] = from_vertex_id(v);
    if (!result.empty() && result.back() == static_cast<VertexId>(local_id))
      continue;
    result.push_back(static_cast<VertexId>(local_id));
  }
  return result;
}

std::ostream&
operator<<(std::ostream& os, const AlignmentGraph& g)
{
  os << "AlignmentGraph(" << g.vertices_.size() << " vertices, "
     << g.arc_count_ << " arcs, "
     << g.layers_ << " layers, "
     << (g.stride_ - 1) << " graph vertices)\n";
  os << "  source=" << g.source_id() << " sink=" << g.sink_id() << "\n";
  return os;
}

} // namespace seq
