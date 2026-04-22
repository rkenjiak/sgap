#include "sequence_graph.h"
#include <algorithm>
#include <cassert>

namespace seq {

SequenceGraph::SequenceGraph(std::size_t vertex_count)
  : vertices_(vertex_count) {}

void
SequenceGraph::add_arc(VertexId u, VertexId v)
{
  assert(u < vertices_.size() && v < vertices_.size());
  vertices_[u].adj.push_back(v);
  ++arc_count_;
}

bool
SequenceGraph::has_arc(VertexId u, VertexId v) const
{
  assert(u < vertices_.size());
  const auto& adj = vertices_[u].adj;
  return std::find(adj.begin(), adj.end(), v) != adj.end();
}

void
SequenceGraph::set_label(VertexId v, Symbol s)
{
  assert(v < vertices_.size());
  vertices_[v].label = s;
}

Symbol
SequenceGraph::label(VertexId v) const
{
  assert(v < vertices_.size());
  return vertices_[v].label;
}

std::span<const VertexId>
SequenceGraph::neighbors(VertexId v) const
{
  assert(v < vertices_.size());
  return vertices_[v].adj;
}

std::size_t
SequenceGraph::vertex_count() const noexcept
{
  return vertices_.size();
}

std::size_t
SequenceGraph::arc_count() const noexcept
{
  return arc_count_;
}

bool
SequenceGraph::is_valid_walk(std::span<const VertexId> walk) const
{
  if (walk.empty()) return true;

  if (walk[0] >= vertices_.size()) return false;

  for (std::size_t i = 1; i < walk.size(); ++i) {
    if (walk[i] >= vertices_.size()) return false;
    if (!has_arc(walk[i - 1], walk[i])) return false;
  }

  return true;
}

bool
SequenceGraph::is_valid_path(std::span<const VertexId> path) const
{
  if (!is_valid_walk(path)) return false;

  std::vector<bool> seen(vertices_.size(), false);
  for (VertexId v : path) {
    if (seen[v]) return false;
    seen[v] = true;
  }

  return true;
}

Sequence
SequenceGraph::induced_sequence(std::span<const VertexId> walk) const
{
  std::string result;
  result.reserve(walk.size());
  for (VertexId v : walk) {
    assert(v < vertices_.size());
    result += vertices_[v].label;
  }
  return Sequence(std::move(result));
}

std::ostream&
operator<<(std::ostream& os, const SequenceGraph& g)
{
  os << "SequenceGraph(" << g.vertices_.size() << " vertices, "
     << g.arc_count_ << " arcs)\n";
  for (std::size_t i = 0; i < g.vertices_.size(); ++i) {
    os << "  " << i << " [" << g.vertices_[i].label << "] ->";
    for (VertexId nb : g.vertices_[i].adj)
      os << " " << nb;
    os << "\n";
  }
  return os;
}

} // namespace seq
