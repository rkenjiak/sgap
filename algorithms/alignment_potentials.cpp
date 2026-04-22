#include "alignment_potentials.h"
#include "alignment_walk.h"
#include <queue>
#include <deque>
#include <limits>
#include <cassert>
#include <iostream>

namespace seq {

namespace {

constexpr int INF = std::numeric_limits<int>::max() / 2;

struct RevAdj {
  struct Arc { VertexId target; int weight; };
  std::vector<std::vector<Arc>> adj;
};

RevAdj build_reverse_adj(const AlignmentGraph& ag) {
  const std::size_t V = ag.vertex_count();
  RevAdj rev;
  rev.adj.resize(V);
  for (VertexId u = 0; u < V; ++u)
    for (const auto& [v, w] : ag.neighbors(u))
      rev.adj[v].push_back({u, w});
  return rev;
}

std::vector<int>
dijkstra_backward(const AlignmentGraph& ag)
{
  const std::size_t V = ag.vertex_count();
  const VertexId snk = ag.sink_id();
  auto rev = build_reverse_adj(ag);

  std::vector<int> dist(V, INF);
  dist[snk] = 0;

  using Entry = std::pair<int, VertexId>;
  std::priority_queue<Entry, std::vector<Entry>, std::greater<>> pq;
  pq.push({0, snk});

  while (!pq.empty()) {
    auto [d, v] = pq.top();
    pq.pop();

    if (d > dist[v])
      continue;

    for (const auto& [u, w] : rev.adj[v]) {
      int nd = d + w;
      if (nd < dist[u]) {
        dist[u] = nd;
        pq.push({nd, u});
      }
    }
  }

  return dist;
}

std::vector<int>
bfs01_backward(const AlignmentGraph& ag)
{
  const std::size_t V = ag.vertex_count();
  const VertexId snk = ag.sink_id();
  auto rev = build_reverse_adj(ag);

  std::vector<int> dist(V, INF);
  std::vector<bool> settled(V, false);
  dist[snk] = 0;

  std::deque<VertexId> dq;
  dq.push_back(snk);

  while (!dq.empty()) {
    VertexId v = dq.front();
    dq.pop_front();

    if (settled[v]) continue;
    settled[v] = true;

    for (const auto& [u, w] : rev.adj[v]) {
      int nd = dist[v] + w;
      if (nd < dist[u]) {
        dist[u] = nd;
        if (w == 0)
          dq.push_front(u);
        else
          dq.push_back(u);
      }
    }
  }

  return dist;
}

std::vector<int>
backward_shortest_path(const AlignmentGraph& ag, ShortestPathMethod method)
{
  switch (method) {
    case ShortestPathMethod::Dijkstra: return dijkstra_backward(ag);
    case ShortestPathMethod::BFS01:    return bfs01_backward(ag);
  }
  return dijkstra_backward(ag);
}

} // namespace

AlignmentPotentials::AlignmentPotentials(const AlignmentGraph& ag,
                                         ShortestPathMethod method)
  : ag_(&ag)
{
  // Forward: shortest path de source (reutiliza walk_alignment_full)
  auto fwd_result = walk_alignment_full(ag, method);
  assert(fwd_result.has_value() && "AG deve ter caminho source->sink");

  walk_cost_ = fwd_result->cost;
  forward_ = std::move(fwd_result->dist);

  // Backward: shortest path reverso de sink
  backward_ = backward_shortest_path(ag, method);

  // Validacao (debug)
  assert(forward_[ag.source_id()] == 0);
  assert(backward_[ag.sink_id()] == 0);
  assert(forward_[ag.sink_id()] == walk_cost_);
  assert(backward_[ag.source_id()] == walk_cost_);

#ifndef NDEBUG
  for (VertexId v = 0; v < ag.vertex_count(); ++v) {
    if (forward_[v] < INF && backward_[v] < INF) {
      assert(forward_[v] + backward_[v] >= walk_cost_
             && "custo reduzido deve ser >= 0");
    }
  }
#endif
}

int
AlignmentPotentials::forward(VertexId ag_id) const
{
  return forward_[ag_id];
}

int
AlignmentPotentials::backward(VertexId ag_id) const
{
  return backward_[ag_id];
}

int
AlignmentPotentials::backward(std::size_t layer, std::size_t local_v) const
{
  return backward_[ag_->to_vertex_id(layer, local_v)];
}

int
AlignmentPotentials::reduced_cost(VertexId ag_id) const
{
  return forward_[ag_id] + backward_[ag_id] - walk_cost_;
}

int
AlignmentPotentials::walk_cost() const noexcept
{
  return walk_cost_;
}

} // namespace seq
