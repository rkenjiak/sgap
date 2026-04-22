#include "alignment_walk.h"
#include <vector>
#include <queue>
#include <deque>
#include <limits>
#include <algorithm>

namespace seq {

namespace {

constexpr int INF = std::numeric_limits<int>::max() / 2;

struct SearchResult {
  std::vector<int> dist;
  std::vector<VertexId> prev;
};

SearchResult
dijkstra_forward(const AlignmentGraph& ag)
{
  const std::size_t V = ag.vertex_count();
  const VertexId src = ag.source_id();

  std::vector<int> dist(V, INF);
  std::vector<VertexId> prev(V, src);

  dist[src] = 0;

  using Entry = std::pair<int, VertexId>;
  std::priority_queue<Entry, std::vector<Entry>, std::greater<>> pq;
  pq.push({0, src});

  while (!pq.empty()) {
    auto [d, v] = pq.top();
    pq.pop();

    if (d > dist[v])
      continue;

    for (const auto& [w, weight] : ag.neighbors(v)) {
      int nd = d + weight;
      if (nd < dist[w]) {
        dist[w] = nd;
        prev[w] = v;
        pq.push({nd, w});
      }
    }
  }

  return {std::move(dist), std::move(prev)};
}

SearchResult
bfs01_forward(const AlignmentGraph& ag)
{
  const std::size_t V = ag.vertex_count();
  const VertexId src = ag.source_id();

  std::vector<int> dist(V, INF);
  std::vector<VertexId> prev(V, src);
  std::vector<bool> settled(V, false);

  dist[src] = 0;

  std::deque<VertexId> dq;
  dq.push_back(src);

  while (!dq.empty()) {
    VertexId v = dq.front();
    dq.pop_front();

    if (settled[v]) continue;
    settled[v] = true;

    for (const auto& [w, weight] : ag.neighbors(v)) {
      int nd = dist[v] + weight;
      if (nd < dist[w]) {
        dist[w] = nd;
        prev[w] = v;
        if (weight == 0)
          dq.push_front(w);
        else
          dq.push_back(w);
      }
    }
  }

  return {std::move(dist), std::move(prev)};
}

SearchResult
shortest_path_forward(const AlignmentGraph& ag, ShortestPathMethod method)
{
  switch (method) {
    case ShortestPathMethod::Dijkstra: return dijkstra_forward(ag);
    case ShortestPathMethod::BFS01:    return bfs01_forward(ag);
  }
  return dijkstra_forward(ag);
}

std::vector<VertexId>
reconstruct_walk(const AlignmentGraph& ag,
                 const std::vector<VertexId>& prev)
{
  const VertexId src = ag.source_id();
  const VertexId snk = ag.sink_id();

  std::vector<VertexId> path;
  for (VertexId v = snk; v != src; v = prev[v])
    path.push_back(v);
  path.push_back(src);

  std::reverse(path.begin(), path.end());
  return ag.to_graph_walk(path);
}

} // namespace

std::optional<AlignmentWalk>
walk_alignment(const AlignmentGraph& ag, ShortestPathMethod method)
{
  auto [dist, prev] = shortest_path_forward(ag, method);

  if (dist[ag.sink_id()] >= INF)
    return std::nullopt;

  auto sg_walk = reconstruct_walk(ag, prev);
  return AlignmentWalk{dist[ag.sink_id()], std::move(sg_walk)};
}

std::optional<AlignmentWalkFull>
walk_alignment_full(const AlignmentGraph& ag, ShortestPathMethod method)
{
  auto [dist, prev] = shortest_path_forward(ag, method);

  if (dist[ag.sink_id()] >= INF)
    return std::nullopt;

  int cost = dist[ag.sink_id()];
  auto sg_walk = reconstruct_walk(ag, prev);
  return AlignmentWalkFull{cost, std::move(sg_walk), std::move(dist)};
}

} // namespace seq
