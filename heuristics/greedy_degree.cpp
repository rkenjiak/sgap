#include "greedy_degree.h"
#include "../algorithms/edit_distance.h"
#include <random>
#include <limits>

namespace seq {

namespace {

// Constroi um caminho a partir de v_start, guiado pelo grau de saida
// do vertice destino. Para quando nao ha vizinhos nao-visitados.
std::vector<VertexId> build_path(
    const SequenceGraph& graph,
    VertexId v_start,
    GreedyDegreeParams::Selection selection,
    std::mt19937& rng)
{
  const std::size_t n = graph.vertex_count();
  std::vector<bool> visited(n, false);
  std::vector<VertexId> path;

  visited[v_start] = true;
  path.push_back(v_start);
  VertexId v = v_start;

  while (true) {
    auto nbrs = graph.neighbors(v);

    // Coletar candidatos nao-visitados
    std::vector<std::size_t> candidates;
    for (std::size_t j = 0; j < nbrs.size(); ++j)
      if (!visited[nbrs[j]])
        candidates.push_back(j);

    if (candidates.empty())
      break;

    std::size_t chosen;

    if (selection == GreedyDegreeParams::Selection::greedy) {
      // Deterministico: maior grau de saida do destino
      chosen = candidates[0];
      std::size_t best_deg = graph.neighbors(nbrs[chosen]).size();
      for (std::size_t k = 1; k < candidates.size(); ++k) {
        std::size_t deg = graph.neighbors(nbrs[candidates[k]]).size();
        if (deg > best_deg) {
          best_deg = deg;
          chosen = candidates[k];
        }
      }
    } else {
      // Roleta proporcional ao grau de saida do destino
      constexpr double EPSILON = 1e-12;
      double total = 0.0;
      for (std::size_t idx : candidates) {
        double deg = static_cast<double>(graph.neighbors(nbrs[idx]).size());
        total += (deg > 0.0) ? deg : EPSILON;
      }

      std::uniform_real_distribution<double> dist(0.0, total);
      double r = dist(rng);
      double acc = 0.0;
      chosen = candidates.back(); // fallback
      for (std::size_t idx : candidates) {
        double deg = static_cast<double>(graph.neighbors(nbrs[idx]).size());
        acc += (deg > 0.0) ? deg : EPSILON;
        if (r <= acc) {
          chosen = idx;
          break;
        }
      }
    }

    v = nbrs[chosen];
    visited[v] = true;
    path.push_back(v);
  }

  return path;
}

} // anonymous namespace

std::vector<VertexId> greedy_degree_path(
    const SequenceGraph& graph,
    GreedyDegreeParams params)
{
  const std::size_t n = graph.vertex_count();
  if (n == 0)
    return {};

  std::mt19937 rng(params.seed);

  if (params.start.has_value()) {
    VertexId v0 = params.start.value();
    if (v0 >= n) return {};
    return build_path(graph, v0, params.selection, rng);
  }

  if (params.selection == GreedyDegreeParams::Selection::roulette) {
    // Roulette sem start definido: usar vertice 0
    return build_path(graph, 0, params.selection, rng);
  }

  // Greedy sem start definido: testar todos, retornar o mais longo
  std::vector<VertexId> best;
  for (std::size_t v = 0; v < n; ++v) {
    auto path = build_path(graph, v, params.selection, rng);
    if (path.size() > best.size())
      best = std::move(path);
  }
  return best;
}

std::optional<GreedyDegreeResult> greedy_degree_best(
    const SequenceGraph& graph,
    const Sequence& query,
    int restarts,
    GreedyDegreeParams base_params)
{
  const std::size_t n = graph.vertex_count();
  if (n == 0 || restarts <= 0)
    return std::nullopt;

  int best_cost = std::numeric_limits<int>::max();
  std::vector<VertexId> best_path;

  for (int i = 0; i < restarts; ++i) {
    GreedyDegreeParams params = base_params;
    params.start = static_cast<VertexId>(i % n);
    params.seed  = base_params.seed + static_cast<unsigned>(i);

    auto path = greedy_degree_path(graph, params);
    if (path.empty())
      continue;

    auto induced = graph.induced_sequence(path);
    int cost = edit_distance(induced, query);

    if (cost < best_cost) {
      best_cost = cost;
      best_path = std::move(path);
    }
  }

  if (best_path.empty())
    return std::nullopt;

  return GreedyDegreeResult{std::move(best_path), best_cost};
}

} // namespace seq
