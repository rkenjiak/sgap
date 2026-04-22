#include "beam_search.h"
#include "../algorithms/edit_distance.h"
#include <algorithm>
#include <limits>
#include <vector>
#include <cstdint>

namespace seq {

namespace {

constexpr int NO_PARENT = -1;

// ============================================================================
//  Bitset compacto de tamanho dinamico
// ============================================================================

struct DynBitset {
  std::vector<std::uint64_t> words;

  DynBitset() = default;
  explicit DynBitset(std::size_t n) : words((n + 63) / 64, 0) {}

  bool test(std::size_t i) const {
    return (words[i / 64] >> (i % 64)) & 1;
  }
  void set(std::size_t i) {
    words[i / 64] |= std::uint64_t(1) << (i % 64);
  }
};

struct State {
  std::size_t sg_vertex = 0;
  int cost_g = 0;
  int cost_f = 0;
  int parent_idx = NO_PARENT;
  DynBitset visited;
};

std::vector<VertexId>
reconstruct_path(const std::vector<std::vector<State>>& beams,
                 std::size_t best_layer, int best_idx,
                 std::size_t dummy)
{
  std::vector<VertexId> path;
  int idx = best_idx;
  for (std::size_t l = best_layer + 1; l-- > 0; ) {
    const auto& st = beams[l][idx];
    if (st.sg_vertex != dummy) {
      if (path.empty() || path.back() != static_cast<VertexId>(st.sg_vertex))
        path.push_back(static_cast<VertexId>(st.sg_vertex));
    }
    idx = st.parent_idx;
    if (idx < 0) break;
  }
  std::reverse(path.begin(), path.end());
  return path;
}

} // namespace

// ============================================================================
//  Beam search
// ============================================================================

std::optional<BeamSearchResult>
beam_search_alignment(const SequenceGraph& graph,
                      const Sequence& query,
                      const AlignmentPotentials& potentials,
                      const AlignmentGraph& ag,
                      BeamSearchParams params)
{
  if (query.empty()) return std::nullopt;

  const int k = params.beam_width;
  const std::size_t m = query.size();
  const std::size_t n = graph.vertex_count();
  const std::size_t dummy = n;

  std::vector<std::vector<State>> beams(m);

  // Camada 0
  {
    std::vector<State> candidates;
    candidates.reserve(n + 1);

    for (std::size_t v = 0; v < n; ++v) {
      int g = (query[0] == graph.label(v)) ? 0 : 1;
      int f = g + potentials.backward(ag.to_vertex_id(0, v));
      DynBitset vis(n);
      vis.set(v);
      candidates.push_back({v, g, f, NO_PARENT, std::move(vis)});
    }
    {
      int g = 1;
      int f = g + potentials.backward(ag.to_vertex_id(0, dummy));
      candidates.push_back({dummy, g, f, NO_PARENT, DynBitset(n)});
    }

    if (static_cast<int>(candidates.size()) > k) {
      std::nth_element(candidates.begin(), candidates.begin() + k,
                       candidates.end(),
                       [](const State& a, const State& b) {
                         return a.cost_f < b.cost_f;
                       });
      candidates.resize(k);
    }
    beams[0] = std::move(candidates);
  }

  // Camadas 1..m-1
  for (std::size_t i = 1; i < m; ++i) {
    std::vector<State> candidates;

    for (int si = 0; si < static_cast<int>(beams[i - 1].size()); ++si) {
      const auto& st = beams[i - 1][si];
      std::size_t u = st.sg_vertex;
      int g = st.cost_g;
      const DynBitset& vis = st.visited;

      if (u == dummy) {
        for (std::size_t v = 0; v < n; ++v) {
          if (vis.test(v)) continue;
          int new_g = g + ((query[i] == graph.label(v)) ? 0 : 1);
          int new_f = new_g + potentials.backward(ag.to_vertex_id(i, v));
          DynBitset new_vis = vis;  // copia O(n/64)
          new_vis.set(v);
          candidates.push_back({v, new_g, new_f, si, std::move(new_vis)});
        }
        {
          int new_g = g + 1;
          int new_f = new_g + potentials.backward(ag.to_vertex_id(i, dummy));
          DynBitset new_vis = vis;
          candidates.push_back({dummy, new_g, new_f, si, std::move(new_vis)});
        }
      } else {
        for (VertexId v : graph.neighbors(u)) {
          if (vis.test(v)) continue;
          int new_g = g + ((query[i] == graph.label(v)) ? 0 : 1);
          int new_f = new_g + potentials.backward(ag.to_vertex_id(i, v));
          DynBitset new_vis = vis;
          new_vis.set(v);
          candidates.push_back({v, new_g, new_f, si, std::move(new_vis)});
        }
        {
          int new_g = g + 1;
          int new_f = new_g + potentials.backward(ag.to_vertex_id(i, u));
          DynBitset new_vis = vis;
          candidates.push_back({u, new_g, new_f, si, std::move(new_vis)});
        }
      }
    }

    if (static_cast<int>(candidates.size()) > k) {
      std::nth_element(candidates.begin(), candidates.begin() + k,
                       candidates.end(),
                       [](const State& a, const State& b) {
                         return a.cost_f < b.cost_f;
                       });
      candidates.resize(k);
    }
    beams[i] = std::move(candidates);
  }

  // Selecionar melhor
  const auto& last = beams[m - 1];
  if (last.empty()) return std::nullopt;

  int best_idx = -1;
  int best_g = std::numeric_limits<int>::max();
  for (int i = 0; i < static_cast<int>(last.size()); ++i) {
    if (last[i].cost_g < best_g) {
      best_g = last[i].cost_g;
      best_idx = i;
    }
  }
  if (best_idx < 0) return std::nullopt;

  auto path = reconstruct_path(beams, m - 1, best_idx, dummy);
  if (path.empty()) return std::nullopt;

  auto induced = graph.induced_sequence(path);
  int cost = edit_distance(induced, query);
  return BeamSearchResult{std::move(path), cost};
}

} // namespace seq
