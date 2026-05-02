#include "warm_start.h"

#include <algorithm>
#include <limits>

namespace seq {

namespace {

inline int
match_cost_local(Symbol q, Symbol s, AlignmentCosts c)
{
  return (q == s) ? 0 : c.substitution;
}

// Localiza j em ag.neighbors(u) tal que arc.target == v e
// arc.weight == expected_weight. Sem isso (peso) arcos paralelos
// ambiguos como (u,u) [self-loop em G via insertion] vs (u,u)
// [delecao do AG] poderiam casar para o mesmo target.
std::optional<std::size_t>
find_arc(const AlignmentGraph& ag, VertexId u, VertexId v, int expected_weight)
{
  auto neigh = ag.neighbors(u);
  for (std::size_t j = 0; j < neigh.size(); ++j) {
    if (neigh[j].target == v && neigh[j].weight == expected_weight)
      return j;
  }
  return std::nullopt;
}

} // namespace

std::optional<WarmStartArcs>
build_warm_start_arcs(const AlignmentGraph&     ag,
                      const SequenceGraph&      g,
                      const Sequence&           query,
                      std::span<const VertexId> walk_g,
                      AlignmentCosts            costs)
{
  const std::size_t m = query.size();
  const std::size_t k = walk_g.size();

  if (m == 0 || k == 0)         return std::nullopt;
  if (!g.is_valid_walk(walk_g)) return std::nullopt;
  if (ag.layer_count() != m)    return std::nullopt;

  constexpr int INF = std::numeric_limits<int>::max() / 2;

  // dp[i][j] = custo minimo para ter consumido i chars de Q e entrado
  // em j vertices de walk_g.
  //   (0, 0)  : source.
  //   (i, 0)  : (i-1, dummy)        (i deletions desde source).
  //   (i, j)  : (i-1, walk_g[j-1])  (alinhamento parcial).
  //   (0, j>0): inviavel (sem arco source->vertice sem consumir Q).
  //
  // Backtrack: 0 = MatchSub, 1 = AgDel (rota i-1->i, mesmo j),
  //            2 = AgIns   (rota j-1->j, mesmo i).
  std::vector<std::vector<int>>  dp(m + 1, std::vector<int>(k + 1, INF));
  std::vector<std::vector<char>> bt(m + 1, std::vector<char>(k + 1, -1));

  dp[0][0] = 0;
  for (std::size_t i = 1; i <= m; ++i) {
    dp[i][0] = static_cast<int>(i) * costs.deletion;
    bt[i][0] = 1;
  }

  for (std::size_t i = 1; i <= m; ++i) {
    for (std::size_t j = 1; j <= k; ++j) {
      Symbol qc = query[i - 1];
      Symbol sc = g.label(walk_g[j - 1]);
      int    mc = match_cost_local(qc, sc, costs);

      int  best = INF;
      char op   = -1;

      if (dp[i - 1][j - 1] < INF) {
        int c = dp[i - 1][j - 1] + mc;
        if (c < best) { best = c; op = 0; }
      }
      if (dp[i - 1][j] < INF) {
        int c = dp[i - 1][j] + costs.deletion;
        if (c < best) { best = c; op = 1; }
      }
      // AgIns so eh valido se o vertice anterior (j-1) ja foi entrado,
      // ou seja, j-1 >= 1. Para j == 1 nao temos arco (i-1, dummy) ->
      // (i-1, walk_g[0]) na build (dummy nao tem insertion arcs).
      if (j >= 2 && dp[i][j - 1] < INF) {
        int c = dp[i][j - 1] + costs.insertion;
        if (c < best) { best = c; op = 2; }
      }

      dp[i][j] = best;
      bt[i][j] = op;
    }
  }

  if (dp[m][k] >= INF) return std::nullopt;

  // Recupera operacoes em ordem.
  std::vector<char> ops;
  ops.reserve(m + k);
  std::size_t i = m, j = k;
  while (i > 0 || j > 0) {
    char o = bt[i][j];
    if (o < 0) return std::nullopt;
    ops.push_back(o);
    if      (o == 0) { --i; --j; }
    else if (o == 1) { --i; }
    else             { --j; }
  }
  std::reverse(ops.begin(), ops.end());

  // Mapeia (i, j) para o vertice corrente no AG.
  const std::size_t dummy_local = ag.graph_vertex_count();
  auto current = [&](std::size_t ci, std::size_t cj) -> VertexId {
    if (ci == 0 && cj == 0) return ag.source_id();
    if (cj == 0)            return ag.to_vertex_id(ci - 1, dummy_local);
    return ag.to_vertex_id(ci - 1, walk_g[cj - 1]);
  };

  WarmStartArcs result;
  result.arcs.reserve(ops.size() + 1);
  result.cost = 0;

  std::size_t ci = 0, cj = 0;
  for (char o : ops) {
    VertexId u = current(ci, cj);

    int         weight;
    std::size_t ni, nj;
    if (o == 0) {
      ni = ci + 1; nj = cj + 1;
      weight = match_cost_local(query[ci], g.label(walk_g[cj]), costs);
    } else if (o == 1) {
      ni = ci + 1; nj = cj;
      weight = costs.deletion;
    } else {
      ni = ci; nj = cj + 1;
      weight = costs.insertion;
    }

    VertexId v = current(ni, nj);
    auto     j_arc = find_arc(ag, u, v, weight);
    if (!j_arc) return std::nullopt;

    result.arcs.emplace_back(u, *j_arc);
    result.cost += weight;

    ci = ni; cj = nj;
  }

  // Arco final (m-1, walk_g[k-1]) -> sink, peso 0.
  VertexId u_last = current(ci, cj);
  auto     j_sink = find_arc(ag, u_last, ag.sink_id(), 0);
  if (!j_sink) return std::nullopt;
  result.arcs.emplace_back(u_last, *j_sink);

  return result;
}

} // namespace seq
