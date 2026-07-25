#include "lagrangian.h"
#include "../io/gurobi_sol.h"  // formato .sol: I/O de texto puro, nao linka Gurobi

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>
#include <queue>
#include <vector>

namespace seq {

namespace {

constexpr double INF = std::numeric_limits<double>::infinity();

// Replica a regra de pli.cpp (linhas 58-69): devolve o local_id do
// vertice do grafo de sequencia que o arco u->target "entra" (e portanto
// conta para in_v), ou -1 quando o arco nao consome nenhum vertice
// (chegada no sink, no dummy, ou arco de delecao que permanece no mesmo
// vertice local).
int
penalized_vertex(const AlignmentGraph& ag, VertexId u, VertexId target)
{
  const std::size_t dummy_local = ag.graph_vertex_count(); // == G

  if (ag.is_sink(target)) return -1;

  auto tgt = ag.from_vertex_id(target);
  if (tgt.local_id == dummy_local) return -1;

  std::size_t src_local;
  if (ag.is_source(u))
    src_local = static_cast<std::size_t>(-1);
  else
    src_local = ag.from_vertex_id(u).local_id;

  if (src_local != tgt.local_id)
    return static_cast<int>(tgt.local_id);
  return -1;
}

struct SubResult {
  std::vector<VertexId> ag_walk;       // source..sink (ids do AG)
  double                modified_cost; // Sum (w + lambda_p) ao longo do caminho
};

// Subproblema Lagrangeano: caminho minimo S->T no grafo de alinhamento
// com custo de arco = w + lambda_p (p = penalized_vertex). Como lambda>=0
// e w>=0, os custos sao nao-negativos => Dijkstra. Reconstroi o caminho
// no proprio AG (precisamos dos ids do AG para SolverResult.ag_walk e
// para recalcular o subgradiente).
std::optional<SubResult>
solve_subproblem(const AlignmentGraph& ag, const std::vector<double>& lambda)
{
  const std::size_t V   = ag.vertex_count();
  const VertexId    src = ag.source_id();
  const VertexId    snk = ag.sink_id();

  std::vector<double>   dist(V, INF);
  std::vector<VertexId> prev(V, src);
  dist[src] = 0.0;

  using Entry = std::pair<double, VertexId>;
  std::priority_queue<Entry, std::vector<Entry>, std::greater<>> pq;
  pq.push({0.0, src});

  while (!pq.empty()) {
    auto [d, v] = pq.top();
    pq.pop();
    if (d > dist[v]) continue;

    for (const auto& a : ag.neighbors(v)) {
      int    p  = penalized_vertex(ag, v, a.target);
      double nd = d + a.weight + (p >= 0 ? lambda[p] : 0.0);
      if (nd < dist[a.target]) {
        dist[a.target] = nd;
        prev[a.target] = v;
        pq.push({nd, a.target});
      }
    }
  }

  if (dist[snk] == INF) return std::nullopt;

  std::vector<VertexId> walk;
  for (VertexId v = snk; v != src; v = prev[v])
    walk.push_back(v);
  walk.push_back(src);
  std::reverse(walk.begin(), walk.end());

  return SubResult{std::move(walk), dist[snk]};
}

// Custo original (sem multiplicadores) de um caminho dado no AG.
int
walk_cost(const AlignmentGraph& ag, const std::vector<VertexId>& walk)
{
  int cost = 0;
  for (std::size_t i = 0; i + 1 < walk.size(); ++i) {
    for (const auto& a : ag.neighbors(walk[i]))
      if (a.target == walk[i + 1]) { cost += a.weight; break; }
  }
  return cost;
}

// in_v(v) do caminho: quantos arcos entram em cada vertice local.
// Tambem indica se o caminho ja eh simples (todos os contadores <= 1).
std::vector<int>
count_in(const AlignmentGraph& ag, const std::vector<VertexId>& walk)
{
  std::vector<int> cnt(ag.graph_vertex_count(), 0);
  for (std::size_t i = 0; i + 1 < walk.size(); ++i) {
    int p = penalized_vertex(ag, walk[i], walk[i + 1]);
    if (p >= 0) ++cnt[p];
  }
  return cnt;
}

bool
is_simple(const std::vector<int>& cnt)
{
  for (int c : cnt) if (c > 1) return false;
  return true;
}

// Arco reverso pre-computado: para cada vertice destino guardamos a
// origem, o peso original e o vertice penalizado, permitindo recomputar
// o custo modificado (w + lambda_p) para qualquer lambda sem reconstruir
// a topologia.
struct RevArc { VertexId from; int w; int p; };
using RevAdj = std::vector<std::vector<RevArc>>;

// Constroi a adjacencia reversa uma unica vez (topologia fixa).
RevAdj
build_reverse_adj(const AlignmentGraph& ag)
{
  const std::size_t V = ag.vertex_count();
  RevAdj rev(V);
  for (VertexId u = 0; u < V; ++u)
    for (const auto& a : ag.neighbors(u))
      rev[a.target].push_back({u, a.weight, penalized_vertex(ag, u, a.target)});
  return rev;
}

// Cost-to-go modificado h[v] = custo minimo de v ate o sink sob os custos
// c_bar = w + lambda_p. Dijkstra para tras a partir do sink (custos
// nao-negativos). Usado como potencial para guiar o reparo por custo
// reduzido.
std::vector<double>
backward_cost_to_go(const AlignmentGraph& ag,
                    const RevAdj&         rev,
                    const std::vector<double>& lambda)
{
  const std::size_t V   = ag.vertex_count();
  const VertexId    snk = ag.sink_id();

  std::vector<double> h(V, INF);
  h[snk] = 0.0;

  using Entry = std::pair<double, VertexId>;
  std::priority_queue<Entry, std::vector<Entry>, std::greater<>> pq;
  pq.push({0.0, snk});

  while (!pq.empty()) {
    auto [d, t] = pq.top();
    pq.pop();
    if (d > h[t]) continue;
    for (const RevArc& r : rev[t]) {       // arco r.from -> t
      double m  = r.w + (r.p >= 0 ? lambda[r.p] : 0.0);
      double nd = d + m;
      if (nd < h[r.from]) { h[r.from] = nd; pq.push({nd, r.from}); }
    }
  }
  return h;
}

// Heuristica de reparo por BEAM guiada pelos custos reduzidos Lagrangeanos.
//
// Para o lambda corrente, computa h[v] (cost-to-go modificado ate o sink)
// e mantem uma frontier de ate B caminhos parciais source->..., expandidos
// em rodadas sincronas. Cada parcial guarda o conjunto de vertices locais
// ja consumidos (used); um arco que entraria num vertice ja usado e
// proibido. Os sucessores sao ranqueados por f = g + c_bar(cur->t) + h[t]
// (g = custo modificado acumulado), equivalente ao custo reduzido
// acumulado; mantem-se as B melhores. Arcos livres (delecao/dummy, p=-1)
// nunca sao bloqueados e garantem alcancar o sink, logo a beam sempre
// produz pelo menos um caminho simples viavel. B=1 reproduz o guloso.
//
// Poda por dominancia relaxada: por rodada mantem-se, no maximo, o melhor f
// por vertice AG (ignorando o conjunto used no dedup), limitando a frontier
// a |V| antes do corte para B. Retorna o walk no AG do melhor candidato
// completo por custo ORIGINAL (walk_cost), nao o modificado.
std::optional<std::vector<VertexId>>
repair_primal(const AlignmentGraph&      ag,
              const RevAdj&              rev,
              const std::vector<double>& lambda,
              int                        beam_width)
{
  const VertexId    src = ag.source_id();
  const VertexId    snk = ag.sink_id();
  const std::size_t V   = ag.vertex_count();
  const std::size_t B   = std::max(1, beam_width);

  const std::vector<double> h = backward_cost_to_go(ag, rev, lambda);
  if (h[src] == INF) return std::nullopt; // source nao alcanca o sink

  struct State {
    VertexId              cur;   // vertice AG atual
    double                g;     // custo modificado acumulado
    double                f;     // g + h[cur] (chave de ordenacao)
    std::vector<char>     used;  // local_ids ja consumidos
    std::vector<VertexId> walk;  // source..cur (ids do AG)
  };

  const std::size_t G = ag.graph_vertex_count();

  std::vector<State> frontier;
  frontier.push_back({src, 0.0, h[src], std::vector<char>(G, 0),
                      std::vector<VertexId>{src}});

  std::optional<State> best_complete;
  std::size_t          rounds = 0;

  while (!frontier.empty() && rounds++ <= V) {
    std::vector<State> succ;
    for (State& s : frontier) {
      if (s.cur == snk) {
        if (!best_complete ||
            walk_cost(ag, s.walk) < walk_cost(ag, best_complete->walk))
          best_complete = s;
        continue;
      }
      for (const auto& a : ag.neighbors(s.cur)) {
        int p = penalized_vertex(ag, s.cur, a.target);
        if (p >= 0 && s.used[p]) continue;   // entraria num vertice ja usado
        if (h[a.target] == INF) continue;    // nao alcanca o sink
        double cbar = a.weight + (p >= 0 ? lambda[p] : 0.0);
        double ng   = s.g + cbar;
        State  ns{a.target, ng, ng + h[a.target], s.used, s.walk};
        if (p >= 0) ns.used[p] = 1;
        ns.walk.push_back(a.target);
        succ.push_back(std::move(ns));
      }
    }
    if (succ.empty()) break;

    // Dominancia relaxada: melhor f por vertice AG (ignora used no dedup).
    std::sort(succ.begin(), succ.end(),
              [](const State& a, const State& b) { return a.f < b.f; });
    std::vector<char>  seen(V, 0);
    std::vector<State> pruned;
    for (State& s : succ) {
      if (seen[s.cur]) continue;
      seen[s.cur] = 1;
      pruned.push_back(std::move(s));
      if (pruned.size() >= B) break;
    }
    frontier = std::move(pruned);
  }

  if (!best_complete) return std::nullopt;
  return best_complete->walk;
}

// Reconstroi um walk no AG a partir de um .sol (arcos com valor >= 0.5).
// Usado para semear um UB inicial. Best-effort: devolve nullopt se o .sol
// nao formar um caminho S->sink consistente.
std::optional<std::vector<VertexId>>
walk_from_sol(const AlignmentGraph& ag, const std::string& sol_path)
{
  std::unordered_map<std::string, double> values;
  try { values = read_gurobi_sol(sol_path); }
  catch (...) { return std::nullopt; }

  const std::size_t V   = ag.vertex_count();
  const VertexId    src = ag.source_id();
  const VertexId    snk = ag.sink_id();

  std::vector<VertexId> walk{src};
  VertexId    cur    = src;
  std::size_t safety = 0;
  while (cur != snk && safety++ <= V) {
    auto   neigh = ag.neighbors(cur);
    bool   moved = false;
    for (std::size_t j = 0; j < neigh.size(); ++j) {
      auto it = values.find(arc_var_name(cur, j));
      if (it != values.end() && it->second >= 0.5) {
        cur = neigh[j].target;
        walk.push_back(cur);
        moved = true;
        break;
      }
    }
    if (!moved) return std::nullopt;
  }
  if (cur != snk) return std::nullopt;
  return walk;
}

// Converte um walk no AG em atribuicoes de variaveis (x_u_j = 1) para
// gravar um .sol compativel com o PLI / warm start.
std::vector<std::pair<std::string, double>>
walk_to_assignments(const AlignmentGraph& ag, const std::vector<VertexId>& walk)
{
  std::vector<std::pair<std::string, double>> assigns;
  for (std::size_t i = 0; i + 1 < walk.size(); ++i) {
    auto neigh = ag.neighbors(walk[i]);
    for (std::size_t j = 0; j < neigh.size(); ++j)
      if (neigh[j].target == walk[i + 1]) {
        assigns.emplace_back(arc_var_name(walk[i], j), 1.0);
        break;
      }
  }
  return assigns;
}

} // namespace

std::optional<SolverResult>
solve_lagrangian(const AlignmentGraph&   ag,
                 const ExperimentConfig& cfg,
                 const std::string&      sol_save_path,
                 const std::string&      warm_start_path)
{
  const std::size_t G       = ag.graph_vertex_count();
  const bool        verbose = cfg.pli_verbose;

  std::vector<double> lambda(G, 0.0);

  // Adjacencia reversa pre-computada uma vez; reaproveitada a cada
  // iteracao pelo reparo guiado por custo reduzido (cost-to-go modificado).
  const RevAdj rev = build_reverse_adj(ag);

  // Largura da beam do reparo primal, usada literalmente (B=1 reproduz o
  // guloso). O clamp para >= 1 fica em repair_primal.
  const int beam_B = cfg.lagr_repair_beam;
  if (verbose) std::cerr << "  [lagr] reparo beam_width=" << beam_B << "\n";

  // ---- UB primal inicial: warm start (se houver) e/ou reparo ----
  double                best_ub = INF;
  std::vector<VertexId> best_walk;

  auto try_ub = [&](const std::vector<VertexId>& walk) {
    if (walk.empty()) return;
    if (!is_simple(count_in(ag, walk))) return; // garante caminho simples
    int c = walk_cost(ag, walk);
    if (c < best_ub) { best_ub = c; best_walk = walk; }
  };

  if (!warm_start_path.empty()) {
    if (auto w = walk_from_sol(ag, warm_start_path)) try_ub(*w);
  }

  // ---- método do subgradiente ----
  double best_lb        = -INF;
  double alpha          = cfg.lagr_alpha;
  int    no_improve     = 0;
  bool   proven_optimal = false;

  const auto   t_start = std::chrono::steady_clock::now();
  const double tl_s    = cfg.lagr_time_limit_s;
  int          iter    = 0;

  for (; iter < cfg.lagr_max_iters; ++iter) {
    if (tl_s > 0.0) {
      double el = std::chrono::duration<double>(
                      std::chrono::steady_clock::now() - t_start).count();
      if (el >= tl_s) break;
    }

    auto sub = solve_subproblem(ag, lambda);
    if (!sub) return std::nullopt; // AG sem caminho S->T (nao deve ocorrer)

    double sum_lambda = 0.0;
    for (double l : lambda) sum_lambda += l;

    const double L = sub->modified_cost - sum_lambda; // limite inferior
    if (L > best_lb) { best_lb = L; no_improve = 0; }
    else             { ++no_improve; }

    // Subgradiente g_v = in_v(v*) - 1, com refinamento de Beasley:
    // zera componentes inativas com multiplicador nulo.
    auto   cnt    = count_in(ag, sub->ag_walk);
    double norm2  = 0.0;
    std::vector<double> g(G);
    for (std::size_t v = 0; v < G; ++v) {
      double gv = static_cast<double>(cnt[v]) - 1.0;
      if (gv < 0.0 && lambda[v] == 0.0) gv = 0.0;
      g[v]   = gv;
      norm2 += gv * gv;
    }

    // Se o passeio do subproblema ja eh simples, ele eh primal viavel.
    if (is_simple(cnt)) try_ub(sub->ag_walk);

    // Reparo guiado pelos custos reduzidos do lambda corrente: gera um
    // candidato a UB viavel a cada iteracao.
    if (auto w = repair_primal(ag, rev, lambda, beam_B)) try_ub(*w);

    // norm2 == 0: solucao do subproblema satisfaz as restricoes ativas
    // com folga complementar => L eh otimo do problema original.
    if (norm2 == 0.0) { proven_optimal = true; break; }

    // UB de referencia para o passo. Se ainda nao ha UB finito, usa uma
    // folga relativa sobre L para nao travar o passo.
    const double ub_ref = std::isfinite(best_ub) ? best_ub
                                                 : L + std::abs(L) + 1.0;
    const double step = alpha * (ub_ref - L) / norm2;

    for (std::size_t v = 0; v < G; ++v)
      lambda[v] = std::max(0.0, lambda[v] + step * g[v]);

    if (no_improve >= cfg.lagr_alpha_patience) {
      alpha *= 0.5;
      no_improve = 0;
    }

    if (verbose)
      std::cerr << "  [lagr] it=" << iter
                << " L=" << L
                << " lb=" << best_lb
                << " ub=" << (std::isfinite(best_ub) ? best_ub : -1.0)
                << " alpha=" << alpha
                << " ||g||^2=" << norm2 << "\n";

    // Gap inteiro fechado: limite inferior alcanca o melhor inteiro viavel.
    if (std::isfinite(best_ub) && std::ceil(best_lb - 1e-9) >= best_ub) {
      proven_optimal = true;
      break;
    }
  }

  // ---- monta o resultado ----
  SolverResult result;
  result.variant  = "lagrangian";
  result.final_lb = std::isfinite(best_lb) ? best_lb : std::nan("");

  const bool have_ub = std::isfinite(best_ub) && !best_walk.empty();

  if (have_ub) {
    result.final_ub = best_ub;
    result.cost     = static_cast<int>(best_ub);
    result.ag_walk  = best_walk;
    result.has_path = true;
  } else {
    result.cost     = -1;
    result.has_path = false;
  }

  // Status proprio do metodo (sem Gurobi); o CSV do sgap-heuristic rotula
  // via lagr_status_label().
  if (proven_optimal && have_ub)
    result.status = LAGR_OPTIMAL;
  else if (tl_s > 0.0 && iter < cfg.lagr_max_iters)
    result.status = LAGR_TIME_LIMIT;
  else
    result.status = LAGR_ITERATION_LIMIT;

  // Grava o melhor caminho viavel como .sol compativel com o PLI.
  if (!sol_save_path.empty() && have_ub) {
    auto assigns = walk_to_assignments(ag, best_walk);
    std::vector<std::string> header = {
        "algorithm=lagrangian",
        "cost=" + std::to_string(static_cast<int>(best_ub)),
    };
    try { write_gurobi_sol(sol_save_path, assigns, header); }
    catch (...) { /* falha de IO no .sol nao invalida o resultado */ }
  }

  return result;
}

} // namespace seq
