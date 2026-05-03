#include "pli.h"
#include "../io/gurobi_sol.h"
#include "gurobi_c++.h"

#include <cmath>
#include <filesystem>
#include <string>
#include <vector>

namespace seq {

std::optional<SolverResult>
solve_pli(const AlignmentGraph&   ag,
          const ExperimentConfig& cfg,
          const std::string&      sol_save_path,
          const std::string&      warm_start_path,
          const std::string&      lp_save_path)
{
  const std::size_t N           = ag.vertex_count();
  const std::size_t G           = ag.graph_vertex_count(); // n (sem dummy)
  const VertexId    source      = ag.source_id();
  const VertexId    sink        = ag.sink_id();
  const std::size_t dummy_local = G;

  GRBEnv env(true);
  const int log_flag = cfg.pli_verbose ? 1 : 0;
  env.set(GRB_IntParam_OutputFlag, log_flag);
  env.set(GRB_IntParam_LogToConsole, log_flag);
  env.start();

  GRBModel model(env);
  if (cfg.pli_time_limit_s > 0.0)
    model.set(GRB_DoubleParam_TimeLimit, cfg.pli_time_limit_s);

  struct ArcRec { VertexId u, v; int w; GRBVar x; };
  std::vector<ArcRec> arcs;
  arcs.reserve(ag.arc_count());

  std::vector<GRBLinExpr> out_expr(N), in_expr(N);
  std::vector<GRBLinExpr> in_v_expr(G);
  GRBLinExpr obj;

  for (VertexId u = 0; u < N; ++u) {
    auto neigh = ag.neighbors(u);
    for (std::size_t j = 0; j < neigh.size(); ++j) {
      const auto& a = neigh[j];

      GRBVar x = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, arc_var_name(u, j));
      arcs.push_back({u, a.target, a.weight, x});

      out_expr[u]       += x;
      in_expr[a.target] += x;
      obj               += a.weight * x;

      if (ag.is_sink(a.target)) continue;
      auto tgt = ag.from_vertex_id(a.target);
      if (tgt.local_id == dummy_local) continue;

      std::size_t src_local;
      if (ag.is_source(u))
        src_local = static_cast<std::size_t>(-1);
      else
        src_local = ag.from_vertex_id(u).local_id;

      if (src_local != tgt.local_id)
        in_v_expr[tgt.local_id] += x;
    }
  }

  model.setObjective(obj, GRB_MINIMIZE);

  // fluxo
  for (VertexId k = 0; k < N; ++k) {
    double rhs = 0.0;
    if      (k == source) rhs =  1.0;
    else if (k == sink)   rhs = -1.0;
    model.addConstr(out_expr[k] - in_expr[k] == rhs);
  }

  // garantir caminho simples
  for (std::size_t v = 0; v < G; ++v)
    model.addConstr(in_v_expr[v] <= 1.0);

  // Warm start: aplica os valores do .sol no atributo Start. Atribuir
  // GRB_UNDEFINED nas demais variaveis nao eh necessario (Gurobi as
  // trata como livres) mas reforcamos a consistencia.
  if (!warm_start_path.empty() && std::filesystem::exists(warm_start_path)) {
    auto values = read_gurobi_sol(warm_start_path);
    for (VertexId u = 0, idx = 0; u < N; ++u) {
      auto neigh = ag.neighbors(u);
      for (std::size_t j = 0; j < neigh.size(); ++j, ++idx) {
        auto it = values.find(arc_var_name(u, j));
        double v = (it != values.end()) ? it->second : GRB_UNDEFINED;
        arcs[idx].x.set(GRB_DoubleAttr_Start, v);
      }
    }
  }

  // Atualiza o modelo para que a escrita do .lp reflita variaveis,
  // restricoes e MIP start adicionados acima.
  if (!lp_save_path.empty()) {
    model.update();
    auto parent = std::filesystem::path(lp_save_path).parent_path();
    if (!parent.empty())
      std::filesystem::create_directories(parent);
    model.write(lp_save_path);
  }

  model.optimize();

  // Aceita otimo ou solucao viavel encontrada antes do time limit.
  const int status = model.get(GRB_IntAttr_Status);
  const bool has_solution = model.get(GRB_IntAttr_SolCount) > 0;
  if (status != GRB_OPTIMAL && !has_solution)
    return std::nullopt;

  if (!sol_save_path.empty()) {
    auto parent = std::filesystem::path(sol_save_path).parent_path();
    if (!parent.empty())
      std::filesystem::create_directories(parent);
    model.write(sol_save_path);
  }

  std::vector<std::vector<std::size_t>> sel_out(N);
  for (std::size_t i = 0; i < arcs.size(); ++i)
    if (arcs[i].x.get(GRB_DoubleAttr_X) >= 0.5)
      sel_out[arcs[i].u].push_back(i);

  std::vector<VertexId> walk;
  walk.reserve(ag.layer_count() + 2);
  walk.push_back(source);

  VertexId    cur    = source;
  std::size_t safety = 0;
  while (cur != sink && safety++ <= N) {
    auto& outs = sel_out[cur];
    if (outs.empty()) return std::nullopt;
    cur = arcs[outs.front()].v;
    walk.push_back(cur);
  }
  if (cur != sink) return std::nullopt;

  int cost = static_cast<int>(std::lround(model.get(GRB_DoubleAttr_ObjVal)));
  return SolverResult{cost, std::move(walk), "path"};
}

} // namespace seq
