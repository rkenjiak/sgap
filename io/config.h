#pragma once

#include "../core/alignment_graph.h"
#include <string>

namespace seq {

struct ExperimentConfig {
  // Algoritmo: "jain" | "greedy" | "beam"
  std::string algorithm = "jain";

  // Jain: "dijkstra" | "bfs01"
  std::string jain_method = "bfs01";

  // Greedy
  int         greedy_restarts  = 1000;
  std::string greedy_selection = "roulette"; // "greedy" | "roulette"
  unsigned    greedy_seed      = 42;

  // Beam Search
  int         beam_width      = 1000;
  std::string beam_width_mode = "constant"; // "constant" | "n" | "density_n"
  std::string beam_potentials = "bfs01";   // "dijkstra" | "bfs01"

  // BRKGA (modo basico, decoder vertex-priority)
  int      brkga_pop_size    = 256;   // p
  double   brkga_pe          = 0.20;  // fracao elite
  double   brkga_pm          = 0.10;  // fracao mutantes
  double   brkga_rhoe        = 0.70;  // prob. heranca elite
  int      brkga_generations = 500;
  unsigned brkga_seed        = 42;

  // Custos de alinhamento
  AlignmentCosts costs = {1, 1, 1};

  // PLI (Gurobi)
  bool      pli_verbose      = false;  // GRB OutputFlag / LogToConsole
  double    pli_time_limit_s = 0.0;    // 0 = sem limite
  long long pli_node_limit   = -1;     // <0 = sem limite; 0 = so no raiz

  // Relaxacao Lagrangeana (algorithm = "lagrangian").
  // Dualiza a restricao de caminho simples; o subproblema vira um
  // caminho minimo S->T resolvido por Dijkstra e o dual eh maximizado
  // pelo metodo do subgradiente.
  int       lagr_max_iters      = 500;  // teto de iteracoes do subgradiente
  double    lagr_time_limit_s   = 0.0;  // 0 = sem limite
  double    lagr_alpha          = 2.0;  // passo inicial do subgradiente
  int       lagr_alpha_patience = 20;   // halva alpha apos N iters sem melhora do LB
  // Largura da beam do reparo primal (guiado por custos reduzidos):
  //   0 = auto (min(densidade*|V|, 64));
  //   1 = guloso puro (comportamento original); N = largura fixa.
  int       lagr_repair_beam    = 0;

  // Persistencia em formato .sol do Gurobi (compartilhado por
  // sgap-solver e sgap-heuristic).
  // gurobi_sol_dir       : se nao vazio, ambos os binarios gravam
  //                        <dir>/<grafo>_<seq>_<config>.sol.
  // pli_warm_start_path  : caminho explicito de um .sol a ser usado
  //                        como MIP start no PLI. Permite escolher
  //                        externamente qual heuristica seeda o solve.
  //                        Ignorado se vazio ou se o arquivo nao
  //                        existir.
  std::string gurobi_sol_dir;
  std::string pli_warm_start_path;

  // gurobi_lp_dir : se nao vazio, sgap-solver grava o modelo do PLI em
  //                 <dir>/<grafo>_<seq>_<config>.lp antes do optimize().
  std::string gurobi_lp_dir;

  // heuristic_out_dir : se nao vazio, sgap-heuristic grava relatorio
  //                     humano-legivel em
  //                     <dir>/<grafo>_<seq>_<config>.out.
  std::string heuristic_out_dir;
};

ExperimentConfig read_config(const std::string& path);

} // namespace seq
