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

  // Custos de alinhamento
  AlignmentCosts costs = {1, 1, 1};

  // PLI (Gurobi)
  bool   pli_verbose      = false;  // GRB OutputFlag / LogToConsole
  double pli_time_limit_s = 0.0;    // 0 = sem limite
};

ExperimentConfig read_config(const std::string& path);

} // namespace seq
