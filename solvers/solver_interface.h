#pragma once

#include "../core/alignment_graph.h"
#include <cmath>
#include <optional>
#include <string>
#include <vector>

namespace seq {

struct SolverResult {
  int                   cost;
  std::vector<VertexId> ag_walk;  // caminho no grafo de alinhamento (inclui S e T)
  std::string           variant;  // ex: "path", "lagrangian_k50", "colgen", "benders"

  // Bounds e status reportados pelo solver. NaN quando nao aplicavel
  // (ex.: sem solucao incumbente -> final_ub=NaN, mip_gap=NaN).
  double final_lb = std::nan("");
  double final_ub = std::nan("");
  double mip_gap  = std::nan("");
  int    status   = -1;            // codigo bruto do solver (ex.: GRB_*)
  bool   has_path = false;         // se ag_walk representa um caminho viavel
};

} // namespace seq
