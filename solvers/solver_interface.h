#pragma once

#include "../core/alignment_graph.h"
#include <optional>
#include <string>
#include <vector>

namespace seq {

struct SolverResult {
  int                   cost;
  std::vector<VertexId> ag_walk;  // caminho no grafo de alinhamento (inclui S e T)
  std::string           variant;  // ex: "path", "lagrangian_k50", "colgen", "benders"
};

} // namespace seq
