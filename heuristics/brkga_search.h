#pragma once

#include "../core/types.h"
#include "../core/sequence_graph.h"
#include "../core/sequence.h"
#include <vector>
#include <optional>

namespace seq {

struct BrkgaParams {
  // Hiperparametros principais (vide BRKGA.h da brkgaAPI).
  unsigned pop_size     = 256;   // p
  double   pe           = 0.20;  // fracao elite
  double   pm           = 0.10;  // fracao mutantes
  double   rhoe         = 0.70;  // prob. de heranca elite

  // Iteracao.
  unsigned generations  = 500;

  unsigned long seed    = 42;
};

struct BrkgaResult {
  std::vector<VertexId> vertices;  // caminho simples em G
  int cost;                        // edit_distance(induced, query)
};

// Roda BRKGA (modo basico, decoder vertex-priority) sobre o SGAP. Retorna a
// melhor solucao encontrada apos generations geracoes, ou nullopt em entrada
// degenerada.
std::optional<BrkgaResult>
brkga_alignment(const SequenceGraph& graph,
                const Sequence&      query,
                BrkgaParams          params = {});

} // namespace seq
