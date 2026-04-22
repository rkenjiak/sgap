#pragma once

#include "../core/alignment_graph.h"
#include "alignment_walk.h"
#include <vector>

namespace seq {

// Potenciais de alinhamento: delta+(i,v) e delta-(i,v) para cada
// vertice do grafo de alinhamento. Permite calcular custo reduzido
// r(i,v) = delta+(i,v) + delta-(i,v) - d*(q,G), que mede quao longe
// o vertice v na camada i esta de participar em um walk otimo.
class AlignmentPotentials {
public:
  // Computa forward e backward (shortest path de source/sink).
  explicit AlignmentPotentials(
      const AlignmentGraph& ag,
      ShortestPathMethod method = ShortestPathMethod::BFS01);

  // Potencial forward: distancia minima de source ate ag_id.
  int forward(VertexId ag_id) const;

  // Potencial backward: distancia minima de ag_id ate sink.
  int backward(VertexId ag_id) const;

  // Potencial backward por (camada, vertice local).
  int backward(std::size_t layer, std::size_t local_v) const;

  // Custo reduzido: forward(v) + backward(v) - walk_cost.
  // Sempre >= 0 para vertices alcancaveis.
  int reduced_cost(VertexId ag_id) const;

  // Custo do walk otimo (source a sink).
  int walk_cost() const noexcept;

private:
  std::vector<int> forward_;
  std::vector<int> backward_;
  int walk_cost_;
  const AlignmentGraph* ag_;
};

} // namespace seq
