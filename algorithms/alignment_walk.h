#pragma once

#include "../core/alignment_graph.h"
#include <optional>
#include <vector>

namespace seq {

enum class ShortestPathMethod {
  Dijkstra,  // O(V + E log V) — priority queue classica
  BFS01,     // O(V + E) — deque, requer pesos em {0, 1}
};

// Encontra o passeio (walk) de custo minimo: caminho mais curto
// de source a sink no grafo de alinhamento.
// Sem restricao de repeticao de vertices do grafo de sequencia.
std::optional<AlignmentWalk> walk_alignment(
    const AlignmentGraph& ag,
    ShortestPathMethod method = ShortestPathMethod::BFS01);

// Mesma coisa, mas retorna tambem o vetor de distancias completo
// (indexado por vertex ID do AG). Usado para computar potenciais forward.
struct AlignmentWalkFull {
  int cost;
  std::vector<VertexId> vertices;  // walk no SG (como AlignmentWalk)
  std::vector<int> dist;           // dist[ag_vertex_id] = delta+(ag_vertex_id)
};

std::optional<AlignmentWalkFull> walk_alignment_full(
    const AlignmentGraph& ag,
    ShortestPathMethod method = ShortestPathMethod::BFS01);

} // namespace seq
