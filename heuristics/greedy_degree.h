#pragma once

#include "../core/types.h"
#include "../core/sequence_graph.h"
#include "../core/sequence.h"
#include <vector>
#include <optional>

namespace seq {

struct GreedyDegreeParams {
  // Vertice inicial. Se nao definido, testa todos e retorna o mais longo.
  std::optional<VertexId> start;

  // Modo de selecao do proximo vizinho:
  //   greedy   = sempre o de maior grau de saida (deterministico)
  //   roulette = probabilistico proporcional ao grau de saida
  enum class Selection { greedy, roulette };
  Selection selection = Selection::greedy;

  // Seed para o modo roulette (ignorado no greedy).
  unsigned seed = 0;
};

// Constroi um caminho (sem repeticao de vertices) guiado pelo grau de
// saida do vertice destino. Retorna vetor de VertexIds do caminho,
// ou vetor vazio se nao conseguir.
std::vector<VertexId> greedy_degree_path(
    const SequenceGraph& graph,
    GreedyDegreeParams params = {});

// Resultado do multi-start: caminho + custo (edit_distance a query).
struct GreedyDegreeResult {
  std::vector<VertexId> vertices;
  int cost; // edit_distance(induced_sequence, query)
};

// Multi-start: constroi N caminhos variando seed e vertice inicial,
// retorna o de menor edit_distance a query.
std::optional<GreedyDegreeResult> greedy_degree_best(
    const SequenceGraph& graph,
    const Sequence& query,
    int restarts,
    GreedyDegreeParams base_params = {});

} // namespace seq
