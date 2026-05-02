#pragma once

#include "../core/alignment_graph.h"
#include "../core/sequence.h"
#include "../core/sequence_graph.h"
#include <cstddef>
#include <optional>
#include <span>
#include <utility>
#include <vector>

namespace seq {

// Resultado de build_warm_start_arcs: lista de arcos do AlignmentGraph
// que correspondem a um caminho source -> sink consistente com um walk
// em G previamente encontrado por uma heuristica.
struct WarmStartArcs {
  // Cada par (u, j) referencia o j-esimo vizinho de u em ag.neighbors(u).
  // Em conjunto formam um caminho source -> sink. O ultimo arco eh o
  // (final_layer_vertex, sink) com peso 0.
  std::vector<std::pair<VertexId, std::size_t>> arcs;

  // Custo total dos arcos escolhidos (= edit_dist consistente com os
  // costs do AG; pode ser <= ao edit_dist reportado pela heuristica
  // caso ela use um alinhamento sub-otimo entre Q e induced(walk_g)).
  int cost = 0;
};

// Dado um walk em G (saida de uma heuristica), monta os arcos do
// AlignmentGraph correspondentes ao melhor alinhamento entre query e
// induced(walk_g), respeitando a topologia do AG (source, sink,
// dummy, layers).
//
// Pre-condicoes:
//   - query nao vazia, walk_g nao vazio.
//   - walk_g eh um walk valido em g (cada par consecutivo eh arco).
//   - ag construido com (g, query, costs). costs deve coincidir.
//
// Retorna nullopt se algo violar essas condicoes ou se nao for possivel
// localizar algum arco esperado em ag.neighbors (ag fora de sincronia
// com os parametros que geraram o walk_g).
std::optional<WarmStartArcs>
build_warm_start_arcs(const AlignmentGraph&     ag,
                      const SequenceGraph&      g,
                      const Sequence&           query,
                      std::span<const VertexId> walk_g,
                      AlignmentCosts            costs);

} // namespace seq
