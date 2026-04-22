#pragma once

#include "../core/types.h"
#include "../core/sequence_graph.h"
#include "../core/sequence.h"
#include "../core/alignment_graph.h"
#include "../algorithms/alignment_potentials.h"
#include <vector>
#include <optional>

namespace seq {

struct BeamSearchParams {
  int beam_width = 1000;
};

struct BeamSearchResult {
  std::vector<VertexId> vertices;  // caminho simples no SG
  int cost;                        // edit_distance(induced_sequence, query)
};

// Beam search no grafo de alinhamento usando potenciais backward como
// heuristica. Constroi caminho simples camada por camada, mantendo
// os k estados mais promissores por camada com bitset visited por estado.
// Retorna o melhor caminho encontrado, ou nullopt se nenhum alcanca o fim.
std::optional<BeamSearchResult> beam_search_alignment(
    const SequenceGraph& graph,
    const Sequence& query,
    const AlignmentPotentials& potentials,
    const AlignmentGraph& ag,
    BeamSearchParams params = {});

} // namespace seq
