#pragma once

#include "../core/types.h"
#include "../core/sequence_graph.h"
#include "../core/sequence.h"
#include <vector>

namespace seq {

// Decoder vertex-priority (modo basico / "vp") para BRKGA aplicado ao SGAP.
//
// Cromossomo = vector<double> de tamanho |V|, cada chave (allele) em [0,1):
//   - start  = argmax_v k[v]
//   - de u, escolhe o vizinho nao-visitado de menor k[v]
//   - para quando nao ha vizinho nao-visitado
//
// Fitness = edit_distance(g.induced_sequence(path), query)  (menor = melhor).
//
// O decoder so le estado const (grafo, query) e aloca buffers locais por
// chamada -- seguro para MAX_THREADS > 1 (sem mutable).
class SGAPDecoder {
public:
  SGAPDecoder(const SequenceGraph& g, const Sequence& query);

  // Tamanho esperado do cromossomo: |V|.
  std::size_t chromosome_length() const;

  // API exigida pela brkgaAPI: retorna fitness (menor = melhor).
  double decode(const std::vector<double>& chromosome) const;

  // Mesmo loop de decode(), mas retorna o caminho em vez do custo.
  std::vector<VertexId> decode_to_path(
      const std::vector<double>& chromosome) const;

private:
  const SequenceGraph& g_;
  const Sequence&      query_;

  // Constroi o caminho vertex-priority; cost_out != nullptr => escreve
  // edit_distance(induced, query).
  std::vector<VertexId> build_path(const std::vector<double>& key,
                                   int* cost_out) const;
};

} // namespace seq
