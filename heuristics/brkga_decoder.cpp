#include "brkga_decoder.h"
#include "../algorithms/edit_distance.h"
#include <limits>

namespace seq {

SGAPDecoder::SGAPDecoder(const SequenceGraph& g, const Sequence& query)
  : g_(g), query_(query)
{}

std::size_t SGAPDecoder::chromosome_length() const
{
  return g_.vertex_count();
}

std::vector<VertexId>
SGAPDecoder::build_path(const std::vector<double>& key, int* cost_out) const
{
  const std::size_t n = g_.vertex_count();
  std::vector<VertexId> path;

  if (n == 0 || key.size() != chromosome_length()) {
    if (cost_out) *cost_out = std::numeric_limits<int>::max();
    return path;
  }

  // start = argmax_v k[v]
  VertexId start = 0;
  double   best_start = -std::numeric_limits<double>::infinity();
  for (std::size_t v = 0; v < n; ++v) {
    if (key[v] > best_start) {
      best_start = key[v];
      start = static_cast<VertexId>(v);
    }
  }

  std::vector<bool> visited(n, false);
  visited[start] = true;
  path.push_back(start);
  VertexId u = start;

  // Estende escolhendo o vizinho nao-visitado de menor chave.
  // Para ao chegar no tamanho da query: induced_sequence() emite 1 simbolo
  // por vertice, entao path.size() == |induced|, e uma induzida maior que a
  // query so adiciona custo (edit_distance >= ||induced| - |query||).
  while (path.size() < query_.size()) {
    auto nbrs = g_.neighbors(u);

    VertexId best = 0;
    double   best_score = std::numeric_limits<double>::infinity();
    bool     found = false;
    for (VertexId v : nbrs) {
      if (visited[v]) continue;
      if (key[v] < best_score) {
        best_score = key[v];
        best   = v;
        found  = true;
      }
    }
    if (!found) break;
    visited[best] = true;
    path.push_back(best);
    u = best;
  }

  if (cost_out) {
    auto induced = g_.induced_sequence(path);
    *cost_out = edit_distance(induced, query_);
  }
  return path;
}

double
SGAPDecoder::decode(const std::vector<double>& chromosome) const
{
  int cost = 0;
  build_path(chromosome, &cost);
  return static_cast<double>(cost);
}

std::vector<VertexId>
SGAPDecoder::decode_to_path(const std::vector<double>& chromosome) const
{
  return build_path(chromosome, nullptr);
}

} // namespace seq
