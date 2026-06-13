#include "brkga_search.h"
#include "brkga_decoder.h"

// brkgaAPI (submodulo em third_party/brkgaAPI/brkgaAPI/)
#include "BRKGA.h"
#include "MTRand.h"

namespace seq {

std::optional<BrkgaResult>
brkga_alignment(const SequenceGraph& graph,
                const Sequence&      query,
                BrkgaParams          params)
{
  const std::size_t n = graph.vertex_count();
  if (n == 0 || query.empty()) return std::nullopt;

  // BRKGA exige p > 0; geracoes >= 2 e par (evolve interno itera de 1).
  if (params.pop_size == 0) return std::nullopt;
  unsigned generations = params.generations;
  if (generations == 0) generations = 2;
  if (generations % 2 != 0) generations += 1;

  SGAPDecoder decoder(graph, query);
  MTRand rng(params.seed);

  using BRKGAType = BRKGA<SGAPDecoder, MTRand>;
  BRKGAType algo(
      static_cast<unsigned>(decoder.chromosome_length()),
      params.pop_size,
      params.pe,
      params.pm,
      params.rhoe,
      decoder,
      rng);   // K = 1, MAX_THREADS = 1 (defaults)

  algo.evolve(generations);

  const auto& best_chromo = algo.getBestChromosome();
  auto path = decoder.decode_to_path(best_chromo);
  if (path.empty()) return std::nullopt;

  int cost = static_cast<int>(algo.getBestFitness());
  return BrkgaResult{std::move(path), cost};
}

} // namespace seq
