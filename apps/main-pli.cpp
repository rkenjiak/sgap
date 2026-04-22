/*
 * Solver exato para alinhamento de sequencia a grafo.
 *
 * Dispatcher: le cfg.algorithm e despacha para o solver adequado.
 * Todos os solvers retornam SolverResult e produzem CSV compativel
 * com sgap-heuristic.
 *
 * Uso:
 *   sgap-solver <grafo.gfa> <sequencia.fastq> <config.ini> [saida.csv]
 */

#include "core/sequence.h"
#include "core/sequence_graph.h"
#include "core/alignment_graph.h"
#include "io/gfa.h"
#include "io/fastq.h"
#include "io/config.h"
#include "solvers/pli.h"
// #include "solvers/lagrangian.h"
// #include "solvers/column_gen.h"
// #include "solvers/benders.h"

#include "gurobi_c++.h"

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

namespace fs = std::filesystem;
using namespace seq;

static const char* csv_header =
    "graph,sequence,vertices,arcs,query_len,"
    "algorithm,variant,edit_dist,induced_len,is_path,time_ms\n";

static void
run(const fs::path& gfa_path,
    const fs::path& fastq_path,
    const fs::path& cfg_path,
    std::ostream&   out)
{
  SequenceGraph    g     = read_gfa(gfa_path);
  Sequence         query = read_fastq(fastq_path);
  ExperimentConfig cfg   = read_config(cfg_path.string());

  std::string graph_name = gfa_path.stem().string();
  std::string seq_name   = fastq_path.stem().string();
  std::string algorithm  = cfg.algorithm;

  std::cerr << graph_name << " x " << seq_name << "  ";

  AlignmentGraph ag(g, query, cfg.costs);

  // ----- despacho por algoritmo -----
  std::optional<SolverResult> result;

  auto t0 = std::chrono::steady_clock::now();

  if (algorithm == "pli") {
    result = solve_pli(ag, cfg);
  }
  // else if (algorithm == "lagrangian") {
  //   result = solve_lagrangian(ag, cfg);
  // }
  // else if (algorithm == "colgen") {
  //   result = solve_colgen(ag, cfg);
  // }
  // else if (algorithm == "benders") {
  //   result = solve_benders(ag, cfg);
  // }
  else {
    std::cerr << "algoritmo desconhecido para solver: " << algorithm << "\n";
    return;
  }

  auto t1 = std::chrono::steady_clock::now();
  double time_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

  int  edit_dist   = -1;
  int  induced_len = 0;
  bool is_path     = false;
  std::string variant = "?";

  if (result) {
    auto walk    = ag.to_graph_walk(result->ag_walk);
    auto induced = g.induced_sequence(walk);
    edit_dist    = result->cost;
    induced_len  = static_cast<int>(induced.size());
    is_path      = g.is_valid_path(walk);
    variant      = result->variant;
  }

  out << graph_name << ","
      << seq_name << ","
      << g.vertex_count() << ","
      << g.arc_count() << ","
      << query.size() << ","
      << algorithm << ","
      << variant << ","
      << edit_dist << ","
      << induced_len << ","
      << (is_path ? 1 : 0) << ","
      << std::fixed << std::setprecision(2)
      << time_ms << "\n"
      << std::defaultfloat;

  std::cerr << algorithm << "(" << variant << ") edit=" << edit_dist << "\n";
}

int
main(int argc, char* argv[])
{
  if (argc < 4) {
    std::cerr << "uso: " << argv[0]
              << " <grafo.gfa> <sequencia.fastq> <config.ini> [saida.csv]\n";
    return 1;
  }

  fs::path gfa_path   = argv[1];
  fs::path fastq_path = argv[2];
  fs::path cfg_path   = argv[3];

  if (!fs::exists(gfa_path)) {
    std::cerr << "arquivo gfa nao encontrado: " << gfa_path << "\n";
    return 1;
  }
  if (!fs::exists(fastq_path)) {
    std::cerr << "arquivo fastq nao encontrado: " << fastq_path << "\n";
    return 1;
  }
  if (!fs::exists(cfg_path)) {
    std::cerr << "arquivo de config nao encontrado: " << cfg_path << "\n";
    return 1;
  }

  try {
    if (argc > 4) {
      fs::path csv_path = argv[4];
      if (csv_path.has_parent_path())
        fs::create_directories(csv_path.parent_path());

      bool file_exists = fs::exists(csv_path);
      std::ofstream csv(csv_path, std::ios::app);
      if (!file_exists) csv << csv_header;

      run(gfa_path, fastq_path, cfg_path, csv);
      csv.close();
      std::cerr << "resultado escrito em " << csv_path << "\n";
    } else {
      std::cout << csv_header;
      run(gfa_path, fastq_path, cfg_path, std::cout);
    }
  }
  catch (const GRBException& e) {
    std::cerr << "erro Gurobi [" << e.getErrorCode() << "]: "
              << e.getMessage() << "\n";
    return 1;
  }
  catch (const std::exception& e) {
    std::cerr << "erro: " << e.what() << "\n";
    return 1;
  }

  return 0;
}
