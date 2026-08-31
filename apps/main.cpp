#include "core/sequence.h"
#include "core/sequence_graph.h"
#include "core/alignment_graph.h"
#include "algorithms/edit_distance.h"
#include "algorithms/alignment_walk.h"
#include "algorithms/alignment_potentials.h"
#include "algorithms/warm_start.h"
#include "heuristics/greedy_degree.h"
#include "heuristics/beam_search.h"
#include "heuristics/brkga_search.h"
#include "solvers/lagrangian.h"
#include "io/gfa.h"
#include "io/fastq.h"
#include "io/config.h"
#include "io/gurobi_sol.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <chrono>
#include <cmath>
#include <string>
#include <filesystem>

namespace fs = std::filesystem;
using namespace seq;

// ============================================================================
//  CSV header (mesmo formato do benchmark original)
// ============================================================================

// Schema unificado com sgap-solver. Colunas que a heuristica nao
// produz ficam vazias (lb, node_limit). status eh sempre "found" ou
// "no_path".
static const char* csv_header =
  "graph,sequence,vertices,arcs,query_len,"
  "algorithm,variant,ub,lb,induced_len,is_path,time_ms,"
  "node_limit,status\n";

// Rotulo do status do Lagrangeano (codigos proprios, sem Gurobi).
static const char*
lagr_status_label(int s)
{
  switch (s) {
    case LAGR_OPTIMAL:         return "optimal";
    case LAGR_TIME_LIMIT:      return "time_limit";
    case LAGR_ITERATION_LIMIT: return "iter_limit";
    default:                   return "other";
  }
}

// ============================================================================
//  Solucao de um algoritmo (vertices + metricas)
// ============================================================================

struct AlgorithmSolution {
  std::string name;
  std::vector<VertexId> vertices;
  Sequence induced;
  int edit_dist = -1;
  bool is_path  = false;
};

// ============================================================================
//  Escrita do arquivo .out (relatorio humano-legivel da heuristica)
// ============================================================================

static void
write_out(const fs::path& out_path,
          const std::string& graph_name,
          const std::string& seq_name,
          const std::string& config_name,
          const SequenceGraph& g,
          const Sequence& query,
          const std::vector<AlgorithmSolution>& solutions)
{
  std::ofstream out(out_path);

  out << "# graph: " << graph_name << "\n"
      << "# sequence: " << seq_name << "\n"
      << "# config: " << config_name << "\n"
      << "# vertices: " << g.vertex_count() << "\n"
      << "# arcs: " << g.arc_count() << "\n"
      << "# query_len: " << query.size() << "\n"
      << "#\n"
      << "# query\n"
      << query.str() << "\n\n";

  for (const auto& s : solutions) {
    if (s.edit_dist < 0) continue;

    out << "## " << s.name << "\n"
        << "# edit_dist: " << s.edit_dist << "\n"
        << "# induced_len: " << s.induced.size() << "\n"
        << "# is_path: " << (s.is_path ? "true" : "false") << "\n"
        << "# induced\n"
        << s.induced.str() << "\n"
        << "# walk (" << s.vertices.size() << " vertices)\n";

    for (std::size_t i = 0; i < s.vertices.size(); ++i) {
      if (i > 0) out << " ";
      out << s.vertices[i];
    }
    out << "\n\n";
  }

  out.close();
}

// ============================================================================
//  Execucao unica: 1 grafo + 1 sequencia + config
// ============================================================================

static void
run(const fs::path& gfa_path,
    const fs::path& fastq_path,
    const fs::path& cfg_path,
    const ExperimentConfig& cfg,
    std::ostream& out)
{
  SequenceGraph g = read_gfa(gfa_path.string());
  Sequence query  = read_fastq(fastq_path.string());

  std::string graph_name  = gfa_path.stem().string();
  std::string seq_name    = fastq_path.stem().string();
  std::string config_name = cfg_path.stem().string();

  std::cerr << graph_name << " x " << seq_name << "  ";

  AlgorithmSolution sol_result;
  std::string algorithm = cfg.algorithm;
  std::string variant;
  int    edit_dist = -1;
  int    induced_len = 0;
  bool   is_path = false;
  double time_ms = 0.0;
  double lb_val  = std::nan("");   // preenchido pelo Lagrangeano e pelo jain
  int    lagr_status = -1;         // status bruto do Lagrangeano
  bool   have_lagr_status = false; // true => usa lagr_status_label no CSV
  bool   report_ub = true;         // false => custo nao eh UB valido (relaxacao)

  if (algorithm == "jain") {
    variant = cfg.jain_method;
    AlignmentGraph ag(g, query, cfg.costs);

    ShortestPathMethod method = (cfg.jain_method == "dijkstra")
        ? ShortestPathMethod::Dijkstra
        : ShortestPathMethod::BFS01;

    auto t0 = std::chrono::steady_clock::now();
    auto result = walk_alignment(ag, method);
    auto t1 = std::chrono::steady_clock::now();
    time_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

    if (result) {
      auto induced = g.induced_sequence(result->vertices);
      edit_dist   = result->cost;
      induced_len = static_cast<int>(induced.size());
      is_path     = g.is_valid_path(result->vertices);

      // O Jain resolve a relaxacao do problema (walk, sem a restricao de
      // caminho simples): o custo eh sempre um LB. So vira UB quando o
      // walk encontrado ja eh um caminho simples, i.e., solucao viavel.
      lb_val     = static_cast<double>(result->cost);
      report_ub  = is_path;

      sol_result = {"jain_" + variant, result->vertices,
                    std::move(induced), edit_dist, is_path};
    }

  } else if (algorithm == "greedy") {
    variant = cfg.greedy_selection;
    GreedyDegreeParams params;
    params.selection = (cfg.greedy_selection == "roulette")
        ? GreedyDegreeParams::Selection::roulette
        : GreedyDegreeParams::Selection::greedy;
    params.seed = cfg.greedy_seed;

    auto t0 = std::chrono::steady_clock::now();
    auto best = greedy_degree_best(g, query, cfg.greedy_restarts, params);
    auto t1 = std::chrono::steady_clock::now();
    time_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

    if (best) {
      auto induced = g.induced_sequence(best->vertices);
      edit_dist   = best->cost;
      induced_len = static_cast<int>(induced.size());
      is_path     = g.is_valid_path(best->vertices);

      sol_result = {"greedy_" + variant, best->vertices,
                    std::move(induced), edit_dist, is_path};
    }

  } else if (algorithm == "beam") {
    variant = cfg.beam_potentials;
    AlignmentGraph ag(g, query, cfg.costs);

    ShortestPathMethod pot_method = (cfg.beam_potentials == "dijkstra")
        ? ShortestPathMethod::Dijkstra
        : ShortestPathMethod::BFS01;

    // Resolver beam_width conforme o modo escolhido
    int resolved_beam_width = cfg.beam_width;
    {
      std::size_t n = g.vertex_count();
      std::size_t a = g.arc_count();
      const auto& mode = cfg.beam_width_mode;
      if (mode == "n") {
        resolved_beam_width = static_cast<int>(n);      
      } else if (mode == "density_n") {
        double density = (n > 0) ? static_cast<double>(a) / n : 1.0;
        resolved_beam_width = static_cast<int>(density * n);
      }
      // mode == "constant" => usa cfg.beam_width como esta
    }
    variant += "_" + cfg.beam_width_mode + "_k" + std::to_string(resolved_beam_width);
    std::cerr << "  beam_width_mode=" << cfg.beam_width_mode
              << " resolved_k=" << resolved_beam_width << "  ";

    BeamSearchParams bp;
    bp.beam_width = resolved_beam_width;

    auto t0 = std::chrono::steady_clock::now();
    AlignmentPotentials pot(ag, pot_method);
    auto result = beam_search_alignment(g, query, pot, ag, bp);
    auto t1 = std::chrono::steady_clock::now();
    time_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

    if (result) {
      auto induced = g.induced_sequence(result->vertices);
      edit_dist   = result->cost;
      induced_len = static_cast<int>(induced.size());
      is_path     = g.is_valid_path(result->vertices);

      sol_result = {"beam_" + variant, result->vertices,
                    std::move(induced), edit_dist, is_path};
    }

  } else if (algorithm == "brkga") {
    BrkgaParams bp;
    bp.pop_size    = static_cast<unsigned>(cfg.brkga_pop_size);
    bp.pe          = cfg.brkga_pe;
    bp.pm          = cfg.brkga_pm;
    bp.rhoe        = cfg.brkga_rhoe;
    bp.generations = static_cast<unsigned>(cfg.brkga_generations);
    bp.seed        = static_cast<unsigned long>(cfg.brkga_seed);

    variant = "vp_p" + std::to_string(bp.pop_size)
            + "_g" + std::to_string(bp.generations)
            + "_s" + std::to_string(bp.seed);

    auto t0 = std::chrono::steady_clock::now();
    auto result = brkga_alignment(g, query, bp);
    auto t1 = std::chrono::steady_clock::now();
    time_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

    if (result) {
      auto induced = g.induced_sequence(result->vertices);
      edit_dist   = result->cost;
      induced_len = static_cast<int>(induced.size());
      is_path     = g.is_valid_path(result->vertices);

      sol_result = {"brkga_" + variant, result->vertices,
                    std::move(induced), edit_dist, is_path};
    }

  } else if (algorithm == "lagrangian") {
    variant = "lagrangian";
    AlignmentGraph ag(g, query, cfg.costs);

    // .sol de saida: passado ao solver, que grava seu proprio caminho
    // exato (walk_to_assignments) como warm start compativel com o PLI.
    std::string sol_path;
    if (!cfg.gurobi_sol_dir.empty()) {
      fs::create_directories(cfg.gurobi_sol_dir);
      sol_path = (fs::path(cfg.gurobi_sol_dir) /
                  (graph_name + "_" + seq_name + "_" + config_name + ".sol")).string();
    }

    // Warm start opcional (ex.: .sol de outra heuristica) como UB inicial.
    std::string warm_path;
    if (!cfg.pli_warm_start_path.empty() && fs::exists(cfg.pli_warm_start_path))
      warm_path = cfg.pli_warm_start_path;

    auto t0 = std::chrono::steady_clock::now();
    auto result = solve_lagrangian(ag, cfg, sol_path, warm_path);
    auto t1 = std::chrono::steady_clock::now();
    time_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

    if (result) {
      lb_val           = result->final_lb;   // LB dual (unico algoritmo que produz)
      lagr_status      = result->status;
      have_lagr_status = true;
      if (result->has_path) {
        auto walk    = ag.to_graph_walk(result->ag_walk);
        auto induced = g.induced_sequence(walk);
        edit_dist    = result->cost;
        induced_len  = static_cast<int>(induced.size());
        is_path      = g.is_valid_path(walk);
        sol_result   = {"lagrangian", walk, std::move(induced), edit_dist, is_path};
      }
    }

  } else {
    std::cerr << "algoritmo desconhecido: " << algorithm << "\n";
    return;
  }

  // -- Escrever linha CSV --
  // Para as heuristicas edit_dist eh um UB, vazio quando o algoritmo
  // falhou em encontrar caminho. O jain resolve a relaxacao (walk):
  // seu custo vai para lb e so aparece em ub se o walk for caminho
  // simples. node_limit fica sempre vazio para manter o schema
  // unificado com o sgap-solver.
  out << graph_name << ","
      << seq_name << ","
      << g.vertex_count() << ","
      << g.arc_count() << ","
      << query.size() << ","
      << algorithm << ","
      << variant << ",";
  if (edit_dist >= 0 && report_ub) out << edit_dist;   // ub (so se viavel)
  out << ",";
  if (std::isfinite(lb_val))               // lb: Lagrangeano e jain (relaxacao)
    out << std::fixed << std::setprecision(6) << lb_val << std::defaultfloat;
  out << ","
      << induced_len << ","
      << (is_path ? 1 : 0) << ","
      << std::fixed << std::setprecision(2) << time_ms
      << std::defaultfloat << ","
      << ","                               // node_limit (vazio)
      << (have_lagr_status ? lagr_status_label(lagr_status)
                           : (edit_dist >= 0 ? "found" : "no_path")) << "\n";

  // -- Escrever arquivo .sol (warm start no formato Gurobi) --
  // Requer: solucao valida + walk eh path (PLI tem restricao de
  // caminho simples: cada vertice de G visitado no maximo uma vez).
  // O Lagrangeano ja gravou seu proprio .sol (caminho exato) dentro de
  // solve_lagrangian; nao reconstruir aqui via build_warm_start_arcs.
  if (algorithm != "lagrangian" &&
      !cfg.gurobi_sol_dir.empty() && edit_dist >= 0 && sol_result.is_path) {
    fs::path sol_path = fs::path(cfg.gurobi_sol_dir) /
                        (graph_name + "_" + seq_name + "_" + config_name + ".sol");

    AlignmentGraph ag(g, query, cfg.costs);
    auto warm = build_warm_start_arcs(ag, g, query, sol_result.vertices, cfg.costs);

    if (warm) {
      std::vector<std::pair<std::string, double>> assignments;
      assignments.reserve(warm->arcs.size());
      for (const auto& [u, j] : warm->arcs)
        assignments.emplace_back(arc_var_name(u, j), 1.0);

      std::vector<std::string> header = {
        "warm start gerado por sgap-heuristic",
        "graph: "           + graph_name,
        "sequence: "        + seq_name,
        "config: "          + config_name,
        "algorithm: "       + sol_result.name,
        "edit_dist: "       + std::to_string(edit_dist),
        "warm_start_cost: " + std::to_string(warm->cost),
      };
      write_gurobi_sol(sol_path.string(), assignments, header);
    } else {
      std::cerr << "[aviso] nao foi possivel construir warm start para "
                << sol_path.filename().string() << "\n";
    }
  }

  // -- Escrever arquivo .out (relatorio humano-legivel) --
  // Mesmo padrao de .sol e .lp: so grava se o diretorio foi
  // explicitamente configurado.
  if (!cfg.heuristic_out_dir.empty() && edit_dist >= 0) {
    fs::path out_dir = cfg.heuristic_out_dir;
    fs::create_directories(out_dir);
    fs::path out_path = out_dir / (graph_name + "_" + seq_name + "_" + config_name + ".out");
    std::vector<AlgorithmSolution> solutions = {std::move(sol_result)};
    write_out(out_path, graph_name, seq_name, config_name, g, query, solutions);
  }

  std::cerr << algorithm << "(" << variant << ") edit=" << edit_dist << "\n";
}

// ============================================================================
//  main
// ============================================================================

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
    ExperimentConfig cfg = read_config(cfg_path.string());

    if (argc > 4) {
      fs::path csv_path = argv[4];

      // Criar diretorio de saida se necessario
      if (csv_path.has_parent_path())
        fs::create_directories(csv_path.parent_path());

      bool file_exists = fs::exists(csv_path);

      std::ofstream csv(csv_path, std::ios::app);
      if (!file_exists)
        csv << csv_header;

      run(gfa_path, fastq_path, cfg_path, cfg, csv);
      csv.close();
      std::cerr << "resultado escrito em " << csv_path << "\n";
    } else {
      std::cout << csv_header;
      run(gfa_path, fastq_path, cfg_path, cfg, std::cout);
    }
  } catch (const std::exception& e) {
    std::cerr << "erro: " << e.what() << "\n";
    return 1;
  }

  return 0;
}
