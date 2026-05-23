#pragma once

#include "solver_interface.h"
#include "../io/config.h"
#include <optional>
#include <string>

namespace seq {

// Resolve o PLI no grafo de alinhamento.
//   sol_save_path   : se nao vazio, escreve o .sol apos optimize().
//   warm_start_path : se nao vazio e o arquivo existir, carrega o .sol
//                     como MIP start antes de optimize() (atributo Start
//                     das variaveis). O Gurobi continua otimizando a
//                     partir dessa solucao.
//
// O mapeamento variavel -> arco depende da ordem de iteracao do
// AlignmentGraph; por isso o .sol so pode ser reusado se o ag for
// reconstruido com os mesmos parametros (gfa, fastq, custos).
// node_limit:
//   < 0  : sem limite (resolucao completa).
//   == 0 : Gurobi processa apenas o no raiz (presolve + LP + cortes +
//          heuristicas de raiz) e para. Util para coletar bounds do no 0.
//   > 0  : limite explicito de nos do branch-and-bound.
//
// Retorna nullopt apenas em erro/inviabilidade. Quando ha apenas bounds
// (sem incumbente, ex.: node_limit=0 sem heuristica achar caminho), o
// SolverResult vem com cost=-1, ag_walk vazio e has_path=false.
std::optional<SolverResult>
solve_pli(const AlignmentGraph&   ag,
          const ExperimentConfig& cfg,
          const std::string&      sol_save_path   = "",
          const std::string&      warm_start_path = "",
          const std::string&      lp_save_path    = "",
          long long               node_limit      = -1);

} // namespace seq
