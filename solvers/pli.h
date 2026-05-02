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
std::optional<SolverResult>
solve_pli(const AlignmentGraph&   ag,
          const ExperimentConfig& cfg,
          const std::string&      sol_save_path   = "",
          const std::string&      warm_start_path = "",
          const std::string&      lp_save_path    = "");

} // namespace seq
