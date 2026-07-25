#pragma once

#include "solver_interface.h"
#include "../io/config.h"
#include <optional>
#include <string>

namespace seq {

// Status reportado pelo solver Lagrangeano. Codigos proprios, independentes
// de qualquer solver externo (o metodo nao usa Gurobi). O CSV do
// sgap-heuristic os rotula via lagr_status_label().
enum LagrangianStatus {
  LAGR_OPTIMAL         = 0,  // gap inteiro fechado ou ||g||^2 = 0 (otimo provado)
  LAGR_TIME_LIMIT      = 1,  // parou por lagr_time_limit_s
  LAGR_ITERATION_LIMIT = 2,  // parou por lagr_max_iters
};

// Relaxacao Lagrangeana do SGAP com restricao de caminho simples.
//
// A unica familia de restricoes que torna o problema NP-dificil eh a de
// caminho simples (in_v(v) <= 1 para cada vertice v do grafo de
// sequencia). Dualizando essa familia com multiplicadores lambda_v >= 0,
// o subproblema que sobra eh apenas conservacao de fluxo S->T, ou seja,
// um caminho de custo minimo no grafo de alinhamento com custos de arco
// modificados (w + lambda_p). Resolvemos esse subproblema por Dijkstra e
// maximizamos o dual pelo metodo do subgradiente.
//
// Observacao teorica: como o subproblema eh um fluxo de custo minimo (com
// a integrality property), o melhor limite Lagrangeano iguala o da
// relaxacao linear; ele nao supera o LP raiz do PLI. O valor pratico esta
// em obter um limite dual sem resolver o MIP e solucoes primais viaveis
// por reparo.
//
//   sol_save_path   : se nao vazio e houver UB viavel, grava o melhor
//                     caminho como .sol (compativel com o PLI / warm start).
//   warm_start_path : se nao vazio e o arquivo existir, carrega um .sol
//                     como solucao primal inicial (UB de partida).
//
// Retorna nullopt apenas em erro grave. Quando nenhum caminho viavel eh
// encontrado, o SolverResult vem com cost=-1, ag_walk vazio e
// has_path=false, mas ainda reporta o limite (final_lb).
std::optional<SolverResult>
solve_lagrangian(const AlignmentGraph&   ag,
                 const ExperimentConfig& cfg,
                 const std::string&      sol_save_path   = "",
                 const std::string&      warm_start_path = "");

} // namespace seq
