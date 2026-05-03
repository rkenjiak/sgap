#pragma once

#include "../core/types.h"
#include <cstddef>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

namespace seq {

// Nome canonico da variavel binaria associada ao j-esimo arco de saida
// do vertice u no AlignmentGraph. Eh o contrato compartilhado entre o
// PLI (que cria a variavel) e qualquer escritor de .sol que queira
// produzir um warm start compativel.
inline std::string
arc_var_name(VertexId u, std::size_t j)
{
  return "x_" + std::to_string(u) + "_" + std::to_string(j);
}

// Le um arquivo .sol do Gurobi e retorna o mapa nome -> valor.
// Lanca std::runtime_error em caso de falha de IO. Linhas iniciadas
// por '#' e linhas vazias sao ignoradas.
std::unordered_map<std::string, double>
read_gurobi_sol(const std::string& path);

// Escreve um arquivo .sol no formato do Gurobi. Cada par
// (nome, valor) vira uma linha "<nome> <valor>". As linhas em
// header_comment (uma por elemento) sao gravadas como comentarios
// no topo do arquivo, permitindo registrar metadados como custo,
// algoritmo de origem, etc.
//
// Cria os diretorios pais conforme necessario.
void
write_gurobi_sol(
    const std::string& path,
    const std::vector<std::pair<std::string, double>>& assignments,
    const std::vector<std::string>& header_comments = {});

} // namespace seq
