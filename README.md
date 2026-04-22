# SGAP — Sequence Graph Alignment Probleam

Projeto de TCC sobre alinhamento de uma sequência a um grafo de sequências com restrição de **caminho simples** (sem repetição de vértices). Problema NP-difícil; o projeto implementa tanto a relaxação polinomial (walk) quanto o problema exato (path) e heurísticas.

## O que está aqui

- `core/` — `Sequence`, `SequenceGraph`, `AlignmentGraph` (DAG em camadas).
- `algorithms/` — Dijkstra / BFS-0-1 no grafo de alinhamento, potenciais forward/backward, edit distance.
- `heuristics/` — beam search (com potencial backward como heurística A\*-like) e greedy-degree multi-start.
- `solvers/` — PLI exata via Gurobi (fluxo + restrição de simplicidade).
- `io/` — leitura de GFA, FASTQ e arquivo de configuração `key=value`.
- `apps/` — dois binários: `sgap-heuristic` (heurísticas) e `sgap-solver` (exato).
- `instances/`, `scripts/` — entrada e utilitários.

## Build

Visual Studio (solução `sgap.sln`) ou CMake (`CMakeLists.txt`). Dependência externa: Gurobi, usada só pelo `sgap-solver`.

## Uso

```
sgap-heuristic  <grafo.gfa> <query.fastq> <config.cfg> <saida.sol>
sgap-solver     <grafo.gfa> <query.fastq> <config.cfg> <saida.sol>
```

Parâmetros relevantes (`config.cfg`): `algorithm`, `beam_width`, `greedy_restarts`, `pli_time_limit_s`, custos de edição. Ver `io/config.h`.
