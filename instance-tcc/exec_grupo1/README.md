# Experimentos — grupo 1

Executa **beam**, **brkga**, **jain** e **pli** sobre as instâncias de `../grupo1`
(tamanhos 25s_50v, 50s_100v, 75s_150v, 100s_200v).

## Como rodar

```
powershell -ExecutionPolicy Bypass -File run_grupo1.ps1
powershell -ExecutionPolicy Bypass -File run_grupo1.ps1 -Algoritmos beam,jain
powershell -ExecutionPolicy Bypass -File run_grupo1.ps1 -Algoritmos brkga -SementesBrkga 1..5
powershell -ExecutionPolicy Bypass -File run_grupo1.ps1 -Algoritmos brkga -SementesBrkga 11,12
```

Sem `-Algoritmos`, roda os quatro na ordem beam → brkga → jain → pli
(o PLI, mais demorado, fica por último). Em cada algoritmo, as instâncias
vão da menor para a maior.

Requisitos: executáveis em `../sgap/x64/Release` (`sgap-heuristic.exe` para
beam/brkga/jain e `sgap-solver.exe` para pli). O Gurobi precisa estar
instalado e licenciado para o PLI.

## Instâncias

Em cada pasta de tamanho, a sequência `<N>s_a4_inst<i>.fastq` é alinhada aos
5 grafos `<M>v_a4_d{1.5,2.0,2.5,3.0,3.5}_inst<i>.gfa` (mesmo `i`).
São 4 tamanhos × 20 instâncias × 5 densidades = **400 execuções por algoritmo**
(no BRKGA, 400 × número de sementes: **4000** com as 10 sementes padrão).

## Várias runs do BRKGA (sementes)

O BRKGA é determinístico para uma semente: rodar duas vezes com a mesma
`brkga_seed` dá o mesmo resultado. Para ter várias runs, o script executa
o BRKGA uma vez por semente, **por padrão as sementes 1 a 10** (`-SementesBrkga`).

- Os parâmetros vêm de `config/brkga.ini`. Para cada semente, o script gera
  `config/brkga_sementes/brkga_s<semente>.ini`, uma cópia com só a `brkga_seed` trocada.
  Esses arquivos são recriados a cada execução, então edite **só** o `brkga.ini`.
- A ordem é semente por semente: todas as 400 instâncias com a semente 1, depois a 2, e assim por diante.
- Todas as runs vão para `resultados/brkga.csv`. A semente aparece no fim da coluna
  `variant` (ex.: `vp_p256_g500_s3`). Na análise, agrupe por `graph,sequence` e
  calcule melhor, média e desvio-padrão do `ub` e do `time_ms`.
- Para adicionar runs depois, passe novas sementes (ex.: `-SementesBrkga 11..20`).
  As sementes que já rodaram são puladas.

## Configurações (`config/`)

| Arquivo     | Algoritmo | Parâmetros |
|-------------|-----------|------------|
| `beam.ini`  | Beam Search | `beam_width_mode = density_n` (k = densidade × vértices = número de arcos), potenciais `bfs01` |
| `brkga.ini` | BRKGA | p = 256, pe = 0.20, pm = 0.10, rhoe = 0.70, 500 gerações; sementes 1..10 (via script) |
| `jain.ini`  | Jain (relaxação walk) | `bfs01` |
| `pli.ini`   | PLI (Gurobi) | sem limite de tempo e sem limite de nós |


## Todas as opções de configuração

Cada `.ini` usa o formato `chave = valor`, com `#` para comentários.
Uma chave ausente fica com o valor padrão abaixo (definido em `sgap/io/config.h`).
Uma chave escrita errado é **ignorada sem aviso**.

### Gerais

| Chave | Padrão | Descrição |
|-------|--------|-----------|
| `algorithm` | `jain` | `beam`, `brkga`, `jain`, `greedy`, `lagrangian` (todos no `sgap-heuristic`) ou `pli` (no `sgap-solver`) |
| `cost_substitution` | `1` | custo de substituição (inteiro) |
| `cost_insertion` | `1` | custo de inserção (inteiro) |
| `cost_deletion` | `1` | custo de remoção (inteiro) |

> Se mudar os custos, mude nos **quatro** `.ini`, senão os resultados não são comparáveis.

### Beam Search (`algorithm = beam`)

| Chave | Padrão | Descrição |
|-------|--------|-----------|
| `beam_width` | `1000` | largura k da beam (usada só com `beam_width_mode = constant`) |
| `beam_width_mode` | `constant` | `constant`: usa `beam_width`; `n`: k = número de vértices; `density_n`: k = (arcos/vértices) × vértices, ou seja, o número de arcos |
| `beam_potentials` | `bfs01` | algoritmo dos potenciais (heurística tipo A\*): `bfs01` ou `dijkstra` |

A coluna `variant` do CSV registra o k usado, ex.: `bfs01_constant_k1000`.

### BRKGA (`algorithm = brkga`)

| Chave | Padrão | Descrição |
|-------|--------|-----------|
| `brkga_pop_size` | `256` | tamanho da população (p) |
| `brkga_pe` | `0.20` | fração da população que é elite |
| `brkga_pm` | `0.10` | fração de mutantes por geração |
| `brkga_rhoe` | `0.70` | probabilidade de herdar o gene do pai elite |
| `brkga_generations` | `500` | número de gerações |
| `brkga_seed` | `42` | semente aleatória. **Nos scripts é sobrescrita** por `-SementesBrkga`; não precisa colocar no `brkga.ini` |

A coluna `variant` registra população, gerações e semente, ex.: `vp_p256_g500_s42`.

### Jain (`algorithm = jain`)

| Chave | Padrão | Descrição |
|-------|--------|-----------|
| `jain_method` | `bfs01` | caminho mínimo no grafo de alinhamento: `bfs01` ou `dijkstra` (mesmo resultado; muda só o tempo) |

### PLI (`algorithm = pli`)

| Chave | Padrão | Descrição |
|-------|--------|-----------|
| `pli_time_limit_s` | `0` | limite de tempo do Gurobi em segundos (`0` = sem limite) |
| `pli_node_limit` | `-1` | `NodeLimit` do Gurobi: `<0` = sem limite; `0` = só o nó raiz; `N` = até N nós |
| `pli_verbose` | `false` | `true` mostra o log do Gurobi no console |
| `pli_warm_start_path` | vazio | `.sol` usado como solução inicial (MIP start). É um único arquivo fixo, então serve para testar uma instância, não para o lote todo |

### Arquivos extras de saída (opcionais)

| Chave | Padrão | Descrição |
|-------|--------|-----------|
| `heuristic_out_dir` | vazio | heurísticas: grava um relatório legível `<grafo>_<seq>_<config>.out` com a sequência induzida e o caminho |
| `gurobi_sol_dir` | vazio | grava a solução como `.sol` do Gurobi (heurísticas: só quando a solução é caminho simples; PLI: a solução final) |
| `gurobi_lp_dir` | vazio | PLI: grava o modelo `.lp` antes de otimizar |

> Caminhos relativos nessas chaves partem da pasta de onde você **executa** o
> script, não da pasta do `.ini`. Prefira caminhos absolutos.

### Outros algoritmos disponíveis (não usados pelos scripts)

- **`greedy`**: `greedy_restarts` (padrão `1000`), `greedy_selection` (`roulette` ou `greedy`), `greedy_seed` (`42`).
- **`lagrangian`**: `lagr_max_iters` (`500`), `lagr_time_limit_s` (`0`), `lagr_alpha` (`2.0`),
  `lagr_alpha_patience` (`20`), `lagr_repair_beam` (`64`).

Para usar um deles, crie `config/<nome>.ini` e adicione o nome em `$Binario`
no script, apontando para `sgap-heuristic.exe`.

## Resultados

Um CSV por algoritmo em `resultados/<algoritmo>.csv`; no BRKGA há uma linha por semente. As colunas são:

```
graph,sequence,vertices,arcs,query_len,algorithm,variant,ub,lb,induced_len,is_path,time_ms,node_limit,status
```

- `ub`: custo da solução viável encontrada (heurísticas e PLI).
- `lb`: limite inferior (Jain, que resolve a relaxação, e o PLI).
- `status`: `found`/`no_path` nas heurísticas; `optimal`, `node_limit` etc. no PLI.

## Retomada

Pares grafo × sequência que já estão no CSV são pulados. No BRKGA, a checagem é feita por semente. Dá para interromper
com Ctrl+C e rodar o script de novo que ele continua de onde parou.
Se você mudar um `.ini`, apague o CSV do algoritmo correspondente para
refazer as execuções.
