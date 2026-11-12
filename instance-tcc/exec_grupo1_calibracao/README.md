# Calibração — grupo 1

Executa **apenas o pli** sobre as instâncias de `../grupo1`
(tamanhos 25s_50v, 50s_100v, 75s_150v, 100s_200v), com a **mesma configuração
do grupo 2**: só o nó raiz (`pli_node_limit = 0`), sem limite de tempo.

## Por que

No grupo 2 o PLI não chega ao ótimo, então não dá para saber diretamente o
quanto o `ub`/`lb` do nó raiz se afastam do ótimo. No grupo 1 o ótimo é
conhecido (`../exec_grupo1/resultados/pli.csv`, sem limite de nós). Rodando o
nó raiz aqui, a comparação dos dois CSVs (por `graph,sequence`) mede:

- **gap do `ub` no raiz** em relação ao ótimo;
- **gap do `lb` no raiz** em relação ao ótimo;
- em quantas instâncias o raiz já resolve (`status = optimal`);
- a fração do tempo do ótimo gasta no raiz.

Isso serve de referência para interpretar os resultados do PLI no grupo 2.

## Como rodar

```
powershell -ExecutionPolicy Bypass -File run_grupo1_calibracao.ps1
```

Requisitos: `../sgap/x64/Release/sgap-solver.exe` e Gurobi instalado e licenciado.
São 4 tamanhos × 20 instâncias × 5 densidades = **400 execuções**.

## Configuração (`config/pli.ini`)

| Chave | Valor | |
|-------|-------|-|
| `pli_node_limit` | `0` | só o nó raiz (igual a `../exec_grupo2/config/pli.ini`) |
| `pli_time_limit_s` | `0` | sem limite de tempo |
| `pli_verbose` | `false` | |

Se mudar o `pli.ini` do grupo 2, mude este também, senão a calibração deixa de valer.
Para as demais opções, veja `../exec_grupo2/README.md`.

## Resultados

`resultados/pli.csv`, com as mesmas colunas dos outros grupos:

```
graph,sequence,vertices,arcs,query_len,algorithm,variant,ub,lb,induced_len,is_path,time_ms,node_limit,status
```

Pares já presentes no CSV são pulados: dá para interromper com Ctrl+C e
retomar. Se mudar o `pli.ini`, apague o CSV para refazer.
