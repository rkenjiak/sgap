import sys
import random
import argparse
import os

# ---------------------------------------------------------------------------
# v2: para cada instancia gera 1 sequencia (FASTQ) e 5 grafos (GFA) com graus
# medios crescentes, tais que G(d1) e subgrafo de G(d2), e assim por diante.
# Os grafos compartilham os mesmos vertices e rotulos; as arestas sao
# acumuladas de forma incremental (o grafo menor e um prefixo do maior).
# ---------------------------------------------------------------------------

def write_sequence(base_name, seq_len, alphabet):
    my_sequence = "".join(random.choice(alphabet) for _ in range(seq_len))
    seq_id = os.path.basename(base_name)

    with open(f"{base_name}.fastq", 'w') as f:
        f.write(f"@{seq_id}\n{my_sequence}\n+\n" + "!" * seq_len + "\n")


def build_nested_edges(nv, edge_targets):
    """Gera uma lista de arestas em ordem de insercao e devolve, para cada
    alvo em edge_targets (crescente), o prefixo da lista com aquele tamanho.
    Assim, o grafo do alvo k e sempre subgrafo do grafo do alvo k+1."""
    edge_list = []
    edge_set = set()
    nodes = list(range(nv))
    random.shuffle(nodes)

    # 1. Garantir conectividade (Spanning Tree = nv - 1 arestas)
    for i in range(1, nv):
        u = nodes[random.randint(0, i - 1)]
        v = nodes[i]
        if random.random() > 0.5: u, v = v, u
        edge_set.add((u, v))
        edge_list.append((u, v))

    # 2. Preencher aleatoriamente ate atingir o maior alvo
    max_possible = nv * (nv - 1)
    final_target = min(edge_targets[-1], max_possible)

    while len(edge_list) < final_target:
        u, v = random.randint(0, nv - 1), random.randint(0, nv - 1)
        if u != v and (u, v) not in edge_set:
            edge_set.add((u, v))
            edge_list.append((u, v))

    # 3. Cada densidade recebe o prefixo correspondente
    return [edge_list[:min(t, max_possible)] for t in edge_targets]


def write_graph(base_name, labels, edges):
    with open(f"{base_name}.gfa", 'w') as f:
        f.write("H\tVN:Z:1.0\n")
        for i, lab in enumerate(labels):
            f.write(f"S\t{i}\t{lab}\n")
        for x, y in edges:
            f.write(f"L\t{x}\t+\t{y}\t+\t0M\n")


def edges_for_degree(d_value, nv):
    # Mesma convencao do graph_gen.py: numero de arcos = d * n (d = grau medio)
    max_possible = nv * (nv - 1)
    if abs(d_value - 1.0) < 1e-9:
        return nv - 1
    return min(int(d_value * nv), max_possible)


def generate_config(outdir, seq_len, nv, degrees, ni, alphabet):
    """Gera ni instancias para uma configuracao (seq_len, nv)."""
    alpha_size = len(alphabet)
    edge_targets = [edges_for_degree(d, nv) for d in degrees]

    cfg_dir = os.path.join(outdir, f"{seq_len}s_{nv}v")
    os.makedirs(cfg_dir, exist_ok=True)

    print(f"  Config: sequencia={seq_len}, vertices={nv}, arcos={edge_targets}  -> {cfg_dir}")

    for i in range(1, ni + 1):
        # Sequencia
        seq_prefix = f"{seq_len}s_a{alpha_size}_inst{i}"
        write_sequence(os.path.join(cfg_dir, seq_prefix), seq_len, alphabet)

        # Grafos aninhados (mesmos rotulos, arestas incrementais)
        labels = [random.choice(alphabet) for _ in range(nv)]
        nested = build_nested_edges(nv, edge_targets)

        for d_value, edges in zip(degrees, nested):
            d_label = f"{d_value:.1f}"
            g_prefix = f"{nv}v_a{alpha_size}_d{d_label}_inst{i}"
            write_graph(os.path.join(cfg_dir, g_prefix), labels, edges)

    print(f"    {ni} sequencias + {ni * len(degrees)} grafos gerados.")
    return ni, ni * len(degrees)


def main():
    parser = argparse.ArgumentParser(
        description="Gerador v2: sequencias FASTQ + grafos GFA aninhados (d1 < d2 < d3) em dois grupos.")

    parser.add_argument("--grupo1", type=int, nargs="+", default=[25, 50, 75, 100],
                        help="Tamanhos de sequencia do grupo 1 (default: 25 50 75 100)")
    parser.add_argument("--grupo2", type=int, nargs="+", default=[125, 250, 375, 500],
                        help="Tamanhos de sequencia do grupo 2 (default: 125 250 375 500)")
    parser.add_argument("--ratio", type=int, default=2,
                        help="Vertices do grafo = ratio * tamanho da sequencia (default: 2)")
    parser.add_argument("-d", "--degrees", type=float, nargs="+", default=[1.5, 2.0, 2.5, 3.0, 3.5],
                        help="Graus medios, em ordem crescente (default: 1.5 2.0 2.5 3.0 3.5)")
    parser.add_argument("--ni", type=int, default=20, help="Instancias por configuracao (default: 20)")

    parser.add_argument("--outdir", type=str, default=".", help="Pasta raiz de saida (default: diretorio atual)")
    parser.add_argument("--seed", type=int, help="Seed randomica")
    parser.add_argument("--alphabet", type=str, default="ACGT", help="Alfabeto (default: ACGT)")

    args = parser.parse_args()

    if args.seed is not None:
        random.seed(args.seed)

    degrees = sorted(args.degrees)
    if degrees != args.degrees:
        print("Aviso: graus reordenados em ordem crescente:", degrees)

    grupos = {
        "grupo1": args.grupo1,
        "grupo2": args.grupo2,
    }

    print("\n" + "=" * 60)
    print("CONFIGURACOES DE GERACAO (v2)")
    print("=" * 60)
    print(f"Grupo 1 (seq)...........: {args.grupo1}")
    print(f"Grupo 2 (seq)...........: {args.grupo2}")
    print(f"Vertices................: {args.ratio} x tamanho da sequencia")
    print(f"Graus medios............: {degrees}")
    print(f"Instancias por config...: {args.ni}")
    print(f"Alfabeto................: {args.alphabet}")
    print(f"Seed....................: {args.seed}")
    print(f"Diretorio de saida......: {os.path.abspath(args.outdir)}")
    print("=" * 60 + "\n")

    total_seqs = total_graphs = 0
    for nome, tamanhos in grupos.items():
        gdir = os.path.join(args.outdir, nome)
        os.makedirs(gdir, exist_ok=True)
        print(f"[{nome}] -> {os.path.abspath(gdir)}")
        for seq_len in tamanhos:
            nv = args.ratio * seq_len
            if nv <= 1:
                print(f"  Erro: numero de vertices ({nv}) deve ser maior que 1.")
                sys.exit(1)
            s, g = generate_config(gdir, seq_len, nv, degrees, args.ni, args.alphabet)
            total_seqs += s
            total_graphs += g
        print("-" * 60)

    print(f"\nConcluido! {total_seqs} sequencias e {total_graphs} grafos gerados.\n")


if __name__ == "__main__":
    main()
