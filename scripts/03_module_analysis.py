"""
03_module_analysis.py - what do the modules built with different alphas
look like? Uses ALL seed genes (nothing hidden).

1. Are the seeds clustered in the network at all? (the paper's premise)
   Size of the largest connected component (LCC) formed by the seeds,
   compared with N_RANDOM_SETS random gene sets of the same size -> z-score.
   (Paper, Fig. 1C-D: z > 1.6 counted as significant.)

2. For each alpha, run DIAMOnD for N_ITERATIONS genes and measure:
   - degree of the added genes          -> does a high alpha let hubs in?
   - fraction that touch a seed directly -> does a high alpha stay local?
   - fraction of seeds joined into one connected module as genes are
     added (paper, Fig. 4F)
   - overlap (Jaccard index) between the gene lists of different alphas

Outputs (results/):
    seed_clustering.csv     LCC size, random mean/sd, z-score
    module_genes.csv        every added gene per alpha, with rank and degree
    module_summary.csv      one row per alpha with the measures above
    seed_integration.csv    fraction of seeds in the module LCC vs iteration
    alpha_overlap.csv       Jaccard index between alphas

Run from the repository root:   python scripts/03_module_analysis.py
"""
import sys
from itertools import combinations
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import config  # noqa: E402
from src.utils import load_diamond, load_network, load_seeds, run_diamond  # noqa: E402


def lcc_size(G, nodes):
    """Number of nodes in the largest connected piece of the subgraph on `nodes`."""
    sub = G.subgraph(nodes)
    if sub.number_of_nodes() == 0:
        return 0
    return max(len(c) for c in nx.connected_components(sub))


def seed_clustering(G, seeds, rng):
    """LCC of the seeds vs random gene sets of the same size (z-score)."""
    observed = lcc_size(G, seeds)
    all_nodes = np.array(list(G.nodes()))
    random_sizes = [lcc_size(G, rng.choice(all_nodes, len(seeds), replace=False))
                    for _ in range(config.N_RANDOM_SETS)]
    mu, sd = np.mean(random_sizes), np.std(random_sizes)
    # empirical p-value: how often a random set is at least as clustered
    p_emp = (np.sum(np.array(random_sizes) >= observed) + 1) / (len(random_sizes) + 1)
    return pd.DataFrame([{
        "n_seeds": len(seeds), "seed_lcc": observed,
        "random_lcc_mean": mu, "random_lcc_sd": sd,
        "z_score": (observed - mu) / sd if sd > 0 else np.nan,
        "empirical_p": p_emp,
    }])


def main():
    G = load_network()
    seeds = load_seeds(G)
    diamond = load_diamond()
    rng = np.random.default_rng(config.RANDOM_SEED)
    config.RESULTS_DIR.mkdir(exist_ok=True)

    # ---- 1. seed clustering --------------------------------------------
    clustering = seed_clustering(G, seeds, rng)
    clustering.to_csv(config.RESULTS_DIR / "seed_clustering.csv", index=False)
    print(clustering.round(3).to_string(index=False))

    # genes that interact directly with at least one seed
    seed_neighbours = set().union(*(set(G[s]) for s in seeds)) - set(seeds)
    degree = dict(G.degree())

    gene_rows, summary_rows, integration_rows, lists = [], [], [], {}
    for alpha in config.ALPHAS:
        print(f"alpha = {alpha}: running DIAMOnD with all {len(seeds)} seeds")
        added = run_diamond(diamond, G, seeds, config.N_ITERATIONS, alpha)
        lists[alpha] = added

        # ---- per-gene table ---------------------------------------------
        for rank, gene in enumerate(added, 1):
            gene_rows.append({"alpha": alpha, "rank": rank, "gene": gene,
                              "degree": degree[gene],
                              "touches_seed": gene in seed_neighbours})

        # ---- how many seeds end up connected, as genes are added ----------
        for i in [0, 10, 25, 50, 100, 150, 200]:
            if i > len(added):
                continue
            module = set(seeds) | set(added[:i])
            biggest = max(nx.connected_components(G.subgraph(module)), key=len)
            integration_rows.append({
                "alpha": alpha, "iteration": i,
                "fraction_seeds_in_lcc": len(biggest & set(seeds)) / len(seeds),
            })

        # ---- summary ------------------------------------------------------
        degs = np.array([degree[g] for g in added])
        summary_rows.append({
            "alpha": alpha,
            "mean_degree_added": degs.mean(),
            "median_degree_added": np.median(degs),
            "fraction_touching_seed": np.mean([g in seed_neighbours for g in added]),
        })

    # reference values to compare the degrees against
    net_deg = np.mean(list(degree.values()))
    seed_deg = np.mean([degree[s] for s in seeds])
    summary = pd.DataFrame(summary_rows)
    summary["mean_degree_network"] = net_deg
    summary["mean_degree_seeds"] = seed_deg

    # ---- overlap between the alpha gene lists -------------------------------
    overlap_rows = []
    for a, b in combinations(config.ALPHAS, 2):
        A, B = set(lists[a]), set(lists[b])
        overlap_rows.append({"alpha_a": a, "alpha_b": b,
                             "jaccard": len(A & B) / len(A | B)})

    pd.DataFrame(gene_rows).to_csv(config.RESULTS_DIR / "module_genes.csv", index=False)
    summary.to_csv(config.RESULTS_DIR / "module_summary.csv", index=False)
    pd.DataFrame(integration_rows).to_csv(
        config.RESULTS_DIR / "seed_integration.csv", index=False)
    pd.DataFrame(overlap_rows).to_csv(config.RESULTS_DIR / "alpha_overlap.csv", index=False)
    print(summary.round(2).to_string(index=False))
    print("wrote module_genes.csv, module_summary.csv, seed_integration.csv, "
          "alpha_overlap.csv")


if __name__ == "__main__":
    main()
