"""
02_leave_out_sweep.py - the main experiment: how does the seed weight alpha
change DIAMOnD's ability to recover hidden disease genes?

Method (the paper's own leave-out test, Fig. 3E-F):
    for each repeat r = 1..N_REPEATS:
        pick a random HIDE_FRACTION of the seed genes and hide them
        for each alpha in ALPHAS:
            run DIAMOnD on the remaining seeds for N_ITERATIONS steps
            recall(i) = fraction of hidden genes among the first i added genes

The same hidden set is used for every alpha within a repeat, so the alphas
are compared on exactly the same task (a paired design).

Output: results/leave_out_recall.csv
        one row per (alpha, repeat, iteration) with the recall at that point,
        plus a row block with alpha = "random" giving the expected recall of
        picking genes at random.

Run from the repository root:   python scripts/02_leave_out_sweep.py
"""
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import config  # noqa: E402
from src.utils import load_diamond, load_network, load_seeds, run_diamond  # noqa: E402

# Each worker process loads the network and the DIAMOnD module once and
# keeps them in these globals (cheaper than sending them with every job).
_G = None
_DIAMOND = None


def _init_worker():
    global _G, _DIAMOND
    _G = load_network()
    _DIAMOND = load_diamond()


def recall_curve(added, hidden):
    """Recall after each iteration.

    In:  added  - genes in the order DIAMOnD added them
         hidden - set of genes that were hidden from DIAMOnD
    Out: array where element i = fraction of hidden genes found among
         the first i+1 added genes
    """
    hits = np.array([gene in hidden for gene in added], dtype=float)
    return np.cumsum(hits) / len(hidden)


def run_job(job):
    """One DIAMOnD run: given kept seeds and alpha, return its recall curve."""
    alpha, repeat, kept, hidden = job
    added = run_diamond(_DIAMOND, _G, kept, config.N_ITERATIONS, alpha)
    recall = recall_curve(added, set(hidden))
    return pd.DataFrame({
        "alpha": alpha,
        "repeat": repeat,
        "iteration": np.arange(1, len(recall) + 1),
        "recall": recall,
    })


def main():
    G = load_network()
    seeds = load_seeds(G)
    n_hide = int(round(config.HIDE_FRACTION * len(seeds)))
    print(f"{len(seeds)} seeds in the network; hiding {n_hide} per repeat, "
          f"{config.N_REPEATS} repeats x {len(config.ALPHAS)} alphas")

    # Draw the hidden sets once, so every alpha sees the same ones.
    rng = np.random.default_rng(config.RANDOM_SEED)
    jobs = []
    for repeat in range(config.N_REPEATS):
        hidden = sorted(rng.choice(seeds, size=n_hide, replace=False))
        kept = sorted(set(seeds) - set(hidden))
        for alpha in config.ALPHAS:
            jobs.append((alpha, repeat, kept, hidden))

    start = time.time()
    with Pool(config.N_JOBS, initializer=_init_worker) as pool:
        frames = []
        for k, frame in enumerate(pool.imap_unordered(run_job, jobs), 1):
            frames.append(frame)
            print(f"  run {k}/{len(jobs)} done ({time.time() - start:.0f} s)")
    results = pd.concat(frames, ignore_index=True)

    # Random baseline: picking i genes at random from the genes that are not
    # kept seeds finds, on average, i * n_hide / (N - n_kept) hidden genes,
    # i.e. recall = i / (N - n_kept).
    n_candidates = G.number_of_nodes() - (len(seeds) - n_hide)
    iters = np.arange(1, config.N_ITERATIONS + 1)
    random_baseline = pd.DataFrame({
        "alpha": "random", "repeat": -1, "iteration": iters,
        "recall": iters / n_candidates,
    })
    results = pd.concat([results, random_baseline], ignore_index=True)

    config.RESULTS_DIR.mkdir(exist_ok=True)
    out = config.RESULTS_DIR / "leave_out_recall.csv"
    results.to_csv(out, index=False)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
