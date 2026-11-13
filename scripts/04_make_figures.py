"""
04_make_figures.py - turn the CSVs from scripts 02 and 03 into figures and
one summary table.

Outputs:
    results/figures/fig1_recall_curves.png    recall vs iteration, one line per alpha
    results/figures/fig2_recall_vs_alpha.png  recall at 50 and 200 genes, per repeat
    results/figures/fig3_degree.png           degree of the added genes per alpha
    results/figures/fig4_locality.png         share of added genes touching a seed
    results/figures/fig5_seed_integration.png share of seeds joined into one module
    results/recall_summary.csv                numbers behind figs 1-2

Run from the repository root:   python scripts/04_make_figures.py
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")                 # draw to files, no window needed
import matplotlib.pyplot as plt       # noqa: E402
import numpy as np                    # noqa: E402
import pandas as pd                   # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import config  # noqa: E402

# One fixed colour per alpha (colour-blind-checked categorical palette),
# always in the same order so alpha=1 is blue in every figure.
PALETTE = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]
COLOURS = {a: PALETTE[i % len(PALETTE)] for i, a in enumerate(config.ALPHAS)}
GREY = "#8a8985"          # random baseline / reference lines
INK = "#52514e"           # axis text

plt.rcParams.update({
    "figure.dpi": 150, "font.size": 9, "axes.edgecolor": "#c9c8c3",
    "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK,
    "axes.grid": True, "grid.color": "#ecebe7", "grid.linewidth": 0.6,
    "axes.spines.top": False, "axes.spines.right": False,
    "lines.linewidth": 2,
})
CHECKPOINTS = [25, 50, 100, 200]


def save(fig, name):
    config.FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(config.FIG_DIR / name, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote figures/{name}")


def recall_figures():
    df = pd.read_csv(config.RESULTS_DIR / "leave_out_recall.csv",
                     dtype={"alpha": str})
    runs = df[df["alpha"] != "random"].copy()
    runs["alpha"] = runs["alpha"].astype(int)
    rand = df[df["alpha"] == "random"]

    # mean and standard deviation over the repeats, per alpha and iteration
    stats = (runs.groupby(["alpha", "iteration"])["recall"]
             .agg(["mean", "std"]).reset_index())

    # ---- Fig 1: recall curves -------------------------------------------
    fig, ax = plt.subplots(figsize=(6.5, 4))
    for alpha in config.ALPHAS:
        s = stats[stats["alpha"] == alpha]
        ax.fill_between(s["iteration"], s["mean"] - s["std"], s["mean"] + s["std"],
                        color=COLOURS[alpha], alpha=0.12, linewidth=0)
        ax.plot(s["iteration"], s["mean"], color=COLOURS[alpha], label=f"α = {alpha}")
    ax.plot(rand["iteration"], rand["recall"], color=GREY, linestyle="--",
            linewidth=1.5, label="random genes")
    ax.set_xlabel("DIAMOnD iteration (number of genes added)")
    ax.set_ylabel("Recall of hidden seed genes")
    ax.set_title(f"Leave-out recall, {int(config.HIDE_FRACTION * 100)}% of seeds "
                 f"hidden (mean ± SD of {config.N_REPEATS} repeats)",
                 fontsize=9, loc="left")
    ax.set_ylim(bottom=0)
    # legend outside the plot so it never covers a curve
    ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1))
    save(fig, "fig1_recall_curves.png")

    # ---- Fig 2: recall at 50 and 200 genes, every repeat shown --------------
    fig, axes = plt.subplots(1, 2, figsize=(7, 3.2), sharey=True)
    for ax, point in zip(axes, [50, 200]):
        at = runs[runs["iteration"] == min(point, runs["iteration"].max())]
        for i, alpha in enumerate(config.ALPHAS):
            vals = at[at["alpha"] == alpha]["recall"]
            jitter = np.random.default_rng(i).uniform(-0.12, 0.12, len(vals))
            ax.scatter(np.full(len(vals), i) + jitter, vals, s=22,
                       color=COLOURS[alpha], edgecolor="white", linewidth=0.8, zorder=3)
            ax.hlines(vals.mean(), i - 0.25, i + 0.25, color=INK, linewidth=2, zorder=4)
        ax.set_xticks(range(len(config.ALPHAS)), [f"α={a}" for a in config.ALPHAS])
        ax.set_title(f"after {point} genes", fontsize=9, loc="left")
    axes[0].set_ylabel("Recall (dots = repeats, bar = mean)")
    save(fig, "fig2_recall_vs_alpha.png")

    # ---- summary table ------------------------------------------------------
    rows = []
    for alpha in config.ALPHAS:
        s = stats[stats["alpha"] == alpha].set_index("iteration")
        final = s["mean"].iloc[-1]
        # plateau = first iteration reaching 95% of the final mean recall
        plateau = int(s.index[s["mean"] >= 0.95 * final][0]) if final > 0 else np.nan
        row = {"alpha": alpha, "plateau_iteration": plateau}
        for c in CHECKPOINTS:
            if c in s.index:
                row[f"recall@{c}_mean"] = round(s.loc[c, "mean"], 3)
                row[f"recall@{c}_sd"] = round(s.loc[c, "std"], 3)
        rows.append(row)
    rand_row = {"alpha": "random", "plateau_iteration": np.nan}
    for c in CHECKPOINTS:
        if c <= rand["iteration"].max():
            rand_row[f"recall@{c}_mean"] = round(
                rand.loc[rand["iteration"] == c, "recall"].iloc[0], 4)
    rows.append(rand_row)
    table = pd.DataFrame(rows)
    table.to_csv(config.RESULTS_DIR / "recall_summary.csv", index=False)
    print(table.to_string(index=False))


def module_figures():
    genes = pd.read_csv(config.RESULTS_DIR / "module_genes.csv")
    summary = pd.read_csv(config.RESULTS_DIR / "module_summary.csv")
    integ = pd.read_csv(config.RESULTS_DIR / "seed_integration.csv")

    # ---- Fig 3: degree of added genes --------------------------------------
    fig, ax = plt.subplots(figsize=(6, 3.4))
    data = [genes.loc[genes["alpha"] == a, "degree"] for a in config.ALPHAS]
    bp = ax.boxplot(data, widths=0.5, patch_artist=True, showfliers=True,
                    medianprops={"color": INK, "linewidth": 2},
                    flierprops={"markersize": 3, "markeredgecolor": GREY})
    for patch, alpha in zip(bp["boxes"], config.ALPHAS):
        patch.set_facecolor(COLOURS[alpha])
        patch.set_alpha(0.55)
        patch.set_edgecolor(COLOURS[alpha])
    ax.set_xticks(range(1, len(config.ALPHAS) + 1), [f"α={a}" for a in config.ALPHAS])
    ax.set_yscale("log")
    # reference lines: average degree of all network genes and of the seeds
    net = summary["mean_degree_network"].iloc[0]
    sd = summary["mean_degree_seeds"].iloc[0]
    ax.axhline(net, color=GREY, linestyle=":", linewidth=1.3,
               label=f"network mean degree ({net:.1f})")
    ax.axhline(sd, color=GREY, linestyle="--", linewidth=1.3,
               label=f"seed mean degree ({sd:.1f})")
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    ax.set_ylabel("Degree of added gene (log scale)")
    ax.set_title(f"Connections of the first {config.N_ITERATIONS} DIAMOnD genes",
                 fontsize=9, loc="left")
    save(fig, "fig3_degree.png")

    # ---- Fig 4: how local is the module? -------------------------------------
    fig, ax = plt.subplots(figsize=(5, 3))
    x = range(len(config.ALPHAS))
    vals = summary.set_index("alpha").loc[config.ALPHAS, "fraction_touching_seed"]
    ax.bar(x, vals, width=0.6, color=[COLOURS[a] for a in config.ALPHAS])
    for i, v in zip(x, vals):
        ax.annotate(f"{v:.0%}", (i, v), xytext=(0, 3), textcoords="offset points",
                    ha="center", fontsize=8, color=INK)
    ax.set_xticks(list(x), [f"α={a}" for a in config.ALPHAS])
    ax.set_ylim(0, 1.08)
    ax.set_ylabel("Share of added genes that\ninteract directly with a seed")
    ax.grid(axis="x", visible=False)
    save(fig, "fig4_locality.png")

    # ---- Fig 5: seeds joined into one module ---------------------------------
    fig, ax = plt.subplots(figsize=(6, 3.4))
    for alpha in config.ALPHAS:
        s = integ[integ["alpha"] == alpha]
        ax.plot(s["iteration"], s["fraction_seeds_in_lcc"], marker="o",
                markersize=4, color=COLOURS[alpha], label=f"α = {alpha}")
    ax.set_xlabel("DIAMOnD iteration")
    ax.set_ylabel("Share of seeds in the\nlargest connected module")
    ax.set_ylim(0, 1.02)
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    save(fig, "fig5_seed_integration.png")


def fmt(mean, sd=None):
    """Format a number (and optional SD) for the README table."""
    if pd.isna(mean):
        return "–"
    return f"{mean:.2f} ± {sd:.2f}" if sd is not None and not pd.isna(sd) else f"{mean:.3f}"


def write_readme_results():
    """Fill the block between the RESULTS markers in README.md with the
    numbers and figures from this run, so nothing has to be copied by hand."""
    recall = pd.read_csv(config.RESULTS_DIR / "recall_summary.csv")
    clust = pd.read_csv(config.RESULTS_DIR / "seed_clustering.csv").iloc[0]
    summ = pd.read_csv(config.RESULTS_DIR / "module_summary.csv")
    integ = pd.read_csv(config.RESULTS_DIR / "seed_integration.csv")
    pct = int(config.HIDE_FRACTION * 100)

    lines = [
        f"Disease: `{config.DISEASE_ID}` · network: {config.NETWORK_LABEL} · "
        f"{int(clust['n_seeds'])} seed genes in the network · "
        f"{pct}% hidden per repeat · {config.N_REPEATS} repeats · "
        f"{config.N_ITERATIONS} DIAMOnD genes per run",
        "",
        "**Are the seeds clustered?** The seeds' largest connected component has "
        f"{int(clust['seed_lcc'])} genes, vs {clust['random_lcc_mean']:.1f} ± "
        f"{clust['random_lcc_sd']:.1f} for {config.N_RANDOM_SETS} random gene sets "
        f"(z = {clust['z_score']:.1f}; the paper counts z > 1.6 as significant).",
        "",
        "**Recall of hidden seed genes** (mean ± SD over repeats)",
        "",
        "| α | recall@50 | recall@200 | plateau iteration |",
        "|---|---|---|---|",
    ]
    for _, r in recall.iterrows():
        plateau = "–" if pd.isna(r["plateau_iteration"]) else int(r["plateau_iteration"])
        lines.append(f"| {r['alpha']} | {fmt(r['recall@50_mean'], r.get('recall@50_sd'))} "
                     f"| {fmt(r['recall@200_mean'], r.get('recall@200_sd'))} | {plateau} |")

    first, last = integ["iteration"].min(), integ["iteration"].max()
    lines += [
        "",
        f"**What the {config.N_ITERATIONS} added genes look like** (all seeds used). "
        f"Mean degree: network {summ['mean_degree_network'].iloc[0]:.1f}, "
        f"seeds {summ['mean_degree_seeds'].iloc[0]:.1f}.",
        "",
        "| α | median degree | mean degree | touch a seed directly | "
        f"seeds in one module (iter {first} → {last}) |",
        "|---|---|---|---|---|",
    ]
    for _, r in summ.iterrows():
        g = integ[integ["alpha"] == r["alpha"]].set_index("iteration")["fraction_seeds_in_lcc"]
        lines.append(f"| {int(r['alpha'])} | {r['median_degree_added']:.0f} | "
                     f"{r['mean_degree_added']:.1f} | {r['fraction_touching_seed']:.0%} | "
                     f"{g.loc[first]:.0%} → {g.loc[last]:.0%} |")

    figs = f"results/{config.NETWORK}/figures"
    lines += ["", f"![Recall curves]({figs}/fig1_recall_curves.png)",
              f"![Recall at 50 and 200 genes]({figs}/fig2_recall_vs_alpha.png)",
              f"![Degree of added genes]({figs}/fig3_degree.png)",
              f"![Share touching a seed]({figs}/fig4_locality.png)",
              f"![Seeds joined into one module]({figs}/fig5_seed_integration.png)"]

    readme = config.ROOT / "README.md"
    text = readme.read_text(encoding="utf-8")
    # each network has its own block in the README
    start = f"<!-- RESULTS:{config.NETWORK}:START -->"
    end = f"<!-- RESULTS:{config.NETWORK}:END -->"
    before, rest = text.split(start, 1)
    _, after = rest.split(end, 1)
    readme.write_text(before + start + "\n" + "\n".join(lines) + "\n" + end + after,
                      encoding="utf-8")
    print(f"  filled the {config.NETWORK} Results block in README.md")


if __name__ == "__main__":
    recall_figures()
    module_figures()
    write_readme_results()
