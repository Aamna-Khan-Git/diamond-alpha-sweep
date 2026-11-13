"""
test_port.py - checks that the Python 3 port computes the same connectivity
p-value as an independent implementation (SciPy's hypergeometric tail).

    p-value = P(X >= kb), X ~ Hypergeometric(N total, s seeds, k draws)

Run from the repository root (after 01_download_data.py):
    python -m pytest tests/      or      python tests/test_port.py
"""
import sys
from pathlib import Path

import numpy as np
from scipy.stats import hypergeom

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.utils import load_diamond  # noqa: E402


def test_pvalue_matches_scipy():
    d = load_diamond()
    N = 1000
    gamma_ln = d.compute_all_gamma_ln(N + 1)
    rng = np.random.default_rng(0)
    for _ in range(200):
        s = int(rng.integers(1, 100))          # seeds in the network
        k = int(rng.integers(1, 60))           # degree of the candidate
        kb = int(rng.integers(1, min(k, s) + 1))  # links to seeds
        ours = float(d.pvalue(kb, k, N, s, gamma_ln))
        ref = hypergeom.sf(kb - 1, N, s, k)    # P(X >= kb)
        assert np.isclose(ours, ref, rtol=1e-6, atol=1e-300), (kb, k, s, ours, ref)


if __name__ == "__main__":
    test_pvalue_matches_scipy()
    print("ported p-value matches scipy.stats.hypergeom")
