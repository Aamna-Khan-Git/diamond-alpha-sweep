"""
utils.py - helpers shared by the analysis scripts.

    load_network()   read the chosen network into a networkx graph
    load_seeds()     read the seed genes that are present in the network
    load_diamond()   import the Python 3 port of the authors' DIAMOnD code
    run_diamond()    run DIAMOnD and return the added genes in order
"""
import importlib.util

import networkx as nx
import pandas as pd

import config


def load_network(path=config.NETWORK_FILE):
    """Read an edge list (first two columns = interacting genes) into a graph.

    In:  tab-separated file, one interaction per line
         (HuRI: Ensembl gene ids; STRING after filtering: gene symbols).
    Out: undirected networkx Graph without self-loops.

    DIAMOnD treats every interaction as undirected and unweighted, so that is
    how the network is built here too.
    """
    edges = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1],
                        names=["a", "b"], dtype=str, comment="#")
    if config.NETWORK == "huri":
        # drop a header line if the file has one (real ids start with "ENSG")
        edges = edges[edges["a"].str.startswith("ENSG")]
    G = nx.from_pandas_edgelist(edges, "a", "b")
    # a protein interacting with itself carries no information about
    # connections to *other* disease proteins, so self-loops are removed
    G.remove_edges_from(nx.selfloop_edges(G))
    return G


def load_seeds(G, path=config.SEEDS_FILE):
    """Return the seed gene ids (sorted list) that are nodes of G.

    Seeds missing from the network can never be found or used by DIAMOnD,
    so they are dropped here (the original code does the same).
    """
    seeds = pd.read_csv(path, sep="\t", dtype=str)
    # HuRI nodes are Ensembl ids, STRING nodes are gene symbols
    in_network = sorted(set(seeds[config.SEED_ID_COLUMN]) & set(G.nodes()))
    return in_network


def load_diamond(path=config.DIAMOND_PY3):
    """Import the ported DIAMOnD module from data/diamond_py3.py."""
    if not path.exists():
        raise FileNotFoundError(
            f"{path} not found - run scripts/01_download_data.py first")
    spec = importlib.util.spec_from_file_location("diamond_py3", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def run_diamond(diamond, G, seeds, n_iterations, alpha):
    """Run the authors' DIAMOnD core loop.

    In:  graph G, iterable of seed ids, number of genes to add, seed weight alpha.
    Out: list of the added gene ids, in the order DIAMOnD added them
         (position 0 = first, most significant, gene).

    diamond_iteration_of_first_X_nodes returns tuples (gene, k, kb, p);
    only the gene id is kept here.
    """
    added = diamond.diamond_iteration_of_first_X_nodes(
        G, set(seeds), n_iterations, alpha)
    return [gene for gene, k, kb, p in added]
