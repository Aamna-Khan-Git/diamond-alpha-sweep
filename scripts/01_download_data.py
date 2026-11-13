"""
01_download_data.py - download everything the analysis needs into data/.
 
    1. The authors' DIAMOnD.py (pinned commit), ported to Python 3
    2. The network chosen in config.NETWORK (STRING physical or HuRI)
    3. Disease-gene associations for config.DISEASE_ID from Open Targets,
       filtered into a seed gene list
 
Run from the repository root:   python scripts/01_download_data.py
Files that already exist are not downloaded again (delete them to refresh).
"""
import sys
from pathlib import Path
 
import pandas as pd
import requests
 
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # repo root
import config  # noqa: E402
from src.port_diamond import port_to_python3  # noqa: E402
from src.utils import load_network  # noqa: E402
 
 
def download(url, dest):
    """Stream a file from url to dest (skips if dest already exists)."""
    if dest.exists():
        print(f"  already have {dest.name}")
        return
    print(f"  downloading {url}")
    with requests.get(url, stream=True, timeout=120) as r:
        r.raise_for_status()
        with open(dest, "wb") as f:
            for chunk in r.iter_content(chunk_size=1 << 20):
                f.write(chunk)
 
 
# ---------------------------------------------------------------------------
# 1. DIAMOnD code
# ---------------------------------------------------------------------------
def get_diamond_code():
    print("[1/3] DIAMOnD code")
    download(config.DIAMOND_URL, config.DIAMOND_ORIGINAL)
    ported = port_to_python3(config.DIAMOND_ORIGINAL.read_text())
    config.DIAMOND_PY3.write_text(ported)
    print(f"  wrote Python 3 port to {config.DIAMOND_PY3.name}")
 
 
# ---------------------------------------------------------------------------
# 2. Network
# ---------------------------------------------------------------------------
def get_huri():
    """Download HuRI (two columns of Ensembl gene ids)."""
    if config.HURI_FILE.exists():
        print(f"  already have {config.HURI_FILE.name}")
        return
    # The HuRI site's security certificate does not always match its
    # address, which makes the download fail. Try a few address variants;
    # if none works, ask for a manual download (one file, done once).
    for url in config.HURI_URLS:
        try:
            download(url, config.HURI_FILE)
            return
        except requests.exceptions.RequestException as err:
            print(f"  could not download from {url} ({type(err).__name__})")
            config.HURI_FILE.unlink(missing_ok=True)
    sys.exit(
        "\n  Automatic HuRI download failed. Download it manually:\n"
        "    1. In a browser, open http://www.interactome-atlas.org/download\n"
        "    2. Under 'HuRI', download the file HuRI.tsv\n"
        f"    3. Save it as: {config.HURI_FILE}\n"
        "    4. Run this script again (it will skip the download).\n")
 
 
def get_string():
    """Download STRING's human physical subnetwork and keep only links with
    strong experimental or curated-database evidence, named by gene symbol."""
    if config.STRING_FILE.exists():
        print(f"  already have {config.STRING_FILE.name}")
        return
    for url, dest in [(config.STRING_LINKS_URL, config.STRING_LINKS_FILE),
                      (config.STRING_INFO_URL, config.STRING_INFO_FILE)]:
        try:
            download(url, dest)
        except requests.exceptions.RequestException as err:
            dest.unlink(missing_ok=True)
            sys.exit(
                f"\n  Could not download {url} ({type(err).__name__}).\n"
                "  Download it manually from https://string-db.org/cgi/download\n"
                "  (choose organism Homo sapiens) and save it as:\n"
                f"    {dest}\n  then run this script again.\n")
 
    # STRING id -> gene symbol, e.g. 9606.ENSP00000257430 -> APC
    info = pd.read_csv(config.STRING_INFO_FILE, sep="\t", usecols=[0, 1])
    info.columns = ["string_id", "symbol"]
    to_symbol = dict(zip(info["string_id"], info["symbol"]))
 
    # columns: protein1 protein2 experimental database textmining combined_score
    links = pd.read_csv(config.STRING_LINKS_FILE, sep=" ")
    print(f"  STRING physical subnetwork: {len(links):,} links (both directions)")
    # find the evidence columns by name (spelled "experimental" in v12,
    # "experiments" in some older files)
    exp_col = next((c for c in links.columns if c.startswith("experiment")), None)
    db_col = next((c for c in links.columns if c.startswith("database")), None)
    if exp_col is None or db_col is None:
        sys.exit(f"  Unexpected STRING columns: {list(links.columns)}")
    evidence = links[[exp_col, db_col]].max(axis=1)
    links = links[evidence >= config.STRING_MIN_SCORE]
 
    edges = pd.DataFrame({"a": links["protein1"].map(to_symbol),
                          "b": links["protein2"].map(to_symbol)}).dropna()
    # STRING lists every link twice (A-B and B-A): put each pair in
    # alphabetical order, then keep one copy
    swap = edges["a"] > edges["b"]
    edges.loc[swap, ["a", "b"]] = edges.loc[swap, ["b", "a"]].values
    edges = edges.drop_duplicates()
    edges = edges[edges["a"] != edges["b"]]
    edges.to_csv(config.STRING_FILE, sep="\t", header=False, index=False)
    print(f"  kept {len(edges):,} links with experiments or database "
          f">= {config.STRING_MIN_SCORE}")
 
 
def get_network():
    print(f"[2/3] Network: {config.NETWORK_LABEL}")
    if config.NETWORK == "huri":
        get_huri()
    else:
        get_string()
    G = load_network()
    print(f"  {G.number_of_nodes():,} proteins, {G.number_of_edges():,} interactions")
 
 
# ---------------------------------------------------------------------------
# 3. Disease genes from Open Targets
# ---------------------------------------------------------------------------
# GraphQL query: one page of genes associated with a disease, with the
# overall score and the score for each evidence type (genetic, literature...).
QUERY = """
query diseaseGenes($id: String!, $index: Int!, $size: Int!) {
  disease(efoId: $id) {
    id
    name
    associatedTargets(page: {index: $index, size: $size}) {
      count
      rows {
        score
        target { id approvedSymbol }
        datatypeScores { id score }
      }
    }
  }
}
"""
 
 
def fetch_associations(disease_id, page_size=500):
    """Download every gene associated with the disease (all pages).
 
    Out: (disease name, DataFrame with ensembl_id, symbol, overall_score,
          and one score_<type> column per evidence type)
    """
    rows, index, name = [], 0, None
    while True:
        r = requests.post(config.OPENTARGETS_API, timeout=120, json={
            "query": QUERY,
            "variables": {"id": disease_id, "index": index, "size": page_size},
        })
        r.raise_for_status()
        payload = r.json()
        if payload.get("errors"):
            raise RuntimeError(payload["errors"])
        disease = payload["data"]["disease"]
        if disease is None:
            raise ValueError(f"Open Targets has no disease with id {disease_id}")
        name = disease["name"]
        page = disease["associatedTargets"]["rows"]
        for row in page:
            record = {
                "ensembl_id": row["target"]["id"],
                "symbol": row["target"]["approvedSymbol"],
                "overall_score": row["score"],
            }
            # one column per evidence type, e.g. score_genetic_association
            for d in row["datatypeScores"]:
                record[f"score_{d['id']}"] = d["score"]
            rows.append(record)
        total = disease["associatedTargets"]["count"]
        index += 1
        if index * page_size >= total or not page:
            break
    return name, pd.DataFrame(rows)
 
 
def get_seeds():
    print(f"[3/3] Disease genes for {config.DISEASE_ID}")
    wanted = [f"score_{t}" for t in config.EVIDENCE_TYPES]
    assoc = None
    if config.ASSOCIATIONS_FILE.exists():
        assoc = pd.read_csv(config.ASSOCIATIONS_FILE, sep="\t")
        if "overall_score" in assoc and any(c.startswith("score_") for c in assoc):
            print(f"  already have {config.ASSOCIATIONS_FILE.name}")
        else:
            assoc = None   # file from an older version of this script: refetch
    if assoc is None:
        name, assoc = fetch_associations(config.DISEASE_ID)
        assoc.to_csv(config.ASSOCIATIONS_FILE, sep="\t", index=False)
        print(f"  {name}: {len(assoc):,} associated genes in total")
 
    # evidence types with no data for this disease simply count as 0
    for col in wanted:
        if col not in assoc:
            assoc[col] = 0.0
    assoc[wanted] = assoc[wanted].fillna(0.0)
    # seed score = best score across the chosen evidence types
    assoc["evidence_score"] = assoc[wanted].max(axis=1)
    assoc["evidence_type"] = assoc[wanted].idxmax(axis=1).str.replace("score_", "")
 
    seeds = (assoc[assoc["evidence_score"] >= config.MIN_EVIDENCE_SCORE]
             .sort_values("evidence_score", ascending=False)
             .head(config.MAX_SEEDS))
 
    # Only seeds that are in the network can be used.
    G = load_network()
    in_net = seeds[seeds[config.SEED_ID_COLUMN].isin(set(G.nodes()))]
    print(f"  {len(seeds)} genes have {' or '.join(config.EVIDENCE_TYPES)} "
          f">= {config.MIN_EVIDENCE_SCORE}; {len(in_net)} of them are in the network")
    print("  seeds in the network by evidence type: "
          + ", ".join(f"{t} {n}" for t, n in in_net["evidence_type"].value_counts().items()))
 
    # Quick check of DIAMOnD's premise: do the seeds interact with each other?
    sub = G.subgraph(in_net[config.SEED_ID_COLUMN])
    print(f"  links between seeds: {sub.number_of_edges()}; seeds with at least "
          f"one seed partner: {sum(1 for n in sub if sub.degree(n) > 0)}")
    if len(in_net) < config.MIN_SEEDS:
        print(f"  WARNING: fewer than {config.MIN_SEEDS} seeds in the network. "
              "Lower MIN_EVIDENCE_SCORE or choose another disease in config.py.")
    seeds.to_csv(config.SEEDS_FILE, sep="\t", index=False)
    print(f"  wrote {config.SEEDS_FILE.name}")
 
 
if __name__ == "__main__":
    config.DATA_DIR.mkdir(exist_ok=True)
    get_diamond_code()
    get_network()
    get_seeds()
    print("Done.")