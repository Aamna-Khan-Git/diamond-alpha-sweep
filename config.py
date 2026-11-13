"""
config.py - every setting used by the pipeline lives here.

Change DISEASE_ID to run the analysis on a different disease. Everything
else can stay as it is.
"""
from pathlib import Path

# ---------------------------------------------------------------------------
# Which network to use
#   "string": STRING v12 physical subnetwork, experimental + curated-database
#             evidence only (main analysis; similar in size to the paper's)
#   "huri":   HuRI binary yeast two-hybrid map (comparison: sparse, but free
#             of literature study bias)
# Results for each network are kept in their own folder: results/<network>/
# ---------------------------------------------------------------------------
NETWORK = "string"

# ---------------------------------------------------------------------------
# Folders
# ---------------------------------------------------------------------------
ROOT = Path(__file__).resolve().parent
DATA_DIR = ROOT / "data"                     # downloaded raw data (not committed)
RESULTS_DIR = ROOT / "results" / NETWORK     # CSV outputs (small, committed)
FIG_DIR = RESULTS_DIR / "figures"            # PNG figures (committed)

# ---------------------------------------------------------------------------
# The authors' DIAMOnD code, pinned to an exact commit so the run is
# reproducible even if the repository changes later.
# ---------------------------------------------------------------------------
DIAMOND_COMMIT = "23e44dad422f78d4af79ade8793f4d56d8d261a8"
DIAMOND_URL = (
    "https://raw.githubusercontent.com/barabasilab/DIAMOnD/"
    f"{DIAMOND_COMMIT}/DIAMOnD.py"
)
DIAMOND_ORIGINAL = DATA_DIR / "DIAMOnD_original.py"   # as downloaded (Python 2)
DIAMOND_PY3 = DATA_DIR / "diamond_py3.py"             # ported to Python 3

# ---------------------------------------------------------------------------
# Network files
# ---------------------------------------------------------------------------
# HuRI (Luck et al. 2020): two columns of Ensembl gene ids. Address variants
# are tried in order, because the site's certificate does not always match
# its address; if all fail, the download script explains the manual download.
HURI_URLS = [
    "https://interactome-atlas.org/data/HuRI.tsv",
    "http://interactome-atlas.org/data/HuRI.tsv",
    "https://www.interactome-atlas.org/data/HuRI.tsv",
]
HURI_FILE = DATA_DIR / "HuRI.tsv"

# STRING v12, human (9606): physical subnetwork with one score per evidence
# channel, plus the protein info table that maps STRING ids to gene symbols.
STRING_LINKS_URL = ("https://stringdb-downloads.org/download/"
                    "protein.physical.links.detailed.v12.0/"
                    "9606.protein.physical.links.detailed.v12.0.txt.gz")
STRING_INFO_URL = ("https://stringdb-downloads.org/download/protein.info.v12.0/"
                   "9606.protein.info.v12.0.txt.gz")
STRING_LINKS_FILE = DATA_DIR / "9606.protein.physical.links.detailed.v12.0.txt.gz"
STRING_INFO_FILE = DATA_DIR / "9606.protein.info.v12.0.txt.gz"
STRING_FILE = DATA_DIR / "string_physical_filtered.tsv"   # built from the two above
# Keep a link only if its experimental OR curated-database score is >= this
# (0-1000). Text-mining-only links are dropped, like the paper, which used
# only interactions with experimental evidence.
STRING_MIN_SCORE = 700

# what the rest of the pipeline uses, depending on NETWORK
NETWORK_FILE = HURI_FILE if NETWORK == "huri" else STRING_FILE
SEED_ID_COLUMN = "ensembl_id" if NETWORK == "huri" else "symbol"
NETWORK_LABEL = "HuRI" if NETWORK == "huri" else "STRING (physical, exp./database)"

# ---------------------------------------------------------------------------
# Disease genes: Open Targets Platform (GraphQL API, no account needed).
# Find a disease ID by searching https://platform.opentargets.org and copying
# the ID from the disease page URL (EFO_..., MONDO_... etc.).
# ---------------------------------------------------------------------------
OPENTARGETS_API = "https://api.platform.opentargets.org/api/v4/graphql"
DISEASE_ID = "MONDO_0005575"    # colorectal cancer (change to your disease)

# Evidence types that can make a gene a seed. A gene's seed score is its
# highest score across these types. Genetic association = inherited variants
# (GWAS, rare-variant studies); somatic mutation = recurrently mutated in
# tumours (e.g. APC, KRAS, TP53). The paper's seeds came from OMIM + GWAS, and
# OMIM's colorectal cancer entry includes somatic driver genes, so both are used.
EVIDENCE_TYPES = ["genetic_association", "somatic_mutation"]
MIN_EVIDENCE_SCORE = 0.5        # Open Targets evidence score, 0-1
MAX_SEEDS = 150                 # keep at most the top-scoring 150 genes
MIN_SEEDS = 20                  # the paper required >= 20 genes per disease

ASSOCIATIONS_FILE = DATA_DIR / f"opentargets_{DISEASE_ID}.tsv"
SEEDS_FILE = DATA_DIR / "seeds.tsv"

# ---------------------------------------------------------------------------
# The experiment
# ---------------------------------------------------------------------------
ALPHAS = [1, 2, 3, 5, 10]   # seed weights to compare (alpha = 1 is the original)
HIDE_FRACTION = 0.2         # fraction of seed genes hidden in each leave-out run
N_REPEATS = 10              # number of different random hidden sets
N_ITERATIONS = 200          # DIAMOnD genes added per run (paper's ~200 cutoff)
N_RANDOM_SETS = 1000        # random gene sets for the seed-clustering z-score
RANDOM_SEED = 42            # makes every random choice reproducible
N_JOBS = None               # parallel processes (None = all CPU cores)
