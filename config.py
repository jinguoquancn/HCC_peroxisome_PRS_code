"""Central configuration for the HCC peroxisome analysis repository.

All scripts resolve their input/output paths through PROJECT_ROOT, so the
repository is relocatable: clone it anywhere, keep the expected layout
(data/processed, scripts/, figures/), and run.
"""
from pathlib import Path

# Repository root (this file lives at <repo>/config.py)
PROJECT_ROOT = Path(__file__).resolve().parent

# Standardized directories
DATA_DIR = PROJECT_ROOT / "data"
PROCESSED_DIR = DATA_DIR / "processed"
GENESETS_DIR = DATA_DIR / "genesets"   # MSigDB gene sets; not bundled (license)
RAW_DIR = DATA_DIR / "raw"             # downloaded raw data; not bundled
FIGURES_DIR = PROJECT_ROOT / "figures"

# Reproducibility
SEED = 42

# Cohort accessions (all publicly available)
SC_ACCESSIONS = ["GSE149614", "GSE151530", "GSE125449", "GSE146115", "GSE98638"]
BULK_ACCESSIONS = ["GSE14520", "GSE76427", "GSE116174", "GSE10186",
                   "GSE25097", "GSE36376", "GSE54236", "GSE87630"]
TCGA_COHORT = "TCGA-LIHC"
ICGC_COHORT = "ICGC-LIRI-JP"
GDSC_RELEASE = "GDSC1"


def processed(*parts):
    """Absolute path inside data/processed/."""
    return str(PROCESSED_DIR.joinpath(*parts))
