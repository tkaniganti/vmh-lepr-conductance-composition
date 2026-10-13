"""config.py - central configuration. Nothing downstream hard-codes a gene list or a path."""
import os

ROOT    = os.environ.get("PROJECT", os.path.dirname(os.path.abspath(__file__)))
DATA    = os.path.join(ROOT, "data")
RESULTS = os.path.join(ROOT, "results")
FIGURES = os.path.join(ROOT, "figures")
for _d in (DATA, RESULTS, FIGURES):
    os.makedirs(_d, exist_ok=True)

# ---------------------------------------------------------------- accessions
ACCESSIONS = {
    "snRNA_VMH":  "GSE172204",   # Affinati 2021 VMH snRNA-seq
    "TRAP_Nr5a1": "GSE172206",   # SF1 VMH TRAP: Bead B1-B3 = 116053-116055, Sup S1-S3 = 116056-116058
    "TRAP_Lepr":  "GSE172200",   # Lepr VMH TRAP: 66853 Bead/PBS, 66854 Sup/PBS, 66855 Bead/Leptin, 66856 Sup/Leptin
    "Rupp_LepRb": "GSE172460",   # Rupp reporter-defined LepR snRNA-seq (NOT GSE172463; that does not exist)
}

# Verified from GSE172206_series_matrix.txt.gz. Do not infer these positionally.
TRAP_NR5A1 = {"bead": ["116053", "116054", "116055"], "sup": ["116056", "116057", "116058"]}
TRAP_LEPR  = {"bead_pbs": "66853", "sup_pbs": "66854",
              "bead_leptin": "66855", "sup_leptin": "66856"}

# ---------------------------------------------------------------- gene panels
PANELS = {
    "HYPER": ["Kcnj11", "Abcc8", "Pde3b", "Kcnj3", "Kcnj6", "Kcnj5", "Kcnj9",
              "Kcnk1", "Kcnk2", "Kcnk3", "Kcnk9", "Kcnq2", "Kcnq3"],
    "DEPOL": ["Trpc1", "Trpc3", "Trpc4", "Trpc5", "Trpc6",
              "Hcn1", "Hcn2", "Hcn3", "Hcn4"],
    "MARKERS": ["Nr5a1", "Slc17a6", "Gad1", "Gad2", "Lepr", "Gck"],
}
KEY_CHANNELS = ["Kcnj11", "Abcc8", "Kcnj3", "Kcnj5", "Kcnj6", "Kcnj9", "Pde3b",
                "Trpc1", "Trpc4", "Trpc5", "Hcn1", "Hcn2"]

# ---------------------------------------------------------------- parameters
NORM_TARGET  = 1e4
CPM_TARGET   = 1e6
N_EXPR_BINS  = 24
N_CTRL_GENES = 100
N_SEEDS      = 20      # control-gene draws; eps2 is reported as mean and range over these
RANDOM_SEED  = 1
N_BOOT       = 5000
N_PERM       = 10000

# Rupp populations with independent ground-truth leptin response sign
RUPP_POPS = {"VMH_SF1": "L3.Nr5a1", "POMC": "L9.Pomc", "AgRP": "L6.Agrp"}

EPS2_BANDS = [(0.01, "negligible"), (0.06, "small"), (0.14, "medium"), (1.01, "large")]

COLORS = {"HYPER": "#4C72B0", "DEPOL": "#DD8452", "LEPR": "#C44E52", "OTHER": "#B0B0B0",
          "snRNA": "#8172B3", "TRAP": "#55A868", "NEUTRAL": "#4C4C4C"}


def eps2_label(e):
    for thr, lab in EPS2_BANDS:
        if e < thr:
            return lab
    return "large"
