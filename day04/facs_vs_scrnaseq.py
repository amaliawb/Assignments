"""
facs_vs_scrnaseq.py
===================

Compare urine cell populations measured two ways:

  1. Flow cytometry / imaging sorting on the BD FACSDiscover S8
     (gated in FlowJo, exported as an Excel table)
  2. 10x Genomics Flex single-cell RNA-seq
     (processed and annotated in a Jupyter notebook, exported as a CSV)

The two methods do NOT have to be run on the same samples. Each sample in the
sample sheet belongs to one method. For example, the scRNA-seq results of a
group of patients can be compared with the FACS results of
  - a POOL of healthy urine from several people (one FACS sample), or
  - SEVERAL patient samples measured by FACS (other people than the scRNA-seq ones).
Because the samples differ, the methods are compared as GROUPS (unpaired),
not sample by sample. Several samples of the same person are first combined
into one (donor_id column), so every person counts once. A FACS pool is ONE
reference value: its proportions are shown next to the Flex proportions and
the reference atlas (plot A4), and the Flex people are tested against it.

FACS panels and gating (choose the panel in SETTINGS, with --panel, or on the website):
  Both panels: Cells -> Single Cells -> Live (Fixable Yellow negative)
               Live  -> CD45+ -> SSC high = Neutrophils / SSC low = Other immune
  Panel "CD45_EpCAM" (CD45 + EpCAM):
               CD45- -> EpCAM high = distal tubule, collecting duct, urothelium
                     -> EpCAM low  = proximal tubule (the CD10+ CD13+ cells)
                                     (+ other EpCAM-negative CD45- cells)
  Panel "CD45_EpCAM_CD10_CD13" (CD45 + EpCAM + CD10 + CD13):
               CD45- -> CD10+ CD13+ = proximal tubule
                     -> NOT (CD10+ CD13+) -> EpCAM high = distal tubule, collecting duct, urothelium
                                          -> EpCAM low  = other CD45- cells

The populations of a panel do not overlap: each live cell falls in exactly one,
so the proportions add up to 100% of live cells. Each population is linked to
the Flex cell-type labels that would fall into its gate, so FACS and Flex are
always compared population by population, whichever panel was used.

What the program produces (in the results folder):
  A. Cell-type PROPORTIONS, scRNA-seq group vs FACS group
       - composition overview: the average proportions of every FACS group
         (e.g. the pool), every scRNA-seq group and the reference atlas side by side
       - every chosen scRNA-seq condition is compared with every chosen FACS
         condition (e.g. "Disease (seq) vs Healthy pool (FACS)")
       - difference and ratio of the medians, with a bootstrap 95% interval
       - groups of >= 3 samples on both sides : Mann-Whitney U test
         one side is a single sample (a pool)  : one-sample Wilcoxon test of the
                                                 other side against that value
       - all p-values corrected for multiple testing (Benjamini-Hochberg FDR)
  B. Differences between CONDITIONS within each method
       - 2 conditions  : Mann-Whitney U test
       - 3+ conditions : Kruskal-Wallis test + pairwise Mann-Whitney U tests
  C. PROTEIN vs RNA per population: typical FACS median intensity of each
     marker vs typical RNA expression of its gene (CD45-PTPRC, EpCAM-EPCAM,
     and with the 4-marker panel also CD10-MME, CD13-ANPEP). Do the same
     populations come out high/low?
  D. S8 IMAGING features per population (FACS samples): differences between
     populations and between conditions
  E. REFERENCE ATLAS (optional): compare FACS and Flex with a published
     single-cell atlas of urine with known cell types
       - proportions of each population vs the atlas (per condition)
       - how different each sample's whole composition is from the atlas
       - marker patterns (atlas RNA vs Flex RNA vs FACS protein)
  F. IMAGING PROFILES per cell type (S8 imaging features)
       - mean +/- 95% confidence interval of each feature per population
       - "morphological fingerprint" of each population (heatmap)
       - PCA map: do the populations separate by imaging alone?
       - outlier samples (robust z-score) and imaging vs marker correlations
       - with a per-cell export (optional): cell distributions, % outlier
         cells per sample (e.g. doublets, clumps) and a cell morphology map

-------------------------------------------------------------------------------
HOW TO RUN
-------------------------------------------------------------------------------
  Try it first with made-up demo data (creates a "demo_data" folder):
      python facs_vs_scrnaseq.py --demo
      python facs_vs_scrnaseq.py --demo --panel CD45_EpCAM      (the 2-marker panel)

  With your real data (after editing the SETTINGS section below):
      python facs_vs_scrnaseq.py

  The program will ask which metadata column holds the conditions, and which
  conditions to include for the scRNA-seq samples and for the FACS samples.
  You can skip the questions by giving them directly:
      python facs_vs_scrnaseq.py --condition-column condition \
          --seq-conditions Disease --facs-conditions "Healthy pool" Disease

  Or use the web interface instead of the terminal (see app.py):
      python app.py        then open http://127.0.0.1:5000

  Needs: pandas, numpy, scipy, matplotlib, openpyxl (all included in Anaconda).

-------------------------------------------------------------------------------
INPUT FILES
-------------------------------------------------------------------------------
1) metadata.csv - one row per sample (FACS samples AND scRNA-seq samples).
   Required columns:
      sample_id   : the sample name (for scRNA-seq: as used in the scRNA-seq file)
      method      : "FACS" or "scRNA-seq" - which method measured this sample
      facs_name   : FACS samples only: the sample name exactly as it appears in
                    the FlowJo export (usually the .fcs file name). Leave it
                    empty for scRNA-seq samples.
   Optional column:
      donor_id    : the PERSON the sample comes from. Several samples of one
                    person (e.g. 12 Flex samples from 8 people) are combined
                    into one per person before the statistics: samples of the
                    same person are not independent. Empty = its own person.
   Plus any number of columns describing the samples (condition, sex,
   flex_probe_barcode, ...). You choose the condition column when you run.
   Example:
      sample_id,method,facs_name,donor_id,condition
      S01,scRNA-seq,,P1,Disease
      S02,scRNA-seq,,P1,Disease       <- a second sample of person P1
      S03,scRNA-seq,,P2,Disease
      POOL1,FACS,healthy_pool.fcs,,Healthy pool
   A pool of urine from several people is ONE FACS sample (one row). Only list
   repeated measurements of the same pool as separate rows if you accept that
   the tests then treat them as independent samples (they are not).
   (An older sheet WITHOUT a method column is still accepted: every row is then
   taken to be measured by both methods.)

2) FlowJo Excel export (Table Editor -> export as .xlsx). One row per FACS sample.
   The first column holds the sample name; the other columns are statistics:
      - the Count of the total gate (e.g. live single cells)
      - the Count of each population gate of the chosen panel
      - the Median of each marker channel of the panel inside each population gate
      - the Median of each S8 imaging feature inside each population gate
   Column names are built from the templates in SETTINGS (COUNT_COLUMN,
   MEDIAN_COLUMN). If a name does not match, the program prints the column
   names it found in your file so you can fix the SETTINGS.

3) scrnaseq_cells.csv - one row per cell, exported from your notebook.
   Columns: cell_barcode, sample_id, cell_type, PTPRC, MME, ANPEP, EPCAM
   (MME and ANPEP are only needed with the CD45_EpCAM_CD10_CD13 panel; cells already demultiplexed by Flex probe barcode -> sample_id, and
   already annotated). Gene values should be log-normalized expression.

   Why this format: it holds everything needed (labels + the 4 marker genes),
   is small, and opens with pandas alone - no scanpy needed here.
   Export it from your Scanpy notebook like this:

       genes = ["PTPRC", "MME", "ANPEP", "EPCAM"]
       # use log-normalized values (not scaled). If adata.X is scaled,
       # use adata.raw.to_adata() or the layer holding log-normalized counts.
       export = adata.obs[["sample_id", "cell_type"]].copy()
       expr = adata[:, genes].X
       export[genes] = expr.toarray() if hasattr(expr, "toarray") else expr
       export.index.name = "cell_barcode"
       export.to_csv("scrnaseq_cells.csv")

4) reference_atlas.csv - OPTIONAL. A reference single-cell atlas of urine with
   known cell types (e.g. a published dataset). If the file is not there,
   section E is simply skipped. Two formats are accepted:
     a) one row per cell (like file 3):
           cell_type, donor_id (optional), PTPRC, MME, ANPEP, EPCAM (optional)
        -> gives proportions AND the marker comparison (plot E3)
     b) one row per cell type (or per cell type x donor), e.g. from a paper's
        supplementary table:
           cell_type, donor_id (optional), n_cells
        -> gives proportions only
   The atlas cell-type names are translated to the populations with
   "atlas_labels" in POPULATIONS (SETTINGS). With a donor_id column, the
   variation between atlas donors is used in the statistics.
   Export from a Scanpy atlas object like this (keep only the columns you have):

       genes = ["PTPRC", "MME", "ANPEP", "EPCAM"]
       export = atlas.obs[["cell_type", "donor_id"]].copy()
       expr = atlas[:, genes].X        # log-normalized values, as for file 3
       export[genes] = expr.toarray() if hasattr(expr, "toarray") else expr
       export.to_csv("reference_atlas.csv", index=False)

5) facs_imaging_cells.csv - OPTIONAL. The single FACS cells with their imaging
   features (section F, cell-level plots). One row per cell:
       facs_name (or sample_id), population, <one column per IMAGING_FEATURE>
   "population" must be one of the population names in POPULATIONS.
   In FlowJo, export each population gate as CSV (scale values) and combine the
   files, e.g. when they are named "<facs_name>__<population>.csv":

       import pandas as pd, pathlib
       parts = []
       for f in pathlib.Path("exports").glob("*.csv"):
           facs_name, population = f.stem.split("__")
           part = pd.read_csv(f)
           part.insert(0, "population", population)
           part.insert(0, "facs_name", facs_name)
           parts.append(part)
       pd.concat(parts).to_csv("facs_imaging_cells.csv", index=False)
"""

# -----------------------------------------------------------------------------
# IMPORTS - the libraries this program uses
# -----------------------------------------------------------------------------
import argparse            # reads options typed on the command line (--demo, --seq-conditions, ...)
import functools           # remembers (caches) results that are needed again
import itertools           # all possible orderings, for the exact Spearman p-value
from pathlib import Path   # handles file and folder paths the same way on Mac and Windows

import matplotlib
matplotlib.use("Agg")  # draw plots straight to files (no window pops up)
import matplotlib.pyplot as plt   # makes the plots
import numpy as np                # fast maths on lists of numbers
import pandas as pd               # tables (read CSV/Excel, group, merge, save)
from scipy import stats           # the statistical tests (Mann-Whitney, Spearman, ...)


# =============================================================================
# SETTINGS - everything you may need to change is in this section
# =============================================================================

# ---- File locations (relative to the folder this script is in) -------------
SCRIPT_DIR = Path(__file__).resolve().parent        # the folder this .py file sits in
DATA_DIR = SCRIPT_DIR / "data"                      # put your real files here
METADATA_FILE = "metadata.csv"                      # sample sheet (one row per sample)
FACS_FILE = "facs_flowjo_export.xlsx"               # FlowJo table export
FACS_SHEET = 0                                      # sheet number or name in the Excel file
SEQ_FILE = "scrnaseq_cells.csv"                     # one row per sequenced cell
OUTPUT_DIR = SCRIPT_DIR / "results"                 # all tables and plots go here

# Column names in the metadata and scRNA-seq files
SAMPLE_COLUMN = "sample_id"      # the sample ID shared by metadata and scRNA-seq file
METHOD_COLUMN = "method"         # which method measured the sample: FACS or scRNA-seq
DONOR_COLUMN = "donor_id"        # optional: the PERSON each sample comes from. Samples of the
                                 # same person (and method and condition) are combined into one
                                 # before the statistics, because they are not independent.
                                 # Empty / missing = every sample is its own person.
FACS_NAME_COLUMN = "facs_name"   # the sample name as FlowJo writes it (FACS samples)
CELL_TYPE_COLUMN = "cell_type"   # the annotation label of each sequenced cell

# The two methods, and the ways they may be written in the method column
# (compared in lower case, without spaces, dashes or underscores)
METHODS = ["scRNA-seq", "FACS"]
METHOD_ALIASES = {
    "facs": "FACS", "flow": "FACS", "flowcytometry": "FACS", "s8": "FACS",
    "scrnaseq": "scRNA-seq", "seq": "scRNA-seq", "rnaseq": "scRNA-seq",
    "flex": "scRNA-seq", "10x": "scRNA-seq", "10xflex": "scRNA-seq",
}

# ---- Reference atlas (optional, section E) -----------------------------------
# Put the atlas file in the data folder (or upload it on the website). If it is
# missing, section E is skipped and everything else runs as before.
ATLAS_FILE = "reference_atlas.csv"
ATLAS_CELL_TYPE_COLUMN = "cell_type"   # column with the atlas' cell-type names
ATLAS_DONOR_COLUMN = "donor_id"        # optional: which atlas donor/sample each row comes from
ATLAS_COUNT_COLUMN = "n_cells"         # optional: number of cells in the row (summary-table format)

# ---- FlowJo column-name templates --------------------------------------------
# FlowJo names a statistic column "<gate path> | <statistic>".
# {gate} is replaced by a gate path below, {parameter} by a channel/feature name.
COUNT_COLUMN = "{gate} | Count"
MEDIAN_COLUMN = "{gate} | Median ({parameter})"

# ---- Which antibodies were used: the FACS PANEL ---------------------------------
# Choose the panel that matches your staining. It decides the gating (which
# populations exist), which Flex labels each population is compared with, and
# which markers are compared with their genes (section C). It can also be chosen
# per run: on the command line with --panel, or on the website.
#   "CD45_EpCAM"           : CD45 + EpCAM (+ Fixable Yellow)
#   "CD45_EpCAM_CD10_CD13" : CD45 + EpCAM + CD10 + CD13 (+ Fixable Yellow)
FACS_PANEL = "CD45_EpCAM_CD10_CD13"

# The viability dye is used in both panels: "Live" = cells NEGATIVE for Fixable
# Yellow (dead cells take up the dye and are bright). Use the gate names exactly as in FlowJo.
LIVE_GATE = "Cells/Single Cells/Live"     # live/dead gate on Fixable Yellow
VIABILITY_DYE = "Fixable Yellow"          # only for documentation; not used in calculations
CD45_POS_GATE = f"{LIVE_GATE}/CD45+"
CD45_NEG_GATE = f"{LIVE_GATE}/CD45-"
NOT_PT_GATE = f"{CD45_NEG_GATE}/NOT CD10+ CD13+"

# Gate that holds ALL good cells (used only if PROPORTION_DENOMINATOR = "total")
FACS_TOTAL_GATE = LIVE_GATE

# ---- The two panels ---------------------------------------------------------------
# For every panel, the populations it can tell apart. For every population:
# "gate"        : the FlowJo gate path of the population
# "seq_labels"  : which cell_type labels from your Flex notebook belong to this
#                 population (list ALL the labels you used, e.g. every immune subtype)
# "atlas_labels": which cell_type names of the REFERENCE ATLAS belong to this
#                 population (an atlas usually names cell types differently from
#                 your notebook). If left out, the seq_labels are used.
# The labels are placeholders - copy the real names from your notebook / atlas.
# The gates of a panel do not overlap and together cover every live cell, so each
# cell is counted exactly once (like in scRNA-seq, where each cell has one label).
# That is why every Flex label must be listed under the population whose GATE
# that cell type would fall into on FACS - e.g. with CD45 + EpCAM only, podocytes
# and fibroblasts (CD45-, EpCAM low) land in the EpCAM-low gate together with the
# proximal tubule cells, so they are listed there too.
#
# Both panels start the same way:
#   Live -> CD45+ -> SSC high : Neutrophils   (granular cells: high side scatter)
#                 -> SSC low  : Other immune   (lymphocytes, monocytes, macrophages)
#        -> CD45-  -> ... (differs per panel, see below)
# EpCAM high = distal tubule, collecting duct and urothelium.
# EpCAM low  = proximal tubule (the same cells as CD10+ CD13+), plus other
#              EpCAM-negative CD45- cells.
IMMUNE_POPULATIONS = {
    "Neutrophils": {
        "gate": f"{CD45_POS_GATE}/SSC high",
        "seq_labels": ["Neutrophil"],
        "atlas_labels": ["Neutrophil"],
    },
    "Other immune": {
        "gate": f"{CD45_POS_GATE}/SSC low",
        "seq_labels": ["T cell", "NK cell", "B cell", "Monocyte", "Macrophage"],
        "atlas_labels": ["T/NK cell", "B cell", "Myeloid cell"],
    },
}
EPCAM_HIGH_LABELS = {
    "seq_labels": ["Distal tubule cell", "Collecting duct cell", "Urothelial cell"],
    "atlas_labels": ["Distal tubule", "Collecting duct", "Urothelial"],
}
# Squamous epithelial cells (from the urethra / vagina) are not listed: add their
# label to whichever gate they fall into in your own data.

PANELS = {
    # CD45 + EpCAM:
    #   CD45- -> EpCAM high : distal tubule, collecting duct, urothelium
    #         -> EpCAM low  : proximal tubule (+ other EpCAM-negative CD45- cells)
    "CD45_EpCAM": {
        "description": "CD45 + EpCAM (+ Fixable Yellow)",
        "markers": ["CD45", "EpCAM"],
        "populations": {
            **IMMUNE_POPULATIONS,
            "EpCAM-low (proximal tubule)": {
                "gate": f"{CD45_NEG_GATE}/EpCAM low",
                "seq_labels": ["Proximal tubule cell", "Podocyte", "Fibroblast"],
                "atlas_labels": ["Proximal tubule", "Podocyte", "Fibroblast"],
            },
            "EpCAM-high epithelium": {
                "gate": f"{CD45_NEG_GATE}/EpCAM high",
                **EPCAM_HIGH_LABELS,
            },
        },
    },
    # CD45 + EpCAM + CD10 + CD13:
    #   CD45- -> CD10+ CD13+         : proximal tubule
    #         -> NOT (CD10+ CD13+) -> EpCAM high : distal tubule, collecting duct, urothelium
    #                              -> EpCAM low  : other CD45- cells (podocytes, fibroblasts, ...)
    # Why CD10/CD13 is gated BEFORE EpCAM: proximal tubule cells can carry some
    # EpCAM. If the EpCAM and CD10+ CD13+ gates both sat directly under CD45-, a
    # cell positive for both would be counted twice. "NOT CD10+ CD13+" is the
    # inverse of the CD10+ CD13+ gate (FlowJo: Boolean NOT gate, or Make Gate -> Invert).
    "CD45_EpCAM_CD10_CD13": {
        "description": "CD45 + EpCAM + CD10 + CD13 (+ Fixable Yellow)",
        "markers": ["CD45", "CD10", "CD13", "EpCAM"],
        "populations": {
            **IMMUNE_POPULATIONS,
            "Proximal tubule": {
                "gate": f"{CD45_NEG_GATE}/CD10+ CD13+",
                "seq_labels": ["Proximal tubule cell"],
                "atlas_labels": ["Proximal tubule"],
            },
            "EpCAM-high epithelium": {
                "gate": f"{NOT_PT_GATE}/EpCAM high",
                **EPCAM_HIGH_LABELS,
            },
            "Other CD45-": {
                "gate": f"{NOT_PT_GATE}/EpCAM low",
                "seq_labels": ["Podocyte", "Fibroblast"],
                "atlas_labels": ["Podocyte", "Fibroblast"],
            },
        },
    },
}

# ---- Antibody markers and their matching genes ------------------------------
# "facs_parameter" is the channel name FlowJo uses. Replace FLUOR1..4 with your
# fluorochromes, e.g. "Comp-BV510-A :: CD45" - copy it exactly from FlowJo.
# "gene" is the column name of the same marker in the scRNA-seq file.
# Only the markers of the chosen panel are used.
ALL_MARKERS = {
    "CD45":  {"gene": "PTPRC", "facs_parameter": "Comp-FLUOR1-A :: CD45"},
    "CD10":  {"gene": "MME",   "facs_parameter": "Comp-FLUOR2-A :: CD10"},
    "CD13":  {"gene": "ANPEP", "facs_parameter": "Comp-FLUOR3-A :: CD13"},
    "EpCAM": {"gene": "EPCAM", "facs_parameter": "Comp-FLUOR4-A :: EpCAM"},
}


def use_panel(name):
    """
    Switch the whole program to one FACS panel: sets POPULATIONS (the gates and
    their Flex / atlas labels) and MARKERS (the antibodies and their genes) that
    every function below uses.
    """
    global ACTIVE_PANEL, POPULATIONS, MARKERS
    if name not in PANELS:
        raise SystemExit(f"Unknown FACS panel '{name}'. Choose one of: {', '.join(PANELS)}")
    ACTIVE_PANEL = name                 # FACS_PANEL itself stays the default from SETTINGS
    POPULATIONS = PANELS[name]["populations"]
    MARKERS = {marker: ALL_MARKERS[marker] for marker in PANELS[name]["markers"]}


use_panel(FACS_PANEL)   # start with the panel chosen above

# ---- S8 imaging features -------------------------------------------------------
# Names of the imaging parameters exactly as they appear in FlowJo.
# These are placeholders - replace with the features you export (add or remove freely).
# An empty list [] switches off section D completely.
IMAGING_FEATURES = [
    "Area (LightLoss)",
    "Eccentricity (LightLoss)",
    "Radial Moment (LightLoss)",
]

# ---- Imaging profiles per cell type (section F) ---------------------------------
# Section F uses the IMAGING_FEATURES above. It always works from the FlowJo table
# (one median per sample and population). If you ALSO export the single cells,
# it adds cell-level plots (distributions, % outlier cells, morphology map).
# Per-cell file: one row per cell, with columns
#   facs_name (or sample_id), population, and one column per IMAGING_FEATURE
# (in FlowJo: export the events of each population gate as CSV and combine them,
#  adding the sample name and the population name as columns).
FACS_CELLS_FILE = "facs_imaging_cells.csv"     # optional per-cell imaging export
FACS_CELLS_POPULATION_COLUMN = "population"    # column with the population of each cell
# A value counts as an OUTLIER when its robust z-score is above this number.
# Robust z = 0.6745 x (value - median) / MAD, where MAD = median absolute deviation.
# 3.5 is the usual cut-off (Iglewicz & Hoaglin). Robust z is not pulled by the
# outliers themselves, unlike the normal mean/SD z-score.
OUTLIER_CUTOFF = 3.5

# ---- Analysis options -----------------------------------------------------------
# How proportions are calculated (same rule for FACS and scRNA-seq):
#   "populations" : population / (sum of all populations of the panel)    <- recommended
#                   both methods then share the same denominator, which removes
#                   differences caused by debris, dead cells, unlabelled cells, etc.
#   "total"       : FACS -> population / FACS_TOTAL_GATE count
#                   seq  -> population / all annotated cells
PROPORTION_DENOMINATOR = "populations"

MIN_SAMPLES_PER_GROUP = 3     # fewer samples than this -> test is skipped (reported as NaN)
                              # (exception: a single FACS pool is compared as ONE value, see section A)
MIN_CELLS_FOR_RNA_MEAN = 10   # a population with fewer cells in a sample gets no RNA mean
ALPHA = 0.05                  # significance threshold applied to FDR-corrected p-values
BOOTSTRAP_ROUNDS = 2000       # resamplings for the 95% interval of the method difference (section A)
MAX_CELLS_PER_VIOLIN = 5000   # violin plots (F6) are drawn from at most this many random cells
                              # per population (keeps big exports fast; the median uses all cells)

# ---- Plot colors --------------------------------------------------------------
# One color per population (same order as POPULATIONS), plus one per method.
POPULATION_COLORS = ["#2a78d6", "#7fb0ea", "#1baf7a", "#eb6834", "#8a63d2"]
METHOD_COLORS = {"FACS": "#52514e", "scRNA-seq": "#b5b3ad"}
ATLAS_COLOR = "#c2477a"       # the reference atlas in the section E plots
PLOT_DPI = 200                # image resolution (dots per inch) of the saved plots


# =============================================================================
# 1. LOADING THE DATA
# =============================================================================
# These functions read the three main input files, check that the expected columns
# are there, and link every sample to its method and condition. (The optional
# reference atlas is read separately, in section 4b.)

def get_column(table, column, file_label):
    """Return one column of a table; if it does not exist, stop and list the real column names."""
    # If the column name is wrong (typo, different gate name in FlowJo, ...)
    # stop the program with a helpful message instead of a cryptic KeyError.
    if column not in table.columns:
        # Put every real column name on its own line so they are easy to compare
        available = "\n    ".join(map(str, table.columns))
        # A missing FACS column often means the wrong panel was chosen
        hint = (f"\n-> Or check that the FACS panel is right: now '{ACTIVE_PANEL}' "
                f"({PANELS[ACTIVE_PANEL]['description']}).") if file_label == "the FACS file" else ""
        raise SystemExit(
            f"\nColumn not found in {file_label}:\n    '{column}'\n"
            f"Columns that ARE in the file:\n    {available}\n"
            "-> Fix the matching name in the SETTINGS section of this script." + hint
        )
    # The column exists: hand it back
    return table[column]


def normalize_method(value):
    """Turn what is written in the method column ('facs', 'Flex', 'scRNA-seq', ...) into FACS or scRNA-seq."""
    key = "".join(ch for ch in str(value).lower() if ch not in " -_")
    return METHOD_ALIASES.get(key)


def load_metadata(path):
    """
    Read the sample sheet. Everything is read as text so IDs like '01' stay intact.
    Returns one row per sample with a column "method" (FACS or scRNA-seq).
    """
    # dtype=str keeps sample IDs as text (otherwise '01' would become the number 1)
    metadata = pd.read_csv(path, dtype=str)
    get_column(metadata, SAMPLE_COLUMN, path.name)

    if METHOD_COLUMN not in metadata.columns:
        # Older sheet without a method column: every sample was measured by BOTH
        # methods, so list each sample once for scRNA-seq and once for FACS
        get_column(metadata, FACS_NAME_COLUMN, path.name)
        print(f"Note: {path.name} has no '{METHOD_COLUMN}' column - every sample is taken to be "
              "measured by both FACS and scRNA-seq.")
        metadata = pd.concat([metadata.assign(method=m) for m in METHODS], ignore_index=True)
    else:
        # Translate the written method names; stop if one is not recognised
        methods = metadata[METHOD_COLUMN].map(normalize_method)
        unknown = sorted(set(metadata.loc[methods.isna(), METHOD_COLUMN].astype(str)))
        if unknown:
            raise SystemExit(f"Unknown value(s) in the '{METHOD_COLUMN}' column of {path.name}: {unknown}\n"
                             "Use 'FACS' or 'scRNA-seq' for every sample.")
        metadata = metadata.drop(columns=METHOD_COLUMN).assign(method=methods.values)

    # FACS samples need their FlowJo name, otherwise they cannot be found in the export
    is_facs = metadata["method"] == "FACS"
    if is_facs.any():
        names = get_column(metadata, FACS_NAME_COLUMN, path.name)
        missing = sorted(metadata.loc[is_facs & names.isna(), SAMPLE_COLUMN])
        if missing:
            raise SystemExit(f"These FACS samples have no '{FACS_NAME_COLUMN}' in {path.name}: {missing}")
    # The same sample ID twice for one method would mix up two samples
    duplicated = metadata.duplicated([SAMPLE_COLUMN, "method"])
    if duplicated.any():
        raise SystemExit(f"Sample IDs listed twice for the same method in {path.name}: "
                         f"{sorted(set(metadata.loc[duplicated, SAMPLE_COLUMN]))}")
    return metadata


def samples_of(metadata, method):
    """The rows of the sample sheet measured by one method."""
    return metadata[metadata["method"] == method]


def ask_for_numbers(prompt, n_options, allow_many):
    """Ask the user to type option number(s); ask again until the answer is valid.
    Returns a list of 0-based indexes, or None if the user pressed Enter (= all)."""
    # Keep asking until we get a valid answer
    while True:
        answer = input(prompt).strip()
        # Empty answer (just Enter) means "take all of them", when that is allowed
        if answer == "" and allow_many:
            return None
        # Turn "1,3" into [0, 2] (the user counts from 1, Python counts from 0)
        try:
            picks = [int(x) - 1 for x in answer.split(",")]
        except ValueError:
            picks = []   # something that is not a number was typed
        # Accept only numbers inside the range, and only one if a single choice is wanted
        if picks and all(0 <= p < n_options for p in picks) and (allow_many or len(picks) == 1):
            return picks
        # Otherwise explain what is expected and loop again
        print(f"  Please type {'numbers' if allow_many else 'one number'} between 1 and {n_options}.")


def condition_candidates(metadata):
    """Every metadata column that could hold the conditions (all except the ID and method columns)."""
    return [c for c in metadata.columns if c not in (SAMPLE_COLUMN, FACS_NAME_COLUMN, "method", DONOR_COLUMN)]


def choose_conditions(metadata, column=None, seq_values=None, facs_values=None):
    """
    Decide which metadata column is the 'condition', and which of its values to
    analyse for the scRNA-seq samples and for the FACS samples (they can differ,
    e.g. scRNA-seq 'Disease' vs FACS 'Healthy pool').
    Asks the user in the terminal for everything not given on the command line.
    Returns (column, {"scRNA-seq": [...], "FACS": [...]}).
    """
    candidates = condition_candidates(metadata)
    if not candidates:
        raise SystemExit("The metadata file has no column describing conditions.")

    # --- which column? ---
    if column is None:
        if len(candidates) == 1:
            # Only one possibility: take it without asking
            column = candidates[0]
        else:
            # Show a numbered list of columns and let the user pick one
            print("\nWhich metadata column holds the conditions?")
            for i, name in enumerate(candidates, start=1):
                print(f"  {i}. {name}")
            pick = ask_for_numbers("Type a number: ", len(candidates), allow_many=False)
            column = candidates[pick[0]]
    # Check the column really exists (matters when it was typed on the command line)
    get_column(metadata, column, "the metadata file")

    # --- which values, separately for each method? ---
    chosen = {}
    for method, values in [("scRNA-seq", seq_values), ("FACS", facs_values)]:
        # Count how many samples of this method each condition has (keeps the order of the file)
        counts = samples_of(metadata, method)[column].value_counts(sort=False)
        available = list(counts.index)
        if not available:
            raise SystemExit(f"The sample sheet has no {method} samples with a value in '{column}'.")
        if values is None:
            # Show the conditions with their sample counts and let the user choose
            print(f"\n{method} samples - conditions found in '{column}':")
            for i, name in enumerate(available, start=1):
                print(f"  {i}. {name}  ({counts[name]} samples)")
            picks = ask_for_numbers(f"Type the numbers of the {method} conditions to include, "
                                    "separated by commas (Enter = all): ", len(available), allow_many=True)
            # None = Enter was pressed = use every condition
            values = available if picks is None else [available[p] for p in picks]
        else:
            # Conditions were given on the command line: check they all exist for this method
            missing = [v for v in values if v not in available]
            if missing:
                raise SystemExit(f"Conditions {missing} not found among the {method} samples in "
                                 f"'{column}'. Available: {available}")
        chosen[method] = list(values)

    print(f"\nUsing column '{column}': scRNA-seq = {', '.join(chosen['scRNA-seq'])} | "
          f"FACS = {', '.join(chosen['FACS'])}")
    return column, chosen


def all_conditions(conditions):
    """The chosen conditions of both methods in one list, without repeats (scRNA-seq first)."""
    return list(dict.fromkeys(conditions["scRNA-seq"] + conditions["FACS"]))


def load_facs(path, metadata):
    """
    Read the FlowJo Excel export and return one row per FACS sample (index = sample_id).
    metadata = the FACS rows of the sample sheet. Rows FlowJo adds itself
    (e.g. 'Mean', 'SD') are dropped because they are not listed in the metadata.
    """
    # Read the Excel sheet into a table
    facs = pd.read_excel(path, sheet_name=FACS_SHEET)
    # FlowJo puts the sample (.fcs) name in the first column: give it a known name
    facs = facs.rename(columns={facs.columns[0]: FACS_NAME_COLUMN})
    facs[FACS_NAME_COLUMN] = facs[FACS_NAME_COLUMN].astype(str)

    # Dictionary: FlowJo sample name -> sample_id (taken from the metadata)
    name_to_id = dict(zip(metadata[FACS_NAME_COLUMN], metadata[SAMPLE_COLUMN]))
    # Warn about samples listed in the metadata that FlowJo did not export
    missing = sorted(set(name_to_id) - set(facs[FACS_NAME_COLUMN]))
    if missing:
        print(f"WARNING: these samples are in the metadata but not in the FACS file: {missing}")

    # Keep only rows that are real samples from the metadata (drops 'Mean', 'SD', ...)
    facs = facs[facs[FACS_NAME_COLUMN].isin(name_to_id)].copy()
    # Add the sample_id and use it as the row label, so FACS and seq can be matched
    facs[SAMPLE_COLUMN] = facs[FACS_NAME_COLUMN].map(name_to_id)
    return facs.set_index(SAMPLE_COLUMN)


def load_seq(path, metadata):
    """
    Read the per-cell scRNA-seq table and add a 'population' column by
    translating each cell_type label using POPULATIONS[...]["seq_labels"].
    Labels not listed anywhere become 'Other'.
    metadata = the scRNA-seq rows of the sample sheet.
    """
    # Read only the columns that are used (the sample, the label and the marker genes):
    # an export with many extra columns then loads much faster and uses less memory.
    # Sample IDs stay text (same as in the metadata).
    needed = {SAMPLE_COLUMN, CELL_TYPE_COLUMN, *(m["gene"] for m in MARKERS.values())}
    cells = pd.read_csv(path, dtype={SAMPLE_COLUMN: str}, usecols=lambda c: c in needed)
    # Check that the sample, cell type and the four gene columns are all there
    get_column(cells, SAMPLE_COLUMN, path.name)
    get_column(cells, CELL_TYPE_COLUMN, path.name)
    for marker in MARKERS.values():
        get_column(cells, marker["gene"], path.name)

    # Dictionary: notebook label -> population, e.g. "T cell" -> "Other immune"
    label_to_population = {
        label: population
        for population, info in POPULATIONS.items()
        for label in info["seq_labels"]
    }
    # Give each cell its population; labels not in SETTINGS become "Other"
    cells["population"] = cells[CELL_TYPE_COLUMN].map(label_to_population).fillna("Other")

    # Tell the user which labels were not assigned (maybe one was forgotten in SETTINGS)
    unmapped = sorted(cells.loc[cells["population"] == "Other", CELL_TYPE_COLUMN].unique())
    if unmapped:
        print(f"Note: these cell_type labels are not in any population (counted as 'Other'): {unmapped}")

    # Warn about samples in the metadata that have no cells in the scRNA-seq file
    missing = sorted(set(metadata[SAMPLE_COLUMN]) - set(cells[SAMPLE_COLUMN]))
    if missing:
        print(f"WARNING: these samples are in the metadata but not in the scRNA-seq file: {missing}")
    # Keep only cells from the samples (and conditions) chosen for this analysis
    return cells[cells[SAMPLE_COLUMN].isin(metadata[SAMPLE_COLUMN])].copy()


# =============================================================================
# 2. BUILDING THE COMPARISON TABLES
# =============================================================================
# These functions turn the raw inputs into "long" tables (one measurement per
# row) that the statistics and plots below can use directly.
# FACS and scRNA-seq measured DIFFERENT samples, so each row belongs to one
# method; FACS and scRNA-seq rows are matched by population (and condition),
# never by sample.

def add_conditions(table, metadata, condition_column):
    """Add the condition of each sample as a column (matched on sample ID AND method)."""
    return table.merge(metadata[[SAMPLE_COLUMN, "method", condition_column]], on=[SAMPLE_COLUMN, "method"])


def build_proportions(facs, cells, metadata, condition_column):
    """
    Long table: one row per sample x population, with the method, the cell count
    and the proportion. FACS samples come from the FlowJo table, scRNA-seq samples
    from the per-cell table. This is the table used for sections A, B and E.
    """
    populations = list(POPULATIONS)

    # FACS counts per sample and population
    # (one column per population, taken from that population's gate "Count" column)
    facs_counts = pd.DataFrame({
        pop: get_column(facs, COUNT_COLUMN.format(gate=info["gate"]), "the FACS file")
        for pop, info in POPULATIONS.items()
    })
    # scRNA-seq counts per sample and population (cells per sample per population)
    seq_counts = pd.crosstab(cells[SAMPLE_COLUMN], cells["population"])
    # Make sure every population has a column, even if it had 0 cells
    seq_counts = seq_counts.reindex(columns=populations + ["Other"], fill_value=0)

    # Choose what each count is divided by (see PROPORTION_DENOMINATOR in SETTINGS)
    if PROPORTION_DENOMINATOR == "populations":
        # sum of all populations of the panel, per sample (= all live cells for FACS)
        facs_total = facs_counts[populations].sum(axis=1)
        seq_total = seq_counts[populations].sum(axis=1)
    elif PROPORTION_DENOMINATOR == "total":
        # all live single cells (FACS) / all cells including "Other" (seq)
        facs_total = get_column(facs, COUNT_COLUMN.format(gate=FACS_TOTAL_GATE), "the FACS file")
        seq_total = seq_counts.sum(axis=1)
    else:
        raise SystemExit("PROPORTION_DENOMINATOR must be 'populations' or 'total'.")

    parts = []
    for method, counts, total in [("scRNA-seq", seq_counts, seq_total), ("FACS", facs_counts, facs_total)]:
        # From a wide table (one column per population) to a long table (one row per population)
        long = (counts[populations].rename_axis(index=SAMPLE_COLUMN, columns=None).reset_index()
                .melt(id_vars=SAMPLE_COLUMN, var_name="population", value_name="count")
                .sort_values(SAMPLE_COLUMN, kind="stable"))
        long.insert(1, "method", method)
        # proportion = count / total of its sample (NaN if the sample has no cells at all)
        long["total"] = long[SAMPLE_COLUMN].map(total)
        long["proportion"] = long["count"] / long["total"].where(long["total"] > 0)
        parts.append(long)
    # Add the condition of each sample (e.g. Healthy pool / Disease) as a column
    return add_conditions(pd.concat(parts, ignore_index=True), metadata, condition_column)


def donor_of_each_row(table, metadata):
    """The person (DONOR_COLUMN) of the sample in each row; a sample without one is its own person."""
    if DONOR_COLUMN not in metadata.columns:
        return table[SAMPLE_COLUMN]
    lookup = metadata.set_index([SAMPLE_COLUMN, "method"])[DONOR_COLUMN]
    keys = pd.MultiIndex.from_frame(table[[SAMPLE_COLUMN, "method"]])
    donors = pd.Series(lookup.reindex(keys).values, index=table.index)
    return donors.fillna(table[SAMPLE_COLUMN])


def combine_per_donor(table, metadata, keys, aggregations):
    """
    Combine the samples of the same person into one row per person (and method
    and keys, e.g. population and condition). Samples of one person are not
    independent, so the statistics must count every person once.
    aggregations = how each column is combined, e.g. {"count": ("count", "sum")}.
    Afterwards the sample_id column holds the PERSON, and n_samples_combined says
    how many samples went into it. Without a donor column nothing changes (1 per row).
    """
    table = table.assign(_sample=table[SAMPLE_COLUMN], **{SAMPLE_COLUMN: donor_of_each_row(table, metadata)})
    return (table.groupby([SAMPLE_COLUMN, "method", *keys], sort=False, dropna=False)
            .agg(**aggregations, n_samples_combined=("_sample", "nunique")).reset_index())


def combine_proportions_per_donor(sample_proportions, metadata, condition_column):
    """
    The proportions per PERSON: the cells of all samples of a person are added up
    (counts and totals), then divided again. A person with two samples therefore
    counts once, like one bigger sample.
    """
    table = combine_per_donor(sample_proportions, metadata, ["population", condition_column],
                              {"count": ("count", "sum"), "total": ("total", "sum")})
    table["proportion"] = table["count"] / table["total"].where(table["total"] > 0)
    for method in METHODS:
        part = sample_proportions[sample_proportions["method"] == method]
        n_samples = part[SAMPLE_COLUMN].nunique()
        n_people = table.loc[table["method"] == method, SAMPLE_COLUMN].nunique()
        if n_people < n_samples:
            print(f"Combined {n_samples} {method} samples into {n_people} people ('{DONOR_COLUMN}' column).")
    return table


def build_facs_markers(facs, metadata, condition_column):
    """
    FACS protein: one row per FACS sample x population x marker, with the median
    intensity of the marker inside the population gate. Used in sections C, E and F.
    """
    parts = []
    for pop, pop_info in POPULATIONS.items():
        for marker, marker_info in MARKERS.items():
            # The "Median (marker channel)" column of this population's gate
            column = MEDIAN_COLUMN.format(gate=pop_info["gate"], parameter=marker_info["facs_parameter"])
            values = get_column(facs, column, "the FACS file")
            parts.append(pd.DataFrame({SAMPLE_COLUMN: facs.index, "method": "FACS", "population": pop,
                                       "marker": marker, "gene": marker_info["gene"],
                                       "protein_median": values.values}))
    return add_conditions(pd.concat(parts, ignore_index=True), metadata, condition_column)


def build_rna_markers(cells, metadata, condition_column):
    """
    scRNA-seq RNA: one row per scRNA-seq sample x population x marker, with the mean
    log-normalized expression of the marker's gene in that population's cells.
    A population with fewer than MIN_CELLS_FOR_RNA_MEAN cells in a sample gets no mean.
    """
    gene_to_marker = {info["gene"]: marker for marker, info in MARKERS.items()}
    known = cells[cells["population"].isin(list(POPULATIONS))]
    grouped = known.groupby([SAMPLE_COLUMN, "population"])
    # Mean expression of each marker gene, per sample and population
    means = grouped[list(gene_to_marker)].mean()
    n_cells = grouped.size()
    # A mean from very few cells is unreliable, so require a minimum
    means.loc[n_cells < MIN_CELLS_FOR_RNA_MEAN] = np.nan
    table = (means.assign(n_cells_rna=n_cells).reset_index()
             .melt(id_vars=[SAMPLE_COLUMN, "population", "n_cells_rna"], var_name="gene", value_name="rna_mean"))
    table.insert(1, "method", "scRNA-seq")
    table.insert(3, "marker", table["gene"].map(gene_to_marker))
    return add_conditions(table, metadata, condition_column)


def build_protein_vs_rna(facs_markers, rna_markers):
    """
    Section C: one row per population x marker. The samples differ between the
    methods, so the TYPICAL level of each population is compared:
      protein = median over the FACS samples of the gate's median intensity
      rna     = median over the scRNA-seq samples of the per-sample mean expression
    plus the spread (25th-75th percentile) and the number of samples behind each.
    """
    def summarise(table, value, name, n_name):
        grouped = table.groupby(["population", "marker"])[value]
        return pd.DataFrame({f"{name}_median": grouped.median(), f"{name}_q25": grouped.quantile(0.25),
                             f"{name}_q75": grouped.quantile(0.75), n_name: grouped.count()})

    # Every population x marker gets a row, even when one method has no value for it
    index = pd.MultiIndex.from_product([list(POPULATIONS), list(MARKERS)], names=["population", "marker"])
    table = pd.concat([summarise(facs_markers, "protein_median", "protein", "n_facs_samples"),
                       summarise(rna_markers, "rna_mean", "rna", "n_seq_samples")], axis=1).reindex(index)
    table[["n_facs_samples", "n_seq_samples"]] = table[["n_facs_samples", "n_seq_samples"]].fillna(0).astype(int)
    table.insert(0, "gene", [MARKERS[m]["gene"] for m in index.get_level_values("marker")])
    return table.reset_index()


def build_imaging(facs, metadata, condition_column):
    """One row per FACS sample x population x imaging feature (FACS median value)."""
    parts = []
    # For every population and every S8 imaging feature ...
    for pop, info in POPULATIONS.items():
        for feature in IMAGING_FEATURES:
            # ... take the "Median (feature)" column of that population's gate, one row per sample
            column = MEDIAN_COLUMN.format(gate=info["gate"], parameter=feature)
            values = get_column(facs, column, "the FACS file")
            parts.append(pd.DataFrame({SAMPLE_COLUMN: facs.index, "population": pop,
                                       "feature": feature, "value": values.values}))
    table = pd.concat(parts, ignore_index=True)
    # Add the condition of each sample as a column (metadata = the FACS samples)
    return table.merge(metadata[[SAMPLE_COLUMN, condition_column]], on=SAMPLE_COLUMN)


# =============================================================================
# 3. STATISTICS
# =============================================================================
# Non-parametric tests are used throughout: urine studies usually have few
# samples, and proportions/intensities are rarely normally distributed.

def fdr_correct(p_values):
    """Benjamini-Hochberg FDR correction. NaN p-values (skipped tests) stay NaN."""
    p = np.asarray(p_values, dtype=float)
    q = np.full_like(p, np.nan)        # start with all NaN
    ok = ~np.isnan(p)                  # which tests actually produced a p-value
    n = ok.sum()                       # number of real tests
    if n == 0:
        return q
    # Sort the p-values from small to large and scale each by n / rank
    order = np.argsort(p[ok])
    ranked = p[ok][order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]   # keep q-values monotonic
    # Put the q-values back in the original order (and cap them at 1)
    corrected = np.empty(n)
    corrected[order] = np.minimum(ranked, 1)
    q[ok] = corrected
    return q


def correlate(x, y):
    """Spearman (rank) and Pearson (linear) correlation between two paired lists."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    # Drop pairs where either value is missing
    keep = ~(np.isnan(x) | np.isnan(y))
    x, y = x[keep], y[keep]
    # Default result: NaN (= test not done)
    result = {"n": len(x), "spearman_r": np.nan, "spearman_p": np.nan,
              "pearson_r": np.nan, "pearson_p": np.nan}
    # Only test if there are enough pairs and the values are not all identical
    if len(x) >= MIN_SAMPLES_PER_GROUP and np.std(x) > 0 and np.std(y) > 0:
        result["spearman_r"], result["spearman_p"] = stats.spearmanr(x, y)
        # scipy's Spearman p-value is an approximation that is far too small for
        # few pairs (e.g. p = 0 for 4 points in perfect order, where the true
        # smallest possible p is 2/24 = 0.083). Use the exact p-value instead.
        if len(x) <= EXACT_SPEARMAN_MAX_N:
            result["spearman_p"] = exact_spearman_p(x, y)
        result["pearson_r"], result["pearson_p"] = stats.pearsonr(x, y)
    return result


EXACT_SPEARMAN_MAX_N = 8   # up to 8 pairs (8! = 40320 orderings) the exact p-value is calculated


@functools.lru_cache(maxsize=None)
def all_orderings(n):
    """Every possible ordering of n items, as rows of a table (cached: made once per n)."""
    return np.array(list(itertools.permutations(range(n))))


def exact_spearman_p(x, y):
    """
    Exact two-sided p-value of the Spearman correlation: the share of ALL possible
    ways of pairing the y values with the x values that give a correlation at least
    as strong as the observed one (a permutation test, done completely).
    """
    rx, ry = stats.rankdata(x), stats.rankdata(y)
    rx, ry = rx - rx.mean(), ry - ry.mean()
    shuffled = ry[all_orderings(len(ry))]                 # one row per possible pairing
    r_all = shuffled @ rx / (np.sqrt(np.sum(rx ** 2)) * np.sqrt(np.sum(ry ** 2)))
    r_observed = np.sum(rx * ry) / (np.sqrt(np.sum(rx ** 2)) * np.sqrt(np.sum(ry ** 2)))
    return float(np.mean(np.abs(r_all) >= abs(r_observed) - 1e-12))


def drop_nan(values):
    """The values as an array of numbers, without the missing ones (NaN)."""
    values = np.asarray(values, dtype=float)
    return values[~np.isnan(values)]


def compare_groups(values_by_group):
    """
    Compare a measurement between groups (conditions or populations).
      2 groups  -> Mann-Whitney U
      3+ groups -> Kruskal-Wallis (overall) + every pair with Mann-Whitney U
    Groups with fewer than MIN_SAMPLES_PER_GROUP values are left out.
    """
    # Remove missing values from every group, then drop groups that are too small to test
    groups = {name: drop_nan(v) for name, v in values_by_group.items()}
    groups = {name: v for name, v in groups.items() if len(v) >= MIN_SAMPLES_PER_GROUP}
    names = list(groups)
    results = []
    # Need at least two groups to compare anything
    if len(names) < 2:
        return results

    # 3 or more groups: one overall test first
    if len(names) > 2:
        statistic, p = stats.kruskal(*groups.values())
        results.append({"test": "Kruskal-Wallis", "comparison": " vs ".join(names),
                        "statistic": statistic, "p_value": p})
    # Then compare every pair of groups (with 2 groups this is the only test)
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            a, b = names[i], names[j]
            statistic, p = stats.mannwhitneyu(groups[a], groups[b], alternative="two-sided")
            results.append({"test": "Mann-Whitney U", "comparison": f"{a} vs {b}",
                            "statistic": statistic, "p_value": p,
                            "median_group1": np.median(groups[a]),
                            "median_group2": np.median(groups[b])})
    return results


def compare_to_reference(values, reference):
    """
    UNPAIRED comparison of a group of values with a reference group (different samples).
      both have >= MIN_SAMPLES_PER_GROUP values -> Mann-Whitney U test
      only `values` has enough                  -> one-sample Wilcoxon test: are the
                                                   values centred on the reference
                                                   (its median, e.g. ONE pooled sample)?
      only `reference` has enough               -> the same the other way round
    Returns (test name, statistic, p-value); ("", NaN, NaN) when no test is possible.
    The test names are generic; the callers translate them for their table.
    """
    values, reference = drop_nan(values), drop_nan(reference)
    enough_values = len(values) >= MIN_SAMPLES_PER_GROUP
    enough_reference = len(reference) >= MIN_SAMPLES_PER_GROUP
    # Every value the same (e.g. a population that is 0 everywhere): nothing to test
    if len(values) and len(reference) and np.ptp(np.r_[values, reference]) == 0:
        return "", np.nan, np.nan
    if enough_values and enough_reference:
        statistic, p = stats.mannwhitneyu(values, reference, alternative="two-sided")
        return "Mann-Whitney U", statistic, p
    # One side is too small for a group test: use it as a single fixed value
    if enough_values and len(reference):
        group, target, test = values, np.median(reference), "Wilcoxon vs reference value"
    elif enough_reference and len(values):
        group, target, test = reference, np.median(values), "Wilcoxon of reference vs value"
    else:
        return "", np.nan, np.nan
    # Every value equal to the target: there is no difference to test
    if not np.any(group != target):
        return "", np.nan, np.nan
    statistic, p = stats.wilcoxon(group - target)
    return test, statistic, p


def bootstrap_median_difference(values, reference, rounds=None):
    """
    95% interval of median(values) - median(reference), by bootstrap: both groups
    are resampled (drawn again with replacement) many times, and the middle 95% of
    the resulting differences is kept. A group of ONE value (e.g. a pool) stays
    fixed. Only rough with few samples. NaN if a group is empty.
    """
    rounds = rounds or BOOTSTRAP_ROUNDS
    values, reference = drop_nan(values), drop_nan(reference)
    if len(values) == 0 or len(reference) == 0 or max(len(values), len(reference)) < 2:
        return np.nan, np.nan
    rng = np.random.default_rng(0)        # fixed seed: the same interval every run

    def resampled_medians(x):
        # rounds x len(x) random picks from x -> one median per round
        return np.median(x[rng.integers(0, len(x), (rounds, len(x)))], axis=1)

    low, high = np.percentile(resampled_medians(values) - resampled_medians(reference), [2.5, 97.5])
    return low, high


# How the generic test names of compare_to_reference() are written in the section A table
METHOD_TEST_NAMES = {
    "Wilcoxon vs reference value": "Wilcoxon: scRNA-seq vs single FACS value",
    "Wilcoxon of reference vs value": "Wilcoxon: FACS vs single scRNA-seq value",
}


def method_comparison_stats(proportions, condition_column, conditions):
    """
    Section A: do scRNA-seq and FACS give the same proportions?
    The two methods measured DIFFERENT samples, so they are compared as groups:
    every chosen scRNA-seq condition against every chosen FACS condition, for each
    population (e.g. "Disease (scRNA-seq) vs Healthy pool (FACS)").
      seq_minus_facs  : median scRNA-seq - median FACS (positive = scRNA-seq finds more)
      seq_to_facs_ratio: median scRNA-seq / median FACS (2 = twice as much in scRNA-seq)
      95% interval of the difference by bootstrap
      test            : see compare_to_reference() (the FACS group is the reference)
    FDR correction over the whole table.
    """
    rows = []
    for seq_group in conditions["scRNA-seq"]:
        for facs_group in conditions["FACS"]:
            for pop in POPULATIONS:
                part = proportions[proportions["population"] == pop]
                # The proportions of this population in the two groups of samples
                seq = drop_nan(part.loc[(part["method"] == "scRNA-seq") &
                                        (part[condition_column] == seq_group), "proportion"])
                facs = drop_nan(part.loc[(part["method"] == "FACS") &
                                         (part[condition_column] == facs_group), "proportion"])
                median_seq = np.median(seq) if len(seq) else np.nan
                median_facs = np.median(facs) if len(facs) else np.nan
                test, statistic, p = compare_to_reference(seq, facs)
                # The interval only where a test was possible (enough samples)
                low, high = bootstrap_median_difference(seq, facs) if test else (np.nan, np.nan)
                rows.append({
                    "population": pop,
                    "comparison": f"{seq_group} (scRNA-seq) vs {facs_group} (FACS)",
                    "seq_group": seq_group, "facs_group": facs_group,
                    "n_seq": len(seq), "n_facs": len(facs),
                    "median_seq": median_seq, "median_facs": median_facs,
                    "seq_minus_facs": median_seq - median_facs,
                    "seq_minus_facs_ci95_low": low, "seq_minus_facs_ci95_high": high,
                    "seq_to_facs_ratio": median_seq / median_facs if median_facs > 0 else np.nan,
                    "test": METHOD_TEST_NAMES.get(test, test), "statistic": statistic, "p_value": p,
                })
    table = pd.DataFrame(rows)
    # Correct for all tests in the table, then flag the significant ones
    table["q_value_fdr"] = fdr_correct(table["p_value"])
    table["significant"] = table["q_value_fdr"] < ALPHA
    return table


def composition_overview(proportions, condition_column, conditions, atlas_props=None):
    """
    Section A: the composition of every group side by side - each FACS group (e.g.
    a pool), each scRNA-seq group and, if given, the reference atlas. One row per
    group x population with the number of people/samples behind it and the mean,
    median, 25th-75th percentile, min and max of the proportion.
    For a single sample (a pool) all of these are simply its own proportion.
    The means of a group add up to 1, so they can be drawn as one stacked bar (A4).
    """
    def describe(values):
        values = drop_nan(values)
        if not len(values):
            return dict.fromkeys(["mean", "median", "q25", "q75", "min", "max"], np.nan)
        return {"mean": values.mean(), "median": np.median(values),
                "q25": np.percentile(values, 25), "q75": np.percentile(values, 75),
                "min": values.min(), "max": values.max()}

    rows = []
    # FACS groups first (the reference, e.g. the pool), then the scRNA-seq groups
    for method in ["FACS", "scRNA-seq"]:
        for group in conditions[method]:
            part = proportions[(proportions["method"] == method) & (proportions[condition_column] == group)]
            for pop in POPULATIONS:
                values = part.loc[part["population"] == pop, "proportion"]
                rows.append({"source": method, "group": group, "label": f"{group} ({method})",
                             "population": pop, "n": part[SAMPLE_COLUMN].nunique(), **describe(values)})
    if atlas_props is not None:
        pooled = atlas_pooled_proportions(atlas_props)
        for pop in POPULATIONS:
            donors = atlas_props.loc[atlas_props["population"] == pop, "proportion"]
            rows.append({"source": "Atlas", "group": "Reference atlas", "label": "Reference atlas",
                         "population": pop, "n": atlas_props["donor"].nunique(), **describe(donors),
                         "atlas_all_cells": pooled[pop]})
    return pd.DataFrame(rows)


def composition_overview_wide(overview):
    """The overview as a compact table: one row per population, one column per group (mean %)."""
    wide = overview.pivot_table(index="population", columns="label", values="mean", sort=False) * 100
    # (order taken from the table itself, so it also works for a table read back from its CSV)
    wide = wide.reindex(index=list(dict.fromkeys(overview["population"])),
                        columns=list(dict.fromkeys(overview["label"])))
    wide.columns = [f"{label} (n={n}), %" for label, n in
                    overview.drop_duplicates("label").set_index("label")["n"].reindex(wide.columns).items()]
    return wide.reset_index()


def condition_stats(table, value_column, split_columns, condition_column, conditions):
    """
    Section B (and D): for every combination of split_columns (e.g. population x
    method) compare value_column between the selected conditions. FDR correction
    is applied over the whole table. Conditions are always compared WITHIN one
    method; a condition without (enough) samples of that method is left out.
    """
    rows = []
    # Go through each group, e.g. ("Neutrophils", "FACS"), ("Neutrophils", "scRNA-seq"), ...
    for keys, part in table.groupby(split_columns, sort=False):
        keys = keys if isinstance(keys, tuple) else (keys,)
        # Collect the values of each condition in this group
        values = {c: part.loc[part[condition_column] == c, value_column].values for c in conditions}
        # Run the tests and label each result with the group it belongs to
        for result in compare_groups(values):
            rows.append({**dict(zip(split_columns, keys)), **result})
    result = pd.DataFrame(rows)
    if not result.empty:
        # Correct for all tests in the table, then flag the significant ones
        result["q_value_fdr"] = fdr_correct(result["p_value"])
        result["significant"] = result["q_value_fdr"] < ALPHA
    return result


def scale_0_1(values):
    """Rescale a series from 0 (lowest value) to 1 (highest); all NaN if every value is equal."""
    low, high = values.min(), values.max()
    return (values - low) / (high - low) if high > low else values * np.nan


def protein_rna_stats(protein_rna):
    """
    Section C: do the populations come out in the same order for protein (FACS)
    and RNA (scRNA-seq)? The samples differ between the methods, so this compares
    the typical level of each POPULATION (from build_protein_vs_rna), not samples.
      one row per marker   : correlation across the populations (few points:
                             read r as "same pattern?", the p-value is weak)
      "All markers (scaled)": every marker first scaled 0-1 across the populations
                             in each method (the units differ), then all
                             population x marker points together - one
                             better-powered answer to "do the patterns match?"
    """
    rows = []
    for marker, part in protein_rna.groupby("marker", sort=False):
        rows.append({"marker": marker, **correlate(part["protein_median"], part["rna_median"])})
    scaled = protein_rna.groupby("marker", sort=False)[["protein_median", "rna_median"]].transform(scale_0_1)
    rows.append({"marker": "All markers (scaled)", **correlate(scaled["protein_median"], scaled["rna_median"])})
    table = pd.DataFrame(rows)
    # Correct the p-values for all correlations in the table
    table["spearman_q_fdr"] = fdr_correct(table["spearman_p"])
    table["pearson_q_fdr"] = fdr_correct(table["pearson_p"])
    return table


def imaging_between_populations(imaging):
    """Section D: do the imaging features differ between the populations?"""
    rows = []
    for feature, part in imaging.groupby("feature", sort=False):
        # For this feature, collect the values of each population
        values = {pop: part.loc[part["population"] == pop, "value"].values for pop in POPULATIONS}
        # Compare the populations (Kruskal-Wallis + pairwise Mann-Whitney)
        for result in compare_groups(values):
            rows.append({"feature": feature, **result})
    table = pd.DataFrame(rows)
    if not table.empty:
        # Correct for all tests, then flag the significant ones
        table["q_value_fdr"] = fdr_correct(table["p_value"])
        table["significant"] = table["q_value_fdr"] < ALPHA
    return table


# =============================================================================
# 4. PLOTS
# =============================================================================
# Each function draws one figure and saves it as a .png in the results folder.

# Shared look for all plots: no top/right frame lines, light grey grid, small font
plt.rcParams.update({
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.edgecolor": "#8a8984", "axes.labelcolor": "#0b0b0b",
    "xtick.color": "#52514e", "ytick.color": "#52514e",
    "axes.grid": True, "grid.color": "#e6e5e0", "grid.linewidth": 0.6,
    "axes.axisbelow": True, "font.size": 9, "legend.frameon": False,
})


def population_color(pop):
    # Each population always gets the same color in every plot
    return POPULATION_COLORS[list(POPULATIONS).index(pop) % len(POPULATION_COLORS)]


def save(fig, name):
    # Fit everything inside the figure, write it to the results folder and free the memory
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / name, dpi=PLOT_DPI)
    plt.close(fig)
    print(f"  saved {name}")


def ordered_samples(table, method, condition_column, conditions):
    """The samples of one method, sorted by condition (in the chosen order), then by sample ID."""
    samples = table.loc[table["method"] == method, [SAMPLE_COLUMN, condition_column]].drop_duplicates()
    rank = samples[condition_column].map({c: i for i, c in enumerate(conditions)})
    return samples.assign(_rank=rank).sort_values(["_rank", SAMPLE_COLUMN]).drop(columns="_rank")


def method_label(method):
    """Short method name for crowded axis labels."""
    return "seq" if method == "scRNA-seq" else method


def plot_stacked_proportions(proportions, condition_column, conditions):
    """Stacked bars, one per person/sample: the scRNA-seq people first, then the FACS samples."""
    # Table: rows = (method, sample), columns = populations, values = proportions
    wide = (proportions.pivot_table(index=["method", SAMPLE_COLUMN], columns="population", values="proportion")
            .reindex(columns=list(POPULATIONS)))
    keys, positions, labels, blocks = [], [], [], []
    x = 0
    for method in METHODS:
        order = ordered_samples(proportions, method, condition_column, conditions[method])
        start = x
        for sample, condition in zip(order[SAMPLE_COLUMN], order[condition_column]):
            keys.append((method, sample))
            positions.append(x)
            labels.append(f"{sample} · {condition}")
            x += 1
        blocks.append((method, start, x - 1))
        x += 1                                  # a gap between the two methods
    values = wide.reindex(pd.MultiIndex.from_tuples(keys)).fillna(0)

    # Figure gets wider when there are more samples
    fig, ax = plt.subplots(figsize=(max(6, len(keys) * 0.45 + 2.5), 4.8))
    bottom = np.zeros(len(keys))
    # One call per population draws that population's block in every bar (stacked on the previous ones)
    for pop in POPULATIONS:
        ax.bar(positions, values[pop].values, bottom=bottom, width=0.85, color=population_color(pop),
               edgecolor="white", linewidth=1, label=pop)
        bottom += values[pop].values
    # Name the two blocks of bars above the plot
    for method, start, end in blocks:
        if end >= start:
            ax.text((start + end) / 2, 1.02, method, transform=ax.get_xaxis_transform(),
                    ha="center", va="bottom", fontweight="bold")
    ax.set_xticks(positions, labels, rotation=90)
    ax.set_ylabel("Proportion of cells")
    ax.set_title("Cell-type proportions per person (scRNA-seq) and per sample (FACS)", loc="left", pad=18)
    ax.grid(axis="x", visible=False)
    ax.legend(loc="upper left", bbox_to_anchor=(1, 1))
    save(fig, "A1_proportions_stacked.png")


def method_groups(conditions):
    """The boxes of plot A2, left to right: (method, condition), scRNA-seq conditions first."""
    return [(method, c) for method in METHODS for c in conditions[method]]


def plot_method_boxes(proportions, condition_column, conditions):
    """
    A2: one panel per population. One box per scRNA-seq condition and per FACS
    condition; every dot is one sample (a pool is a single dot).
    """
    groups = method_groups(conditions)
    fig, axes = plt.subplots(1, len(POPULATIONS), figsize=(max(3.4, 0.9 * len(groups) + 1.2) * len(POPULATIONS), 4))
    for ax, pop in zip(np.atleast_1d(axes), POPULATIONS):
        part = proportions[proportions["population"] == pop]
        for i, (method, condition) in enumerate(groups):
            values = part.loc[(part["method"] == method) & (part[condition_column] == condition), "proportion"]
            boxplot_with_points(ax, [values], [i], METHOD_COLORS[method], width=0.5)
        ax.set_xticks(range(len(groups)), [f"{c}\n({method_label(m)})" for m, c in groups])
        ax.set_title(pop, loc="left")
        ax.set_ylabel("Proportion of cells")
    # Legend entries for the two method colors
    last = np.atleast_1d(axes)[-1]
    for method in METHODS:
        last.scatter([], [], color=METHOD_COLORS[method], s=22, label=method)
    last.legend(loc="upper left", bbox_to_anchor=(1, 1))
    save(fig, "A2_proportions_by_method.png")


def plot_method_difference(comparison):
    """
    A3: one panel per comparison (scRNA-seq group vs FACS group). Each row is a
    population: the dot is median scRNA-seq - median FACS, the line its bootstrap
    95% interval. Dots right of 0 = scRNA-seq finds more of that population.
    """
    comparisons = list(dict.fromkeys(comparison["comparison"]))
    pops = list(POPULATIONS)
    fig, axes = plt.subplots(1, len(comparisons), figsize=(4.4 * len(comparisons) + 1.2, 0.55 * len(pops) + 1.8),
                             squeeze=False, sharex=True)
    for k, (ax, name) in enumerate(zip(axes[0], comparisons)):
        part = comparison[comparison["comparison"] == name].set_index("population").reindex(pops)
        for i, pop in enumerate(pops):
            row, y = part.loc[pop], len(pops) - 1 - i          # first population at the top
            if np.isnan(row["seq_minus_facs"]):
                continue
            color = population_color(pop)
            if np.isfinite(row["seq_minus_facs_ci95_low"]):
                ax.plot([row["seq_minus_facs_ci95_low"], row["seq_minus_facs_ci95_high"]], [y, y],
                        color=color, linewidth=2, solid_capstyle="round")
            ax.scatter(row["seq_minus_facs"], y, s=50, color=color, edgecolor="white", linewidth=1.2, zorder=3)
            # The FDR q-value just outside the right edge of the panel (never on top of a line)
            q = row["q_value_fdr"]
            ax.text(1.02, y, "no test" if np.isnan(q) else f"q = {q:.2g}", transform=ax.get_yaxis_transform(),
                    ha="left", va="center", fontsize=7, color="#52514e")
        ax.axvline(0, color="#8a8984", linewidth=1, linestyle="--")
        ax.set_ylim(-0.6, len(pops) - 0.4)
        ax.set_yticks(range(len(pops)), pops[::-1] if k == 0 else [])
        ax.set_title(name.replace(" vs ", "\nvs "), loc="left")
        ax.set_xlabel("Median scRNA-seq − median FACS")
    save(fig, "A3_method_difference.png")


def plot_composition_overview(overview):
    """
    A4: one stacked bar per group - the FACS groups (e.g. the pool), the scRNA-seq
    groups and the reference atlas - showing the average proportion of each
    population, with the percentage written in each block.
    """
    labels = list(dict.fromkeys(overview["label"]))
    n_of = overview.drop_duplicates("label").set_index("label")["n"]
    source_of = overview.drop_duplicates("label").set_index("label")["source"]
    unit = {"FACS": "sample", "scRNA-seq": "people", "Atlas": "donors"}
    fig, ax = plt.subplots(figsize=(max(5, 1.5 * len(labels) + 2.5), 4.8))
    bottom = np.zeros(len(labels))
    for pop in POPULATIONS:
        part = overview[overview["population"] == pop].set_index("label")["mean"].reindex(labels).fillna(0).values
        ax.bar(range(len(labels)), part, bottom=bottom, width=0.7, color=population_color(pop),
               edgecolor="white", linewidth=1, label=pop)
        # The percentage inside the block (only when there is room for it)
        for x, (value, start) in enumerate(zip(part, bottom)):
            if value >= 0.04:
                ax.text(x, start + value / 2, f"{100 * value:.0f}%", ha="center", va="center",
                        fontsize=8, color="white", fontweight="bold")
        bottom += part
    names = []
    for label in labels:
        n = int(n_of[label])
        what = "sample" if source_of[label] == "FACS" and n == 1 else ("samples" if source_of[label] == "FACS"
                                                                       else unit[source_of[label]])
        names.append(f"{label.replace(' (', chr(10) + '(')}\nn = {n} {what}")
    ax.set_xticks(range(len(labels)), names)
    ax.set_ylabel("Average proportion of cells")
    ax.set_ylim(0, 1)
    ax.set_title("Composition overview: FACS, scRNA-seq and reference atlas", loc="left")
    ax.grid(axis="x", visible=False)
    ax.legend(loc="upper left", bbox_to_anchor=(1, 1))
    save(fig, "A4_composition_overview.png")


def boxplot_with_points(ax, data_by_x, positions, color, width=0.35):
    """Draw thin box plots with every sample shown as a dot on top (empty groups are skipped)."""
    # Fixed random seed so the dot jitter looks the same every run
    rng = np.random.default_rng(0)
    # Remove missing values; skip empty groups
    groups = [(drop_nan(values), x) for values, x in zip(data_by_x, positions)]
    groups = [(values, x) for values, x in groups if len(values)]
    if not groups:
        return
    # All boxes in ONE call (median and quartiles, no separate outlier markers):
    # much faster than one call per box when a figure has many boxes
    ax.boxplot([values for values, _ in groups], positions=[x for _, x in groups], widths=width,
               showfliers=False, manage_ticks=False,
               boxprops={"color": color, "linewidth": 1.2},
               whiskerprops={"color": color}, capprops={"color": color},
               medianprops={"color": "#0b0b0b", "linewidth": 1.5})
    # The individual samples, shifted a little sideways so they don't overlap (one call for all dots)
    xs = np.concatenate([x + rng.uniform(-width / 4, width / 4, len(values)) for values, x in groups])
    ys = np.concatenate([values for values, _ in groups])
    ax.scatter(xs, ys, s=22, color=color, edgecolor="white", linewidth=1, zorder=3)


def plot_proportions_by_condition(proportions, condition_column, conditions):
    """
    B1: one panel per population: proportions per condition, scRNA-seq and FACS
    side by side. conditions = all chosen conditions of both methods (one list);
    a method without samples in a condition simply has no box there.
    """
    fig, axes = plt.subplots(1, len(POPULATIONS),
                             figsize=(max(4, 1.6 * len(conditions)) * len(POPULATIONS) / 1.3, 4))
    for ax, pop in zip(np.atleast_1d(axes), POPULATIONS):
        part = proportions[proportions["population"] == pop]
        for k, method in enumerate(METHODS):
            # One list of proportions per condition, for this method
            data = [part[(part["method"] == method) & (part[condition_column] == c)]["proportion"]
                    for c in conditions]
            # scRNA-seq box slightly left, FACS box slightly right of each condition
            positions = [i + (k - 0.5) * 0.4 for i in range(len(conditions))]
            boxplot_with_points(ax, data, positions, METHOD_COLORS[method])
            ax.scatter([], [], color=METHOD_COLORS[method], s=22, label=method)  # legend entry
        ax.set_xticks(range(len(conditions)), conditions)
        ax.set_title(pop, loc="left")
        ax.set_ylabel("Proportion of cells")
    np.atleast_1d(axes)[-1].legend(loc="upper left", bbox_to_anchor=(1, 1))
    save(fig, "B1_proportions_by_condition.png")


def plot_protein_vs_rna(protein_rna, rna_stats):
    """
    C1: one panel per marker. Each dot is one population: x = typical RNA of the
    gene (scRNA-seq samples), y = typical protein (FACS samples); the lines show
    the 25th-75th percentile over the samples of each method.
    """
    fig, axes = plt.subplots(1, len(MARKERS), figsize=(3.6 * len(MARKERS), 3.8))
    for ax, (marker, info) in zip(np.atleast_1d(axes), MARKERS.items()):
        part = protein_rna[protein_rna["marker"] == marker]
        for pop in POPULATIONS:
            row = part[part["population"] == pop].iloc[0]
            if np.isnan(row["rna_median"]) or np.isnan(row["protein_median"]):
                continue
            color = population_color(pop)
            # Spread over the samples of each method (25th-75th percentile)
            ax.plot([row["rna_q25"], row["rna_q75"]], [row["protein_median"]] * 2, color=color, linewidth=1.2)
            ax.plot([row["rna_median"]] * 2, [row["protein_q25"], row["protein_q75"]], color=color, linewidth=1.2)
            ax.scatter(row["rna_median"], row["protein_median"], s=50, color=color,
                       edgecolor="white", linewidth=1.5, label=pop, zorder=3)
        # The correlation across populations from section C in the title
        stat = rna_stats[rna_stats["marker"] == marker].iloc[0]
        ax.set_title(f"{marker} protein vs {info['gene']} RNA\nSpearman r = {stat['spearman_r']:.2f}",
                     loc="left")
        ax.set_xlabel(f"{info['gene']} expression (scRNA-seq)")
        ax.set_ylabel(f"{marker} median intensity (FACS)")
    np.atleast_1d(axes)[-1].legend(loc="upper left", bbox_to_anchor=(1, 1))
    save(fig, "C1_protein_vs_rna.png")


def plot_imaging(imaging, condition_column, conditions):
    """Grid: one row per imaging feature, one column per population, x = FACS condition."""
    n_rows, n_cols = len(IMAGING_FEATURES), len(POPULATIONS)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(3.6 * n_cols, 2.8 * n_rows), squeeze=False)
    for r, feature in enumerate(IMAGING_FEATURES):
        for c, pop in enumerate(POPULATIONS):
            ax = axes[r][c]
            # Values of this feature in this population, one list per condition
            part = imaging[(imaging["feature"] == feature) & (imaging["population"] == pop)]
            data = [part[part[condition_column] == cond]["value"] for cond in conditions]
            boxplot_with_points(ax, data, range(len(conditions)), population_color(pop), width=0.5)
            ax.set_xticks(range(len(conditions)), conditions)
            # Population name only on the top row, feature name only on the left column
            if r == 0:
                ax.set_title(pop, loc="left")
            if c == 0:
                ax.set_ylabel(feature)
    save(fig, "D1_imaging_features.png")


# =============================================================================
# 4b. REFERENCE ATLAS COMPARISON (section E - optional)
# =============================================================================
# Everything for section E is in this block: reading the atlas, building its
# tables, the statistics and the three E plots. It only runs when an atlas file
# is given; the sections A-D above do not use anything from here.
#
# The idea: a reference atlas is an independent, well-annotated single-cell
# dataset of urine. If our FACS and Flex results look like the atlas, that
# supports our gating and annotation; where they differ, it shows either a
# biological effect (e.g. a disease condition) or a technical bias of a method.
# This is also useful when FACS and scRNA-seq measured different people: the
# atlas is a third, shared reference for both.

ATLAS_SOURCES = ["Atlas RNA", "Flex RNA", "FACS protein"]   # the three panels of plot E3


def load_atlas(path):
    """
    Read the reference atlas and give every row a 'population', using
    POPULATIONS[...]["atlas_labels"] (or the seq_labels if none are given).
    Adds a "cells_in_row" column: 1 per row for a per-cell file, or the number from
    ATLAS_COUNT_COLUMN for a summary table. Labels not listed become 'Other'.
    """
    # Read only the columns that are used (a published atlas can have many more);
    # cell types and donors stay text
    needed = {ATLAS_CELL_TYPE_COLUMN, ATLAS_DONOR_COLUMN, ATLAS_COUNT_COLUMN,
              *(m["gene"] for m in MARKERS.values())}
    atlas = pd.read_csv(path, dtype={ATLAS_CELL_TYPE_COLUMN: str, ATLAS_DONOR_COLUMN: str},
                        usecols=lambda c: c in needed)
    get_column(atlas, ATLAS_CELL_TYPE_COLUMN, path.name)

    # Dictionary: atlas label -> population, e.g. "Myeloid cell" -> "Other immune"
    label_to_population = {
        label: population
        for population, info in POPULATIONS.items()
        for label in info.get("atlas_labels", info["seq_labels"])
    }
    atlas["population"] = atlas[ATLAS_CELL_TYPE_COLUMN].map(label_to_population).fillna("Other")

    # How many cells each row stands for: 1 (per-cell file) or the count column (summary table)
    if ATLAS_COUNT_COLUMN in atlas.columns:
        atlas["cells_in_row"] = pd.to_numeric(atlas[ATLAS_COUNT_COLUMN], errors="coerce").fillna(0)
    else:
        atlas["cells_in_row"] = 1

    # No donor column: treat the whole atlas as a single reference group
    if ATLAS_DONOR_COLUMN not in atlas.columns:
        atlas[ATLAS_DONOR_COLUMN] = "atlas"

    # Tell the user which atlas labels were not assigned to a population
    unmapped = sorted(atlas.loc[atlas["population"] == "Other", ATLAS_CELL_TYPE_COLUMN].unique())
    if unmapped:
        print(f"Note: these atlas cell types are not in any population (counted as 'Other'): {unmapped}")
    print(f"Loaded reference atlas: {int(atlas['cells_in_row'].sum())} cells "
          f"from {atlas[ATLAS_DONOR_COLUMN].nunique()} donor(s).")
    return atlas


def build_atlas_proportions(atlas):
    """
    One row per atlas donor x population: cell count, the total it is divided
    by, and the proportion. Uses the same PROPORTION_DENOMINATOR rule as our own
    samples, so the numbers can be compared directly.
    """
    populations = list(POPULATIONS)
    # Count table: rows = donors, columns = populations (+ Other), values = number of cells
    counts = atlas.pivot_table(index=ATLAS_DONOR_COLUMN, columns="population",
                               values="cells_in_row", aggfunc="sum", fill_value=0)
    counts = counts.reindex(columns=populations + ["Other"], fill_value=0)
    # Same denominator rule as build_proportions()
    if PROPORTION_DENOMINATOR == "populations":
        total = counts[populations].sum(axis=1)
    else:
        total = counts.sum(axis=1)

    rows = []
    for donor in counts.index:
        for pop in populations:
            rows.append({
                "donor": donor,
                "population": pop,
                "count": counts.loc[donor, pop],
                "total": total[donor],
                "proportion": counts.loc[donor, pop] / total[donor] if total[donor] > 0 else np.nan,
            })
    return pd.DataFrame(rows)


def atlas_pooled_proportions(atlas_props):
    """The atlas proportion of each population with all atlas cells pooled together."""
    # Add up counts and totals over all donors, then divide
    sums = atlas_props.groupby("population", sort=False)[["count", "total"]].sum()
    return (sums["count"] / sums["total"].where(sums["total"] > 0)).reindex(list(POPULATIONS))


# How the generic test names of compare_to_reference() are written in the section E table
ATLAS_TEST_NAMES = {
    "Mann-Whitney U": "Mann-Whitney U vs atlas donors",
    "Wilcoxon vs reference value": "Wilcoxon vs atlas proportion",
    "Wilcoxon of reference vs value": "Wilcoxon of atlas donors vs sample value",
}


def atlas_comparison_stats(proportions, atlas_props, condition_column, conditions):
    """
    Section E: for every population x method x condition (the conditions chosen
    for that method), compare our sample proportions with the atlas
    (compare_to_reference, with the atlas as the reference):
      - >= MIN_SAMPLES_PER_GROUP samples AND atlas donors:
            Mann-Whitney U test, our samples vs the atlas donors
            (takes into account that healthy people also differ from each other)
      - enough samples, but an atlas without donors (or too few):
            one-sample Wilcoxon test: are our samples centred on the atlas proportion?
      - a single sample (e.g. a FACS pool) and enough atlas donors:
            one-sample Wilcoxon test: are the atlas donors centred on that sample's value?
    FDR correction is applied over the whole table.
    conditions = {"scRNA-seq": [...], "FACS": [...]}
    """
    pooled = atlas_pooled_proportions(atlas_props)
    rows = []
    for pop in POPULATIONS:
        # The proportion of this population in every atlas donor
        donors = atlas_props.loc[atlas_props["population"] == pop, "proportion"].dropna().values
        # The reference: the donors when there are enough of them, otherwise the pooled atlas value
        reference = donors if len(donors) >= MIN_SAMPLES_PER_GROUP else [pooled[pop]]
        for method in METHODS:
            for group in conditions[method]:
                # Our samples of this method and condition
                values = proportions.loc[(proportions["population"] == pop) &
                                         (proportions["method"] == method) &
                                         (proportions[condition_column] == group), "proportion"].dropna().values
                median = np.median(values) if len(values) else np.nan
                test, statistic, p = compare_to_reference(values, reference)
                rows.append({"population": pop, "method": method, "group": group,
                             "n_samples": len(values), "median_samples": median,
                             "atlas_proportion": pooled[pop],
                             "atlas_donor_median": np.median(donors) if len(donors) else np.nan,
                             "difference": median - pooled[pop],       # positive = more than in the atlas
                             "n_atlas_donors": len(donors),
                             "test": ATLAS_TEST_NAMES.get(test, test), "statistic": statistic, "p_value": p})
    table = pd.DataFrame(rows)
    # Correct for all tests in the table, then flag the significant ones
    table["q_value_fdr"] = fdr_correct(table["p_value"])
    table["significant"] = table["q_value_fdr"] < ALPHA
    return table


def composition_difference(proportions, atlas_props, condition_column):
    """
    Section E: how different is each sample's WHOLE composition from the atlas?
      composition difference = half the sum of |sample proportion - atlas proportion|
                               over all populations
    This is the fraction of cells that would have to change population for the
    sample to look exactly like the atlas (0 = identical, 1 = nothing in common).
    The atlas donors are scored as well (each donor vs the pooled atlas), to show
    how much reference individuals differ among themselves.
    """
    pooled = atlas_pooled_proportions(atlas_props)
    rows = []
    # Our samples: one score per sample and method
    for (sample, method, group), part in proportions.groupby([SAMPLE_COLUMN, "method", condition_column], sort=False):
        p = part.set_index("population")["proportion"].reindex(pooled.index)
        rows.append({"sample": sample, "method": method, "group": group,
                     "composition_difference": 0.5 * np.abs(p - pooled).sum(min_count=1)})
    # The atlas donors (only meaningful if there is more than one)
    if atlas_props["donor"].nunique() > 1:
        for donor, part in atlas_props.groupby("donor", sort=False):
            p = part.set_index("population")["proportion"].reindex(pooled.index)
            rows.append({"sample": donor, "method": "Atlas donor", "group": "Atlas",
                         "composition_difference": 0.5 * np.abs(p - pooled).sum(min_count=1)})
    return pd.DataFrame(rows)


def build_marker_profiles(atlas, cells, facs_markers):
    """
    Section E: the typical level of each marker in each population, from three sources:
      Atlas RNA    : mean expression of the gene in the atlas cells
      Flex RNA     : mean expression of the gene in our sequenced cells
      FACS protein : median (over FACS samples) of the antibody's median intensity
    The three sources use different units, so inside each source every marker is
    scaled 0-1 across the populations (0 = lowest population, 1 = highest).
    Only the PATTERN is compared: which population is high or low for a marker.
    Returns None when the atlas has no per-cell gene columns.
    """
    genes = [m["gene"] for m in MARKERS.values()]
    # Needs a per-cell atlas with all four gene columns
    if ATLAS_COUNT_COLUMN in atlas.columns or not all(g in atlas.columns for g in genes):
        print("Note: the atlas has no per-cell marker gene columns - the marker comparison (E3) is skipped.")
        return None
    populations = list(POPULATIONS)
    # Average gene expression per population in the atlas and in our Flex data
    atlas_means = atlas[atlas["population"].isin(populations)].groupby("population")[genes].mean()
    seq_means = cells[cells["population"].isin(populations)].groupby("population")[genes].mean()
    # Typical FACS intensity per population and marker (median over FACS samples)
    facs_medians = facs_markers.groupby(["population", "marker"])["protein_median"].median()

    rows = []
    for marker, info in MARKERS.items():
        for pop in populations:
            for source, value in [("Atlas RNA", atlas_means[info["gene"]].get(pop, np.nan)),
                                  ("Flex RNA", seq_means[info["gene"]].get(pop, np.nan)),
                                  ("FACS protein", facs_medians.get((pop, marker), np.nan))]:
                rows.append({"source": source, "marker": marker, "gene": info["gene"],
                             "population": pop, "value": value})
    table = pd.DataFrame(rows)
    # Scale 0-1 inside each source x marker: (value - lowest) / (highest - lowest)
    grouped = table.groupby(["source", "marker"])["value"]
    low, high = grouped.transform("min"), grouped.transform("max")
    table["scaled"] = ((table["value"] - low) / (high - low)).where(high > low)
    return table


def plot_atlas_proportions(proportions, atlas_props, condition_column, conditions):
    """
    E1: one panel per population. Our samples per condition (scRNA-seq and FACS
    boxes, as in B1), the atlas donors as an extra box, and a dashed line at the
    pooled atlas proportion. conditions = all chosen conditions of both methods.
    """
    pooled = atlas_pooled_proportions(atlas_props)
    has_donors = atlas_props["donor"].nunique() > 1
    groups = list(conditions) + (["Atlas"] if has_donors else [])
    fig, axes = plt.subplots(1, len(POPULATIONS),
                             figsize=(max(4, 1.6 * len(groups)) * len(POPULATIONS) / 1.3, 4))
    for ax, pop in zip(np.atleast_1d(axes), POPULATIONS):
        part = proportions[proportions["population"] == pop]
        # Our samples: scRNA-seq box slightly left, FACS box slightly right of each condition
        for k, method in enumerate(METHODS):
            data = [part[(part["method"] == method) & (part[condition_column] == c)]["proportion"]
                    for c in conditions]
            positions = [i + (k - 0.5) * 0.4 for i in range(len(conditions))]
            boxplot_with_points(ax, data, positions, METHOD_COLORS[method])
            ax.scatter([], [], color=METHOD_COLORS[method], s=22, label=method)  # legend entry
        # The atlas donors as their own box at the right
        if has_donors:
            donors = atlas_props.loc[atlas_props["population"] == pop, "proportion"]
            boxplot_with_points(ax, [donors], [len(conditions)], ATLAS_COLOR)
            ax.scatter([], [], color=ATLAS_COLOR, s=22, label="Atlas donors")
        # The pooled atlas proportion as a reference line across the panel
        ax.axhline(pooled[pop], color=ATLAS_COLOR, linewidth=1.2, linestyle="--", label="Atlas (all cells)")
        ax.set_xticks(range(len(groups)), groups)
        ax.set_title(pop, loc="left")
        ax.set_ylabel("Proportion of cells")
    np.atleast_1d(axes)[-1].legend(loc="upper left", bbox_to_anchor=(1, 1))
    save(fig, "E1_atlas_proportions.png")


def plot_atlas_difference(differences, conditions):
    """E2: composition difference from the atlas per sample (in %), per condition and method."""
    fig, ax = plt.subplots(figsize=(max(5, 1.8 * (len(conditions) + 1)), 4))
    # Our samples, scRNA-seq and FACS side by side for each condition
    for k, method in enumerate(METHODS):
        data = [differences[(differences["method"] == method) &
                            (differences["group"] == c)]["composition_difference"] * 100
                for c in conditions]
        positions = [i + (k - 0.5) * 0.4 for i in range(len(conditions))]
        boxplot_with_points(ax, data, positions, METHOD_COLORS[method])
        ax.scatter([], [], color=METHOD_COLORS[method], s=22, label=method)  # legend entry
    # The atlas donors: how much reference individuals differ from the pooled atlas
    groups = list(conditions)
    donors = differences[differences["method"] == "Atlas donor"]["composition_difference"] * 100
    if len(donors) > 0:
        boxplot_with_points(ax, [donors], [len(conditions)], ATLAS_COLOR)
        ax.scatter([], [], color=ATLAS_COLOR, s=22, label="Atlas donors")
        groups.append("Atlas donors")
    ax.set_xticks(range(len(groups)), groups)
    ax.set_ylim(bottom=0)
    ax.set_ylabel("Composition difference from atlas (%)")
    ax.set_title("How different is each sample from the reference atlas?", loc="left")
    ax.legend(loc="upper left", bbox_to_anchor=(1, 1))
    save(fig, "E2_atlas_composition_difference.png")


def plot_marker_profiles(profiles):
    """
    E3: three heatmaps side by side (atlas RNA, Flex RNA, FACS protein).
    Rows = populations, columns = markers, color = scaled level (0-1).
    The titles give the Spearman correlation of each pattern with the atlas pattern.
    """
    markers, pops = list(MARKERS), list(POPULATIONS)
    # Three heatmap panels + one narrow panel for the color scale
    fig, axes = plt.subplots(1, 4, figsize=(3.2 * 3 + 1.2, 0.55 * len(pops) + 2.2),
                             gridspec_kw={"width_ratios": [1, 1, 1, 0.07]})
    atlas_scaled = (profiles[profiles["source"] == "Atlas RNA"]
                    .set_index(["population", "marker"])["scaled"])
    for i, (ax, source) in enumerate(zip(axes[:3], ATLAS_SOURCES)):
        part = profiles[profiles["source"] == source]
        grid = part.pivot(index="population", columns="marker", values="scaled").reindex(index=pops, columns=markers)
        image = ax.imshow(grid.values, cmap="RdPu", vmin=0, vmax=1, aspect="auto")
        # Write the value in every cell (dark text on light cells, white on dark cells)
        for r in range(len(pops)):
            for c in range(len(markers)):
                value = grid.values[r, c]
                ax.text(c, r, "–" if np.isnan(value) else f"{value:.2f}", ha="center", va="center",
                        fontsize=8, color="white" if value > 0.6 else "#2e2228")
        ax.set_xticks(range(len(markers)), markers)
        ax.set_yticks(range(len(pops)), pops if i == 0 else [])
        ax.grid(False)
        # Compare this pattern with the atlas pattern (same population x marker cells)
        title = source
        if source != "Atlas RNA":
            this = part.set_index(["population", "marker"])["scaled"].reindex(atlas_scaled.index)
            r_value = correlate(atlas_scaled.values, this.values)["spearman_r"]
            title += f"\nSpearman r vs atlas = {r_value:.2f}"
        ax.set_title(title, loc="left")
    fig.colorbar(image, cax=axes[3], label="Scaled level (0 = lowest, 1 = highest)")
    save(fig, "E3_atlas_marker_profiles.png")


def run_atlas_comparison(atlas, atlas_props, proportions, cells, facs_markers, condition_column, conditions):
    """
    Section E from start to end (the atlas is loaded in run_analysis, because the
    composition overview of section A also uses it): run the statistics, save the
    tables and draw the plots. Returns the test table.
    """
    tests = atlas_comparison_stats(proportions, atlas_props, condition_column, conditions)
    differences = composition_difference(proportions, atlas_props, condition_column)
    profiles = build_marker_profiles(atlas, cells, facs_markers)

    # Tables
    atlas_props.to_csv(OUTPUT_DIR / "E_atlas_proportions.csv", index=False)
    tests.to_csv(OUTPUT_DIR / "E_atlas_proportion_tests.csv", index=False)
    differences.to_csv(OUTPUT_DIR / "E_atlas_composition_difference.csv", index=False)
    if profiles is not None:
        profiles.to_csv(OUTPUT_DIR / "E_atlas_marker_profiles.csv", index=False)

    # Plots
    plot_atlas_proportions(proportions, atlas_props, condition_column, all_conditions(conditions))
    plot_atlas_difference(differences, all_conditions(conditions))
    if profiles is not None:
        plot_marker_profiles(profiles)
    return tests


# =============================================================================
# 4c. IMAGING PROFILES PER CELL TYPE (section F)
# =============================================================================
# Everything for section F is in this block. It describes each population with
# all the S8 imaging features at once: typical values, a "fingerprint", a PCA
# map, outliers, and links between shape and marker intensity. It only uses the
# FACS samples (imaging comes from the S8). Sections A-E do
# not use anything from here.
#
# Two levels of data:
#   - SAMPLE level (always): the FlowJo table gives, per sample and population,
#     the MEDIAN of each feature (the "imaging" table built for section D).
#   - CELL level (optional): a per-cell export gives every single cell, so we
#     can look at whole distributions and at unusual cells inside a sample.

CONDITION_MARKERS = ["o", "s", "^", "D", "v", "P", "X"]   # point shapes for conditions in plot F3


def robust_z(values):
    """
    Robust z-score (also called modified z-score):
        0.6745 x (value - median) / MAD,   MAD = median of |value - median|
    Like a normal z-score, but based on the median, so a few extreme values
    cannot hide themselves by inflating the mean and SD. Values above
    OUTLIER_CUTOFF (3.5) are called outliers. NaN if all values are identical.
    """
    x = np.asarray(values, dtype=float)
    median = np.nanmedian(x) if np.any(~np.isnan(x)) else np.nan
    mad = np.nanmedian(np.abs(x - median)) if not np.isnan(median) else np.nan
    if not mad > 0:          # no spread (or no data): cannot judge outliers
        return np.full_like(x, np.nan)
    return 0.6745 * (x - median) / mad


def imaging_outlier_samples(imaging, condition_column):
    """
    Sample-level outliers: for every population x feature, the robust z-score of
    each sample's median compared with all samples of that population.
    One row per sample x population x feature, with 'outlier' True/False.
    """
    table = imaging[[SAMPLE_COLUMN, condition_column, "population", "feature", "value"]].copy()
    # robust z inside each population x feature (samples are compared with each other)
    table["robust_z"] = table.groupby(["population", "feature"])["value"].transform(robust_z)
    table["outlier"] = table["robust_z"].abs() > OUTLIER_CUTOFF
    return table


def imaging_population_summary(imaging, outliers):
    """
    One row per population x feature with the typical value over samples:
    n, mean, SD, 95% confidence interval of the mean, median, CV%, min, max,
    and how many samples are outliers.
    """
    rows = []
    for pop in POPULATIONS:
        for feature in IMAGING_FEATURES:
            keep = (imaging["population"] == pop) & (imaging["feature"] == feature)
            values = imaging.loc[keep, "value"].dropna().values
            n = len(values)
            mean = values.mean() if n else np.nan
            sd = values.std(ddof=1) if n > 1 else np.nan
            # 95% confidence interval of the mean (t distribution, suits small n)
            half_width = stats.t.ppf(0.975, n - 1) * sd / np.sqrt(n) if n > 1 else np.nan
            flagged = outliers[(outliers["population"] == pop) & (outliers["feature"] == feature)]
            rows.append({
                "population": pop, "feature": feature, "n_samples": n,
                "mean": mean, "sd": sd,
                "ci95_low": mean - half_width, "ci95_high": mean + half_width,
                "median": np.median(values) if n else np.nan,
                "cv_percent": 100 * sd / mean if n > 1 and mean != 0 else np.nan,   # spread relative to the mean
                "min": values.min() if n else np.nan, "max": values.max() if n else np.nan,
                "n_outlier_samples": int(flagged["outlier"].sum()),
            })
    return pd.DataFrame(rows)


def imaging_pca(imaging, condition_column):
    """
    PCA (principal component analysis) of the imaging features. Each point is
    one sample x population, described by all its feature medians. PCA finds
    the two directions (PC1, PC2) that capture most of the differences, so a
    many-feature description can be drawn on a flat map.
    Features are standardized first (mean 0, SD 1) so each counts equally.
    Returns (scores table, loadings table), or (None, None) if there are fewer
    than 2 usable features or 3 points.
    """
    # Table: rows = sample x population, columns = features
    wide = (imaging.pivot_table(index=[SAMPLE_COLUMN, "population", condition_column],
                                columns="feature", values="value")
            .reindex(columns=IMAGING_FEATURES).dropna())
    # Features without any variation carry no information: leave them out
    wide = wide.loc[:, wide.std(ddof=1) > 0]
    if wide.shape[1] < 2 or len(wide) < 3:
        return None, None
    z = (wide - wide.mean()) / wide.std(ddof=1)                 # standardize each feature
    u, s, vt = np.linalg.svd(z.values, full_matrices=False)    # the PCA itself
    explained = s ** 2 / np.sum(s ** 2)                        # share of variation per PC

    scores = wide.index.to_frame(index=False)
    scores["PC1"], scores["PC2"] = u[:, 0] * s[0], u[:, 1] * s[1]
    scores["PC1_explained"], scores["PC2_explained"] = explained[0], explained[1]
    # Loadings: how much each feature contributes to PC1 and PC2 (the arrows in F3)
    loadings = pd.DataFrame({"feature": wide.columns, "PC1": vt[0], "PC2": vt[1]})
    return scores, loadings


def imaging_marker_correlations(imaging, facs_markers):
    """
    Inside each population: across FACS samples, does an imaging feature go together
    with a marker's FACS intensity? (e.g. do larger proximal tubule cells carry
    more CD13?) Spearman correlation for every population x feature x marker,
    FDR-corrected over the whole table. (Both come from the same FACS samples.)
    """
    # Pair the two measurements of the same sample and population, once for the whole table
    both = imaging.merge(facs_markers[[SAMPLE_COLUMN, "population", "marker", "protein_median"]],
                         on=[SAMPLE_COLUMN, "population"])
    pairs = dict(iter(both.groupby(["population", "feature", "marker"], sort=False)))
    no_pairs = both.iloc[:0]
    rows = []
    for pop in POPULATIONS:
        for feature in IMAGING_FEATURES:
            for marker in MARKERS:
                part = pairs.get((pop, feature, marker), no_pairs)
                rows.append({"population": pop, "feature": feature, "marker": marker,
                             **correlate(part["value"], part["protein_median"])})
    table = pd.DataFrame(rows)
    table["spearman_q_fdr"] = fdr_correct(table["spearman_p"])
    return table


def load_facs_cells(path, metadata, condition_column):
    """
    Read the optional per-cell imaging export. Each cell gets its sample_id
    (directly, or translated from facs_name) and the condition of its sample.
    Only cells of known populations and of the chosen samples are kept.
    metadata = the FACS rows of the sample sheet.
    """
    # Read only the columns that are used (a FlowJo export has many more channels)
    needed = {SAMPLE_COLUMN, FACS_NAME_COLUMN, FACS_CELLS_POPULATION_COLUMN, *IMAGING_FEATURES}
    cells_img = pd.read_csv(path, dtype={SAMPLE_COLUMN: str, FACS_NAME_COLUMN: str},
                            usecols=lambda c: c in needed)
    # Which sample is each cell from? Use sample_id, or translate facs_name via the metadata
    if SAMPLE_COLUMN not in cells_img.columns:
        if FACS_NAME_COLUMN not in cells_img.columns:
            raise SystemExit(f"{path.name} needs a '{SAMPLE_COLUMN}' or a '{FACS_NAME_COLUMN}' column "
                             "saying which sample each cell comes from.")
        name_to_id = dict(zip(metadata[FACS_NAME_COLUMN], metadata[SAMPLE_COLUMN]))
        cells_img[SAMPLE_COLUMN] = cells_img[FACS_NAME_COLUMN].map(name_to_id)
    # Check the population column and every imaging feature column
    get_column(cells_img, FACS_CELLS_POPULATION_COLUMN, path.name)
    for feature in IMAGING_FEATURES:
        get_column(cells_img, feature, path.name)
    cells_img = cells_img.rename(columns={FACS_CELLS_POPULATION_COLUMN: "population"})

    # Warn about population names that are not in SETTINGS (those cells are left out)
    unknown = sorted(set(cells_img["population"].astype(str)) - set(POPULATIONS))
    if unknown:
        print(f"Note: these populations in {path.name} are not in POPULATIONS and are left out: {unknown}")
    keep = cells_img["population"].isin(list(POPULATIONS)) & cells_img[SAMPLE_COLUMN].isin(metadata[SAMPLE_COLUMN])
    cells_img = cells_img[keep].merge(metadata[[SAMPLE_COLUMN, condition_column]], on=SAMPLE_COLUMN)
    print(f"Loaded {len(cells_img)} single FACS cells with imaging features "
          f"from {cells_img[SAMPLE_COLUMN].nunique()} samples.")
    return cells_img


def cell_outlier_flags(cells_img):
    """
    Cell-level robust z-scores: every cell is compared with ALL cells of its
    population (all chosen samples together). Returns a table of True/False,
    one column per feature, plus 'Any feature' (outlier in at least one feature).
    """
    flags = pd.DataFrame(index=cells_img.index)
    for feature in IMAGING_FEATURES:
        z = cells_img.groupby("population")[feature].transform(robust_z)
        flags[feature] = z.abs() > OUTLIER_CUTOFF
    flags["Any feature"] = flags[IMAGING_FEATURES].any(axis=1)
    return flags


def cell_level_summary(cells_img, flags, condition_column):
    """
    One row per sample x population x feature, from the single cells:
    number of cells, mean, median, SD, CV%, 5th and 95th percentile, and the
    % of outlier cells. Extra rows with feature 'Any feature' give the % of
    cells that are an outlier in at least one feature.
    Calculated per feature for all groups at once (fast also for millions of cells).
    """
    keys = [SAMPLE_COLUMN, "population", condition_column]
    grouped = cells_img.groupby(keys, sort=False)
    # The outlier flags grouped the same way (flags has the same rows as cells_img)
    flag_groups = flags.groupby([cells_img[k] for k in keys], sort=False)
    parts = []
    for feature in IMAGING_FEATURES:
        values = grouped[feature]
        part = pd.DataFrame({"n_cells": values.count(), "mean": values.mean(), "median": values.median(),
                             "sd": values.std(ddof=1),
                             "p5": values.quantile(0.05), "p95": values.quantile(0.95),   # the middle 90% of cells
                             "pct_outlier_cells": 100 * flag_groups[feature].mean()})
        part["cv_percent"] = 100 * part["sd"] / part["mean"].where(part["mean"] != 0)
        parts.append(part.assign(feature=feature))
    # Cells that are unusual in at least one feature
    parts.append(pd.DataFrame({"n_cells": grouped.size(),
                               "pct_outlier_cells": 100 * flag_groups["Any feature"].mean()})
                 .assign(feature="Any feature"))
    columns = keys + ["feature", "n_cells", "mean", "median", "sd", "cv_percent", "p5", "p95", "pct_outlier_cells"]
    return (pd.concat(parts).reset_index().reindex(columns=columns)
            .sort_values(keys[:2], kind="stable", ignore_index=True))


def feature_grid(n_panels, width=4.2, height=3.6):
    """A grid of panels (at most 3 per row) for 'one panel per feature' plots; unused panels hidden."""
    n_cols = min(3, n_panels)
    n_rows = int(np.ceil(n_panels / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(width * n_cols, height * n_rows), squeeze=False)
    for ax in axes.flat[n_panels:]:
        ax.set_visible(False)
    return fig, list(axes.flat[:n_panels])


def plot_imaging_means(summary, outliers):
    """
    F1: one panel per feature. For each population: a bar at the mean over
    samples, the 95% confidence interval as an error bar, every sample as a dot.
    Outlier samples get a red ring and their name.
    """
    pops = list(POPULATIONS)
    fig, axes = feature_grid(len(IMAGING_FEATURES))
    rng = np.random.default_rng(0)            # same dot jitter every run
    for ax, feature in zip(axes, IMAGING_FEATURES):
        for i, pop in enumerate(pops):
            color = population_color(pop)
            row = summary[(summary["population"] == pop) & (summary["feature"] == feature)].iloc[0]
            if np.isnan(row["mean"]):
                continue
            # Bar = mean, error bar = 95% confidence interval of the mean
            ax.bar(i, row["mean"], width=0.6, color=color, alpha=0.25, edgecolor=color)
            if np.isfinite(row["ci95_low"]):
                ax.errorbar(i, row["mean"], yerr=[[row["mean"] - row["ci95_low"]], [row["ci95_high"] - row["mean"]]],
                            color=color, capsize=4, linewidth=1.5)
            # Dots = samples; outliers with a red ring and the sample name
            part = outliers[(outliers["population"] == pop) & (outliers["feature"] == feature)].dropna(subset=["value"])
            x = i + rng.uniform(-0.15, 0.15, len(part))
            ax.scatter(x, part["value"], s=18, color=color, edgecolor="white", linewidth=0.8, zorder=3)
            flagged = part["outlier"].values
            ax.scatter(x[flagged], part["value"].values[flagged], s=70, facecolor="none",
                       edgecolor="#d62728", linewidth=1.5, zorder=4)
            for xi, (_, out) in zip(x[flagged], part[flagged].iterrows()):
                ax.annotate(out[SAMPLE_COLUMN], (xi, out["value"]), xytext=(5, 3),
                            textcoords="offset points", fontsize=7, color="#d62728")
        ax.set_xticks(range(len(pops)), pops, rotation=20, ha="right")
        ax.set_title(feature, loc="left")
        ax.set_ylabel("Mean of sample medians")
    # Legend explaining the red rings
    axes[0].scatter([], [], s=70, facecolor="none", edgecolor="#d62728", linewidth=1.5,
                    label=f"Outlier sample (robust z > {OUTLIER_CUTOFF})")
    axes[0].legend(loc="upper left", fontsize=7)
    save(fig, "F1_imaging_feature_means.png")


def plot_imaging_fingerprint(summary):
    """
    F2: heatmap, rows = populations, columns = features. The color is the
    population's mean compared with the other populations (z-score per feature:
    red = higher than the other populations, blue = lower); the number written
    in each cell is the actual mean value.
    """
    pops, features = list(POPULATIONS), IMAGING_FEATURES
    means = summary.pivot(index="population", columns="feature", values="mean").reindex(index=pops, columns=features)
    # z-score per feature across the populations (each column centred on its own average)
    z = (means - means.mean()) / means.std(ddof=0).where(lambda s: s > 0)
    fig, (ax, cax) = plt.subplots(1, 2, figsize=(1.7 * len(features) + 3, 0.6 * len(pops) + 2),
                                  gridspec_kw={"width_ratios": [1, 0.05]})
    image = ax.imshow(z.values, cmap="RdBu_r", vmin=-2, vmax=2, aspect="auto")
    for r in range(len(pops)):
        for c in range(len(features)):
            ax.text(c, r, f"{means.values[r, c]:.3g}", ha="center", va="center", fontsize=8,
                    color="white" if abs(np.nan_to_num(z.values[r, c])) > 1.3 else "#2e2228")
    ax.set_xticks(range(len(features)), features, rotation=20, ha="right")
    ax.set_yticks(range(len(pops)), pops)
    ax.grid(False)
    ax.set_title("Imaging fingerprint of each population (numbers = mean value)", loc="left")
    fig.colorbar(image, cax=cax, label="Compared with other populations (z-score)")
    save(fig, "F2_imaging_fingerprint.png")


def plot_imaging_pca(scores, loadings, condition_column, conditions):
    """
    F3: PCA map. Each point = one sample x population; color = population,
    shape = condition. Arrows show which feature pulls points in which direction.
    """
    # Left: the map. Right: a narrow empty panel that holds the legend (so it is never cut off)
    fig, (ax, legend_ax) = plt.subplots(1, 2, figsize=(8, 5), gridspec_kw={"width_ratios": [1, 0.32]})
    legend_ax.axis("off")
    for pop in POPULATIONS:
        for k, cond in enumerate(conditions):
            part = scores[(scores["population"] == pop) & (scores[condition_column] == cond)]
            ax.scatter(part["PC1"], part["PC2"], s=45, color=population_color(pop),
                       marker=CONDITION_MARKERS[k % len(CONDITION_MARKERS)],
                       edgecolor="white", linewidth=1, zorder=3)
    # Feature arrows (loadings), scaled so that they stay inside the cloud of points
    # in both directions (one scale for x and y, so the arrow angles stay correct)
    scale = 0.8 * min(np.nanmax(np.abs(scores["PC1"])) / max(np.max(np.abs(loadings["PC1"])), 1e-9),
                      np.nanmax(np.abs(scores["PC2"])) / max(np.max(np.abs(loadings["PC2"])), 1e-9))
    for _, row in loadings.iterrows():
        ax.annotate("", xy=(row["PC1"] * scale, row["PC2"] * scale), xytext=(0, 0),
                    arrowprops={"arrowstyle": "->", "color": "#52514e", "linewidth": 1})
        ax.text(row["PC1"] * scale * 1.08, row["PC2"] * scale * 1.08, row["feature"], fontsize=7, color="#52514e")
    # Invisible points just beyond the arrow tips, so the axes always include every arrow and label
    ax.scatter(loadings["PC1"] * scale * 1.25, loadings["PC2"] * scale * 1.25, s=0)
    # Legend: colors = populations, shapes = conditions
    for pop in POPULATIONS:
        ax.scatter([], [], color=population_color(pop), s=40, label=pop)
    for k, cond in enumerate(conditions):
        ax.scatter([], [], color="#8a8984", marker=CONDITION_MARKERS[k % len(CONDITION_MARKERS)], s=40, label=cond)
    legend_ax.legend(*ax.get_legend_handles_labels(), loc="upper left")
    ax.axhline(0, color="#e6e5e0", linewidth=1, zorder=1)
    ax.axvline(0, color="#e6e5e0", linewidth=1, zorder=1)
    ax.set_xlabel(f"PC1 ({scores['PC1_explained'].iloc[0]:.0%} of variation)")
    ax.set_ylabel(f"PC2 ({scores['PC2_explained'].iloc[0]:.0%} of variation)")
    ax.set_title("Imaging map: samples x populations (PCA of all features)", loc="left")
    save(fig, "F3_imaging_pca.png")


def plot_imaging_outliers(outliers, condition_column, conditions):
    """
    F4: heatmap of robust z-scores. Rows = samples (grouped by condition),
    columns = population x feature. Outliers (|z| > OUTLIER_CUTOFF) are labelled
    with their z-value.
    """
    # Order the samples by condition, then by name
    order = (outliers[[SAMPLE_COLUMN, condition_column]].drop_duplicates()
             .assign(_c=lambda d: d[condition_column].map(list(conditions).index))
             .sort_values(["_c", SAMPLE_COLUMN]))
    columns = [(pop, feature) for pop in POPULATIONS for feature in IMAGING_FEATURES]
    grid = (outliers.pivot_table(index=SAMPLE_COLUMN, columns=["population", "feature"], values="robust_z")
            .reindex(index=order[SAMPLE_COLUMN], columns=pd.MultiIndex.from_tuples(columns)))
    fig, (ax, cax) = plt.subplots(1, 2, figsize=(0.5 * len(columns) + 3.5, 0.32 * len(grid) + 3),
                                  gridspec_kw={"width_ratios": [1, 0.03]})
    image = ax.imshow(np.clip(grid.values, -6, 6), cmap="RdBu_r", vmin=-6, vmax=6, aspect="auto")
    # Write the z-value only where it is an outlier
    for r in range(grid.shape[0]):
        for c in range(grid.shape[1]):
            value = grid.values[r, c]
            if np.isfinite(value) and abs(value) > OUTLIER_CUTOFF:
                ax.text(c, r, f"{value:.1f}", ha="center", va="center", fontsize=7,
                        fontweight="bold", color="white")
    ax.set_xticks(range(len(columns)), [f"{p} · {f}" for p, f in columns], rotation=90, fontsize=7)
    ax.set_yticks(range(len(grid)), [f"{s} ({c})" for s, c in zip(order[SAMPLE_COLUMN], order[condition_column])],
                  fontsize=7)
    ax.grid(False)
    ax.set_title(f"Outlier samples: robust z-score (labelled when |z| > {OUTLIER_CUTOFF})", loc="left")
    fig.colorbar(image, cax=cax, label="Robust z-score")
    save(fig, "F4_imaging_outlier_samples.png")


def plot_imaging_marker_correlations(correlations):
    """
    F5: one heatmap per population. Rows = imaging features, columns = markers,
    color = Spearman r across samples. * = significant after FDR correction.
    """
    pops, markers = list(POPULATIONS), list(MARKERS)
    fig, axes = plt.subplots(1, len(pops) + 1,
                             figsize=(2.2 * len(pops) * len(markers) / 3 + 2.5, 0.5 * len(IMAGING_FEATURES) + 2.4),
                             gridspec_kw={"width_ratios": [1] * len(pops) + [0.07]})
    for i, (ax, pop) in enumerate(zip(axes[:-1], pops)):
        part = correlations[correlations["population"] == pop]
        r = part.pivot(index="feature", columns="marker", values="spearman_r").reindex(index=IMAGING_FEATURES, columns=markers)
        q = part.pivot(index="feature", columns="marker", values="spearman_q_fdr").reindex(index=IMAGING_FEATURES, columns=markers)
        image = ax.imshow(r.values, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
        for a in range(r.shape[0]):
            for b in range(r.shape[1]):
                value = r.values[a, b]
                star = "*" if q.values[a, b] < ALPHA else ""
                ax.text(b, a, "–" if np.isnan(value) else f"{value:.2f}{star}", ha="center", va="center",
                        fontsize=7, color="white" if abs(np.nan_to_num(value)) > 0.6 else "#2e2228")
        ax.set_xticks(range(len(markers)), markers, fontsize=8)
        ax.set_yticks(range(len(IMAGING_FEATURES)), IMAGING_FEATURES if i == 0 else [], fontsize=8)
        ax.grid(False)
        ax.set_title(pop, loc="left")
    fig.colorbar(image, cax=axes[-1], label="Spearman r (across samples)")
    save(fig, "F5_imaging_vs_markers.png")


def plot_cell_distributions(cells_img):
    """
    F6: one panel per feature; violin plot of the single cells of each population.
    The violin shape is estimated from at most MAX_CELLS_PER_VIOLIN random cells
    (the shape barely changes, but drawing gets much faster with large exports);
    the median line uses ALL cells.
    """
    pops = list(POPULATIONS)
    rng = np.random.default_rng(0)            # same random cells every run
    fig, axes = feature_grid(len(IMAGING_FEATURES))
    for ax, feature in zip(axes, IMAGING_FEATURES):
        for i, pop in enumerate(pops):
            values = cells_img.loc[cells_img["population"] == pop, feature].dropna().values
            if len(values) < 2:
                continue
            shown = values
            if len(values) > MAX_CELLS_PER_VIOLIN:
                shown = rng.choice(values, MAX_CELLS_PER_VIOLIN, replace=False)
            # The violin's width shows how many cells have each value
            part = ax.violinplot(shown, positions=[i], widths=0.8, showextrema=False)
            for body in part["bodies"]:
                body.set_facecolor(population_color(pop))
                body.set_edgecolor(population_color(pop))
                body.set_alpha(0.45)
            # The median of all cells as a black line across the violin
            ax.plot([i - 0.25, i + 0.25], [np.median(values)] * 2, color="#0b0b0b", linewidth=1.5)
        ax.set_xticks(range(len(pops)), pops, rotation=20, ha="right")
        ax.set_title(feature, loc="left")
        ax.set_ylabel("Value per cell")
    save(fig, "F6_imaging_cell_distributions.png")


def plot_cell_outliers(cell_summary):
    """
    F7: one panel per feature (+ 'Any feature'). Each dot = one sample: the % of
    its cells that are outliers compared with all cells of the same population.
    """
    pops = list(POPULATIONS)
    panels = IMAGING_FEATURES + ["Any feature"]
    fig, axes = feature_grid(len(panels), width=4, height=3.4)
    for ax, feature in zip(axes, panels):
        part = cell_summary[cell_summary["feature"] == feature]
        for i, pop in enumerate(pops):
            values = part.loc[part["population"] == pop, "pct_outlier_cells"]
            boxplot_with_points(ax, [values], [i], population_color(pop), width=0.5)
        ax.set_xticks(range(len(pops)), pops, rotation=20, ha="right")
        ax.set_ylim(bottom=0)
        ax.set_title(feature, loc="left")
        ax.set_ylabel("% outlier cells per sample")
    save(fig, "F7_imaging_cell_outliers.png")


def plot_morphology_map(cells_img, flags, condition_column, conditions):
    """
    F8: single cells on a map of the first two imaging features, one panel per
    condition, colored by population. Outlier cells (any feature) are drawn as
    black crosses. At most 4000 cells per panel are drawn, picked at random.
    """
    x_feature, y_feature = IMAGING_FEATURES[0], IMAGING_FEATURES[1]
    rng = np.random.default_rng(0)
    fig, axes = plt.subplots(1, len(conditions), figsize=(4.6 * len(conditions) + 1.5, 4.2), squeeze=False)
    for ax, cond in zip(axes[0], conditions):
        part = cells_img[cells_img[condition_column] == cond]
        if len(part) > 4000:
            part = part.loc[rng.choice(part.index, 4000, replace=False)]
        for pop in POPULATIONS:
            sub = part[part["population"] == pop]
            ax.scatter(sub[x_feature], sub[y_feature], s=4, alpha=0.4, color=population_color(pop),
                       linewidths=0, label=pop)
        odd = part[flags.loc[part.index, "Any feature"]]
        ax.scatter(odd[x_feature], odd[y_feature], s=10, marker="x", color="#0b0b0b", linewidths=0.6,
                   label="Outlier cell")
        ax.set_title(cond, loc="left")
        ax.set_xlabel(x_feature)
        ax.set_ylabel(y_feature)
    axes[0][-1].legend(loc="upper left", bbox_to_anchor=(1, 1), markerscale=3)
    save(fig, "F8_imaging_morphology_map.png")


def run_imaging_profiles(imaging, facs_markers, data_dir, metadata, condition_column, conditions):
    """
    Section F from start to end. Sample-level tables and plots F1-F5 always;
    cell-level table and plots F6-F8 only if the per-cell export is there.
    Uses only the FACS samples: metadata = the FACS rows of the sample sheet,
    conditions = the chosen FACS conditions. Returns the population summary table.
    """
    # --- Sample level (from the FlowJo table) ---
    outliers = imaging_outlier_samples(imaging, condition_column)
    summary = imaging_population_summary(imaging, outliers)
    scores, loadings = imaging_pca(imaging, condition_column)
    correlations = imaging_marker_correlations(imaging, facs_markers)

    summary.to_csv(OUTPUT_DIR / "F_imaging_population_summary.csv", index=False)
    outliers.to_csv(OUTPUT_DIR / "F_imaging_outlier_samples.csv", index=False)
    correlations.to_csv(OUTPUT_DIR / "F_imaging_marker_correlations.csv", index=False)
    if scores is not None:
        scores.to_csv(OUTPUT_DIR / "F_imaging_pca.csv", index=False)
        loadings.to_csv(OUTPUT_DIR / "F_imaging_pca_loadings.csv", index=False)

    plot_imaging_means(summary, outliers)
    plot_imaging_fingerprint(summary)
    if scores is not None:
        plot_imaging_pca(scores, loadings, condition_column, conditions)
    else:
        print("Note: the imaging PCA (F3) needs at least 2 imaging features with variation - skipped.")
    plot_imaging_outliers(outliers, condition_column, conditions)
    plot_imaging_marker_correlations(correlations)

    # --- Cell level (optional per-cell export) ---
    cells_path = data_dir / FACS_CELLS_FILE
    if cells_path.exists():
        cells_img = load_facs_cells(cells_path, metadata, condition_column)
        flags = cell_outlier_flags(cells_img)
        cell_summary = cell_level_summary(cells_img, flags, condition_column)
        cell_summary.to_csv(OUTPUT_DIR / "F_imaging_cell_summary.csv", index=False)
        plot_cell_distributions(cells_img)
        plot_cell_outliers(cell_summary)
        if len(IMAGING_FEATURES) >= 2:
            plot_morphology_map(cells_img, flags, condition_column, conditions)
    else:
        print(f"No per-cell imaging export ({FACS_CELLS_FILE}) found - cell-level plots F6-F8 skipped.")
    return summary


# =============================================================================
# 5. DEMO DATA (only used with --demo)
# =============================================================================
# Makes realistic-looking fake input files so the whole program can be tested
# before real data exists. The numbers are random but follow fixed patterns
# (e.g. more neutrophils in "Disease", fewer neutrophils and urothelial cells
# captured by scRNA-seq) so the statistics have something to find.
# The demo is built from individual CELL TYPES (the Flex labels); each FACS gate
# then gets the cells of the cell types listed under its population. The FlowJo
# export contains the gates of EVERY panel, so the demo works with either panel.

# Per cell type: atlas name, typical share (Healthy, Disease), how well scRNA-seq
# captures it, marker levels on FACS, and imaging features (area, eccentricity, radial moment)
DEMO_CELL_TYPES = {
    #  label                   atlas name         Healthy Disease capture  CD45  CD10  CD13  EpCAM  imaging
    "Neutrophil":           ("Neutrophil",       0.12, 0.28, 0.5, {"CD45": 7000, "CD10": 1500, "CD13": 3000, "EpCAM": 300},  [110, 0.35, 5.0]),
    "T cell":               ("T/NK cell",        0.08, 0.10, 1.2, {"CD45": 9500, "CD10": 300,  "CD13": 300,  "EpCAM": 300},  [70, 0.25, 3.5]),
    "Macrophage":           ("Myeloid cell",     0.05, 0.07, 1.0, {"CD45": 8000, "CD10": 400,  "CD13": 3500, "EpCAM": 300},  [150, 0.40, 5.5]),
    "Proximal tubule cell": ("Proximal tubule",  0.18, 0.13, 1.0, {"CD45": 250,  "CD10": 6000, "CD13": 8000, "EpCAM": 1200}, [210, 0.50, 7.0]),
    "Distal tubule cell":   ("Distal tubule",    0.12, 0.09, 1.0, {"CD45": 250,  "CD10": 600,  "CD13": 800,  "EpCAM": 6000}, [180, 0.45, 6.5]),
    "Collecting duct cell": ("Collecting duct",  0.08, 0.06, 1.0, {"CD45": 250,  "CD10": 500,  "CD13": 700,  "EpCAM": 7000}, [190, 0.45, 6.5]),
    "Urothelial cell":      ("Urothelial",       0.30, 0.20, 0.7, {"CD45": 250,  "CD10": 800,  "CD13": 900,  "EpCAM": 7000}, [340, 0.60, 9.5]),
    "Podocyte":             ("Podocyte",         0.03, 0.03, 1.0, {"CD45": 250,  "CD10": 500,  "CD13": 700,  "EpCAM": 300},  [160, 0.45, 5.5]),
    "Fibroblast":           ("Fibroblast",       0.04, 0.04, 1.0, {"CD45": 250,  "CD10": 700,  "CD13": 900,  "EpCAM": 300},  [170, 0.70, 6.0]),
}


def demo_gates():
    """Every population gate of every panel: gate path -> the demo cell types in it (one entry per gate)."""
    gates = {}
    for panel in PANELS.values():
        for info in panel["populations"].values():
            gates[info["gate"]] = [t for t in DEMO_CELL_TYPES if t in info["seq_labels"]]
    return gates


def make_demo_data(folder):
    """
    Write small, made-up input files in the exact format the program expects,
    so you can test the program before you have real data. As in a real study,
    the two methods measured DIFFERENT people:
      scRNA-seq : 12 samples from 8 people (donor_id P1-P8): P1-P4 'Healthy',
                  P5-P8 'Disease'; P1, P2, P5 and P6 gave two samples each
      FACS      : POOL = one pool of healthy urine from many donors ('Healthy pool')
                  + FACS01-FACS06 = 6 other patients ('Disease')
    Plus a reference atlas of 8 healthy donors and a per-cell imaging export.
    The FlowJo and per-cell files hold the populations of EVERY panel, so the
    demo can be analysed with either panel.
    """
    rng = np.random.default_rng(1)          # fixed seed -> the same demo data every time
    folder.mkdir(parents=True, exist_ok=True)
    types = list(DEMO_CELL_TYPES)
    share = {condition: np.array([DEMO_CELL_TYPES[t][1 + k] for t in types])
             for k, condition in enumerate(["Healthy", "Disease"])}
    capture = np.array([DEMO_CELL_TYPES[t][3] for t in types])
    protein = {t: DEMO_CELL_TYPES[t][4] for t in types}
    imaging = {t: DEMO_CELL_TYPES[t][5] for t in types}

    seq_samples = [f"SEQ{i:02d}" for i in range(1, 13)]
    seq_donors = ["P1", "P1", "P2", "P2", "P3", "P4", "P5", "P5", "P6", "P6", "P7", "P8"]
    seq_conditions = ["Healthy"] * 6 + ["Disease"] * 6
    facs_samples = ["POOL"] + [f"FACS{i:02d}" for i in range(1, 7)]
    facs_conditions = ["Healthy pool"] + ["Disease"] * 6
    # The sample sheet: one row per sample, with the method that measured it
    metadata = pd.concat([
        pd.DataFrame({SAMPLE_COLUMN: seq_samples, METHOD_COLUMN: "scRNA-seq", FACS_NAME_COLUMN: "",
                      DONOR_COLUMN: seq_donors, "condition": seq_conditions,
                      "flex_probe_barcode": [f"BC{i:03d}" for i in range(1, 13)]}),
        pd.DataFrame({SAMPLE_COLUMN: facs_samples, METHOD_COLUMN: "FACS",
                      FACS_NAME_COLUMN: [f"{s}_urine.fcs" for s in facs_samples], DONOR_COLUMN: "",
                      "condition": facs_conditions, "flex_probe_barcode": ""}),
    ])

    # --- FACS: one row per FACS sample, with the gates of every panel ---
    gates = demo_gates()
    facs_rows = []
    for sample, condition in zip(facs_samples, facs_conditions):
        biology = "Healthy" if condition.startswith("Healthy") else "Disease"
        # A pool averages many people, so it varies less around the typical values
        concentration = 300 if sample == "POOL" else 60
        fraction = dict(zip(types, rng.dirichlet(share[biology] * concentration)))
        live = int(rng.integers(20000, 60000))       # number of live single cells
        row = {"Name": f"{sample}_urine.fcs", COUNT_COLUMN.format(gate=FACS_TOTAL_GATE): live}
        marker_factor = {m: rng.lognormal(0, 0.25) for m in ALL_MARKERS}   # sample-to-sample variation
        for gate, members in gates.items():
            weights = np.array([fraction[t] for t in members])
            # Gate count (85% of live cells fall in one of the population gates)
            row[COUNT_COLUMN.format(gate=gate)] = int(live * weights.sum() * 0.85)
            weights = weights / weights.sum()
            # Marker medians: the mix of the cell types in the gate, with some random noise
            for marker, info in ALL_MARKERS.items():
                level = sum(w * protein[t][marker] for w, t in zip(weights, members))
                value = level * marker_factor[marker] * rng.lognormal(0, 0.1)
                row[MEDIAN_COLUMN.format(gate=gate, parameter=info["facs_parameter"])] = round(value, 1)
            # Imaging medians; proximal tubule cells are 15% larger in Disease
            for f, feature in enumerate(IMAGING_FEATURES):
                level = sum(w * imaging[t][f % 3] * (1.15 if biology == "Disease" and t == "Proximal tubule cell" else 1)
                            for w, t in zip(weights, members))
                row[MEDIAN_COLUMN.format(gate=gate, parameter=feature)] = round(level * rng.lognormal(0, 0.08), 3)
        facs_rows.append(row)

    # --- scRNA-seq: one table of cells per scRNA-seq sample ---
    # Every person has their own composition; their samples vary a little around it
    person_fraction = {}
    for donor, condition in dict(zip(seq_donors, seq_conditions)).items():
        expected = share[condition] * capture          # as scRNA-seq captures the cell types
        person_fraction[donor] = rng.dirichlet(np.append(expected / expected.sum(), 0.05) * 60)
    cell_tables = []
    for sample, donor in zip(seq_samples, seq_donors):
        # Cell types of this sample (plus a few 'Other' cells)
        fraction = rng.dirichlet(person_fraction[donor] * 300)
        n_cells = int(rng.integers(800, 2500))
        marker_factor = {m: rng.lognormal(0, 0.25) for m in ALL_MARKERS}
        labels = rng.choice(types + ["Other"], size=n_cells, p=fraction)
        labels = np.where(labels == "Other", rng.choice(["Unknown", "Doublet"], size=n_cells), labels)
        cells = pd.DataFrame({"cell_barcode": [f"{sample}_{i:05d}" for i in range(n_cells)],
                              SAMPLE_COLUMN: sample, CELL_TYPE_COLUMN: labels})
        # Gene expression that follows the protein level of the cell type (all four genes)
        for marker, info in ALL_MARKERS.items():
            level = {t: np.log1p(protein[t][marker] / 1500 * marker_factor[marker]) for t in types}
            mean = np.array([level.get(label, 0.3) for label in labels])
            cells[info["gene"]] = np.clip(rng.normal(mean, 0.5), 0, None).round(3)
        cell_tables.append(cells)

    # --- Reference atlas: 8 healthy donors with known cell types (per-cell format) ---
    atlas_tables = []
    for d in range(1, 9):
        # Healthy-like composition with donor-to-donor variation, plus a few unassigned cells
        fraction = rng.dirichlet(np.append(share["Healthy"], 0.03) * 60)
        n_atlas = int(rng.integers(1500, 4000))
        atlas_type = rng.choice(types + ["Other"], size=n_atlas, p=fraction)
        # Atlas names for the cell types (they differ from the notebook labels)
        labels = [DEMO_CELL_TYPES[t][0] if t != "Other" else "Unknown" for t in atlas_type]
        atlas = pd.DataFrame({ATLAS_DONOR_COLUMN: f"D{d:02d}", ATLAS_CELL_TYPE_COLUMN: labels})
        for marker, info in ALL_MARKERS.items():
            level = {t: np.log1p(protein[t][marker] / 1500) for t in types}
            mean = np.array([level.get(t, 0.3) for t in atlas_type])
            atlas[info["gene"]] = np.clip(rng.normal(mean, 0.5), 0, None).round(3)
        atlas_tables.append(atlas)

    # --- Per-cell imaging export: the single cells of every population gate (section F) ---
    # (the populations of all panels; a real export only has those of your panel)
    all_populations = {pop: info for panel in PANELS.values() for pop, info in panel["populations"].items()}
    cell_rng = np.random.default_rng(2)
    imaging_cells = []
    for row in facs_rows:
        for pop, info in all_populations.items():
            n_cells_gate = int(cell_rng.integers(150, 300))
            part = {FACS_NAME_COLUMN: row["Name"], FACS_CELLS_POPULATION_COLUMN: pop}
            # About 2% odd cells (e.g. doublets or clumps) that sit far from the rest
            odd = cell_rng.random(n_cells_gate) < 0.02
            for feature in IMAGING_FEATURES:
                # Cells scatter around the sample's median from the FlowJo table
                median = row[MEDIAN_COLUMN.format(gate=info["gate"], parameter=feature)]
                values = median * cell_rng.lognormal(0, 0.12, n_cells_gate)
                values[odd] *= 1.9
                if "Eccentricity" in feature:          # eccentricity is between 0 and 1
                    values = np.clip(values, 0, 0.99)
                part[feature] = values.round(3)
            imaging_cells.append(pd.DataFrame(part))

    # Write the files in the same formats as the real inputs
    metadata.to_csv(folder / METADATA_FILE, index=False)
    pd.DataFrame(facs_rows).to_excel(folder / FACS_FILE, index=False)
    pd.concat(cell_tables).to_csv(folder / SEQ_FILE, index=False)
    pd.concat(atlas_tables).to_csv(folder / ATLAS_FILE, index=False)
    pd.concat(imaging_cells).to_csv(folder / FACS_CELLS_FILE, index=False)
    print(f"Demo data written to {folder}")


# =============================================================================
# 6. MAIN PROGRAM - runs the steps above in order
# =============================================================================
# run_analysis() does the actual work and is shared by the terminal (main) and
# the web app (app.py). main() only handles the command-line options.

def run_analysis(data_dir, output_dir, condition_column=None, seq_conditions=None, facs_conditions=None,
                 panel=None):
    """
    Run steps 1-5 on the input files in data_dir, plus section E (if an atlas
    file is there) and section F (if imaging features are set), and write every
    table and plot to output_dir. Used by the command line (main, below) and by the web
    interface (app.py). Anything not given (condition column, the scRNA-seq
    conditions, the FACS conditions) is asked in the terminal.
    panel = which FACS panel was used (a key of PANELS); None = FACS_PANEL from SETTINGS.
    Returns the summary tables so the caller can show them.
    """
    # OUTPUT_DIR is changed for the whole script so save() writes to the right folder
    global OUTPUT_DIR
    data_dir = Path(data_dir)
    OUTPUT_DIR = Path(output_dir)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)   # create the results folder if needed
    # Switch the gating, the Flex labels and the markers to the chosen panel
    use_panel(panel or FACS_PANEL)
    print(f"FACS panel: {PANELS[ACTIVE_PANEL]['description']} -> populations: {', '.join(POPULATIONS)}")

    # --- Step 1: load files and choose conditions ---
    metadata = load_metadata(data_dir / METADATA_FILE)
    condition_column, conditions = choose_conditions(metadata, condition_column, seq_conditions, facs_conditions)
    # Keep only the samples of the chosen conditions, separately for each method
    chosen = pd.concat([samples_of(metadata, m).loc[lambda d: d[condition_column].isin(conditions[m])]
                        for m in METHODS], ignore_index=True)
    seq_metadata, facs_metadata = samples_of(chosen, "scRNA-seq"), samples_of(chosen, "FACS")
    facs = load_facs(data_dir / FACS_FILE, facs_metadata)
    cells = load_seq(data_dir / SEQ_FILE, seq_metadata)
    if facs.empty:
        raise SystemExit(f"None of the chosen FACS samples were found in {FACS_FILE}. "
                         f"Check the '{FACS_NAME_COLUMN}' column of the sample sheet.")
    if cells.empty:
        raise SystemExit(f"None of the chosen scRNA-seq samples have cells in {SEQ_FILE}. "
                         f"Check the '{SAMPLE_COLUMN}' column of both files.")
    print(f"Loaded {len(facs)} FACS samples and {len(cells)} cells "
          f"from {cells[SAMPLE_COLUMN].nunique()} scRNA-seq samples.")

    # --- Step 2: build tables ---
    # Proportions per sample, then per PERSON (samples of the same person combined):
    # every statistic and plot below counts each person once
    sample_proportions = build_proportions(facs, cells, chosen, condition_column)
    proportions = combine_proportions_per_donor(sample_proportions, chosen, condition_column)
    # Marker levels per sample (FACS protein is also used per sample in section F) and per person (section C)
    facs_markers = build_facs_markers(facs, facs_metadata, condition_column)
    rna_markers = build_rna_markers(cells, seq_metadata, condition_column)
    keys = ["population", "marker", "gene", condition_column]
    protein_rna = build_protein_vs_rna(
        combine_per_donor(facs_markers, facs_metadata, keys, {"protein_median": ("protein_median", "mean")}),
        combine_per_donor(rna_markers, seq_metadata, keys, {"rna_mean": ("rna_mean", "mean")}))
    # Imaging table only if imaging features are listed in SETTINGS
    imaging = build_imaging(facs, facs_metadata, condition_column) if IMAGING_FEATURES else None
    # The reference atlas (optional): loaded here because section A's overview shows it too
    atlas = atlas_props = None
    atlas_path = data_dir / ATLAS_FILE
    if atlas_path.exists():
        atlas = load_atlas(atlas_path)
        atlas_props = build_atlas_proportions(atlas)

    # --- Step 3: statistics ---
    overview = composition_overview(proportions, condition_column, conditions, atlas_props)  # section A
    comparison = method_comparison_stats(proportions, condition_column, conditions)       # section A
    by_condition = condition_stats(proportions, "proportion", ["population", "method"],
                                   condition_column, all_conditions(conditions))           # section B
    rna_stats = protein_rna_stats(protein_rna)                                              # section C

    # --- Step 4: save tables ---
    print("\nSaving results to", OUTPUT_DIR)
    sample_proportions.to_csv(OUTPUT_DIR / "A_proportions_per_sample.csv", index=False)
    proportions.to_csv(OUTPUT_DIR / "A_proportions.csv", index=False)
    overview.to_csv(OUTPUT_DIR / "A_composition_overview.csv", index=False)
    comparison.to_csv(OUTPUT_DIR / "A_method_comparison_stats.csv", index=False)
    by_condition.to_csv(OUTPUT_DIR / "B_condition_tests_proportions.csv", index=False)
    facs_markers.to_csv(OUTPUT_DIR / "C_protein_per_sample.csv", index=False)
    rna_markers.to_csv(OUTPUT_DIR / "C_rna_per_sample.csv", index=False)
    protein_rna.to_csv(OUTPUT_DIR / "C_protein_vs_rna.csv", index=False)
    rna_stats.to_csv(OUTPUT_DIR / "C_protein_vs_rna_stats.csv", index=False)

    # --- Step 5: plots ---
    plot_stacked_proportions(proportions, condition_column, conditions)
    plot_method_boxes(proportions, condition_column, conditions)
    plot_method_difference(comparison)
    plot_composition_overview(overview)
    # Comparing conditions only makes sense when a method has at least two of them
    if any(len(conditions[m]) > 1 for m in METHODS):
        plot_proportions_by_condition(proportions, condition_column, all_conditions(conditions))
    plot_protein_vs_rna(protein_rna, rna_stats)

    # --- Section D: imaging statistics, tables and plot (FACS only; skipped if no features) ---
    imaging_pop = imaging_cond = None
    if imaging is not None:
        imaging_pop = imaging_between_populations(imaging)
        imaging_cond = condition_stats(imaging, "value", ["feature", "population"],
                                       condition_column, conditions["FACS"])
        imaging.to_csv(OUTPUT_DIR / "D_imaging_features.csv", index=False)
        imaging_pop.to_csv(OUTPUT_DIR / "D_imaging_tests_between_populations.csv", index=False)
        imaging_cond.to_csv(OUTPUT_DIR / "D_imaging_tests_between_conditions.csv", index=False)
        plot_imaging(imaging, condition_column, conditions["FACS"])

    # --- Section E: comparison with the reference atlas (skipped if no atlas file) ---
    atlas_tests = None
    if atlas is not None:
        print("\nComparing with the reference atlas ...")
        atlas_tests = run_atlas_comparison(atlas, atlas_props, proportions, cells, facs_markers,
                                           condition_column, conditions)
    else:
        print(f"\nNo reference atlas ({ATLAS_FILE}) found - section E skipped.")

    # --- Section F: imaging profiles per cell type (FACS only; skipped if no imaging features) ---
    imaging_summary = None
    if imaging is not None:
        print("\nImaging profiles per cell type ...")
        imaging_summary = run_imaging_profiles(imaging, facs_markers, data_dir, facs_metadata,
                                               condition_column, conditions["FACS"])

    # Hand the main results back to whoever called this (terminal or web app)
    return {
        "imaging_summary": imaging_summary,
        "atlas_tests": atlas_tests,
        "panel": ACTIVE_PANEL,
        "condition_column": condition_column,
        "conditions": conditions,
        "n_facs_samples": len(facs),
        "n_cells": len(cells),
        "n_seq_samples": cells[SAMPLE_COLUMN].nunique(),
        "n_seq_people": proportions.loc[proportions["method"] == "scRNA-seq", SAMPLE_COLUMN].nunique(),
        "overview": overview,
        "comparison": comparison,
        "by_condition": by_condition,
        "rna_stats": rna_stats,
        "imaging_pop": imaging_pop,
        "imaging_cond": imaging_cond,
    }


def main():
    # Define the command-line options the program understands
    parser = argparse.ArgumentParser(description="Compare FACS (BD S8) and 10x Flex scRNA-seq cell populations.")
    parser.add_argument("--demo", action="store_true", help="create and analyse made-up demo data")
    parser.add_argument("--condition-column", help="metadata column that holds the conditions")
    parser.add_argument("--seq-conditions", nargs="+", help="conditions of the scRNA-seq samples to include (default: ask)")
    parser.add_argument("--facs-conditions", nargs="+", help="conditions of the FACS samples to include (default: ask)")
    parser.add_argument("--panel", choices=list(PANELS),
                        help=f"antibodies used on FACS (default: FACS_PANEL in SETTINGS = {FACS_PANEL})")
    args = parser.parse_args()

    # Use the real data folders, or create demo data and use demo folders instead
    data_dir, output_dir = DATA_DIR, OUTPUT_DIR
    if args.demo:
        data_dir = SCRIPT_DIR / "demo_data"
        output_dir = SCRIPT_DIR / "results_demo"
        make_demo_data(data_dir)

    # Run the whole analysis (steps 1-5, plus E and F)
    results = run_analysis(data_dir, output_dir, args.condition_column, args.seq_conditions, args.facs_conditions,
                           args.panel)
    comparison, by_condition = results["comparison"], results["by_condition"]

    # --- Step 6: short summary in the terminal ---
    # How the scRNA-seq group(s) compare with the FACS group(s), per population
    print("\nscRNA-seq vs FACS (proportions, median per group):")
    print(comparison[["comparison", "population", "n_seq", "n_facs", "median_seq", "median_facs",
                      "seq_minus_facs", "test", "q_value_fdr"]].round(3).to_string(index=False))
    # Which condition differences are significant after FDR correction
    if not by_condition.empty:
        hits = by_condition[by_condition["significant"]]
        print(f"\nSignificant condition differences in proportions (FDR < {ALPHA}): {len(hits)}")
        if not hits.empty:
            print(hits[["population", "method", "test", "comparison", "q_value_fdr"]]
                  .round(4).to_string(index=False))
    # Comparison with the reference atlas (section E), if an atlas was given
    atlas_tests = results["atlas_tests"]
    if atlas_tests is not None:
        print("\nProportions compared with the reference atlas (median of samples vs atlas):")
        print(atlas_tests[["population", "method", "group", "median_samples", "atlas_proportion",
                           "q_value_fdr", "significant"]].round(3).to_string(index=False))
    print("\nDone.")


# Run main() only when this file is started directly (not when app.py imports it)
if __name__ == "__main__":
    main()
