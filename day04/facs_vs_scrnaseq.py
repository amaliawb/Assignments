"""
facs_vs_scrnaseq.py
===================

Compare urine cell populations measured two ways on the SAME samples:

  1. Flow cytometry / imaging sorting on the BD FACSDiscover S8
     (gated in FlowJo, exported as an Excel table)
  2. 10x Genomics Flex single-cell RNA-seq
     (processed and annotated in a Jupyter notebook, exported as a CSV)

FACS gating:
  Cells -> Single Cells -> Live (Fixable Yellow negative)
    Live  -> CD45+ / CD45-
    CD45- -> CD10+ CD13+ / NOT (CD10+ CD13+)
    NOT (CD10+ CD13+) -> EpCAM+ / EpCAM-

Populations (they do not overlap: each live cell falls in exactly one,
so the proportions add up to 100% of live cells):
  - Immune cells           : CD45+
  - Proximal tubule cells  : CD45- CD10+ CD13+          (EpCAM+ or EpCAM-)
  - Epithelial cells       : CD45- NOT(CD10+ CD13+) EpCAM+
  - Other CD45- cells      : CD45- NOT(CD10+ CD13+) EpCAM-

What the program produces (in the results folder):
  A. Cell-type PROPORTIONS per sample, FACS vs scRNA-seq
       - agreement between methods: Spearman + Pearson correlation,
         paired Wilcoxon signed-rank test, Bland-Altman bias / limits of agreement
  B. Differences between CONDITIONS (chosen when you run the program)
       - 2 conditions  : Mann-Whitney U test
       - 3+ conditions : Kruskal-Wallis test + pairwise Mann-Whitney U tests
       - all p-values corrected for multiple testing (Benjamini-Hochberg FDR)
  C. PROTEIN vs RNA: FACS median intensity of each marker vs. mean RNA
     expression of its gene (CD45-PTPRC, CD10-MME, CD13-ANPEP, EpCAM-EPCAM)
  D. S8 IMAGING features per population: differences between populations
     and between conditions
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

  With your real data (after editing the SETTINGS section below):
      python facs_vs_scrnaseq.py

  The program will ask which metadata column holds the conditions and which
  conditions to include. You can skip the questions by giving them directly:
      python facs_vs_scrnaseq.py --condition-column condition --conditions Healthy Disease

  Or use the web interface instead of the terminal (see app.py):
      python app.py        then open http://127.0.0.1:5000

  Needs: pandas, numpy, scipy, matplotlib, openpyxl (all included in Anaconda).

-------------------------------------------------------------------------------
INPUT FILES
-------------------------------------------------------------------------------
1) metadata.csv - one row per sample. Required columns:
      sample_id   : the sample name used in the scRNA-seq file
      facs_name   : the sample name exactly as it appears in the FlowJo export
                    (usually the .fcs file name)
   Plus any number of columns describing the samples (condition, donor, sex,
   flex_probe_barcode, ...). You choose the condition column when you run.

2) FlowJo Excel export (Table Editor -> export as .xlsx). One row per sample.
   The first column holds the sample name; the other columns are statistics:
      - the Count of the total gate (e.g. live single cells)
      - the Count of each population gate
      - the Median of each marker channel inside each population gate
      - the Median of each S8 imaging feature inside each population gate
   Column names are built from the templates in SETTINGS (COUNT_COLUMN,
   MEDIAN_COLUMN). If a name does not match, the program prints the column
   names it found in your file so you can fix the SETTINGS.

3) scrnaseq_cells.csv - one row per cell, exported from your notebook.
   Columns: cell_barcode, sample_id, cell_type, PTPRC, MME, ANPEP, EPCAM
   (cells already demultiplexed by Flex probe barcode -> sample_id, and
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
import argparse            # reads options typed on the command line (--demo, --conditions, ...)
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
FACS_NAME_COLUMN = "facs_name"   # the sample name as FlowJo writes it
CELL_TYPE_COLUMN = "cell_type"   # the annotation label of each sequenced cell

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

# ---- Gating hierarchy ------------------------------------------------------------
#   Cells -> Single Cells -> Live -> CD45+                                (Immune)
#                                 -> CD45- -> CD10+ CD13+                  (Proximal tubule)
#                                          -> NOT CD10+ CD13+ -> EpCAM+    (Epithelial)
#                                                             -> EpCAM-    (Other CD45-)
# "Live" = cells NEGATIVE for the Fixable Yellow live/dead dye (dead cells take
# up the dye and are bright).
# Why CD10/CD13 is gated BEFORE EpCAM: proximal tubule cells can express EpCAM.
# If the EpCAM+ and CD10+ CD13+ gates both sat directly under CD45-, a cell
# positive for both would be counted twice and the proportions would be wrong.
# "NOT CD10+ CD13+" is the inverse of the CD10+ CD13+ gate (FlowJo: Boolean
# NOT gate, or Make Gate -> Invert). Use the gate names exactly as in FlowJo.
LIVE_GATE = "Cells/Single Cells/Live"     # live/dead gate on Fixable Yellow
VIABILITY_DYE = "Fixable Yellow"          # only for documentation; not used in calculations
CD45_NEG_GATE = f"{LIVE_GATE}/CD45-"
NOT_PT_GATE = f"{CD45_NEG_GATE}/NOT CD10+ CD13+"

# Gate that holds ALL good cells (used only if PROPORTION_DENOMINATOR = "total")
FACS_TOTAL_GATE = LIVE_GATE

# ---- The four populations -------------------------------------------------------
# "gate"       : the FlowJo gate path of the population
# "seq_labels" : which cell_type labels from your notebook belong to this
#                population (list ALL the labels you used, e.g. every immune subtype)
# "atlas_labels": which cell_type names of the REFERENCE ATLAS belong to this
#                population (an atlas usually names cell types differently from
#                your notebook). If left out, the seq_labels are used.
#                These are placeholders - copy the real names from your atlas.
# The four gates do not overlap and together cover every live cell, so each
# cell is counted exactly once (like in scRNA-seq, where each cell has one label).
POPULATIONS = {
    "Immune": {
        "gate": f"{LIVE_GATE}/CD45+",
        "seq_labels": ["T cell", "Macrophage", "Neutrophil"],
        "atlas_labels": ["T/NK cell", "Myeloid cell", "Neutrophil"],
    },
    "Proximal tubule": {
        "gate": f"{CD45_NEG_GATE}/CD10+ CD13+",
        "seq_labels": ["Proximal tubule cell"],
        "atlas_labels": ["Proximal tubule"],
    },
    "Epithelial": {
        "gate": f"{NOT_PT_GATE}/EpCAM+",
        "seq_labels": ["Urothelial cell", "Squamous epithelial cell"],
        "atlas_labels": ["Urothelial", "Squamous"],
    },
    "Other CD45-": {
        "gate": f"{NOT_PT_GATE}/EpCAM-",
        # Placeholder labels: every notebook label that is not immune, not proximal
        # tubule and not EpCAM+ epithelial (e.g. podocytes, fibroblasts)
        "seq_labels": ["Podocyte", "Fibroblast"],
        "atlas_labels": ["Podocyte", "Fibroblast"],
    },
}

# ---- Antibody markers and their matching genes ------------------------------
# "facs_parameter" is the channel name FlowJo uses. Replace FLUOR1..4 with your
# fluorochromes, e.g. "Comp-BV510-A :: CD45" - copy it exactly from FlowJo.
# "gene" is the column name of the same marker in the scRNA-seq file.
MARKERS = {
    "CD45":  {"gene": "PTPRC", "facs_parameter": "Comp-FLUOR1-A :: CD45"},
    "CD10":  {"gene": "MME",   "facs_parameter": "Comp-FLUOR2-A :: CD10"},
    "CD13":  {"gene": "ANPEP", "facs_parameter": "Comp-FLUOR3-A :: CD13"},
    "EpCAM": {"gene": "EPCAM", "facs_parameter": "Comp-FLUOR4-A :: EpCAM"},
}

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
#   "populations" : population / (sum of the four populations)    <- recommended
#                   both methods then share the same denominator, which removes
#                   differences caused by debris, dead cells, unlabelled cells, etc.
#   "total"       : FACS -> population / FACS_TOTAL_GATE count
#                   seq  -> population / all annotated cells
PROPORTION_DENOMINATOR = "populations"

MIN_SAMPLES_PER_GROUP = 3     # fewer samples than this -> test is skipped (reported as NaN)
MIN_CELLS_FOR_RNA_MEAN = 10   # a population with fewer cells in a sample gets no RNA mean
ALPHA = 0.05                  # significance threshold applied to FDR-corrected p-values

# ---- Plot colors --------------------------------------------------------------
# One color per population (same order as POPULATIONS), plus one per method.
POPULATION_COLORS = ["#2a78d6", "#1baf7a", "#eb6834", "#8a63d2"]
METHOD_COLORS = {"FACS": "#52514e", "scRNA-seq": "#b5b3ad"}
ATLAS_COLOR = "#c2477a"       # the reference atlas in the section E plots
PLOT_DPI = 200                # image resolution (dots per inch) of the saved plots


# =============================================================================
# 1. LOADING THE DATA
# =============================================================================
# These functions read the three main input files, check that the expected columns
# are there, and match the samples between the files. (The optional reference
# atlas is read separately, in section 4b.)

def get_column(table, column, file_label):
    """Return one column of a table; if it does not exist, stop and list the real column names."""
    # If the column name is wrong (typo, different gate name in FlowJo, ...)
    # stop the program with a helpful message instead of a cryptic KeyError.
    if column not in table.columns:
        # Put every real column name on its own line so they are easy to compare
        available = "\n    ".join(map(str, table.columns))
        raise SystemExit(
            f"\nColumn not found in {file_label}:\n    '{column}'\n"
            f"Columns that ARE in the file:\n    {available}\n"
            "-> Fix the matching name in the SETTINGS section of this script."
        )
    # The column exists: hand it back
    return table[column]


def load_metadata(path):
    """Read the sample sheet. Everything is read as text so IDs like '01' stay intact."""
    # dtype=str keeps sample IDs as text (otherwise '01' would become the number 1)
    metadata = pd.read_csv(path, dtype=str)
    # Make sure the two required columns exist (stops with a message if not)
    get_column(metadata, SAMPLE_COLUMN, path.name)
    get_column(metadata, FACS_NAME_COLUMN, path.name)
    return metadata


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


def choose_conditions(metadata, column=None, values=None):
    """
    Decide which metadata column is the 'condition' and which of its values to analyse.
    Asks the user in the terminal unless both were given on the command line.
    """
    # Every metadata column except the two ID columns could hold the conditions
    candidates = [c for c in metadata.columns if c not in (SAMPLE_COLUMN, FACS_NAME_COLUMN)]
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

    # --- which values? ---
    # Count how many samples each condition has (keeps the order of the file)
    counts = metadata[column].value_counts(sort=False)
    available = list(counts.index)
    if values is None:
        # Show the conditions with their sample counts and let the user choose
        print(f"\nConditions found in '{column}':")
        for i, name in enumerate(available, start=1):
            print(f"  {i}. {name}  ({counts[name]} samples)")
        picks = ask_for_numbers("Type the numbers to include separated by commas (Enter = all): ",
                                len(available), allow_many=True)
        # None = Enter was pressed = use every condition
        values = available if picks is None else [available[p] for p in picks]
    else:
        # Conditions were given on the command line: check they all exist
        missing = [v for v in values if v not in available]
        if missing:
            raise SystemExit(f"Conditions {missing} not found in '{column}'. Available: {available}")

    print(f"\nUsing column '{column}' with conditions: {', '.join(values)}")
    return column, list(values)


def load_facs(path, metadata):
    """
    Read the FlowJo Excel export and return one row per sample (index = sample_id).
    Rows FlowJo adds itself (e.g. 'Mean', 'SD') are dropped because they are not
    listed in the metadata.
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
    """
    # Read the per-cell table; sample IDs stay text (same as in the metadata)
    cells = pd.read_csv(path, dtype={SAMPLE_COLUMN: str})
    # Check that the sample, cell type and the four gene columns are all there
    get_column(cells, SAMPLE_COLUMN, path.name)
    get_column(cells, CELL_TYPE_COLUMN, path.name)
    for marker in MARKERS.values():
        get_column(cells, marker["gene"], path.name)

    # Dictionary: notebook label -> population, e.g. "T cell" -> "Immune"
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
    return cells[cells[SAMPLE_COLUMN].isin(metadata[SAMPLE_COLUMN])]


# =============================================================================
# 2. BUILDING THE COMPARISON TABLES
# =============================================================================
# These functions turn the raw inputs into "long" tables (one measurement per
# row) that the statistics and plots below can use directly.

def build_proportions(facs, cells, metadata, condition_column):
    """
    Long table: one row per sample x population x method, with the cell count
    and the proportion. This is the table used for sections A and B.
    """
    populations = list(POPULATIONS)

    # FACS counts per population
    # (one column per population, taken from that population's gate "Count" column)
    facs_counts = pd.DataFrame({
        pop: get_column(facs, COUNT_COLUMN.format(gate=info["gate"]), "the FACS file")
        for pop, info in POPULATIONS.items()
    })
    # scRNA-seq counts per population (cells per sample per population)
    seq_counts = pd.crosstab(cells[SAMPLE_COLUMN], cells["population"])
    # Make sure every population has a column, even if it had 0 cells
    seq_counts = seq_counts.reindex(columns=populations + ["Other"], fill_value=0)

    # Choose what each count is divided by (see PROPORTION_DENOMINATOR in SETTINGS)
    if PROPORTION_DENOMINATOR == "populations":
        # sum of the four populations, per sample (= all live cells for FACS)
        facs_total = facs_counts[populations].sum(axis=1)
        seq_total = seq_counts[populations].sum(axis=1)
    elif PROPORTION_DENOMINATOR == "total":
        # all live single cells (FACS) / all cells including "Other" (seq)
        facs_total = get_column(facs, COUNT_COLUMN.format(gate=FACS_TOTAL_GATE), "the FACS file")
        seq_total = seq_counts.sum(axis=1)
    else:
        raise SystemExit("PROPORTION_DENOMINATOR must be 'populations' or 'total'.")

    # Build one row for every sample x population, first for FACS, then for scRNA-seq
    rows = []
    for method, counts, total in [("FACS", facs_counts, facs_total),
                                  ("scRNA-seq", seq_counts, seq_total)]:
        for sample in counts.index:
            for pop in populations:
                rows.append({
                    SAMPLE_COLUMN: sample,
                    "population": pop,
                    "method": method,
                    "count": counts.loc[sample, pop],
                    # proportion = count / total (NaN if the sample has no cells at all)
                    "proportion": counts.loc[sample, pop] / total[sample] if total[sample] > 0 else np.nan,
                })
    table = pd.DataFrame(rows)
    # Add the condition of each sample (e.g. Healthy / Disease) as a column
    return table.merge(metadata[[SAMPLE_COLUMN, condition_column]], on=SAMPLE_COLUMN)


def build_protein_vs_rna(facs, cells, metadata, condition_column):
    """
    One row per sample x population x marker:
      protein = FACS median intensity of the marker inside the population gate
      rna     = mean log-normalized expression of the gene in that population's cells
    """
    # Mean expression of each of the four genes, per sample and population
    rna_means = cells.groupby([SAMPLE_COLUMN, "population"])[
        [m["gene"] for m in MARKERS.values()]].mean()
    # Number of sequenced cells per sample and population
    cell_numbers = cells.groupby([SAMPLE_COLUMN, "population"]).size()

    rows = []
    # Go through every population and every marker
    for pop, pop_info in POPULATIONS.items():
        for marker, marker_info in MARKERS.items():
            # FACS: the "Median (marker channel)" column of this population's gate
            column = MEDIAN_COLUMN.format(gate=pop_info["gate"], parameter=marker_info["facs_parameter"])
            protein = get_column(facs, column, "the FACS file")
            for sample in facs.index:
                key = (sample, pop)
                # A mean from very few cells is unreliable, so require a minimum
                enough_cells = cell_numbers.get(key, 0) >= MIN_CELLS_FOR_RNA_MEAN
                rows.append({
                    SAMPLE_COLUMN: sample,
                    "population": pop,
                    "marker": marker,
                    "gene": marker_info["gene"],
                    "protein_median": protein[sample],
                    "rna_mean": rna_means.loc[key, marker_info["gene"]] if enough_cells else np.nan,
                    "n_cells_rna": cell_numbers.get(key, 0),
                })
    table = pd.DataFrame(rows)
    # Add the condition of each sample as a column
    return table.merge(metadata[[SAMPLE_COLUMN, condition_column]], on=SAMPLE_COLUMN)


def build_imaging(facs, metadata, condition_column):
    """One row per sample x population x imaging feature (FACS median value)."""
    rows = []
    # For every population and every S8 imaging feature ...
    for pop, info in POPULATIONS.items():
        for feature in IMAGING_FEATURES:
            # ... take the "Median (feature)" column of that population's gate
            column = MEDIAN_COLUMN.format(gate=info["gate"], parameter=feature)
            values = get_column(facs, column, "the FACS file")
            # ... and add one row per sample
            for sample in facs.index:
                rows.append({SAMPLE_COLUMN: sample, "population": pop,
                             "feature": feature, "value": values[sample]})
    table = pd.DataFrame(rows)
    # Add the condition of each sample as a column
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
        result["pearson_r"], result["pearson_p"] = stats.pearsonr(x, y)
    return result


def compare_groups(values_by_group):
    """
    Compare a measurement between groups (conditions or populations).
      2 groups  -> Mann-Whitney U
      3+ groups -> Kruskal-Wallis (overall) + every pair with Mann-Whitney U
    Groups with fewer than MIN_SAMPLES_PER_GROUP values are left out.
    """
    # Remove missing values from every group
    groups = {name: np.asarray(v, float)[~np.isnan(np.asarray(v, float))]
              for name, v in values_by_group.items()}
    # Drop groups that are too small to test
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


def method_agreement_stats(proportions):
    """
    Section A: do FACS and scRNA-seq give the same proportions?
      - correlation : do samples rank / scale the same way in both methods?
      - Wilcoxon    : is one method systematically higher (paired by sample)?
      - Bland-Altman: average difference (bias) and the range where 95% of
                      sample differences fall (limits of agreement)
    """
    # Reshape so each row has the FACS and the scRNA-seq proportion side by side
    wide = proportions.pivot_table(index=[SAMPLE_COLUMN, "population"],
                                   columns="method", values="proportion").reset_index()
    rows = []
    for pop in POPULATIONS:
        # Only samples that have a value from BOTH methods can be paired
        part = wide[wide["population"] == pop].dropna(subset=["FACS", "scRNA-seq"])
        facs_values, seq_values = part["FACS"].values, part["scRNA-seq"].values
        difference = seq_values - facs_values          # positive = seq gives more
        row = {"population": pop, **correlate(facs_values, seq_values),
               "mean_FACS": facs_values.mean(), "mean_scRNA-seq": seq_values.mean(),
               # Bland-Altman: mean difference and mean +/- 1.96 SD
               "bland_altman_bias": difference.mean(),
               "bland_altman_lower": difference.mean() - 1.96 * difference.std(ddof=1),
               "bland_altman_upper": difference.mean() + 1.96 * difference.std(ddof=1),
               "wilcoxon_stat": np.nan, "wilcoxon_p": np.nan}
        # Paired Wilcoxon test (needs enough samples and at least one non-zero difference)
        if len(part) >= MIN_SAMPLES_PER_GROUP and np.any(difference != 0):
            row["wilcoxon_stat"], row["wilcoxon_p"] = stats.wilcoxon(facs_values, seq_values)
        rows.append(row)

    table = pd.DataFrame(rows)
    # Correct each kind of p-value for testing several populations
    for p_col in ["spearman_p", "pearson_p", "wilcoxon_p"]:
        table[p_col.replace("_p", "_q_fdr")] = fdr_correct(table[p_col])
    return table


def condition_stats(table, value_column, split_columns, condition_column, conditions):
    """
    Section B (and D): for every combination of split_columns (e.g. population x
    method) compare value_column between the selected conditions. FDR correction
    is applied over the whole table.
    """
    rows = []
    # Go through each group, e.g. ("Immune", "FACS"), ("Immune", "scRNA-seq"), ...
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


def protein_rna_stats(protein_rna):
    """
    Section C: correlation between protein (FACS) and RNA (scRNA-seq) for each marker,
    once with all populations pooled and once inside each population.
    """
    rows = []
    for marker, part in protein_rna.groupby("marker", sort=False):
        # All populations together (mostly shows that the marker separates populations)
        rows.append({"marker": marker, "population": "All populations",
                     **correlate(part["protein_median"], part["rna_mean"])})
        # Inside each population (does protein follow RNA between samples?)
        for pop, sub in part.groupby("population", sort=False):
            rows.append({"marker": marker, "population": pop,
                         **correlate(sub["protein_median"], sub["rna_mean"])})
    table = pd.DataFrame(rows)
    # Correct the p-values for all correlations in the table
    table["spearman_q_fdr"] = fdr_correct(table["spearman_p"])
    table["pearson_q_fdr"] = fdr_correct(table["pearson_p"])
    return table


def imaging_between_populations(imaging):
    """Section D: do the imaging features differ between the four populations?"""
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


def plot_stacked_proportions(proportions, condition_column, conditions):
    """Stacked bars: each sample gets a FACS bar and a scRNA-seq bar side by side."""
    # Order samples by condition (in the order chosen), then by sample ID
    order = (proportions[[SAMPLE_COLUMN, condition_column]].drop_duplicates()
             .assign(_c=lambda d: d[condition_column].map(conditions.index))
             .sort_values(["_c", SAMPLE_COLUMN]))
    samples = list(order[SAMPLE_COLUMN])
    # Table: rows = (sample, method), columns = populations, values = proportions
    wide = proportions.pivot_table(index=[SAMPLE_COLUMN, "method"], columns="population",
                                   values="proportion")

    # Figure gets wider when there are more samples
    fig, ax = plt.subplots(figsize=(max(6, len(samples) * 0.9), 4.5))
    positions, labels = [], []
    for i, sample in enumerate(samples):
        for j, method in enumerate(["FACS", "scRNA-seq"]):
            x = i * 3 + j          # two bars per sample, then a gap
            bottom = 0
            # Stack the populations on top of each other
            for pop in POPULATIONS:
                value = wide.loc[(sample, method), pop] if (sample, method) in wide.index else 0
                ax.bar(x, value, bottom=bottom, width=0.85, color=population_color(pop),
                       edgecolor="white", linewidth=1, label=pop if (i == 0 and j == 0) else None)
                bottom += value
            positions.append(x)
            labels.append(f"{sample} · {'FACS' if method == 'FACS' else 'seq'}")
    # Axis labels, title and legend
    ax.set_xticks(positions, labels, rotation=90)
    ax.set_ylabel("Proportion of cells")
    ax.set_title("Cell-type proportions per sample (FACS vs scRNA-seq)", loc="left")
    ax.grid(axis="x", visible=False)
    ax.legend(loc="upper left", bbox_to_anchor=(1, 1))
    save(fig, "A1_proportions_stacked.png")


def plot_method_scatter(proportions, agreement):
    """One panel per population: FACS proportion (x) vs scRNA-seq proportion (y)."""
    # FACS and scRNA-seq proportions side by side for each sample and population
    wide = proportions.pivot_table(index=[SAMPLE_COLUMN, "population"],
                                   columns="method", values="proportion").reset_index()
    fig, axes = plt.subplots(1, len(POPULATIONS), figsize=(4 * len(POPULATIONS), 4))
    for ax, pop in zip(np.atleast_1d(axes), POPULATIONS):
        part = wide[wide["population"] == pop]
        # One dot per sample
        ax.scatter(part["FACS"], part["scRNA-seq"], s=45, color=population_color(pop),
                   edgecolor="white", linewidth=1.5, zorder=3)
        # Dashed diagonal: where the dots would be if both methods agreed perfectly
        limit = np.nanmax(part[["FACS", "scRNA-seq"]].values) * 1.1
        ax.plot([0, limit], [0, limit], color="#8a8984", linewidth=1, linestyle="--", label="identity (x = y)")
        ax.set_xlim(0, limit)
        ax.set_ylim(0, limit)
        # Write the correlation from section A in the title
        row = agreement[agreement["population"] == pop].iloc[0]
        ax.set_title(f"{pop}\nSpearman r = {row['spearman_r']:.2f}, q = {row['spearman_q_fdr']:.3g}",
                     loc="left")
        ax.set_xlabel("FACS proportion")
        ax.set_ylabel("scRNA-seq proportion")
    np.atleast_1d(axes)[0].legend(loc="upper left")
    save(fig, "A2_proportions_FACS_vs_seq.png")


def plot_bland_altman(proportions, agreement):
    """
    Bland-Altman: x = mean of the two methods, y = difference (seq - FACS).
    Points around 0 = agreement; the dashed lines = 95% limits of agreement.
    """
    # FACS and scRNA-seq proportions side by side for each sample and population
    wide = proportions.pivot_table(index=[SAMPLE_COLUMN, "population"],
                                   columns="method", values="proportion").reset_index()
    fig, axes = plt.subplots(1, len(POPULATIONS), figsize=(4 * len(POPULATIONS), 3.8))
    for ax, pop in zip(np.atleast_1d(axes), POPULATIONS):
        part = wide[wide["population"] == pop]
        mean = (part["FACS"] + part["scRNA-seq"]) / 2     # x: average of the two methods
        diff = part["scRNA-seq"] - part["FACS"]           # y: how much they differ
        ax.scatter(mean, diff, s=45, color=population_color(pop), edgecolor="white", linewidth=1.5, zorder=3)
        row = agreement[agreement["population"] == pop].iloc[0]
        ax.axhline(0, color="#8a8984", linewidth=1)                               # zero line
        ax.axhline(row["bland_altman_bias"], color="#0b0b0b", linewidth=1.2)      # bias line
        for limit in ["bland_altman_lower", "bland_altman_upper"]:                # limits of agreement
            ax.axhline(row[limit], color="#52514e", linewidth=1, linestyle="--")
        ax.set_title(f"{pop}\nbias = {row['bland_altman_bias']:+.3f}", loc="left")
        ax.set_xlabel("Mean proportion (FACS + seq) / 2")
        ax.set_ylabel("scRNA-seq − FACS")
    save(fig, "A3_bland_altman.png")


def boxplot_with_points(ax, data_by_x, positions, color, width=0.35):
    """Draw thin box plots with every sample shown as a dot on top."""
    # Fixed random seed so the dot jitter looks the same every run
    rng = np.random.default_rng(0)
    for values, x in zip(data_by_x, positions):
        # Remove missing values; skip empty groups
        values = np.asarray(values, float)
        values = values[~np.isnan(values)]
        if len(values) == 0:
            continue
        # The box (median and quartiles), without separate outlier markers
        ax.boxplot(values, positions=[x], widths=width, showfliers=False,
                   boxprops={"color": color, "linewidth": 1.2},
                   whiskerprops={"color": color}, capprops={"color": color},
                   medianprops={"color": "#0b0b0b", "linewidth": 1.5})
        # The individual samples, shifted a little sideways so they don't overlap
        jitter = rng.uniform(-width / 4, width / 4, len(values))
        ax.scatter(x + jitter, values, s=22, color=color, edgecolor="white", linewidth=1, zorder=3)


def plot_proportions_by_condition(proportions, condition_column, conditions):
    """One panel per population: proportions per condition, FACS and scRNA-seq side by side."""
    fig, axes = plt.subplots(1, len(POPULATIONS),
                             figsize=(max(4, 1.6 * len(conditions)) * len(POPULATIONS) / 1.3, 4))
    for ax, pop in zip(np.atleast_1d(axes), POPULATIONS):
        part = proportions[proportions["population"] == pop]
        for k, method in enumerate(["FACS", "scRNA-seq"]):
            # One list of proportions per condition, for this method
            data = [part[(part["method"] == method) & (part[condition_column] == c)]["proportion"]
                    for c in conditions]
            # FACS box slightly left, scRNA-seq box slightly right of each condition
            positions = [i + (k - 0.5) * 0.4 for i in range(len(conditions))]
            boxplot_with_points(ax, data, positions, METHOD_COLORS[method])
            ax.scatter([], [], color=METHOD_COLORS[method], s=22, label=method)  # legend entry
        ax.set_xticks(range(len(conditions)), conditions)
        ax.set_title(pop, loc="left")
        ax.set_ylabel("Proportion of cells")
    np.atleast_1d(axes)[-1].legend(loc="upper left", bbox_to_anchor=(1, 1))
    save(fig, "B1_proportions_by_condition.png")


def plot_protein_vs_rna(protein_rna, rna_stats):
    """One panel per marker: FACS median intensity vs mean RNA, colored by population."""
    fig, axes = plt.subplots(1, len(MARKERS), figsize=(3.6 * len(MARKERS), 3.8))
    for ax, (marker, info) in zip(np.atleast_1d(axes), MARKERS.items()):
        part = protein_rna[protein_rna["marker"] == marker]
        # One dot per sample, colored by population
        for pop in POPULATIONS:
            sub = part[part["population"] == pop]
            ax.scatter(sub["rna_mean"], sub["protein_median"], s=40, color=population_color(pop),
                       edgecolor="white", linewidth=1.5, label=pop, zorder=3)
        # Show the pooled ("All populations") correlation from section C in the title
        row = rna_stats[(rna_stats["marker"] == marker) &
                        (rna_stats["population"] == "All populations")].iloc[0]
        ax.set_title(f"{marker} protein vs {info['gene']} RNA\nSpearman r = {row['spearman_r']:.2f}",
                     loc="left")
        ax.set_xlabel(f"{info['gene']} mean expression (log-norm)")
        ax.set_ylabel(f"{marker} median intensity (FACS)")
    np.atleast_1d(axes)[-1].legend(loc="upper left", bbox_to_anchor=(1, 1))
    save(fig, "C1_protein_vs_rna.png")


def plot_imaging(imaging, condition_column, conditions):
    """Grid: one row per imaging feature, one column per population, x = condition."""
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

ATLAS_SOURCES = ["Atlas RNA", "Flex RNA", "FACS protein"]   # the three panels of plot E3


def load_atlas(path):
    """
    Read the reference atlas and give every row a 'population', using
    POPULATIONS[...]["atlas_labels"] (or the seq_labels if none are given).
    Adds a "cells_in_row" column: 1 per row for a per-cell file, or the number from
    ATLAS_COUNT_COLUMN for a summary table. Labels not listed become 'Other'.
    """
    # Read the file; cell types and donors stay text
    atlas = pd.read_csv(path, dtype={ATLAS_CELL_TYPE_COLUMN: str, ATLAS_DONOR_COLUMN: str})
    get_column(atlas, ATLAS_CELL_TYPE_COLUMN, path.name)

    # Dictionary: atlas label -> population, e.g. "Myeloid cell" -> "Immune"
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


def atlas_comparison_stats(proportions, atlas_props, condition_column, conditions):
    """
    Section E: for every population x method x condition, compare our sample
    proportions with the atlas.
      - atlas with >= MIN_SAMPLES_PER_GROUP donors:
            Mann-Whitney U test, our samples vs the atlas donors
            (takes into account that healthy people also differ from each other)
      - otherwise (atlas without donors, or too few):
            one-sample Wilcoxon test: are our samples centred on the atlas proportion?
    FDR correction is applied over the whole table.
    """
    pooled = atlas_pooled_proportions(atlas_props)
    rows = []
    for pop in POPULATIONS:
        # The proportion of this population in every atlas donor
        donors = atlas_props.loc[atlas_props["population"] == pop, "proportion"].dropna().values
        for method in ["FACS", "scRNA-seq"]:
            for group in conditions:
                # Our samples of this method and condition
                values = proportions.loc[(proportions["population"] == pop) &
                                         (proportions["method"] == method) &
                                         (proportions[condition_column] == group), "proportion"].dropna().values
                median = np.median(values) if len(values) else np.nan
                row = {"population": pop, "method": method, "group": group,
                       "n_samples": len(values), "median_samples": median,
                       "atlas_proportion": pooled[pop],
                       "atlas_donor_median": np.median(donors) if len(donors) else np.nan,
                       "difference": median - pooled[pop],       # positive = more than in the atlas
                       "n_atlas_donors": len(donors),
                       "test": "", "statistic": np.nan, "p_value": np.nan}
                # Choose the test depending on how many atlas donors there are
                if len(donors) >= MIN_SAMPLES_PER_GROUP and len(values) >= MIN_SAMPLES_PER_GROUP:
                    row["test"] = "Mann-Whitney U vs atlas donors"
                    row["statistic"], row["p_value"] = stats.mannwhitneyu(values, donors, alternative="two-sided")
                elif (len(values) >= MIN_SAMPLES_PER_GROUP and not np.isnan(pooled[pop])
                      and np.any(values != pooled[pop])):
                    row["test"] = "Wilcoxon vs atlas proportion"
                    row["statistic"], row["p_value"] = stats.wilcoxon(values - pooled[pop])
                rows.append(row)
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


def build_marker_profiles(atlas, cells, protein_rna):
    """
    Section E: the typical level of each marker in each population, from three sources:
      Atlas RNA    : mean expression of the gene in the atlas cells
      Flex RNA     : mean expression of the gene in our sequenced cells
      FACS protein : median (over samples) of the antibody's median intensity
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
    # Typical FACS intensity per population and marker (median over samples)
    facs_medians = protein_rna.groupby(["population", "marker"])["protein_median"].median()

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
    E1: one panel per population. Our samples per condition (FACS and scRNA-seq
    boxes, as in B1), the atlas donors as an extra box, and a dashed line at the
    pooled atlas proportion.
    """
    pooled = atlas_pooled_proportions(atlas_props)
    has_donors = atlas_props["donor"].nunique() > 1
    groups = list(conditions) + (["Atlas"] if has_donors else [])
    fig, axes = plt.subplots(1, len(POPULATIONS),
                             figsize=(max(4, 1.6 * len(groups)) * len(POPULATIONS) / 1.3, 4))
    for ax, pop in zip(np.atleast_1d(axes), POPULATIONS):
        part = proportions[proportions["population"] == pop]
        # Our samples: FACS box slightly left, scRNA-seq box slightly right of each condition
        for k, method in enumerate(["FACS", "scRNA-seq"]):
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
    # Our samples, FACS and scRNA-seq side by side for each condition
    for k, method in enumerate(["FACS", "scRNA-seq"]):
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


def run_atlas_comparison(atlas_path, proportions, cells, protein_rna, condition_column, conditions):
    """
    Section E from start to end: load the atlas, build the tables, run the
    statistics, save the tables and draw the plots. Returns the test table.
    """
    atlas = load_atlas(atlas_path)
    atlas_props = build_atlas_proportions(atlas)
    tests = atlas_comparison_stats(proportions, atlas_props, condition_column, conditions)
    differences = composition_difference(proportions, atlas_props, condition_column)
    profiles = build_marker_profiles(atlas, cells, protein_rna)

    # Tables
    atlas_props.to_csv(OUTPUT_DIR / "E_atlas_proportions.csv", index=False)
    tests.to_csv(OUTPUT_DIR / "E_atlas_proportion_tests.csv", index=False)
    differences.to_csv(OUTPUT_DIR / "E_atlas_composition_difference.csv", index=False)
    if profiles is not None:
        profiles.to_csv(OUTPUT_DIR / "E_atlas_marker_profiles.csv", index=False)

    # Plots
    plot_atlas_proportions(proportions, atlas_props, condition_column, conditions)
    plot_atlas_difference(differences, conditions)
    if profiles is not None:
        plot_marker_profiles(profiles)
    return tests


# =============================================================================
# 4c. IMAGING PROFILES PER CELL TYPE (section F)
# =============================================================================
# Everything for section F is in this block. It describes each population with
# all the S8 imaging features at once: typical values, a "fingerprint", a PCA
# map, outliers, and links between shape and marker intensity. Sections A-E do
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


def imaging_marker_correlations(imaging, protein_rna):
    """
    Inside each population: across samples, does an imaging feature go together
    with a marker's FACS intensity? (e.g. do larger proximal tubule cells carry
    more CD13?) Spearman correlation for every population x feature x marker,
    FDR-corrected over the whole table.
    """
    protein = protein_rna[[SAMPLE_COLUMN, "population", "marker", "protein_median"]]
    rows = []
    for pop in POPULATIONS:
        for feature in IMAGING_FEATURES:
            image_values = imaging[(imaging["population"] == pop) &
                                   (imaging["feature"] == feature)][[SAMPLE_COLUMN, "value"]]
            for marker in MARKERS:
                marker_values = protein[(protein["population"] == pop) &
                                        (protein["marker"] == marker)][[SAMPLE_COLUMN, "protein_median"]]
                # Pair the two measurements of the same sample
                both = image_values.merge(marker_values, on=SAMPLE_COLUMN)
                rows.append({"population": pop, "feature": feature, "marker": marker,
                             **correlate(both["value"], both["protein_median"])})
    table = pd.DataFrame(rows)
    table["spearman_q_fdr"] = fdr_correct(table["spearman_p"])
    return table


def load_facs_cells(path, metadata, condition_column):
    """
    Read the optional per-cell imaging export. Each cell gets its sample_id
    (directly, or translated from facs_name) and the condition of its sample.
    Only cells of known populations and of the chosen samples are kept.
    """
    cells_img = pd.read_csv(path, dtype={SAMPLE_COLUMN: str, FACS_NAME_COLUMN: str})
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
    """
    rows = []
    for (sample, pop, group), part in cells_img.groupby([SAMPLE_COLUMN, "population", condition_column], sort=False):
        for feature in IMAGING_FEATURES:
            values = part[feature].dropna()
            mean = values.mean() if len(values) else np.nan
            sd = values.std(ddof=1) if len(values) > 1 else np.nan
            rows.append({
                SAMPLE_COLUMN: sample, "population": pop, condition_column: group, "feature": feature,
                "n_cells": len(values), "mean": mean, "median": values.median(), "sd": sd,
                "cv_percent": 100 * sd / mean if mean else np.nan,
                "p5": values.quantile(0.05), "p95": values.quantile(0.95),   # the middle 90% of cells
                "pct_outlier_cells": 100 * flags.loc[part.index, feature].mean(),
            })
        # Cells that are unusual in at least one feature
        rows.append({SAMPLE_COLUMN: sample, "population": pop, condition_column: group,
                     "feature": "Any feature", "n_cells": len(part),
                     "pct_outlier_cells": 100 * flags.loc[part.index, "Any feature"].mean()})
    return pd.DataFrame(rows)


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
    """F6: one panel per feature; violin plot of ALL single cells of each population."""
    pops = list(POPULATIONS)
    fig, axes = feature_grid(len(IMAGING_FEATURES))
    for ax, feature in zip(axes, IMAGING_FEATURES):
        for i, pop in enumerate(pops):
            values = cells_img.loc[cells_img["population"] == pop, feature].dropna().values
            if len(values) < 2:
                continue
            # The violin's width shows how many cells have each value
            part = ax.violinplot(values, positions=[i], widths=0.8, showmedians=True, showextrema=False)
            for body in part["bodies"]:
                body.set_facecolor(population_color(pop))
                body.set_edgecolor(population_color(pop))
                body.set_alpha(0.45)
            part["cmedians"].set_color("#0b0b0b")
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


def run_imaging_profiles(imaging, protein_rna, data_dir, metadata, condition_column, conditions):
    """
    Section F from start to end. Sample-level tables and plots F1-F5 always;
    cell-level table and plots F6-F8 only if the per-cell export is there.
    Returns the population summary table.
    """
    # --- Sample level (from the FlowJo table) ---
    outliers = imaging_outlier_samples(imaging, condition_column)
    summary = imaging_population_summary(imaging, outliers)
    scores, loadings = imaging_pca(imaging, condition_column)
    correlations = imaging_marker_correlations(imaging, protein_rna)

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
# (e.g. more immune cells in "Disease") so the statistics have something to find.

def make_demo_data(folder):
    """
    Write small, made-up input files in the exact format the program expects,
    so you can test the program before you have real data.
    10 samples: 5 'Healthy', 5 'Disease'.
    """
    rng = np.random.default_rng(1)          # fixed seed -> the same demo data every time
    folder.mkdir(exist_ok=True)
    pops = list(POPULATIONS)
    # 10 samples S01..S10, the first 5 Healthy and the last 5 Disease
    samples = [f"S{i:02d}" for i in range(1, 11)]
    conditions = ["Healthy"] * 5 + ["Disease"] * 5
    metadata = pd.DataFrame({
        SAMPLE_COLUMN: samples,
        FACS_NAME_COLUMN: [f"{s}_urine.fcs" for s in samples],
        "condition": conditions,
        "flex_probe_barcode": [f"BC{i:03d}" for i in range(1, 11)],
    })

    # Typical share of Immune / Proximal tubule / Epithelial / Other CD45- cells per condition
    typical_fraction = {"Healthy": [0.25, 0.18, 0.50, 0.07], "Disease": [0.45, 0.13, 0.35, 0.07]}
    protein_level = {  # typical median intensity per population and marker
        "Immune":          {"CD45": 9000, "CD10": 600,  "CD13": 1500, "EpCAM": 300},
        "Proximal tubule": {"CD45": 250,  "CD10": 6000, "CD13": 8000, "EpCAM": 2000},
        "Epithelial":      {"CD45": 250,  "CD10": 800,  "CD13": 900,  "EpCAM": 7000},
        "Other CD45-":     {"CD45": 250,  "CD10": 500,  "CD13": 700,  "EpCAM": 300},
    }
    # Typical value of each imaging feature per population
    imaging_level = {"Immune": [90, 0.30, 4.0], "Proximal tubule": [210, 0.50, 7.0],
                     "Epithelial": [320, 0.60, 9.0], "Other CD45-": [150, 0.45, 5.5]}

    facs_rows, cell_tables = [], []
    for sample, condition in zip(samples, conditions):
        # --- FACS row for this sample ---
        # Random population fractions around the typical values of its condition
        true_fraction = rng.dirichlet(np.array(typical_fraction[condition]) * 40)
        live = int(rng.integers(20000, 60000))       # number of live single cells
        row = {"Name": f"{sample}_urine.fcs", COUNT_COLUMN.format(gate=FACS_TOTAL_GATE): live}
        marker_factor = {m: rng.lognormal(0, 0.25) for m in MARKERS}  # sample-to-sample variation
        for pop, fraction, img in zip(pops, true_fraction, imaging_level.values()):
            gate = POPULATIONS[pop]["gate"]
            # Gate count (85% of live cells fall in one of the four gates)
            row[COUNT_COLUMN.format(gate=gate)] = int(live * fraction * 0.85)
            # Marker medians with some random noise
            for marker, info in MARKERS.items():
                value = protein_level[pop][marker] * marker_factor[marker] * rng.lognormal(0, 0.1)
                row[MEDIAN_COLUMN.format(gate=gate, parameter=info["facs_parameter"])] = round(value, 1)
            # Imaging medians; proximal tubule cells are 15% larger in Disease
            shift = 1.15 if condition == "Disease" and pop == "Proximal tubule" else 1.0
            for feature, base in zip(IMAGING_FEATURES, img):
                row[MEDIAN_COLUMN.format(gate=gate, parameter=feature)] = round(base * shift * rng.lognormal(0, 0.08), 3)
        facs_rows.append(row)

        # scRNA-seq: similar fractions with method-specific noise, plus a few 'Other' cells
        seq_fraction = rng.dirichlet(np.append(true_fraction, 0.05) * 60)
        n_cells = int(rng.integers(800, 2500))
        # Give every cell a population, then a matching notebook label
        cell_pop = rng.choice(pops + ["Other"], size=n_cells, p=seq_fraction)
        labels = [rng.choice(POPULATIONS[p]["seq_labels"]) if p != "Other"
                  else rng.choice(["Unknown", "Doublet"]) for p in cell_pop]
        cells = pd.DataFrame({"cell_barcode": [f"{sample}_{i:05d}" for i in range(n_cells)],
                              SAMPLE_COLUMN: sample, CELL_TYPE_COLUMN: labels})
        # Gene expression that follows the protein level of the cell's population
        for marker, info in MARKERS.items():
            mean = np.array([np.log1p(protein_level[p][marker] / 1500 * marker_factor[marker])
                             if p != "Other" else 0.3 for p in cell_pop])
            cells[info["gene"]] = np.clip(rng.normal(mean, 0.5), 0, None).round(3)
        cell_tables.append(cells)

    # --- Reference atlas: 8 healthy donors with known cell types (per-cell format) ---
    # Made after the sample data, so the three files above are the same as before.
    atlas_tables = []
    for d in range(1, 9):
        # Healthy-like composition with donor-to-donor variation, plus a few unassigned cells
        fraction = rng.dirichlet(np.append(typical_fraction["Healthy"], 0.03) * 50)
        n_atlas = int(rng.integers(1500, 4000))
        atlas_pop = rng.choice(pops + ["Other"], size=n_atlas, p=fraction)
        # Atlas names for the cell types (they differ from the notebook labels)
        labels = [rng.choice(POPULATIONS[p].get("atlas_labels", POPULATIONS[p]["seq_labels"]))
                  if p != "Other" else "Unknown" for p in atlas_pop]
        atlas = pd.DataFrame({ATLAS_DONOR_COLUMN: f"D{d:02d}", ATLAS_CELL_TYPE_COLUMN: labels})
        # Marker genes with the same population pattern as in our Flex data
        for marker, info in MARKERS.items():
            mean = np.array([np.log1p(protein_level[p][marker] / 1500) if p != "Other" else 0.3
                             for p in atlas_pop])
            atlas[info["gene"]] = np.clip(rng.normal(mean, 0.5), 0, None).round(3)
        atlas_tables.append(atlas)

    # --- Per-cell imaging export: the single cells of every gate (section F) ---
    # Uses its own random numbers, so all the other demo files stay exactly the same.
    cell_rng = np.random.default_rng(2)
    imaging_cells = []
    for row in facs_rows:
        for pop in pops:
            gate = POPULATIONS[pop]["gate"]
            n_cells_gate = int(cell_rng.integers(150, 300))
            part = {FACS_NAME_COLUMN: row["Name"], FACS_CELLS_POPULATION_COLUMN: pop}
            # About 2% odd cells (e.g. doublets or clumps) that sit far from the rest
            odd = cell_rng.random(n_cells_gate) < 0.02
            for feature in IMAGING_FEATURES:
                # Cells scatter around the sample's median from the FlowJo table
                median = row[MEDIAN_COLUMN.format(gate=gate, parameter=feature)]
                values = median * cell_rng.lognormal(0, 0.12, n_cells_gate)
                values[odd] *= 1.9
                if "Eccentricity" in feature:          # eccentricity is between 0 and 1
                    values = np.clip(values, 0, 0.99)
                part[feature] = values.round(3)
            imaging_cells.append(pd.DataFrame(part))

    # Write the files in the same formats as the real inputs
    metadata.to_csv(folder / METADATA_FILE, index=False)
    pd.concat(imaging_cells).to_csv(folder / FACS_CELLS_FILE, index=False)
    pd.DataFrame(facs_rows).to_excel(folder / FACS_FILE, index=False)
    pd.concat(cell_tables).to_csv(folder / SEQ_FILE, index=False)
    pd.concat(atlas_tables).to_csv(folder / ATLAS_FILE, index=False)
    print(f"Demo data written to {folder}")


# =============================================================================
# 6. MAIN PROGRAM - runs the steps above in order
# =============================================================================
# run_analysis() does the actual work and is shared by the terminal (main) and
# the web app (app.py). main() only handles the command-line options.

def run_analysis(data_dir, output_dir, condition_column=None, conditions=None):
    """
    Run steps 1-5 on the input files in data_dir, plus section E (if an atlas
    file is there) and section F (if imaging features are set), and write every
    table and plot to output_dir. Used by the command line (main, below) and by the web
    interface (app.py). If condition_column or conditions is None, the user is
    asked in the terminal. Returns the summary tables so the caller can show them.
    """
    # OUTPUT_DIR is changed for the whole script so save() writes to the right folder
    global OUTPUT_DIR
    data_dir = Path(data_dir)
    OUTPUT_DIR = Path(output_dir)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)   # create the results folder if needed

    # --- Step 1: load files and choose conditions ---
    metadata = load_metadata(data_dir / METADATA_FILE)
    condition_column, conditions = choose_conditions(metadata, condition_column, conditions)
    # Keep only the samples that belong to the chosen conditions
    metadata = metadata[metadata[condition_column].isin(conditions)]
    facs = load_facs(data_dir / FACS_FILE, metadata)
    cells = load_seq(data_dir / SEQ_FILE, metadata)
    print(f"Loaded {len(facs)} FACS samples and {len(cells)} cells "
          f"from {cells[SAMPLE_COLUMN].nunique()} scRNA-seq samples.")

    # --- Step 2: build tables ---
    proportions = build_proportions(facs, cells, metadata, condition_column)
    protein_rna = build_protein_vs_rna(facs, cells, metadata, condition_column)
    # Imaging table only if imaging features are listed in SETTINGS
    imaging = build_imaging(facs, metadata, condition_column) if IMAGING_FEATURES else None

    # --- Step 3: statistics ---
    agreement = method_agreement_stats(proportions)                        # section A
    by_condition = condition_stats(proportions, "proportion", ["population", "method"],
                                   condition_column, conditions)          # section B
    rna_stats = protein_rna_stats(protein_rna)                             # section C

    # --- Step 4: save tables ---
    print("\nSaving results to", OUTPUT_DIR)
    proportions.to_csv(OUTPUT_DIR / "A_proportions.csv", index=False)
    agreement.to_csv(OUTPUT_DIR / "A_method_agreement_stats.csv", index=False)
    by_condition.to_csv(OUTPUT_DIR / "B_condition_tests_proportions.csv", index=False)
    protein_rna.to_csv(OUTPUT_DIR / "C_protein_vs_rna.csv", index=False)
    rna_stats.to_csv(OUTPUT_DIR / "C_protein_vs_rna_stats.csv", index=False)

    # --- Step 5: plots ---
    plot_stacked_proportions(proportions, condition_column, conditions)
    plot_method_scatter(proportions, agreement)
    plot_bland_altman(proportions, agreement)
    # Comparing conditions only makes sense with at least two of them
    if len(conditions) > 1:
        plot_proportions_by_condition(proportions, condition_column, conditions)
    plot_protein_vs_rna(protein_rna, rna_stats)

    # --- Section D: imaging statistics, tables and plot (skipped if no features) ---
    imaging_pop = imaging_cond = None
    if imaging is not None:
        imaging_pop = imaging_between_populations(imaging)
        imaging_cond = condition_stats(imaging, "value", ["feature", "population"],
                                       condition_column, conditions)
        imaging.to_csv(OUTPUT_DIR / "D_imaging_features.csv", index=False)
        imaging_pop.to_csv(OUTPUT_DIR / "D_imaging_tests_between_populations.csv", index=False)
        imaging_cond.to_csv(OUTPUT_DIR / "D_imaging_tests_between_conditions.csv", index=False)
        plot_imaging(imaging, condition_column, conditions)

    # --- Section E: comparison with the reference atlas (skipped if no atlas file) ---
    atlas_tests = None
    atlas_path = data_dir / ATLAS_FILE
    if atlas_path.exists():
        print("\nComparing with the reference atlas ...")
        atlas_tests = run_atlas_comparison(atlas_path, proportions, cells, protein_rna,
                                           condition_column, conditions)
    else:
        print(f"\nNo reference atlas ({ATLAS_FILE}) found - section E skipped.")

    # --- Section F: imaging profiles per cell type (skipped if no imaging features) ---
    imaging_summary = None
    if imaging is not None:
        print("\nImaging profiles per cell type ...")
        imaging_summary = run_imaging_profiles(imaging, protein_rna, data_dir, metadata,
                                               condition_column, conditions)

    # Hand the main results back to whoever called this (terminal or web app)
    return {
        "imaging_summary": imaging_summary,
        "atlas_tests": atlas_tests,
        "condition_column": condition_column,
        "conditions": conditions,
        "n_facs_samples": len(facs),
        "n_cells": len(cells),
        "n_seq_samples": cells[SAMPLE_COLUMN].nunique(),
        "agreement": agreement,
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
    parser.add_argument("--conditions", nargs="+", help="conditions to include (default: ask)")
    args = parser.parse_args()

    # Use the real data folders, or create demo data and use demo folders instead
    data_dir, output_dir = DATA_DIR, OUTPUT_DIR
    if args.demo:
        data_dir = SCRIPT_DIR / "demo_data"
        output_dir = SCRIPT_DIR / "results_demo"
        make_demo_data(data_dir)

    # Run the whole analysis (steps 1-5)
    results = run_analysis(data_dir, output_dir, args.condition_column, args.conditions)
    agreement, by_condition = results["agreement"], results["by_condition"]

    # --- Step 6: short summary in the terminal ---
    # How well FACS and scRNA-seq agree, per population
    print("\nFACS vs scRNA-seq agreement (proportions):")
    print(agreement[["population", "n", "spearman_r", "spearman_q_fdr",
                     "bland_altman_bias", "wilcoxon_q_fdr"]].round(3).to_string(index=False))
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
