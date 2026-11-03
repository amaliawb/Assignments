"""
facs_vs_scrnaseq.py
===================

Compare urine cell populations measured two ways on the SAME samples:

  1. Flow cytometry / imaging sorting on the BD FACSDiscover S8
     (gated in FlowJo, exported as an Excel table)
  2. 10x Genomics Flex single-cell RNA-seq
     (processed and annotated in a Jupyter notebook, exported as a CSV)

Populations (FACS gate logic):
  - Immune cells           : CD45+
  - Epithelial cells       : CD45- EpCAM+
  - Proximal tubule cells  : CD45- CD10+ CD13+

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
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # draw plots straight to files (no window pops up)
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


# =============================================================================
# SETTINGS - everything you may need to change is in this section
# =============================================================================

# ---- File locations (relative to the folder this script is in) -------------
SCRIPT_DIR = Path(__file__).resolve().parent
DATA_DIR = SCRIPT_DIR / "data"                      # put your real files here
METADATA_FILE = "metadata.csv"
FACS_FILE = "facs_flowjo_export.xlsx"
FACS_SHEET = 0                                      # sheet number or name in the Excel file
SEQ_FILE = "scrnaseq_cells.csv"
OUTPUT_DIR = SCRIPT_DIR / "results"                 # all tables and plots go here

# Column names in the metadata and scRNA-seq files
SAMPLE_COLUMN = "sample_id"
FACS_NAME_COLUMN = "facs_name"
CELL_TYPE_COLUMN = "cell_type"

# ---- FlowJo column-name templates --------------------------------------------
# FlowJo names a statistic column "<gate path> | <statistic>".
# {gate} is replaced by a gate path below, {parameter} by a channel/feature name.
COUNT_COLUMN = "{gate} | Count"
MEDIAN_COLUMN = "{gate} | Median ({parameter})"

# Gate that holds ALL good cells (used only if PROPORTION_DENOMINATOR = "total")
FACS_TOTAL_GATE = "Cells/Single Cells/Live"

# ---- The three populations ----------------------------------------------------
# "gate"       : the FlowJo gate path of the population
# "seq_labels" : which cell_type labels from your notebook belong to this
#                population (list ALL the labels you used, e.g. every immune subtype)
POPULATIONS = {
    "Immune": {
        "gate": "Cells/Single Cells/Live/CD45+",
        "seq_labels": ["T cell", "Macrophage", "Neutrophil"],
    },
    "Epithelial": {
        "gate": "Cells/Single Cells/Live/CD45-/EpCAM+",
        "seq_labels": ["Urothelial cell", "Squamous epithelial cell"],
    },
    "Proximal tubule": {
        "gate": "Cells/Single Cells/Live/CD45-/CD10+ CD13+",
        "seq_labels": ["Proximal tubule cell"],
    },
}

# ---- Antibody markers and their matching genes ------------------------------
# "facs_parameter" is the channel name FlowJo uses. Replace FLUOR1..4 with your
# fluorochromes, e.g. "Comp-BV510-A :: CD45" - copy it exactly from FlowJo.
MARKERS = {
    "CD45":  {"gene": "PTPRC", "facs_parameter": "Comp-FLUOR1-A :: CD45"},
    "CD10":  {"gene": "MME",   "facs_parameter": "Comp-FLUOR2-A :: CD10"},
    "CD13":  {"gene": "ANPEP", "facs_parameter": "Comp-FLUOR3-A :: CD13"},
    "EpCAM": {"gene": "EPCAM", "facs_parameter": "Comp-FLUOR4-A :: EpCAM"},
}

# ---- S8 imaging features -------------------------------------------------------
# Names of the imaging parameters exactly as they appear in FlowJo.
# These are placeholders - replace with the features you export (add or remove freely).
IMAGING_FEATURES = [
    "Area (LightLoss)",
    "Eccentricity (LightLoss)",
    "Radial Moment (LightLoss)",
]

# ---- Analysis options -----------------------------------------------------------
# How proportions are calculated (same rule for FACS and scRNA-seq):
#   "populations" : population / (sum of the three populations)   <- recommended
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
POPULATION_COLORS = ["#2a78d6", "#eb6834", "#1baf7a"]
METHOD_COLORS = {"FACS": "#52514e", "scRNA-seq": "#b5b3ad"}
PLOT_DPI = 200


# =============================================================================
# 1. LOADING THE DATA
# =============================================================================

def get_column(table, column, file_label):
    """Return one column of a table; if it does not exist, stop and list the real column names."""
    if column not in table.columns:
        available = "\n    ".join(map(str, table.columns))
        raise SystemExit(
            f"\nColumn not found in {file_label}:\n    '{column}'\n"
            f"Columns that ARE in the file:\n    {available}\n"
            "-> Fix the matching name in the SETTINGS section of this script."
        )
    return table[column]


def load_metadata(path):
    """Read the sample sheet. Everything is read as text so IDs like '01' stay intact."""
    metadata = pd.read_csv(path, dtype=str)
    get_column(metadata, SAMPLE_COLUMN, path.name)
    get_column(metadata, FACS_NAME_COLUMN, path.name)
    return metadata


def ask_for_numbers(prompt, n_options, allow_many):
    """Ask the user to type option number(s); ask again until the answer is valid.
    Returns a list of 0-based indexes, or None if the user pressed Enter (= all)."""
    while True:
        answer = input(prompt).strip()
        if answer == "" and allow_many:
            return None
        try:
            picks = [int(x) - 1 for x in answer.split(",")]
        except ValueError:
            picks = []
        if picks and all(0 <= p < n_options for p in picks) and (allow_many or len(picks) == 1):
            return picks
        print(f"  Please type {'numbers' if allow_many else 'one number'} between 1 and {n_options}.")


def choose_conditions(metadata, column=None, values=None):
    """
    Decide which metadata column is the 'condition' and which of its values to analyse.
    Asks the user in the terminal unless both were given on the command line.
    """
    candidates = [c for c in metadata.columns if c not in (SAMPLE_COLUMN, FACS_NAME_COLUMN)]
    if not candidates:
        raise SystemExit("The metadata file has no column describing conditions.")

    # --- which column? ---
    if column is None:
        if len(candidates) == 1:
            column = candidates[0]
        else:
            print("\nWhich metadata column holds the conditions?")
            for i, name in enumerate(candidates, start=1):
                print(f"  {i}. {name}")
            pick = ask_for_numbers("Type a number: ", len(candidates), allow_many=False)
            column = candidates[pick[0]]
    get_column(metadata, column, "the metadata file")

    # --- which values? ---
    counts = metadata[column].value_counts(sort=False)
    available = list(counts.index)
    if values is None:
        print(f"\nConditions found in '{column}':")
        for i, name in enumerate(available, start=1):
            print(f"  {i}. {name}  ({counts[name]} samples)")
        picks = ask_for_numbers("Type the numbers to include separated by commas (Enter = all): ",
                                len(available), allow_many=True)
        values = available if picks is None else [available[p] for p in picks]
    else:
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
    facs = pd.read_excel(path, sheet_name=FACS_SHEET)
    facs = facs.rename(columns={facs.columns[0]: FACS_NAME_COLUMN})
    facs[FACS_NAME_COLUMN] = facs[FACS_NAME_COLUMN].astype(str)

    name_to_id = dict(zip(metadata[FACS_NAME_COLUMN], metadata[SAMPLE_COLUMN]))
    missing = sorted(set(name_to_id) - set(facs[FACS_NAME_COLUMN]))
    if missing:
        print(f"WARNING: these samples are in the metadata but not in the FACS file: {missing}")

    facs = facs[facs[FACS_NAME_COLUMN].isin(name_to_id)].copy()
    facs[SAMPLE_COLUMN] = facs[FACS_NAME_COLUMN].map(name_to_id)
    return facs.set_index(SAMPLE_COLUMN)


def load_seq(path, metadata):
    """
    Read the per-cell scRNA-seq table and add a 'population' column by
    translating each cell_type label using POPULATIONS[...]["seq_labels"].
    Labels not listed anywhere become 'Other'.
    """
    cells = pd.read_csv(path, dtype={SAMPLE_COLUMN: str})
    get_column(cells, SAMPLE_COLUMN, path.name)
    get_column(cells, CELL_TYPE_COLUMN, path.name)
    for marker in MARKERS.values():
        get_column(cells, marker["gene"], path.name)

    label_to_population = {
        label: population
        for population, info in POPULATIONS.items()
        for label in info["seq_labels"]
    }
    cells["population"] = cells[CELL_TYPE_COLUMN].map(label_to_population).fillna("Other")

    unmapped = sorted(cells.loc[cells["population"] == "Other", CELL_TYPE_COLUMN].unique())
    if unmapped:
        print(f"Note: these cell_type labels are not in any population (counted as 'Other'): {unmapped}")

    missing = sorted(set(metadata[SAMPLE_COLUMN]) - set(cells[SAMPLE_COLUMN]))
    if missing:
        print(f"WARNING: these samples are in the metadata but not in the scRNA-seq file: {missing}")
    return cells[cells[SAMPLE_COLUMN].isin(metadata[SAMPLE_COLUMN])]


# =============================================================================
# 2. BUILDING THE COMPARISON TABLES
# =============================================================================

def build_proportions(facs, cells, metadata, condition_column):
    """
    Long table: one row per sample x population x method, with the cell count
    and the proportion. This is the table used for sections A and B.
    """
    populations = list(POPULATIONS)

    # FACS counts per population
    facs_counts = pd.DataFrame({
        pop: get_column(facs, COUNT_COLUMN.format(gate=info["gate"]), "the FACS file")
        for pop, info in POPULATIONS.items()
    })
    # scRNA-seq counts per population (cells per sample per population)
    seq_counts = pd.crosstab(cells[SAMPLE_COLUMN], cells["population"])
    seq_counts = seq_counts.reindex(columns=populations + ["Other"], fill_value=0)

    if PROPORTION_DENOMINATOR == "populations":
        facs_total = facs_counts[populations].sum(axis=1)
        seq_total = seq_counts[populations].sum(axis=1)
    elif PROPORTION_DENOMINATOR == "total":
        facs_total = get_column(facs, COUNT_COLUMN.format(gate=FACS_TOTAL_GATE), "the FACS file")
        seq_total = seq_counts.sum(axis=1)
    else:
        raise SystemExit("PROPORTION_DENOMINATOR must be 'populations' or 'total'.")

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
                    "proportion": counts.loc[sample, pop] / total[sample] if total[sample] > 0 else np.nan,
                })
    table = pd.DataFrame(rows)
    return table.merge(metadata[[SAMPLE_COLUMN, condition_column]], on=SAMPLE_COLUMN)


def build_protein_vs_rna(facs, cells, metadata, condition_column):
    """
    One row per sample x population x marker:
      protein = FACS median intensity of the marker inside the population gate
      rna     = mean log-normalized expression of the gene in that population's cells
    """
    rna_means = cells.groupby([SAMPLE_COLUMN, "population"])[
        [m["gene"] for m in MARKERS.values()]].mean()
    cell_numbers = cells.groupby([SAMPLE_COLUMN, "population"]).size()

    rows = []
    for pop, pop_info in POPULATIONS.items():
        for marker, marker_info in MARKERS.items():
            column = MEDIAN_COLUMN.format(gate=pop_info["gate"], parameter=marker_info["facs_parameter"])
            protein = get_column(facs, column, "the FACS file")
            for sample in facs.index:
                key = (sample, pop)
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
    return table.merge(metadata[[SAMPLE_COLUMN, condition_column]], on=SAMPLE_COLUMN)


def build_imaging(facs, metadata, condition_column):
    """One row per sample x population x imaging feature (FACS median value)."""
    rows = []
    for pop, info in POPULATIONS.items():
        for feature in IMAGING_FEATURES:
            column = MEDIAN_COLUMN.format(gate=info["gate"], parameter=feature)
            values = get_column(facs, column, "the FACS file")
            for sample in facs.index:
                rows.append({SAMPLE_COLUMN: sample, "population": pop,
                             "feature": feature, "value": values[sample]})
    table = pd.DataFrame(rows)
    return table.merge(metadata[[SAMPLE_COLUMN, condition_column]], on=SAMPLE_COLUMN)


# =============================================================================
# 3. STATISTICS
# =============================================================================
# Non-parametric tests are used throughout: urine studies usually have few
# samples, and proportions/intensities are rarely normally distributed.

def fdr_correct(p_values):
    """Benjamini-Hochberg FDR correction. NaN p-values (skipped tests) stay NaN."""
    p = np.asarray(p_values, dtype=float)
    q = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    n = ok.sum()
    if n == 0:
        return q
    order = np.argsort(p[ok])
    ranked = p[ok][order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]   # keep q-values monotonic
    corrected = np.empty(n)
    corrected[order] = np.minimum(ranked, 1)
    q[ok] = corrected
    return q


def correlate(x, y):
    """Spearman (rank) and Pearson (linear) correlation between two paired lists."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    keep = ~(np.isnan(x) | np.isnan(y))
    x, y = x[keep], y[keep]
    result = {"n": len(x), "spearman_r": np.nan, "spearman_p": np.nan,
              "pearson_r": np.nan, "pearson_p": np.nan}
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
    groups = {name: np.asarray(v, float)[~np.isnan(np.asarray(v, float))]
              for name, v in values_by_group.items()}
    groups = {name: v for name, v in groups.items() if len(v) >= MIN_SAMPLES_PER_GROUP}
    names = list(groups)
    results = []
    if len(names) < 2:
        return results

    if len(names) > 2:
        statistic, p = stats.kruskal(*groups.values())
        results.append({"test": "Kruskal-Wallis", "comparison": " vs ".join(names),
                        "statistic": statistic, "p_value": p})
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
    wide = proportions.pivot_table(index=[SAMPLE_COLUMN, "population"],
                                   columns="method", values="proportion").reset_index()
    rows = []
    for pop in POPULATIONS:
        part = wide[wide["population"] == pop].dropna(subset=["FACS", "scRNA-seq"])
        facs_values, seq_values = part["FACS"].values, part["scRNA-seq"].values
        difference = seq_values - facs_values
        row = {"population": pop, **correlate(facs_values, seq_values),
               "mean_FACS": facs_values.mean(), "mean_scRNA-seq": seq_values.mean(),
               "bland_altman_bias": difference.mean(),
               "bland_altman_lower": difference.mean() - 1.96 * difference.std(ddof=1),
               "bland_altman_upper": difference.mean() + 1.96 * difference.std(ddof=1),
               "wilcoxon_stat": np.nan, "wilcoxon_p": np.nan}
        if len(part) >= MIN_SAMPLES_PER_GROUP and np.any(difference != 0):
            row["wilcoxon_stat"], row["wilcoxon_p"] = stats.wilcoxon(facs_values, seq_values)
        rows.append(row)

    table = pd.DataFrame(rows)
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
    for keys, part in table.groupby(split_columns, sort=False):
        keys = keys if isinstance(keys, tuple) else (keys,)
        values = {c: part.loc[part[condition_column] == c, value_column].values for c in conditions}
        for result in compare_groups(values):
            rows.append({**dict(zip(split_columns, keys)), **result})
    result = pd.DataFrame(rows)
    if not result.empty:
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
        rows.append({"marker": marker, "population": "All populations",
                     **correlate(part["protein_median"], part["rna_mean"])})
        for pop, sub in part.groupby("population", sort=False):
            rows.append({"marker": marker, "population": pop,
                         **correlate(sub["protein_median"], sub["rna_mean"])})
    table = pd.DataFrame(rows)
    table["spearman_q_fdr"] = fdr_correct(table["spearman_p"])
    table["pearson_q_fdr"] = fdr_correct(table["pearson_p"])
    return table


def imaging_between_populations(imaging):
    """Section D: do the imaging features differ between the three populations?"""
    rows = []
    for feature, part in imaging.groupby("feature", sort=False):
        values = {pop: part.loc[part["population"] == pop, "value"].values for pop in POPULATIONS}
        for result in compare_groups(values):
            rows.append({"feature": feature, **result})
    table = pd.DataFrame(rows)
    if not table.empty:
        table["q_value_fdr"] = fdr_correct(table["p_value"])
        table["significant"] = table["q_value_fdr"] < ALPHA
    return table


# =============================================================================
# 4. PLOTS
# =============================================================================

plt.rcParams.update({
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.edgecolor": "#8a8984", "axes.labelcolor": "#0b0b0b",
    "xtick.color": "#52514e", "ytick.color": "#52514e",
    "axes.grid": True, "grid.color": "#e6e5e0", "grid.linewidth": 0.6,
    "axes.axisbelow": True, "font.size": 9, "legend.frameon": False,
})


def population_color(pop):
    return POPULATION_COLORS[list(POPULATIONS).index(pop) % len(POPULATION_COLORS)]


def save(fig, name):
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / name, dpi=PLOT_DPI)
    plt.close(fig)
    print(f"  saved {name}")


def plot_stacked_proportions(proportions, condition_column, conditions):
    """Stacked bars: each sample gets a FACS bar and a scRNA-seq bar side by side."""
    order = (proportions[[SAMPLE_COLUMN, condition_column]].drop_duplicates()
             .assign(_c=lambda d: d[condition_column].map(conditions.index))
             .sort_values(["_c", SAMPLE_COLUMN]))
    samples = list(order[SAMPLE_COLUMN])
    wide = proportions.pivot_table(index=[SAMPLE_COLUMN, "method"], columns="population",
                                   values="proportion")

    fig, ax = plt.subplots(figsize=(max(6, len(samples) * 0.9), 4.5))
    positions, labels = [], []
    for i, sample in enumerate(samples):
        for j, method in enumerate(["FACS", "scRNA-seq"]):
            x = i * 3 + j
            bottom = 0
            for pop in POPULATIONS:
                value = wide.loc[(sample, method), pop] if (sample, method) in wide.index else 0
                ax.bar(x, value, bottom=bottom, width=0.85, color=population_color(pop),
                       edgecolor="white", linewidth=1, label=pop if (i == 0 and j == 0) else None)
                bottom += value
            positions.append(x)
            labels.append(f"{sample} · {'FACS' if method == 'FACS' else 'seq'}")
    ax.set_xticks(positions, labels, rotation=90)
    ax.set_ylabel("Proportion of cells")
    ax.set_title("Cell-type proportions per sample (FACS vs scRNA-seq)", loc="left")
    ax.grid(axis="x", visible=False)
    ax.legend(loc="upper left", bbox_to_anchor=(1, 1))
    save(fig, "A1_proportions_stacked.png")


def plot_method_scatter(proportions, agreement):
    """One panel per population: FACS proportion (x) vs scRNA-seq proportion (y)."""
    wide = proportions.pivot_table(index=[SAMPLE_COLUMN, "population"],
                                   columns="method", values="proportion").reset_index()
    fig, axes = plt.subplots(1, len(POPULATIONS), figsize=(4 * len(POPULATIONS), 4))
    for ax, pop in zip(np.atleast_1d(axes), POPULATIONS):
        part = wide[wide["population"] == pop]
        ax.scatter(part["FACS"], part["scRNA-seq"], s=45, color=population_color(pop),
                   edgecolor="white", linewidth=1.5, zorder=3)
        limit = np.nanmax(part[["FACS", "scRNA-seq"]].values) * 1.1
        ax.plot([0, limit], [0, limit], color="#8a8984", linewidth=1, linestyle="--", label="identity (x = y)")
        ax.set_xlim(0, limit)
        ax.set_ylim(0, limit)
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
    wide = proportions.pivot_table(index=[SAMPLE_COLUMN, "population"],
                                   columns="method", values="proportion").reset_index()
    fig, axes = plt.subplots(1, len(POPULATIONS), figsize=(4 * len(POPULATIONS), 3.8))
    for ax, pop in zip(np.atleast_1d(axes), POPULATIONS):
        part = wide[wide["population"] == pop]
        mean = (part["FACS"] + part["scRNA-seq"]) / 2
        diff = part["scRNA-seq"] - part["FACS"]
        ax.scatter(mean, diff, s=45, color=population_color(pop), edgecolor="white", linewidth=1.5, zorder=3)
        row = agreement[agreement["population"] == pop].iloc[0]
        ax.axhline(0, color="#8a8984", linewidth=1)
        ax.axhline(row["bland_altman_bias"], color="#0b0b0b", linewidth=1.2)
        for limit in ["bland_altman_lower", "bland_altman_upper"]:
            ax.axhline(row[limit], color="#52514e", linewidth=1, linestyle="--")
        ax.set_title(f"{pop}\nbias = {row['bland_altman_bias']:+.3f}", loc="left")
        ax.set_xlabel("Mean proportion (FACS + seq) / 2")
        ax.set_ylabel("scRNA-seq − FACS")
    save(fig, "A3_bland_altman.png")


def boxplot_with_points(ax, data_by_x, positions, color, width=0.35):
    """Draw thin box plots with every sample shown as a dot on top."""
    rng = np.random.default_rng(0)
    for values, x in zip(data_by_x, positions):
        values = np.asarray(values, float)
        values = values[~np.isnan(values)]
        if len(values) == 0:
            continue
        ax.boxplot(values, positions=[x], widths=width, showfliers=False,
                   boxprops={"color": color, "linewidth": 1.2},
                   whiskerprops={"color": color}, capprops={"color": color},
                   medianprops={"color": "#0b0b0b", "linewidth": 1.5})
        jitter = rng.uniform(-width / 4, width / 4, len(values))
        ax.scatter(x + jitter, values, s=22, color=color, edgecolor="white", linewidth=1, zorder=3)


def plot_proportions_by_condition(proportions, condition_column, conditions):
    """One panel per population: proportions per condition, FACS and scRNA-seq side by side."""
    fig, axes = plt.subplots(1, len(POPULATIONS),
                             figsize=(max(4, 1.6 * len(conditions)) * len(POPULATIONS) / 1.3, 4))
    for ax, pop in zip(np.atleast_1d(axes), POPULATIONS):
        part = proportions[proportions["population"] == pop]
        for k, method in enumerate(["FACS", "scRNA-seq"]):
            data = [part[(part["method"] == method) & (part[condition_column] == c)]["proportion"]
                    for c in conditions]
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
        for pop in POPULATIONS:
            sub = part[part["population"] == pop]
            ax.scatter(sub["rna_mean"], sub["protein_median"], s=40, color=population_color(pop),
                       edgecolor="white", linewidth=1.5, label=pop, zorder=3)
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
            part = imaging[(imaging["feature"] == feature) & (imaging["population"] == pop)]
            data = [part[part[condition_column] == cond]["value"] for cond in conditions]
            boxplot_with_points(ax, data, range(len(conditions)), population_color(pop), width=0.5)
            ax.set_xticks(range(len(conditions)), conditions)
            if r == 0:
                ax.set_title(pop, loc="left")
            if c == 0:
                ax.set_ylabel(feature)
    save(fig, "D1_imaging_features.png")


# =============================================================================
# 5. DEMO DATA (only used with --demo)
# =============================================================================

def make_demo_data(folder):
    """
    Write small, made-up input files in the exact format the program expects,
    so you can test the program before you have real data.
    10 samples: 5 'Healthy', 5 'Disease'.
    """
    rng = np.random.default_rng(1)
    folder.mkdir(exist_ok=True)
    pops = list(POPULATIONS)
    samples = [f"S{i:02d}" for i in range(1, 11)]
    conditions = ["Healthy"] * 5 + ["Disease"] * 5
    metadata = pd.DataFrame({
        SAMPLE_COLUMN: samples,
        FACS_NAME_COLUMN: [f"{s}_urine.fcs" for s in samples],
        "condition": conditions,
        "flex_probe_barcode": [f"BC{i:03d}" for i in range(1, 11)],
    })

    typical_fraction = {"Healthy": [0.25, 0.55, 0.20], "Disease": [0.45, 0.40, 0.15]}
    protein_level = {  # typical median intensity per population and marker
        "Immune":          {"CD45": 9000, "CD10": 600,  "CD13": 1500, "EpCAM": 300},
        "Epithelial":      {"CD45": 250,  "CD10": 800,  "CD13": 900,  "EpCAM": 7000},
        "Proximal tubule": {"CD45": 250,  "CD10": 6000, "CD13": 8000, "EpCAM": 2000},
    }
    imaging_level = {"Immune": [90, 0.30, 4.0], "Epithelial": [320, 0.60, 9.0],
                     "Proximal tubule": [210, 0.50, 7.0]}

    facs_rows, cell_tables = [], []
    for sample, condition in zip(samples, conditions):
        true_fraction = rng.dirichlet(np.array(typical_fraction[condition]) * 40)
        live = int(rng.integers(20000, 60000))
        row = {"Name": f"{sample}_urine.fcs", COUNT_COLUMN.format(gate=FACS_TOTAL_GATE): live}
        marker_factor = {m: rng.lognormal(0, 0.25) for m in MARKERS}  # sample-to-sample variation
        for pop, fraction, img in zip(pops, true_fraction, imaging_level.values()):
            gate = POPULATIONS[pop]["gate"]
            row[COUNT_COLUMN.format(gate=gate)] = int(live * fraction * 0.85)
            for marker, info in MARKERS.items():
                value = protein_level[pop][marker] * marker_factor[marker] * rng.lognormal(0, 0.1)
                row[MEDIAN_COLUMN.format(gate=gate, parameter=info["facs_parameter"])] = round(value, 1)
            shift = 1.15 if condition == "Disease" and pop == "Proximal tubule" else 1.0
            for feature, base in zip(IMAGING_FEATURES, img):
                row[MEDIAN_COLUMN.format(gate=gate, parameter=feature)] = round(base * shift * rng.lognormal(0, 0.08), 3)
        facs_rows.append(row)

        # scRNA-seq: similar fractions with method-specific noise, plus a few 'Other' cells
        seq_fraction = rng.dirichlet(np.append(true_fraction, 0.05) * 60)
        n_cells = int(rng.integers(800, 2500))
        cell_pop = rng.choice(pops + ["Other"], size=n_cells, p=seq_fraction)
        labels = [rng.choice(POPULATIONS[p]["seq_labels"]) if p != "Other"
                  else rng.choice(["Podocyte", "Unknown"]) for p in cell_pop]
        cells = pd.DataFrame({"cell_barcode": [f"{sample}_{i:05d}" for i in range(n_cells)],
                              SAMPLE_COLUMN: sample, CELL_TYPE_COLUMN: labels})
        for marker, info in MARKERS.items():
            mean = np.array([np.log1p(protein_level[p][marker] / 1500 * marker_factor[marker])
                             if p != "Other" else 0.3 for p in cell_pop])
            cells[info["gene"]] = np.clip(rng.normal(mean, 0.5), 0, None).round(3)
        cell_tables.append(cells)

    metadata.to_csv(folder / METADATA_FILE, index=False)
    pd.DataFrame(facs_rows).to_excel(folder / FACS_FILE, index=False)
    pd.concat(cell_tables).to_csv(folder / SEQ_FILE, index=False)
    print(f"Demo data written to {folder}")


# =============================================================================
# 6. MAIN PROGRAM - runs the steps above in order
# =============================================================================

def main():
    global OUTPUT_DIR

    parser = argparse.ArgumentParser(description="Compare FACS (BD S8) and 10x Flex scRNA-seq cell populations.")
    parser.add_argument("--demo", action="store_true", help="create and analyse made-up demo data")
    parser.add_argument("--condition-column", help="metadata column that holds the conditions")
    parser.add_argument("--conditions", nargs="+", help="conditions to include (default: ask)")
    args = parser.parse_args()

    data_dir = DATA_DIR
    if args.demo:
        data_dir = SCRIPT_DIR / "demo_data"
        OUTPUT_DIR = SCRIPT_DIR / "results_demo"
        make_demo_data(data_dir)
    OUTPUT_DIR.mkdir(exist_ok=True)

    # --- Step 1: load files and choose conditions ---
    metadata = load_metadata(data_dir / METADATA_FILE)
    condition_column, conditions = choose_conditions(metadata, args.condition_column, args.conditions)
    metadata = metadata[metadata[condition_column].isin(conditions)]
    facs = load_facs(data_dir / FACS_FILE, metadata)
    cells = load_seq(data_dir / SEQ_FILE, metadata)
    print(f"Loaded {len(facs)} FACS samples and {len(cells)} cells "
          f"from {cells[SAMPLE_COLUMN].nunique()} scRNA-seq samples.")

    # --- Step 2: build tables ---
    proportions = build_proportions(facs, cells, metadata, condition_column)
    protein_rna = build_protein_vs_rna(facs, cells, metadata, condition_column)
    imaging = build_imaging(facs, metadata, condition_column) if IMAGING_FEATURES else None

    # --- Step 3: statistics ---
    agreement = method_agreement_stats(proportions)
    by_condition = condition_stats(proportions, "proportion", ["population", "method"],
                                   condition_column, conditions)
    rna_stats = protein_rna_stats(protein_rna)

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
    if len(conditions) > 1:
        plot_proportions_by_condition(proportions, condition_column, conditions)
    plot_protein_vs_rna(protein_rna, rna_stats)

    if imaging is not None:
        imaging_pop = imaging_between_populations(imaging)
        imaging_cond = condition_stats(imaging, "value", ["feature", "population"],
                                       condition_column, conditions)
        imaging.to_csv(OUTPUT_DIR / "D_imaging_features.csv", index=False)
        imaging_pop.to_csv(OUTPUT_DIR / "D_imaging_tests_between_populations.csv", index=False)
        imaging_cond.to_csv(OUTPUT_DIR / "D_imaging_tests_between_conditions.csv", index=False)
        plot_imaging(imaging, condition_column, conditions)

    # --- Step 6: short summary in the terminal ---
    print("\nFACS vs scRNA-seq agreement (proportions):")
    print(agreement[["population", "n", "spearman_r", "spearman_q_fdr",
                     "bland_altman_bias", "wilcoxon_q_fdr"]].round(3).to_string(index=False))
    if not by_condition.empty:
        hits = by_condition[by_condition["significant"]]
        print(f"\nSignificant condition differences in proportions (FDR < {ALPHA}): {len(hits)}")
        if not hits.empty:
            print(hits[["population", "method", "test", "comparison", "q_value_fdr"]]
                  .round(4).to_string(index=False))
    print("\nDone.")


if __name__ == "__main__":
    main()
