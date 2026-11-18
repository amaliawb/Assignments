"""
test_business_logic.py - unit tests for the CALCULATIONS ("business logic")
of facs_vs_scrnaseq.py.

test_app.py checks the whole program and the website from the outside (upload
files, does a results page appear?). This file checks the inside: each test calls
ONE function with a tiny, hand-made table whose correct answer we know in
advance, and compares the function's answer with it.

Run from the day04 folder:
    python -m pytest test_business_logic.py -v

No plots are drawn and nothing in day04 is changed (files are only written to
pytest's temporary folder, tmp_path). See TESTS_README.md for an explanation.
"""

import numpy as np
import pandas as pd
import pytest
from scipy import stats

import facs_vs_scrnaseq as analysis

# Short names for settings used everywhere below.
# The tests read the population / feature names from the SETTINGS of the script,
# so they keep working if you rename a population there.
S = analysis.SAMPLE_COLUMN            # "sample_id"
COND = "condition"                    # name of the condition column in our test tables
POPS = list(analysis.POPULATIONS)     # e.g. ["Immune", "Proximal tubule", ...]
FEATURES = analysis.IMAGING_FEATURES  # e.g. ["Area (LightLoss)", ...]
CUTOFF = analysis.OUTLIER_CUTOFF      # 3.5

@pytest.fixture(autouse=True)
def default_panel():
    """Every test starts and ends with the default FACS panel (some tests switch panels)."""
    analysis.use_panel(analysis.FACS_PANEL)
    yield
    analysis.use_panel(analysis.FACS_PANEL)


# Imaging tests are skipped if IMAGING_FEATURES = [] in the settings
needs_imaging = pytest.mark.skipif(not FEATURES, reason="no IMAGING_FEATURES in SETTINGS")
needs_two_features = pytest.mark.skipif(len(FEATURES) < 2, reason="PCA needs 2+ imaging features")


# =============================================================================
# HELPERS: build tiny input tables by hand
# =============================================================================

def count_column(pop):
    """FlowJo column name of a population's cell count, e.g. 'Cells/.../CD45+ | Count'."""
    return analysis.COUNT_COLUMN.format(gate=analysis.POPULATIONS[pop]["gate"])


def median_column(pop, parameter):
    """FlowJo column name of a population's median, e.g. '... | Median (Comp-FLUOR1-A :: CD45)'."""
    return analysis.MEDIAN_COLUMN.format(gate=analysis.POPULATIONS[pop]["gate"], parameter=parameter)


def make_facs(counts, live=None):
    """
    A FACS table like load_facs() returns (one row per sample).
    counts = {sample: [count of POPS[0], count of POPS[1], ...]}
    live   = {sample: count of the live gate}; default = sum of the populations
    """
    facs = pd.DataFrame({count_column(p): [counts[s][i] for s in counts] for i, p in enumerate(POPS)},
                        index=pd.Index(list(counts), name=S))
    live_column = analysis.COUNT_COLUMN.format(gate=analysis.FACS_TOTAL_GATE)
    facs[live_column] = [live[s] if live else sum(counts[s]) for s in counts]
    return facs


def make_cells(counts):
    """
    A scRNA-seq table like load_seq() returns: one row per cell.
    counts = {sample: {population: number of cells}}
    """
    rows = [{S: sample, "population": pop}
            for sample, per_pop in counts.items()
            for pop, n in per_pop.items()
            for _ in range(n)]
    return pd.DataFrame(rows)


def make_metadata(conditions, methods=("scRNA-seq", "FACS")):
    """
    A sample sheet like load_metadata() returns (with its "method" column).
    conditions = {sample: condition}; every sample is listed once for each method
    in `methods` (so the same tiny sample can be used for FACS and scRNA-seq).
    """
    return pd.concat([pd.DataFrame({S: list(conditions),
                                    analysis.FACS_NAME_COLUMN: [f"{s}.fcs" for s in conditions],
                                    COND: list(conditions.values()), "method": method})
                      for method in methods], ignore_index=True)


def make_proportions(values, condition="Healthy", methods=("FACS", "scRNA-seq")):
    """
    A proportions table like build_proportions() returns.
    values = {sample: [proportion of POPS[0], of POPS[1], ...]} (same for every method)
    """
    # (each sample gets a row for every method; real FACS and scRNA-seq samples differ,
    #  but the functions never pair them by sample, so that does not matter here)
    rows = [{S: sample, "population": pop, "method": method, "proportion": props[i], COND: condition}
            for sample, props in values.items()
            for method in methods
            for i, pop in enumerate(POPS)]
    return pd.DataFrame(rows)


def make_atlas(cells_per_donor):
    """
    An atlas table like load_atlas() returns.
    cells_per_donor = {donor: {population: number of cells}}
    """
    rows = [{analysis.ATLAS_DONOR_COLUMN: donor, "population": pop, "cells_in_row": n}
            for donor, per_pop in cells_per_donor.items()
            for pop, n in per_pop.items()]
    return pd.DataFrame(rows)


def write_csv(path, text):
    """Write a small CSV file (text with one line per row) and return its path."""
    path.write_text(text.strip() + "\n")
    return path


# =============================================================================
# 0. THE GATING SETTINGS THEMSELVES
# =============================================================================
# The proportions are only correct if every cell is counted exactly once.

# Every check runs for BOTH panels (CD45 + EpCAM, and CD45 + EpCAM + CD10 + CD13)
every_panel = pytest.mark.parametrize("panel", list(analysis.PANELS))


@every_panel
def test_no_population_gate_sits_inside_another(panel):
    # If gate B were a child of gate A (path "A/..."), B's cells would be counted twice
    gates = [info["gate"] for info in analysis.PANELS[panel]["populations"].values()]
    for a in gates:
        for b in gates:
            assert not b.startswith(a + "/"), f"{b} is inside {a}: its cells would be counted twice"


@every_panel
def test_every_population_gate_is_below_the_live_gate(panel):
    # Both panels use the Fixable Yellow live/dead gate
    for info in analysis.PANELS[panel]["populations"].values():
        assert info["gate"].startswith(analysis.LIVE_GATE + "/")


@every_panel
def test_no_cell_type_label_belongs_to_two_populations(panel):
    # A scRNA-seq or atlas label listed under two populations would be counted twice
    for key in ["seq_labels", "atlas_labels"]:
        labels = [label for info in analysis.PANELS[panel]["populations"].values() for label in info.get(key, [])]
        assert len(labels) == len(set(labels)), f"a label appears twice in {key}"


@every_panel
def test_every_demo_cell_type_falls_into_exactly_one_population(panel):
    # Every cell type must land in one gate of the panel, like every live cell on FACS
    populations = analysis.PANELS[panel]["populations"].values()
    for cell_type in analysis.DEMO_CELL_TYPES:
        assert sum(cell_type in info["seq_labels"] for info in populations) == 1, cell_type


@every_panel
def test_neutrophils_are_the_ssc_high_cd45_positive_cells(panel):
    neutrophils = analysis.PANELS[panel]["populations"]["Neutrophils"]
    assert neutrophils["gate"] == f"{analysis.LIVE_GATE}/CD45+/SSC high"
    assert "Neutrophil" in neutrophils["seq_labels"]


def test_use_panel_switches_populations_and_markers():
    analysis.use_panel("CD45_EpCAM")
    assert list(analysis.MARKERS) == ["CD45", "EpCAM"]
    epcam_low = analysis.POPULATIONS["EpCAM-low (proximal tubule)"]
    # Without CD10/CD13, proximal tubule cells are found as CD45- EpCAM low
    assert epcam_low["gate"].endswith("/CD45-/EpCAM low") and "Proximal tubule cell" in epcam_low["seq_labels"]
    analysis.use_panel("CD45_EpCAM_CD10_CD13")
    assert list(analysis.MARKERS) == ["CD45", "CD10", "CD13", "EpCAM"]
    assert analysis.POPULATIONS["Proximal tubule"]["gate"].endswith("/CD45-/CD10+ CD13+")
    with pytest.raises(SystemExit, match="Unknown FACS panel"):
        analysis.use_panel("CD3_CD19")


def test_two_marker_panel_counts_podocytes_with_the_epcam_low_cells(tmp_path):
    # CD45 + EpCAM cannot separate podocytes from proximal tubule cells (both EpCAM low),
    # so their Flex labels belong to the same population
    analysis.use_panel("CD45_EpCAM")
    genes = [m["gene"] for m in analysis.MARKERS.values()]
    lines = [f"{S},{analysis.CELL_TYPE_COLUMN}," + ",".join(genes),
             "01,Podocyte," + ",".join("1" for _ in genes),
             "01,Proximal tubule cell," + ",".join("1" for _ in genes)]
    cells = analysis.load_seq(write_csv(tmp_path / "cells.csv", "\n".join(lines)),
                              make_metadata({"01": "A"}, methods=["scRNA-seq"]))
    assert set(cells["population"]) == {"EpCAM-low (proximal tubule)"}


def test_missing_facs_column_message_names_the_panel():
    with pytest.raises(SystemExit, match="FACS panel is right"):
        analysis.get_column(pd.DataFrame({"a": [1]}), "missing", "the FACS file")


# =============================================================================
# 1. LOADING THE DATA
# =============================================================================

def test_get_column_returns_the_column():
    table = pd.DataFrame({"a": [1, 2], "b": [3, 4]})
    assert list(analysis.get_column(table, "b", "test file")) == [3, 4]


def test_get_column_missing_stops_and_lists_real_columns():
    table = pd.DataFrame({"Real column": [1]})
    with pytest.raises(SystemExit) as error:
        analysis.get_column(table, "Typo column", "test file")
    # The message names the missing column AND the columns that do exist
    assert "Typo column" in str(error.value)
    assert "Real column" in str(error.value)


METHOD = analysis.METHOD_COLUMN       # "method"
NAME = analysis.FACS_NAME_COLUMN      # "facs_name"


def test_load_metadata_keeps_sample_ids_as_text(tmp_path):
    path = write_csv(tmp_path / "metadata.csv", f"{S},{METHOD},{NAME},{COND}\n01,FACS,a.fcs,Healthy")
    metadata = analysis.load_metadata(path)
    assert metadata.loc[0, S] == "01"        # not the number 1


def test_load_metadata_reads_the_method_of_each_sample(tmp_path):
    # Different samples per method; method names written in different ways are understood
    path = write_csv(tmp_path / "metadata.csv", "\n".join([
        f"{S},{METHOD},{NAME},{COND}",
        "P01,scRNA-seq,,Disease",
        "P02,Flex,,Disease",
        "POOL,facs,pool.fcs,Healthy pool",
    ]))
    metadata = analysis.load_metadata(path)
    assert list(metadata["method"]) == ["scRNA-seq", "scRNA-seq", "FACS"]
    assert list(analysis.samples_of(metadata, "FACS")[S]) == ["POOL"]


def test_load_metadata_unknown_method_stops(tmp_path):
    path = write_csv(tmp_path / "metadata.csv", f"{S},{METHOD},{NAME},{COND}\n01,microscopy,a.fcs,Healthy")
    with pytest.raises(SystemExit, match="microscopy"):
        analysis.load_metadata(path)


def test_load_metadata_facs_sample_without_flowjo_name_stops(tmp_path):
    # A FACS sample without facs_name cannot be found in the FlowJo export
    path = write_csv(tmp_path / "metadata.csv", f"{S},{METHOD},{NAME},{COND}\nPOOL,FACS,,Healthy")
    with pytest.raises(SystemExit, match="POOL"):
        analysis.load_metadata(path)


def test_load_metadata_same_sample_twice_for_one_method_stops(tmp_path):
    path = write_csv(tmp_path / "metadata.csv", f"{S},{METHOD},{NAME},{COND}\n01,FACS,a.fcs,A\n01,FACS,b.fcs,B")
    with pytest.raises(SystemExit, match="twice"):
        analysis.load_metadata(path)


def test_load_metadata_old_sheet_without_method_means_both_methods(tmp_path):
    # An older sheet (every sample measured by both methods) is still accepted
    path = write_csv(tmp_path / "metadata.csv", f"{S},{NAME},{COND}\n01,a.fcs,Healthy\n02,b.fcs,Disease")
    metadata = analysis.load_metadata(path)
    assert sorted(zip(metadata[S], metadata["method"])) == [
        ("01", "FACS"), ("01", "scRNA-seq"), ("02", "FACS"), ("02", "scRNA-seq")]


def test_load_metadata_old_sheet_without_facs_name_column_stops(tmp_path):
    path = write_csv(tmp_path / "metadata.csv", f"{S},{COND}\n01,Healthy")
    with pytest.raises(SystemExit):
        analysis.load_metadata(path)


def test_choose_conditions_separately_for_each_method():
    metadata = pd.concat([make_metadata({"p1": "Disease", "p2": "Disease"}, methods=["scRNA-seq"]),
                          make_metadata({"pool": "Healthy pool", "p3": "Disease"}, methods=["FACS"])])
    column, conditions = analysis.choose_conditions(metadata, COND, ["Disease"], ["Healthy pool"])
    assert column == COND
    assert conditions == {"scRNA-seq": ["Disease"], "FACS": ["Healthy pool"]}
    # "Healthy pool" exists only among the FACS samples, so it cannot be chosen for scRNA-seq
    with pytest.raises(SystemExit, match="not found among the scRNA-seq samples"):
        analysis.choose_conditions(metadata, COND, ["Healthy pool"], ["Disease"])


def test_load_facs_drops_flowjo_summary_rows_and_uses_sample_ids(tmp_path):
    # FlowJo adds 'Mean' and 'SD' rows under the samples: they must be ignored
    export = pd.DataFrame({"Name": ["s1.fcs", "s2.fcs", "Mean", "SD"],
                           count_column(POPS[0]): [100, 200, 150, 70]})
    path = tmp_path / "facs.xlsx"
    export.to_excel(path, index=False)
    facs = analysis.load_facs(path, make_metadata({"s1": "Healthy", "s2": "Disease"}, methods=["FACS"]))
    assert list(facs.index) == ["s1", "s2"]
    assert list(facs[count_column(POPS[0])]) == [100, 200]


def test_load_seq_translates_labels_to_populations(tmp_path):
    genes = [m["gene"] for m in analysis.MARKERS.values()]
    # One cell for the first label of every population, plus one unknown label
    labels = [info["seq_labels"][0] for info in analysis.POPULATIONS.values()] + ["Mystery cell"]
    lines = [f"{S},{analysis.CELL_TYPE_COLUMN}," + ",".join(genes)]
    lines += [f"01,{label}," + ",".join("1" for _ in genes) for label in labels]
    lines += [f"99,{labels[0]}," + ",".join("1" for _ in genes)]   # sample not in the metadata
    path = write_csv(tmp_path / "cells.csv", "\n".join(lines))

    cells = analysis.load_seq(path, make_metadata({"01": "Healthy"}, methods=["scRNA-seq"]))
    assert list(cells["population"]) == POPS + ["Other"]   # unknown label -> "Other"
    assert set(cells[S]) == {"01"}                          # cell of sample 99 removed


# =============================================================================
# 2. BUILDING THE TABLES
# =============================================================================

def test_build_proportions_divides_each_count_by_the_sample_total():
    facs_counts = [10 * (i + 1) for i in range(len(POPS))]       # 10, 20, 30, 40
    facs = make_facs({"s1": facs_counts})
    # scRNA-seq: same counts plus 100 "Other" cells, which must NOT be counted
    seq_counts = {pop: facs_counts[i] for i, pop in enumerate(POPS)}
    cells = make_cells({"s1": {**seq_counts, "Other": 100}})

    table = analysis.build_proportions(facs, cells, make_metadata({"s1": "Healthy"}), COND)
    for method in ["FACS", "scRNA-seq"]:
        part = table[table["method"] == method].set_index("population")
        expected = np.array(facs_counts) / sum(facs_counts)
        assert part.loc[POPS, "proportion"].values == pytest.approx(expected)


def test_build_proportions_sum_to_one_per_sample_and_method():
    facs = make_facs({"s1": [5, 0, 7, 3][:len(POPS)] + [1] * (len(POPS) - 4),
                      "s2": [1] * len(POPS)})
    cells = make_cells({"s1": {POPS[0]: 3, POPS[1]: 9}, "s2": {POPS[-1]: 4, "Other": 2}})
    table = analysis.build_proportions(facs, cells, make_metadata({"s1": "A", "s2": "B"}), COND)
    sums = table.groupby([S, "method"])["proportion"].sum()
    assert sums.values == pytest.approx(np.ones(4))


def test_build_proportions_total_denominator_includes_other_cells(monkeypatch):
    # With "total", FACS divides by the live gate and scRNA-seq by ALL cells incl. "Other"
    monkeypatch.setattr(analysis, "PROPORTION_DENOMINATOR", "total")
    counts = [10] * len(POPS)
    facs = make_facs({"s1": counts}, live={"s1": 200})
    cells = make_cells({"s1": {**{pop: 10 for pop in POPS}, "Other": 200 - 10 * len(POPS)}})
    table = analysis.build_proportions(facs, cells, make_metadata({"s1": "Healthy"}), COND)
    assert table["proportion"].values == pytest.approx(np.full(2 * len(POPS), 10 / 200))


def test_build_proportions_sample_without_cells_gives_nan_not_crash():
    facs = make_facs({"s1": [0] * len(POPS), "s2": [1] * len(POPS)})
    cells = make_cells({"s2": {POPS[0]: 5}})
    table = analysis.build_proportions(facs, cells, make_metadata({"s1": "A", "s2": "A"}), COND)
    assert table.loc[table[S] == "s1", "proportion"].isna().all()


def test_build_proportions_unknown_denominator_stops(monkeypatch):
    monkeypatch.setattr(analysis, "PROPORTION_DENOMINATOR", "something else")
    facs = make_facs({"s1": [1] * len(POPS)})
    with pytest.raises(SystemExit):
        analysis.build_proportions(facs, make_cells({"s1": {POPS[0]: 1}}), make_metadata({"s1": "A"}), COND)


def test_build_proportions_methods_can_have_different_samples():
    # FACS: one pool; scRNA-seq: two other people. Each keeps its own samples and condition.
    facs = make_facs({"pool": [10] * len(POPS)})
    cells = make_cells({"p1": {POPS[0]: 2, POPS[1]: 2}, "p2": {POPS[0]: 4}})
    metadata = pd.concat([make_metadata({"pool": "Healthy pool"}, methods=["FACS"]),
                          make_metadata({"p1": "Disease", "p2": "Disease"}, methods=["scRNA-seq"])])
    table = analysis.build_proportions(facs, cells, metadata, COND)
    assert set(zip(table[S], table["method"], table[COND])) == {
        ("pool", "FACS", "Healthy pool"), ("p1", "scRNA-seq", "Disease"), ("p2", "scRNA-seq", "Disease")}
    p2 = table[table[S] == "p2"].set_index("population")["proportion"]
    assert p2[POPS[0]] == pytest.approx(1.0)


def test_build_rna_markers_needs_enough_cells_for_an_rna_mean():
    # 12 cells of POPS[0] (values 1..12, mean 6.5) and only 5 cells of POPS[1]
    many = analysis.MIN_CELLS_FOR_RNA_MEAN + 2
    cells = make_cells({"s1": {POPS[0]: many, POPS[1]: 5}})
    for marker in analysis.MARKERS.values():
        cells[marker["gene"]] = np.r_[np.arange(1, many + 1), np.zeros(5)]

    table = analysis.build_rna_markers(cells, make_metadata({"s1": "A"}, methods=["scRNA-seq"]), COND)
    first_marker = list(analysis.MARKERS)[0]
    row = table[(table["marker"] == first_marker)].set_index("population")
    assert row.loc[POPS[0], "rna_mean"] == pytest.approx((many + 1) / 2)   # mean of 1..12
    assert np.isnan(row.loc[POPS[1], "rna_mean"])                          # 5 cells: too few
    assert row.loc[POPS[1], "n_cells_rna"] == 5


def test_build_facs_markers_reads_the_gate_medians():
    facs = make_facs({"s1": [1] * len(POPS), "s2": [1] * len(POPS)})
    for pop in POPS:
        for marker in analysis.MARKERS.values():
            facs[median_column(pop, marker["facs_parameter"])] = [500.0, 700.0]
    table = analysis.build_facs_markers(facs, make_metadata({"s1": "A", "s2": "A"}, methods=["FACS"]), COND)
    assert len(table) == 2 * len(POPS) * len(analysis.MARKERS)
    assert sorted(set(table["protein_median"])) == [500.0, 700.0]


def test_build_protein_vs_rna_compares_typical_levels_of_different_samples():
    marker = list(analysis.MARKERS)[0]
    # FACS: 3 samples with protein 100, 200, 900 -> median 200; scRNA-seq: 2 other samples, RNA 1 and 3 -> median 2
    facs_markers = pd.DataFrame({S: ["f1", "f2", "f3"], "population": POPS[0], "marker": marker,
                                 "protein_median": [100.0, 200.0, 900.0]})
    rna_markers = pd.DataFrame({S: ["r1", "r2"], "population": POPS[0], "marker": marker,
                                "rna_mean": [1.0, 3.0]})
    table = analysis.build_protein_vs_rna(facs_markers, rna_markers).set_index(["population", "marker"])
    row = table.loc[(POPS[0], marker)]
    assert (row["protein_median"], row["rna_median"]) == (200.0, 2.0)
    assert (row["n_facs_samples"], row["n_seq_samples"]) == (3, 2)
    # Every population x marker has a row, also without data (n = 0, no value)
    assert len(table) == len(POPS) * len(analysis.MARKERS)
    assert table.loc[(POPS[1], marker), "n_facs_samples"] == 0


def test_combine_per_donor_counts_every_person_once():
    # Person P1 gave two samples (3 of 10 and 7 of 30 cells), person P2 one sample
    sample_props = pd.DataFrame({
        S: ["s1", "s2", "s3"], "method": "scRNA-seq", "population": POPS[0],
        "count": [3, 7, 5], "total": [10, 30, 20], COND: "Disease"})
    sample_props["proportion"] = sample_props["count"] / sample_props["total"]
    metadata = make_metadata({"s1": "Disease", "s2": "Disease", "s3": "Disease"}, methods=["scRNA-seq"])
    metadata[analysis.DONOR_COLUMN] = ["P1", "P1", "P2"]
    table = analysis.combine_proportions_per_donor(sample_props, metadata, COND).set_index(S)
    assert list(table.index) == ["P1", "P2"]
    # P1: the cells of both samples added up: (3 + 7) / (10 + 30)
    assert table.loc["P1", "proportion"] == pytest.approx(10 / 40)
    assert table.loc["P1", "n_samples_combined"] == 2
    assert table.loc["P2", "proportion"] == pytest.approx(5 / 20)


def test_combine_per_donor_without_donor_column_changes_nothing():
    sample_props = pd.DataFrame({S: ["s1", "s2"], "method": "FACS", "population": POPS[0],
                                 "count": [1, 2], "total": [10, 10], COND: "A"})
    table = analysis.combine_proportions_per_donor(sample_props, make_metadata({"s1": "A", "s2": "A"}), COND)
    assert list(table[S]) == ["s1", "s2"] and (table["n_samples_combined"] == 1).all()


def test_composition_overview_puts_pool_flex_and_atlas_side_by_side():
    rest = [0.0] * (len(POPS) - 2)
    proportions = pd.concat([
        make_proportions({"pool": [0.4, 0.6] + rest}, "Healthy pool", ["FACS"]),
        make_proportions({"P1": [0.2, 0.8] + rest, "P2": [0.4, 0.6] + rest}, "Healthy", ["scRNA-seq"]),
    ], ignore_index=True)
    atlas_props = analysis.build_atlas_proportions(make_atlas({"d1": {POPS[0]: 1, POPS[1]: 1}}))
    overview = analysis.composition_overview(proportions, COND, {"FACS": ["Healthy pool"],
                                                                 "scRNA-seq": ["Healthy"]}, atlas_props)
    first = overview[overview["population"] == POPS[0]].set_index("label")
    assert first.loc["Healthy pool (FACS)", "mean"] == pytest.approx(0.4)        # the pool: its own value
    assert first.loc["Healthy (scRNA-seq)", "mean"] == pytest.approx(0.3)        # mean of 0.2 and 0.4
    assert first.loc["Healthy (scRNA-seq)", "n"] == 2
    assert first.loc["Reference atlas", "mean"] == pytest.approx(0.5)
    # The means of each group add up to 1 (they are drawn as one stacked bar)
    assert overview.groupby("label")["mean"].sum().values == pytest.approx([1, 1, 1])
    wide = analysis.composition_overview_wide(overview)
    assert list(wide.columns) == ["population", "Healthy pool (FACS) (n=1), %",
                                  "Healthy (scRNA-seq) (n=2), %", "Reference atlas (n=1), %"]


# =============================================================================
# 3. STATISTICS
# =============================================================================

def test_fdr_correct_matches_hand_calculation():
    # Benjamini-Hochberg by hand: sort p, multiply by n/rank, keep monotonic
    #   sorted p : 0.005  0.01  0.03  0.04
    #   x 4/rank : 0.02   0.02  0.04  0.04
    q = analysis.fdr_correct([0.01, 0.04, 0.03, 0.005])
    assert q == pytest.approx([0.02, 0.04, 0.04, 0.02])


def test_fdr_correct_is_never_smaller_than_p_and_never_above_one():
    p = np.array([0.001, 0.2, 0.5, 0.9, 0.95, 0.99])
    q = analysis.fdr_correct(p)
    # (1e-12 allows for computer rounding: 0.99 * 6 / 6 can come out as 0.98999999...)
    assert np.all(q >= p - 1e-12)
    assert np.all(q <= 1)


def test_fdr_correct_keeps_skipped_tests_as_nan():
    q = analysis.fdr_correct([0.01, np.nan, 0.02])
    assert np.isnan(q[1])
    # The NaN is not counted as a test: n = 2, not 3
    assert q[0] == pytest.approx(0.02) and q[2] == pytest.approx(0.02)
    assert np.isnan(analysis.fdr_correct([np.nan, np.nan])).all()


def test_correlate_perfect_line():
    result = analysis.correlate([1, 2, 3, 4, 5], [3, 5, 7, 9, 11])     # y = 2x + 1
    assert result["spearman_r"] == pytest.approx(1)
    assert result["pearson_r"] == pytest.approx(1)


def test_correlate_curve_spearman_one_pearson_below_one():
    # Always going up, but not in a straight line: rank correlation is perfect, linear is not
    result = analysis.correlate([1, 2, 3, 4, 5], [1, 8, 27, 64, 125])
    assert result["spearman_r"] == pytest.approx(1)
    assert result["pearson_r"] < 1


def test_correlate_spearman_p_is_exact_for_few_points():
    # 4 points in perfect order: 2 of the 24 possible orderings are this extreme
    # (perfectly rising or perfectly falling) -> exact p = 2/24, not 0
    result = analysis.correlate([1, 2, 3, 4], [10, 20, 30, 40])
    assert result["spearman_p"] == pytest.approx(2 / 24)
    # Many points: scipy's own p-value is used and agrees with the rank correlation
    x = np.arange(20)
    assert analysis.correlate(x, x ** 2)["spearman_p"] < 1e-6


def test_correlate_drops_missing_pairs_and_skips_impossible_tests():
    result = analysis.correlate([1, 2, np.nan, 4, 5], [2, 4, 6, np.nan, 10])
    assert result["n"] == 3                                              # 2 pairs dropped
    assert np.isnan(analysis.correlate([1, 2], [1, 2])["spearman_r"])    # too few pairs
    assert np.isnan(analysis.correlate([1, 2, 3], [5, 5, 5])["pearson_r"])   # no variation


def test_compare_groups_two_groups_uses_mann_whitney():
    results = analysis.compare_groups({"A": [1, 2, 3, 4, 5], "B": [10, 11, 12, 13, 14]})
    assert len(results) == 1
    assert results[0]["test"] == "Mann-Whitney U"
    assert results[0]["median_group1"] == 3 and results[0]["median_group2"] == 12
    # Completely separated groups of 5: exact p = 2 / 252
    assert results[0]["p_value"] == pytest.approx(2 / 252)


def test_compare_groups_three_groups_kruskal_then_every_pair():
    results = analysis.compare_groups({"A": [1, 2, 3], "B": [4, 5, 6], "C": [7, 8, 9]})
    assert [r["test"] for r in results] == ["Kruskal-Wallis"] + ["Mann-Whitney U"] * 3
    assert [r["comparison"] for r in results[1:]] == ["A vs B", "A vs C", "B vs C"]


def test_compare_groups_leaves_out_too_small_groups():
    # C has only 2 values, and A has only 2 after removing NaN -> nothing to compare
    assert analysis.compare_groups({"A": [1, 2, np.nan], "C": [1, 2]}) == []
    results = analysis.compare_groups({"A": [1, 2, 3], "B": [4, 5, 6], "C": [7, 8]})
    assert [r["comparison"] for r in results] == ["A vs B"]


def test_compare_to_reference_two_groups_uses_mann_whitney():
    test, _, p = analysis.compare_to_reference([10, 11, 12, 13, 14], [1, 2, 3, 4, 5])
    assert test == "Mann-Whitney U"
    assert p == pytest.approx(2 / 252)          # completely separated groups of 5


def test_compare_to_reference_against_a_single_pool_uses_wilcoxon():
    # 5 samples all above one pooled value: smallest possible two-sided Wilcoxon p = 2/32
    test, _, p = analysis.compare_to_reference([0.3, 0.31, 0.32, 0.33, 0.34], [0.1])
    assert test == "Wilcoxon vs reference value" and p == pytest.approx(2 / 32)
    # The other way round: one value against a group
    test, _, _ = analysis.compare_to_reference([0.1], [0.3, 0.31, 0.32, 0.33, 0.34])
    assert test == "Wilcoxon of reference vs value"


def test_compare_to_reference_too_few_samples_or_no_difference_gives_no_test():
    assert analysis.compare_to_reference([0.1, 0.2], [0.3])[0] == ""        # 2 vs 1: too few
    assert analysis.compare_to_reference([0.2, 0.2, 0.2], [0.2])[0] == ""   # all equal to the pool
    assert analysis.compare_to_reference([0.1, 0.2, 0.3], [])[0] == ""      # nothing to compare with


def test_bootstrap_interval_keeps_a_single_pool_fixed():
    # All scRNA-seq samples are 0.3 and the pool is 0.1: every resampling gives 0.2
    low, high = analysis.bootstrap_median_difference([0.3, 0.3, 0.3], [0.1])
    assert (low, high) == (pytest.approx(0.2), pytest.approx(0.2))
    # The interval contains the observed difference and widens with spread
    low, high = analysis.bootstrap_median_difference([0.2, 0.3, 0.4, 0.5, 0.6], [0.1])
    assert low <= 0.3 <= high and high - low > 0
    assert np.isnan(analysis.bootstrap_median_difference([0.3], [0.1])[0])   # 1 vs 1: no interval


def make_method_groups(seq, facs, seq_condition="Disease", facs_condition="Healthy pool"):
    """Proportions of separate scRNA-seq and FACS samples; seq/facs = list of POPS[0] proportions."""
    rest = [0.0] * (len(POPS) - 1)
    return pd.concat([
        make_proportions({f"seq{i}": [v] + rest for i, v in enumerate(seq)}, seq_condition, ["scRNA-seq"]),
        make_proportions({f"facs{i}": [v] + rest for i, v in enumerate(facs)}, facs_condition, ["FACS"]),
    ], ignore_index=True)


def test_method_comparison_against_a_pool():
    proportions = make_method_groups(seq=[0.3, 0.4, 0.5, 0.6, 0.7], facs=[0.25])
    table = analysis.method_comparison_stats(proportions, COND, {"scRNA-seq": ["Disease"],
                                                                  "FACS": ["Healthy pool"]})
    row = table[table["population"] == POPS[0]].iloc[0]
    assert row["comparison"] == "Disease (scRNA-seq) vs Healthy pool (FACS)"
    assert (row["n_seq"], row["n_facs"]) == (5, 1)
    assert row["seq_minus_facs"] == pytest.approx(0.5 - 0.25)
    assert row["seq_to_facs_ratio"] == pytest.approx(2.0)
    assert row["test"] == "Wilcoxon: scRNA-seq vs single FACS value"
    assert row["p_value"] == pytest.approx(2 / 32)
    assert row["seq_minus_facs_ci95_low"] <= 0.25 <= row["seq_minus_facs_ci95_high"]


def test_method_comparison_two_groups_and_every_combination():
    proportions = pd.concat([
        make_method_groups(seq=[0.5, 0.6, 0.7, 0.8, 0.9], facs=[0.1, 0.2, 0.3, 0.4, 0.45], facs_condition="Disease"),
        make_proportions({"pool": [0.2] + [0.0] * (len(POPS) - 1)}, "Healthy pool", ["FACS"]),
    ], ignore_index=True)
    table = analysis.method_comparison_stats(proportions, COND, {"scRNA-seq": ["Disease"],
                                                                  "FACS": ["Disease", "Healthy pool"]})
    # 1 scRNA-seq condition x 2 FACS conditions x every population
    assert len(table) == 2 * len(POPS)
    row = table[(table["population"] == POPS[0]) & (table["facs_group"] == "Disease")].iloc[0]
    assert row["test"] == "Mann-Whitney U" and row["p_value"] == pytest.approx(2 / 252)
    # Populations that are 0 everywhere: no difference, no test, no ratio (cannot divide by 0)
    empty = table[table["population"] == POPS[1]]
    assert (empty["test"] == "").all() and empty["p_value"].isna().all()
    assert empty["seq_to_facs_ratio"].isna().all()
    # Only POPS[0] was tested: once against the FACS Disease group, once against the pool
    assert table["q_value_fdr"].notna().sum() == 2


def test_condition_stats_compares_conditions_and_flags_significance():
    table = pd.DataFrame({
        "population": [POPS[0]] * 8,
        "value": [1, 2, 3, 4, 10, 11, 12, 13],
        COND: ["Healthy"] * 4 + ["Disease"] * 4,
    })
    result = analysis.condition_stats(table, "value", ["population"], COND, ["Healthy", "Disease"])
    assert len(result) == 1
    assert result.loc[0, "p_value"] == pytest.approx(2 / 70)    # exact Mann-Whitney, 4 vs 4
    assert result.loc[0, "q_value_fdr"] == pytest.approx(2 / 70)  # only one test: q = p
    assert bool(result.loc[0, "significant"])


def test_condition_stats_with_one_condition_gives_empty_table():
    table = pd.DataFrame({"population": [POPS[0]] * 3, "value": [1, 2, 3], COND: ["Healthy"] * 3})
    assert analysis.condition_stats(table, "value", ["population"], COND, ["Healthy"]).empty


def test_condition_stats_compares_conditions_only_within_one_method():
    # scRNA-seq has Healthy and Disease; FACS has Disease and a single pool -> only scRNA-seq is tested
    rows = ([("scRNA-seq", "Healthy", v) for v in [1, 2, 3]] + [("scRNA-seq", "Disease", v) for v in [7, 8, 9]] +
            [("FACS", "Disease", v) for v in [4, 5, 6]] + [("FACS", "Healthy pool", 1)])
    table = pd.DataFrame(rows, columns=["method", COND, "value"]).assign(population=POPS[0])
    result = analysis.condition_stats(table, "value", ["population", "method"], COND,
                                      ["Healthy", "Disease", "Healthy pool"])
    assert list(result["method"]) == ["scRNA-seq"]
    assert result.loc[0, "comparison"] == "Healthy vs Disease"


def test_protein_rna_stats_same_order_of_populations_gives_r_one():
    # Every marker: protein and RNA rise together over the populations (different units)
    rows = [{"population": pop, "marker": m, "protein_median": 1000.0 * (i + 1), "rna_median": 0.5 * (i + 1)}
            for m in analysis.MARKERS for i, pop in enumerate(POPS)]
    table = analysis.protein_rna_stats(pd.DataFrame(rows)).set_index("marker")
    assert table["spearman_r"].values == pytest.approx(np.ones(len(analysis.MARKERS) + 1))
    # One row per marker (n = number of populations) + all markers together after 0-1 scaling
    assert table.loc[list(analysis.MARKERS)[0], "n"] == len(POPS)
    assert table.loc["All markers (scaled)", "n"] == len(POPS) * len(analysis.MARKERS)


def test_scale_0_1():
    assert list(analysis.scale_0_1(pd.Series([10.0, 20.0, 30.0]))) == [0, 0.5, 1]
    assert analysis.scale_0_1(pd.Series([5.0, 5.0])).isna().all()        # no spread: no scale


# =============================================================================
# 4b. SECTION E - REFERENCE ATLAS
# =============================================================================

def test_load_atlas_per_cell_file(tmp_path):
    labels = [info.get("atlas_labels", info["seq_labels"])[0] for info in analysis.POPULATIONS.values()]
    lines = [f"{analysis.ATLAS_CELL_TYPE_COLUMN},{analysis.ATLAS_DONOR_COLUMN}"]
    lines += [f"{label},d1" for label in labels] + ["Mystery cell,d2"]
    atlas = analysis.load_atlas(write_csv(tmp_path / "atlas.csv", "\n".join(lines)))
    assert list(atlas["population"]) == POPS + ["Other"]
    assert (atlas["cells_in_row"] == 1).all()                 # one row = one cell


def test_load_atlas_summary_table_without_donors(tmp_path):
    label = analysis.POPULATIONS[POPS[0]].get("atlas_labels", analysis.POPULATIONS[POPS[0]]["seq_labels"])[0]
    text = f"{analysis.ATLAS_CELL_TYPE_COLUMN},{analysis.ATLAS_COUNT_COLUMN}\n{label},250"
    atlas = analysis.load_atlas(write_csv(tmp_path / "atlas.csv", text))
    assert atlas.loc[0, "cells_in_row"] == 250                # one row = 250 cells
    assert atlas.loc[0, analysis.ATLAS_DONOR_COLUMN] == "atlas"   # whole atlas = one group


def test_build_atlas_proportions_ignores_other_cells():
    atlas = make_atlas({"d1": {POPS[0]: 30, POPS[1]: 70, "Other": 50}})
    props = analysis.build_atlas_proportions(atlas).set_index("population")["proportion"]
    assert props[POPS[0]] == pytest.approx(0.3)
    assert props[POPS[1]] == pytest.approx(0.7)
    assert props.sum() == pytest.approx(1)


def test_atlas_pooled_proportions_weights_donors_by_cell_number():
    # Small donor: 1 of 10 cells; large donor: 90 of 100 cells.
    # Pooled = 91 / 110, NOT the plain average of 0.1 and 0.9 (= 0.5)
    atlas = make_atlas({"small": {POPS[0]: 1, POPS[1]: 9}, "large": {POPS[0]: 90, POPS[1]: 10}})
    pooled = analysis.atlas_pooled_proportions(analysis.build_atlas_proportions(atlas))
    assert pooled[POPS[0]] == pytest.approx(91 / 110)


def test_composition_difference_zero_half_and_one():
    atlas_props = analysis.build_atlas_proportions(make_atlas({"d1": {POPS[0]: 100}}))  # atlas = only POPS[0]
    zeros = [0.0] * (len(POPS) - 2)
    proportions = make_proportions({
        "same": [1.0, 0.0] + zeros,     # identical to the atlas         -> 0
        "half": [0.5, 0.5] + zeros,     # half of the cells differ        -> 0.5
        "other": [0.0, 1.0] + zeros,    # nothing in common with the atlas -> 1
    }, methods=["FACS"])
    result = analysis.composition_difference(proportions, atlas_props, COND).set_index("sample")
    assert result.loc["same", "composition_difference"] == pytest.approx(0)
    assert result.loc["half", "composition_difference"] == pytest.approx(0.5)
    assert result.loc["other", "composition_difference"] == pytest.approx(1)
    assert "Atlas donor" not in set(result["method"])     # a single donor is not scored


def test_composition_difference_scores_atlas_donors_when_there_are_several():
    atlas_props = analysis.build_atlas_proportions(
        make_atlas({"d1": {POPS[0]: 50, POPS[1]: 50}, "d2": {POPS[0]: 50, POPS[1]: 50}}))
    proportions = make_proportions({"s1": [0.5, 0.5] + [0.0] * (len(POPS) - 2)}, methods=["FACS"])
    result = analysis.composition_difference(proportions, atlas_props, COND)
    donors = result[result["method"] == "Atlas donor"]
    assert len(donors) == 2
    assert donors["composition_difference"].values == pytest.approx([0, 0])


def test_atlas_comparison_uses_mann_whitney_with_enough_donors_else_wilcoxon():
    proportions = make_proportions({f"s{i}": [0.1 * i] + [0.0] * (len(POPS) - 1) for i in range(1, 5)})
    one_donor = analysis.build_atlas_proportions(make_atlas({"d1": {POPS[0]: 9, POPS[1]: 1}}))
    many_donors = analysis.build_atlas_proportions(
        make_atlas({f"d{i}": {POPS[0]: 9, POPS[1]: i} for i in range(1, 5)}))

    conditions = {"scRNA-seq": ["Healthy"], "FACS": ["Healthy"]}

    table = analysis.atlas_comparison_stats(proportions, one_donor, COND, conditions)
    first = table[table["population"] == POPS[0]].iloc[0]
    assert first["test"] == "Wilcoxon vs atlas proportion"
    assert first["atlas_proportion"] == pytest.approx(0.9)
    # difference = median of our samples (0.25) minus the atlas (0.9)
    assert first["difference"] == pytest.approx(0.25 - 0.9)

    table = analysis.atlas_comparison_stats(proportions, many_donors, COND, conditions)
    assert table.loc[table["population"] == POPS[0], "test"].iloc[0] == "Mann-Whitney U vs atlas donors"


def test_atlas_comparison_single_pool_against_the_donors():
    # One FACS pool: the atlas donors are tested against the pool's value
    proportions = make_proportions({"pool": [0.5] + [0.0] * (len(POPS) - 1)}, "Healthy pool", ["FACS"])
    donors = analysis.build_atlas_proportions(
        make_atlas({f"d{i}": {POPS[0]: 9, POPS[1]: i} for i in range(1, 6)}))
    table = analysis.atlas_comparison_stats(proportions, donors, COND,
                                            {"scRNA-seq": [], "FACS": ["Healthy pool"]})
    row = table[table["population"] == POPS[0]].iloc[0]
    assert (row["method"], row["n_samples"]) == ("FACS", 1)
    assert row["test"] == "Wilcoxon of atlas donors vs sample value"
    assert row["p_value"] == pytest.approx(2 / 32)     # all 5 donors above 0.5


def test_build_marker_profiles_scales_each_marker_from_0_to_1():
    genes = [m["gene"] for m in analysis.MARKERS.values()]
    # Every source: population number i has expression i (POPS[0] lowest, POPS[-1] highest)
    atlas = pd.DataFrame([{"population": pop, **{g: i for g in genes}} for i, pop in enumerate(POPS)])
    cells = pd.DataFrame([{"population": pop, **{g: 10 * i for g in genes}} for i, pop in enumerate(POPS)])
    protein_rna = pd.DataFrame([{"population": pop, "marker": m, "protein_median": 1000 * i}
                                for i, pop in enumerate(POPS) for m in analysis.MARKERS])
    table = analysis.build_marker_profiles(atlas, cells, protein_rna)
    # Different units in each source, but after scaling all three agree:
    # lowest population -> 0, highest -> 1, the others evenly in between
    for i, pop in enumerate(POPS):
        assert table.loc[table["population"] == pop, "scaled"].values == pytest.approx(
            np.full(3 * len(genes), i / (len(POPS) - 1)))


def test_build_marker_profiles_skipped_for_summary_atlas():
    atlas = pd.DataFrame({"population": [POPS[0]], analysis.ATLAS_COUNT_COLUMN: [10]})
    assert analysis.build_marker_profiles(atlas, None, None) is None


# =============================================================================
# 4c. SECTION F - IMAGING PROFILES
# =============================================================================

def test_robust_z_hand_calculation():
    # median = 3, distances from median = 2,1,0,1,97 -> MAD = 1
    z = analysis.robust_z([1, 2, 3, 4, 100])
    assert z == pytest.approx(0.6745 * np.array([-2, -1, 0, 1, 97]))
    assert abs(z[-1]) > CUTOFF and np.all(np.abs(z[:-1]) < CUTOFF)


def test_robust_z_no_spread_gives_nan_and_ignores_missing_values():
    assert np.isnan(analysis.robust_z([5, 5, 5, 5])).all()
    z = analysis.robust_z([1, 2, np.nan, 3, 4, 100])
    assert np.isnan(z[2]) and z[-1] > CUTOFF


def make_imaging(values, pop=None, feature=None):
    """Imaging table like build_imaging(): one sample per value, one population and feature."""
    return pd.DataFrame({S: [f"s{i}" for i in range(len(values))], COND: "Healthy",
                         "population": pop or POPS[0], "feature": feature or FEATURES[0],
                         "value": values})


@needs_imaging
def test_imaging_outlier_samples_flags_only_the_odd_sample():
    imaging = make_imaging([10, 11, 10, 12, 11, 50])
    table = analysis.imaging_outlier_samples(imaging, COND)
    assert list(table["outlier"]) == [False] * 5 + [True]


@needs_imaging
def test_imaging_population_summary_mean_sd_ci_cv():
    imaging = make_imaging([10, 12, 14])
    summary = analysis.imaging_population_summary(imaging, analysis.imaging_outlier_samples(imaging, COND))
    row = summary[(summary["population"] == POPS[0]) & (summary["feature"] == FEATURES[0])].iloc[0]
    half_width = stats.t.ppf(0.975, 2) * 2 / np.sqrt(3)       # t-based 95% CI, n = 3, SD = 2
    assert (row["n_samples"], row["mean"], row["sd"], row["median"]) == (3, 12, 2, 12)
    assert row["ci95_low"] == pytest.approx(12 - half_width)
    assert row["ci95_high"] == pytest.approx(12 + half_width)
    assert row["cv_percent"] == pytest.approx(100 * 2 / 12)
    # Populations without data are listed with n = 0 instead of crashing
    assert (summary.loc[summary["population"] != POPS[0], "n_samples"] == 0).all()


@needs_two_features
def test_imaging_pca_basic_properties():
    rng = np.random.default_rng(0)
    imaging = pd.concat([make_imaging(rng.normal(100, 10, 8), feature=f) for f in FEATURES])
    scores, loadings = analysis.imaging_pca(imaging, COND)
    assert len(scores) == 8 and list(loadings["feature"]) == FEATURES
    # PC1 explains the most, and the two together cannot explain more than 100%
    assert scores["PC1_explained"].iloc[0] >= scores["PC2_explained"].iloc[0]
    assert scores["PC1_explained"].iloc[0] + scores["PC2_explained"].iloc[0] <= 1 + 1e-9
    # Scores are centred on 0; each PC direction has length 1
    assert scores["PC1"].mean() == pytest.approx(0, abs=1e-9)
    assert np.sum(loadings["PC1"] ** 2) == pytest.approx(1)


@needs_two_features
def test_imaging_pca_needs_enough_points_and_varying_features():
    two_samples = pd.concat([make_imaging([1.0, 2.0], feature=f) for f in FEATURES])
    assert analysis.imaging_pca(two_samples, COND) == (None, None)
    # Only one feature varies -> nothing for a 2D map
    one_varies = pd.concat([make_imaging([1.0, 2.0, 3.0, 4.0], feature=FEATURES[0])] +
                           [make_imaging([5.0] * 4, feature=f) for f in FEATURES[1:]])
    assert analysis.imaging_pca(one_varies, COND) == (None, None)


@needs_imaging
def test_imaging_marker_correlation_finds_a_perfect_relation():
    imaging = make_imaging([1, 2, 3, 4, 5])
    marker = list(analysis.MARKERS)[0]
    protein_rna = pd.DataFrame({S: [f"s{i}" for i in range(5)], "population": POPS[0],
                                "marker": marker, "protein_median": [10, 20, 30, 40, 50]})
    table = analysis.imaging_marker_correlations(imaging, protein_rna)
    row = table[(table["population"] == POPS[0]) & (table["feature"] == FEATURES[0]) &
                (table["marker"] == marker)].iloc[0]
    assert row["spearman_r"] == pytest.approx(1)


def make_cells_img(n_per_sample, odd_value=None):
    """
    Per-cell imaging table for two samples (s1, s2) of POPS[0]. Every feature has
    the values 100, 101, 102, 103, 104 repeated (no outliers); optionally the first
    cell of s1 gets odd_value in FEATURES[0].
    """
    rows = []
    for sample in ["s1", "s2"]:
        for i in range(n_per_sample):
            rows.append({S: sample, COND: "Healthy", "population": POPS[0],
                         **{f: 100 + i % 5 for f in FEATURES}})
    cells_img = pd.DataFrame(rows)
    if odd_value is not None:
        cells_img.loc[0, FEATURES[0]] = odd_value
    return cells_img


@needs_imaging
def test_cell_outlier_flags_only_the_giant_cell():
    cells_img = make_cells_img(10, odd_value=1000)
    flags = analysis.cell_outlier_flags(cells_img)
    assert list(flags.index[flags["Any feature"]]) == [0]
    assert flags.loc[0, FEATURES[0]]


@needs_imaging
def test_cell_level_summary_percent_outlier_cells():
    cells_img = make_cells_img(10, odd_value=1000)
    summary = analysis.cell_level_summary(cells_img, analysis.cell_outlier_flags(cells_img), COND)
    first = summary[summary["feature"] == FEATURES[0]].set_index(S)
    assert first.loc["s1", "n_cells"] == 10
    assert first.loc["s1", "pct_outlier_cells"] == pytest.approx(10)   # 1 of 10 cells
    assert first.loc["s2", "pct_outlier_cells"] == pytest.approx(0)
    # Middle 90% of s2's cells: between the 5th and 95th percentile of 100..104
    assert first.loc["s2", "p5"] == pytest.approx(np.percentile([100 + i % 5 for i in range(10)], 5))


# =============================================================================
# 5. THE DEFAULT FLEX DATA (.h5ad -> cells table + sample sheet)
# =============================================================================

GENES = [m["gene"] for m in analysis.ALL_MARKERS.values()]     # PTPRC, MME, ANPEP, EPCAM


def write_fake_h5ad(path, cells, extra_genes=("GAPDH", "ACTB")):
    """
    Write a tiny .h5ad like the real Flex file. cells = a table with one row per
    cell: sample_id, leiden_atlas, patient_id, role, timepoint and one column per
    marker gene. Two extra genes are added (with values) so the reader must pick
    the right columns out of the sparse matrix.
    """
    anndata = pytest.importorskip("anndata")
    from scipy import sparse
    genes = [*extra_genes[:1], *GENES, *extra_genes[1:]]          # marker genes in the middle
    values = np.column_stack([np.full(len(cells), 9.0) if g in extra_genes else cells[g].to_numpy(float)
                              for g in genes])
    obs = cells[[S, "leiden_atlas", "patient_id", "role", "timepoint"]].astype("category")
    obs.index = [f"cell{i}" for i in range(len(cells))]
    adata = anndata.AnnData(X=sparse.csr_matrix(values * 0 + 1), obs=obs,
                            var=pd.DataFrame(index=genes))
    adata.layers["log1p_cp10k"] = sparse.csr_matrix(values)
    adata.write_h5ad(path)
    return path


def fake_flex_cells():
    """5 cells: sample A (person P1, healthy) and samples B1, B2 (both person P2, recipient)."""
    return pd.DataFrame({
        S: ["A", "A", "B1", "B1", "B2"],
        "leiden_atlas": ["4: T_NK", "6: myeloid", "6: myeloid", "3: ascending_thick_LOH_kidney ?", "6: myeloid"],
        "patient_id": ["P1", "P1", "P2", "P2", "P2"],
        "role": ["healthy_control", "healthy_control", "recipient", "recipient", "recipient"],
        "timepoint": ["control", "control", "T00", "T00", "T07"],
        "PTPRC": [2.0, 1.5, 0.0, 0.0, 3.0],
        "MME": [0.0, 0.0, 0.5, 0.0, 0.0],
        "ANPEP": [0.0, 0.25, 0.0, 0.0, 0.0],
        "EPCAM": [0.0, 0.0, 0.0, 1.75, 0.0],
    })


@pytest.fixture
def fake_flex(tmp_path, monkeypatch):
    """Point the default Flex settings at a tiny fake .h5ad and a temporary cache folder."""
    source = write_fake_h5ad(tmp_path / "flex.h5ad", fake_flex_cells())
    monkeypatch.setattr(analysis, "DEFAULT_FLEX_H5AD", source)
    monkeypatch.setattr(analysis, "DEFAULT_FLEX_CACHE_DIR", tmp_path / "cache")
    return source


def test_read_h5ad_cells_takes_the_right_genes_from_the_sparse_matrix(fake_flex):
    table = analysis.read_h5ad_cells(fake_flex, GENES, "log1p_cp10k", [S, "leiden_atlas"])
    expected = fake_flex_cells()
    assert list(table[S]) == list(expected[S])
    for gene in GENES:            # the values of each marker gene, not of the extra genes (9.0)
        assert table[gene].tolist() == pytest.approx(expected[gene].tolist())


def test_read_h5ad_cells_reports_a_missing_gene_or_layer(fake_flex):
    with pytest.raises(SystemExit, match="NOTAGENE"):
        analysis.read_h5ad_cells(fake_flex, ["NOTAGENE"], "log1p_cp10k", [S])
    with pytest.raises(SystemExit, match="no layer 'counts'"):
        analysis.read_h5ad_cells(fake_flex, GENES, "counts", [S])


def test_default_flex_tables_make_one_sheet_row_per_sample(fake_flex):
    cells, sheet = analysis.default_flex_tables()
    assert list(cells.columns) == [S, analysis.CELL_TYPE_COLUMN, *GENES]
    assert len(cells) == 5
    sheet = sheet.set_index(S)
    assert list(sheet.index) == ["A", "B1", "B2"]
    assert set(sheet[analysis.METHOD_COLUMN]) == {"scRNA-seq"}
    # B1 and B2 come from the same person, so they are combined later
    assert sheet.loc["B1", analysis.DONOR_COLUMN] == sheet.loc["B2", analysis.DONOR_COLUMN] == "P2"
    assert sheet.loc["A", "condition"] == "healthy_control" and sheet.loc["B2", "timepoint"] == "T07"


def test_default_flex_cache_is_used_when_the_drive_is_missing(fake_flex, tmp_path, monkeypatch):
    analysis.default_flex_tables()                                   # builds the cache
    monkeypatch.setattr(analysis, "DEFAULT_FLEX_H5AD", tmp_path / "not_connected.h5ad")
    cells, sheet = analysis.default_flex_tables()                    # from the cache only
    assert len(cells) == 5 and len(sheet) == 3
    # no cache and no file: a clear message
    monkeypatch.setattr(analysis, "DEFAULT_FLEX_CACHE_DIR", tmp_path / "empty_cache")
    with pytest.raises(SystemExit, match="default Flex data was not found"):
        analysis.default_flex_tables()


def test_add_default_flex_adds_the_flex_samples_to_a_facs_only_sheet(fake_flex, tmp_path):
    data = tmp_path / "data"
    data.mkdir()
    pd.DataFrame({S: ["POOL"], "method": ["FACS"], "facs_name": ["pool.fcs"],
                  "condition": ["Healthy pool"]}).to_csv(data / analysis.METADATA_FILE, index=False)
    assert analysis.add_default_flex(data) is True
    metadata = analysis.load_metadata(data / analysis.METADATA_FILE)
    assert len(analysis.samples_of(metadata, "FACS")) == 1
    assert sorted(analysis.samples_of(metadata, "scRNA-seq")[S]) == ["A", "B1", "B2"]
    assert len(pd.read_csv(data / analysis.SEQ_FILE)) == 5


def test_add_default_flex_keeps_your_own_scrnaseq_rows_and_files(fake_flex, tmp_path):
    data = tmp_path / "data"
    data.mkdir()
    sheet = pd.DataFrame({S: ["POOL", "A"], "method": ["FACS", "Flex"], "facs_name": ["pool.fcs", None],
                          "condition": ["Healthy pool", "my own label"]})
    sheet.to_csv(data / analysis.METADATA_FILE, index=False)
    # The sheet already lists scRNA-seq samples: they are used as written, nothing is added
    assert analysis.add_default_flex(data) is True
    metadata = analysis.load_metadata(data / analysis.METADATA_FILE)
    assert list(analysis.samples_of(metadata, "scRNA-seq")["condition"]) == ["my own label"]
    # An uploaded scRNA-seq file is never replaced
    (data / analysis.SEQ_FILE).write_text("my own file")
    assert analysis.add_default_flex(data) is False
    assert (data / analysis.SEQ_FILE).read_text() == "my own file"


@pytest.mark.parametrize("panel", list(analysis.PANELS))
def test_every_flex_cluster_label_is_in_exactly_one_population(panel):
    analysis.use_panel(panel)
    flex_labels = ["0: epithelium_prostate_urethra ?", "1: epithelium_prostate_urethra ?",
                   "2: epithelium_prostate_urethra ?", "3: ascending_thick_LOH_kidney ?", "4: T_NK",
                   "5: myeloid", "6: myeloid", "7: myeloid", "8: myeloid", "9: myeloid", "10: myeloid",
                   "11: myeloid ?", "12: B_plasma"]
    for label in flex_labels:
        homes = [name for name, info in analysis.POPULATIONS.items() if label in info["seq_labels"]]
        assert len(homes) == 1, (label, homes)
    # Neutrophil clusters (CSF3R/FCGR3B high) are neutrophils; macrophage clusters are other immune
    assert "6: myeloid" in analysis.POPULATIONS["Neutrophils"]["seq_labels"]
    assert "9: myeloid" in analysis.POPULATIONS["Other immune"]["seq_labels"]


def test_condition_candidates_leave_out_ids_barcodes_and_empty_columns():
    metadata = pd.DataFrame({S: ["a", "b"], "method": ["FACS", "scRNA-seq"], "facs_name": ["a.fcs", None],
                             "donor_id": ["P1", "P2"], "condition": ["Healthy pool", "Disease"],
                             "flex_probe_barcode": [None, "BC001"], "raw_file": ["x", "y"],
                             "batch_id": ["1", "2"], "notes": [None, None], "timepoint": [None, "T00"]})
    assert analysis.condition_candidates(metadata) == ["condition", "timepoint"]
