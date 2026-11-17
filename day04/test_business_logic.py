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


def make_metadata(conditions):
    """A sample sheet. conditions = {sample: condition}"""
    return pd.DataFrame({S: list(conditions),
                         analysis.FACS_NAME_COLUMN: [f"{s}.fcs" for s in conditions],
                         COND: list(conditions.values())})


def make_proportions(values, condition="Healthy", methods=("FACS", "scRNA-seq")):
    """
    A proportions table like build_proportions() returns.
    values = {sample: [proportion of POPS[0], of POPS[1], ...]} (same for every method)
    """
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

def test_no_population_gate_sits_inside_another():
    # If gate B were a child of gate A (path "A/..."), B's cells would be counted twice
    gates = [info["gate"] for info in analysis.POPULATIONS.values()]
    for a in gates:
        for b in gates:
            assert not b.startswith(a + "/"), f"{b} is inside {a}: its cells would be counted twice"


def test_every_population_gate_is_below_the_live_gate():
    for info in analysis.POPULATIONS.values():
        assert info["gate"].startswith(analysis.LIVE_GATE + "/")


def test_no_cell_type_label_belongs_to_two_populations():
    # A scRNA-seq or atlas label listed under two populations would be counted twice
    for key in ["seq_labels", "atlas_labels"]:
        labels = [label for info in analysis.POPULATIONS.values() for label in info.get(key, [])]
        assert len(labels) == len(set(labels)), f"a label appears twice in {key}"


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


def test_load_metadata_keeps_sample_ids_as_text(tmp_path):
    path = write_csv(tmp_path / "metadata.csv", f"{S},{analysis.FACS_NAME_COLUMN},{COND}\n01,a.fcs,Healthy")
    metadata = analysis.load_metadata(path)
    assert metadata.loc[0, S] == "01"        # not the number 1


def test_load_metadata_without_facs_name_column_stops(tmp_path):
    path = write_csv(tmp_path / "metadata.csv", f"{S},{COND}\n01,Healthy")
    with pytest.raises(SystemExit):
        analysis.load_metadata(path)


def test_load_facs_drops_flowjo_summary_rows_and_uses_sample_ids(tmp_path):
    # FlowJo adds 'Mean' and 'SD' rows under the samples: they must be ignored
    export = pd.DataFrame({"Name": ["s1.fcs", "s2.fcs", "Mean", "SD"],
                           count_column(POPS[0]): [100, 200, 150, 70]})
    path = tmp_path / "facs.xlsx"
    export.to_excel(path, index=False)
    facs = analysis.load_facs(path, make_metadata({"s1": "Healthy", "s2": "Disease"}))
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

    cells = analysis.load_seq(path, make_metadata({"01": "Healthy"}))
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


def test_build_protein_vs_rna_needs_enough_cells_for_an_rna_mean():
    # FACS table: every population x marker median = 500
    facs = make_facs({"s1": [1] * len(POPS)})
    for pop in POPS:
        for marker in analysis.MARKERS.values():
            facs[median_column(pop, marker["facs_parameter"])] = 500.0
    # 12 cells of POPS[0] (values 1..12, mean 6.5) and only 5 cells of POPS[1]
    many = analysis.MIN_CELLS_FOR_RNA_MEAN + 2
    cells = make_cells({"s1": {POPS[0]: many, POPS[1]: 5}})
    for marker in analysis.MARKERS.values():
        cells[marker["gene"]] = np.r_[np.arange(1, many + 1), np.zeros(5)]

    table = analysis.build_protein_vs_rna(facs, cells, make_metadata({"s1": "A"}), COND)
    first_marker = list(analysis.MARKERS)[0]
    row = table[(table["marker"] == first_marker)].set_index("population")
    assert row.loc[POPS[0], "rna_mean"] == pytest.approx((many + 1) / 2)   # mean of 1..12
    assert np.isnan(row.loc[POPS[1], "rna_mean"])                          # 5 cells: too few
    assert row.loc[POPS[1], "n_cells_rna"] == 5
    assert (table["protein_median"] == 500.0).all()


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


def test_method_agreement_identical_methods_have_no_bias():
    values = {f"s{i}": [0.1 * i] * len(POPS) for i in range(1, 6)}
    table = analysis.method_agreement_stats(make_proportions(values))
    assert table["bland_altman_bias"].values == pytest.approx(np.zeros(len(POPS)))
    assert table["spearman_r"].values == pytest.approx(np.ones(len(POPS)))
    # All differences are 0, so the paired Wilcoxon test is skipped
    assert table["wilcoxon_p"].isna().all()


def test_method_agreement_finds_a_constant_shift():
    # scRNA-seq always gives 0.05 more than FACS
    facs = make_proportions({f"s{i}": [0.1 * i] * len(POPS) for i in range(1, 6)}, methods=["FACS"])
    seq = facs.assign(method="scRNA-seq", proportion=facs["proportion"] + 0.05)
    table = analysis.method_agreement_stats(pd.concat([facs, seq]))
    assert table["bland_altman_bias"].values == pytest.approx(np.full(len(POPS), 0.05))
    assert table["bland_altman_lower"].values == pytest.approx(np.full(len(POPS), 0.05))
    # 5 samples, all shifted the same way: smallest possible two-sided Wilcoxon p = 2/32
    assert table["wilcoxon_p"].values == pytest.approx(np.full(len(POPS), 2 / 32))


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

    table = analysis.atlas_comparison_stats(proportions, one_donor, COND, ["Healthy"])
    first = table[table["population"] == POPS[0]].iloc[0]
    assert first["test"] == "Wilcoxon vs atlas proportion"
    assert first["atlas_proportion"] == pytest.approx(0.9)
    # difference = median of our samples (0.25) minus the atlas (0.9)
    assert first["difference"] == pytest.approx(0.25 - 0.9)

    table = analysis.atlas_comparison_stats(proportions, many_donors, COND, ["Healthy"])
    assert table.loc[table["population"] == POPS[0], "test"].iloc[0] == "Mann-Whitney U vs atlas donors"


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
