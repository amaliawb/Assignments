# Tests for the business logic of `facs_vs_scrnaseq.py`

This file explains the unit tests in [`test_business_logic.py`](test_business_logic.py): how they were written and what each one checks.

**Run them** from the `day04` folder:

```bash
cd day04
python -m pytest test_business_logic.py -v      # only the business-logic tests (~5 seconds)
python -m pytest -v                             # all tests: these + test_app.py (~1.5 minutes)
```

All **77 tests pass**, along with the 20 tests in `test_app.py` (97 in total).

---

## What is "business logic"?

The **business logic** is the part of the program that does the real work: the calculations and rules. In `facs_vs_scrnaseq.py` that means:

- turning FlowJo gate counts and scRNA-seq cell labels into **population proportions**;
- the **statistics**: Mann-Whitney U, Kruskal-Wallis, one-sample Wilcoxon against a pool, the bootstrap interval, correlations (with exact Spearman p-values), FDR correction;
- the **atlas comparison**: pooled atlas proportions and the composition difference;
- the **imaging calculations**: robust z-scores, outlier flags, mean/SD/CI/CV, PCA.

The rest is "plumbing": drawing the plots, the website, and saving files. That part is already checked by `test_app.py`.

## How these tests differ from `test_app.py`

| | `test_app.py` (already existed) | `test_business_logic.py` (new) |
|---|---|---|
| **What it checks** | The whole program and the website, from the outside | One function at a time, from the inside |
| **Input** | The demo data (12 Flex samples from 8 people, 7 FACS samples, thousands of cells) | Tiny tables made by hand inside each test |
| **Question it answers** | "Does everything run, and do all 40 result files appear?" | "Is this number **correct**?" |
| **Speed** | About 1.5 minutes | About 5 seconds |
| **Lab analogy** | Running the whole protocol on a control sample | Checking each reagent and each machine separately |

`test_app.py` would still pass if the FDR correction gave slightly wrong numbers, because the files would still appear. `test_business_logic.py` would catch that.

---

## How the tests were written

### 1. The idea: a known input with a known answer

Each test follows the same three steps (often called **Arrange – Act – Assert**):

1. **Arrange:** build a tiny input table where the right answer can be worked out by hand. Example: a sample with 10, 20, 30, 40, ... cells in the populations.
2. **Act:** call one function of the program on it, for example `build_proportions(...)`.
3. **Assert:** compare the result with the answer worked out by hand. Here that is 0.1, 0.2, 0.3 and 0.4. `assert` means "this must be true, otherwise the test fails."

This is like running a **standard of known concentration** through an assay: if the assay doesn't return the known value, something is wrong with the assay.

### 2. Where the expected answers come from

The expected answers were **not** copied from the program's own output; that would only test that the program agrees with itself. They come from:

- **hand calculations**, written as comments in the test. For example, the Benjamini-Hochberg q-values for p = 0.01, 0.04, 0.03 and 0.005 are worked out step by step;
- **exact values from statistics**. For example, two completely separated groups of 5 give an exact Mann-Whitney p-value of 2/252, 5 samples all above one pooled value give a one-sample Wilcoxon p-value of 2/32, and 4 points in perfect order give an exact Spearman p-value of 2/24;
- **logical rules that must always hold**. For example, proportions of one sample must add up to 1, and a q-value can never be smaller than its p-value;
- **extreme cases**. For example, a sample identical to the atlas must have a composition difference of 0, and a sample with nothing in common with it must have 1.

### 3. Helper functions build the small tables

At the top of the file there are small helpers, such as `make_facs`, `make_cells`, `make_metadata`, `make_proportions`, `make_atlas` and `make_imaging`. They build tables **in exactly the shape the real loading functions return**, so a test can be written in a few lines. For example:

```python
facs = make_facs({"s1": [10, 20, 30, 40]})     # FlowJo counts for one sample
```

### 4. The tests read the names from SETTINGS

The tests don't type population names like `"Neutrophils"` themselves. They read them from the `POPULATIONS`, `MARKERS` and `IMAGING_FEATURES` settings of the program (for the default FACS panel; the gating checks run for every panel). So if you rename a population, change a gate path or add an imaging feature, the tests still work. If `IMAGING_FEATURES` is empty, the imaging tests are **skipped** (shown as "s"), not failed.

### 5. Nothing on disk is changed

Tests that need a real file, such as a CSV for `load_seq` or an Excel file for `load_facs`, write it into pytest's **temporary folder** (`tmp_path`), which is deleted afterwards. Tests that change a setting, for example `PROPORTION_DENOMINATOR = "total"`, use pytest's `monkeypatch`, which puts the original value back after the test.

### 6. Checking that the tests can actually fail

A test that can never fail is useless. To check this, deliberate mistakes were put into a **scratch copy** of the program. The real program was not changed. All the new tests were then run against the broken copy:

| Deliberate mistake | Caught? |
|---|---|
| Composition difference without the "half" (× 1 instead of × 0.5) | ✅ 1 test failed |
| FDR correction that forgets to multiply by n / rank | ✅ 2 tests failed |
| Proportions divided by the largest population instead of the sum | ✅ 2 tests failed |
| Robust z-score without the 0.6745 factor | ✅ 1 test failed |

Each mistake was caught. With the correct program, all tests pass again.

### 7. One lesson learned while writing them

One test first failed: "a q-value is never smaller than its p-value." The reason was **computer rounding**: 0.99 × 6 / 6 comes out as 0.98999999… instead of 0.99. It wasn't a bug in the program, so the test now allows a difference of 0.000000000001. Comparing decimal numbers with `pytest.approx(...)` instead of `==` avoids this problem everywhere else in the file.

---

## What each test checks

### 0. The gating settings (no double counting)

The first five tests run once for **each FACS panel** (CD45 + EpCAM, and CD45 + EpCAM + CD10 + CD13), so they appear twice in the results.

| Test | What it checks |
|---|---|
| `test_no_population_gate_sits_inside_another` | No population gate is a child of another population gate. If one were, its cells would be counted in both populations, and the proportions would be wrong. |
| `test_every_population_gate_is_below_the_live_gate` | Every population comes from live cells (Fixable Yellow negative). |
| `test_no_cell_type_label_belongs_to_two_populations` | A scRNA-seq or atlas label isn't listed under two populations, which would also count cells twice. |
| `test_every_demo_cell_type_falls_into_exactly_one_population` | Every cell type (neutrophil, T cell, proximal tubule, ...) lands in exactly one gate of the panel, like every live cell on FACS. |
| `test_neutrophils_are_the_ssc_high_cd45_positive_cells` | In both panels, neutrophils are the CD45⁺ SSC-high gate and are matched with the Flex label "Neutrophil". |
| `test_use_panel_switches_populations_and_markers` | Switching panels changes the populations and markers. CD45 + EpCAM finds proximal tubule as CD45⁻ EpCAM-low, while the 4-marker panel uses CD10⁺CD13⁺. An unknown panel stops with a message. |
| `test_two_marker_panel_counts_podocytes_with_the_epcam_low_cells` | With CD45 + EpCAM, podocytes and proximal tubule cells (both EpCAM-low) are counted in the same population, as they are on FACS. |
| `test_missing_facs_column_message_names_the_panel` | A missing FlowJo column also asks whether the right panel was chosen. |

### 1. Loading the data

| Test | What it checks |
|---|---|
| `test_get_column_returns_the_column` | Asking for an existing column returns it. |
| `test_get_column_missing_stops_and_lists_real_columns` | A wrong column or gate name stops the program with a message that shows both the wrong name and the names that do exist. |
| `test_load_metadata_keeps_sample_ids_as_text` | Sample ID `01` stays `"01"` and doesn't become the number 1. Otherwise it wouldn't match the other files. |
| `test_load_metadata_reads_the_method_of_each_sample` | Each sample gets its method (FACS or scRNA-seq), and ways of writing it such as `facs` or `Flex` are understood. scRNA-seq samples don't need a FlowJo name. |
| `test_load_metadata_unknown_method_stops` | A method the program doesn't know (e.g. "microscopy") stops with a message naming it. |
| `test_load_metadata_facs_sample_without_flowjo_name_stops` | A FACS sample without `facs_name` is rejected, because it couldn't be found in the FlowJo export. |
| `test_load_metadata_same_sample_twice_for_one_method_stops` | The same sample ID twice for one method is rejected, because it would mix up two samples. |
| `test_load_metadata_old_sheet_without_method_means_both_methods` | An older sheet without a `method` column still works: every sample is listed once for each method. |
| `test_load_metadata_old_sheet_without_facs_name_column_stops` | An old-style sheet missing the `facs_name` column is rejected. |
| `test_choose_conditions_separately_for_each_method` | Conditions are chosen per method. A FACS-only condition such as "Healthy pool" can't be chosen for scRNA-seq. |
| `test_load_facs_drops_flowjo_summary_rows_and_uses_sample_ids` | The "Mean" and "SD" rows that FlowJo adds under the samples are ignored, and `.fcs` names are translated to sample IDs. |
| `test_load_seq_translates_labels_to_populations` | Each scRNA-seq cell label gets the right population. Unknown labels become "Other", and cells from samples not in the sample sheet are removed. |

### 2. Building the tables

| Test | What it checks |
|---|---|
| `test_build_proportions_divides_each_count_by_the_sample_total` | Counts 10/20/30/40 give proportions 0.1/0.2/0.3/0.4, for both FACS and scRNA-seq. scRNA-seq "Other" cells are not counted. |
| `test_build_proportions_sum_to_one_per_sample_and_method` | The four proportions of every sample add up to exactly 1, for both methods. |
| `test_build_proportions_total_denominator_includes_other_cells` | With `PROPORTION_DENOMINATOR = "total"`, FACS divides by the live gate and scRNA-seq by all cells including "Other". |
| `test_build_proportions_sample_without_cells_gives_nan_not_crash` | A sample with 0 cells gives "no value" (NaN) instead of a division-by-zero crash. |
| `test_build_proportions_unknown_denominator_stops` | A typo in `PROPORTION_DENOMINATOR` stops the program with a message. |
| `test_build_proportions_methods_can_have_different_samples` | A FACS pool and two other scRNA-seq people each keep their own samples, method and condition in the proportions table. |
| `test_combine_per_donor_counts_every_person_once` | Two samples of one person (3 of 10 and 7 of 30 cells) become one person with (3 + 7) / (10 + 30) = 0.25, so every person counts once. |
| `test_combine_per_donor_without_donor_column_changes_nothing` | Without a `donor_id` column, every sample stays its own person. |
| `test_composition_overview_puts_pool_flex_and_atlas_side_by_side` | The overview shows the pool's own value (0.4), the Flex average over people (0.3) and the atlas (0.5). Each group adds up to 100%. |
| `test_build_rna_markers_needs_enough_cells_for_an_rna_mean` | The RNA mean is calculated correctly (1…12 → 6.5) when there are enough cells. With fewer than `MIN_CELLS_FOR_RNA_MEAN` (10) cells it is left empty, because a mean of 5 cells is unreliable. |
| `test_build_facs_markers_reads_the_gate_medians` | The FACS marker medians are read for every sample × population × marker. |
| `test_build_protein_vs_rna_compares_typical_levels_of_different_samples` | Protein (median of 100, 200, 900 = 200) and RNA (median of 1 and 3 = 2) come from different samples and are summarised per population, with the number of samples. Populations without data still get a row. |

### 3. Statistics

| Test | What it checks |
|---|---|
| `test_fdr_correct_matches_hand_calculation` | Benjamini-Hochberg q-values match a calculation done by hand. |
| `test_fdr_correct_is_never_smaller_than_p_and_never_above_one` | Correction can only make p-values larger, and never above 1. |
| `test_fdr_correct_keeps_skipped_tests_as_nan` | Skipped tests stay empty and aren't counted as tests. |
| `test_correlate_perfect_line` | y = 2x + 1 gives Spearman r = 1 and Pearson r = 1. |
| `test_correlate_curve_spearman_one_pearson_below_one` | y = x³ (always rising, but curved) gives Spearman = 1 and Pearson < 1. This shows the difference between the two correlations. |
| `test_correlate_spearman_p_is_exact_for_few_points` | 4 points in perfect order give the exact p = 2/24. scipy's approximation would give 0. |
| `test_correlate_drops_missing_pairs_and_skips_impossible_tests` | Pairs with a missing value are dropped. With fewer than 3 pairs, or with no variation, no correlation is computed. |
| `test_compare_groups_two_groups_uses_mann_whitney` | Two groups → Mann-Whitney U, with the exact p-value (2/252) and the correct medians. |
| `test_compare_groups_three_groups_kruskal_then_every_pair` | Three groups → Kruskal-Wallis first, then every pair (A vs B, A vs C, B vs C). |
| `test_compare_groups_leaves_out_too_small_groups` | Groups with fewer than 3 values, counted after removing missing values, are left out. |
| `test_compare_to_reference_two_groups_uses_mann_whitney` | Two groups of different samples → Mann-Whitney U, exact p = 2/252. |
| `test_compare_to_reference_against_a_single_pool_uses_wilcoxon` | 5 samples against ONE pooled value → one-sample Wilcoxon, exact p = 2/32 (also the other way round). |
| `test_compare_to_reference_too_few_samples_or_no_difference_gives_no_test` | 2 samples vs 1, all values equal to the pool, or nothing to compare with → no test. |
| `test_bootstrap_interval_keeps_a_single_pool_fixed` | The bootstrap interval of the difference: exactly 0.2 when there is no spread, contains the observed difference otherwise, and is empty for 1 vs 1. |
| `test_method_comparison_against_a_pool` | scRNA-seq group vs a FACS pool: medians, difference, ratio, test name, p = 2/32 and an interval that contains the difference. |
| `test_method_comparison_two_groups_and_every_combination` | Every scRNA-seq condition is compared with every FACS condition. Two groups use Mann-Whitney. Populations that are 0 everywhere get no test and no ratio. |
| `test_condition_stats_compares_conditions_and_flags_significance` | Healthy vs Disease is compared correctly (p = 2/70), the q-value is added, and the result is marked significant. |
| `test_condition_stats_with_one_condition_gives_empty_table` | With only one condition there's nothing to compare, so the table is empty instead of the program crashing. |
| `test_condition_stats_compares_conditions_only_within_one_method` | Conditions are only compared within a method. A FACS side with one big group and a single pool gives no test. |
| `test_protein_rna_stats_same_order_of_populations_gives_r_one` | Protein and RNA that put the populations in the same order give r = 1, per marker (n = 4 populations) and for all markers together (n = 16). |
| `test_scale_0_1` | Values are rescaled from 0 (lowest) to 1 (highest). Without spread there is no scale. |

### 4b. Section E – reference atlas

| Test | What it checks |
|---|---|
| `test_load_atlas_per_cell_file` | Atlas labels are translated to the populations through `atlas_labels`. Each row is one cell. |
| `test_load_atlas_summary_table_without_donors` | In a summary table, a row with `n_cells = 250` counts as 250 cells. Without a donor column, the whole atlas is treated as one group. |
| `test_build_atlas_proportions_ignores_other_cells` | Atlas proportions use the same rule as our own samples, so they can be compared fairly. |
| `test_atlas_pooled_proportions_weights_donors_by_cell_number` | A donor with 10 cells and a donor with 100 cells are pooled by **cell number**: 91/110, not the simple average 0.5. |
| `test_composition_difference_zero_half_and_one` | Identical to the atlas → 0; half the cells different → 0.5; nothing in common → 1. |
| `test_composition_difference_scores_atlas_donors_when_there_are_several` | With several atlas donors, each donor is also scored against the pooled atlas. |
| `test_atlas_comparison_uses_mann_whitney_with_enough_donors_else_wilcoxon` | With 1 atlas donor, the program uses a one-sample Wilcoxon test against the atlas proportion. With 3 or more donors, it uses Mann-Whitney U against the donors. It also checks that "difference" = our median − atlas proportion. |
| `test_atlas_comparison_single_pool_against_the_donors` | A single FACS pool is compared with the atlas by testing the 5 donors against the pool's value (p = 2/32). |
| `test_build_marker_profiles_scales_each_marker_from_0_to_1` | The three sources (atlas RNA, Flex RNA, FACS protein) use different units, but after scaling the lowest population is 0, the highest is 1, and the three agree. |
| `test_build_marker_profiles_skipped_for_summary_atlas` | A summary-table atlas has no gene values, so the marker comparison is skipped. |

### 4c. Section F – imaging profiles

| Test | What it checks |
|---|---|
| `test_robust_z_hand_calculation` | Robust z for 1, 2, 3, 4 and 100 matches the hand calculation (median 3, MAD 1). Only 100 is above the 3.5 cutoff. |
| `test_robust_z_no_spread_gives_nan_and_ignores_missing_values` | If all values are equal, outliers can't be judged, so the result is empty. Missing values are ignored. |
| `test_imaging_outlier_samples_flags_only_the_odd_sample` | Of six samples (10, 11, 10, 12, 11, 50), only the 50 is flagged. |
| `test_imaging_population_summary_mean_sd_ci_cv` | For 10, 12, 14: mean 12, SD 2, median 12, CV 16.7%, and the t-based 95% confidence interval. Populations without data show n = 0 instead of crashing. |
| `test_imaging_pca_basic_properties` | PC1 explains at least as much as PC2, together they explain at most 100%, the scores are centred on 0, and each PC direction has length 1. |
| `test_imaging_pca_needs_enough_points_and_varying_features` | With too few samples, or only one feature that varies, no PCA is attempted. |
| `test_imaging_marker_correlation_finds_a_perfect_relation` | An imaging feature that rises exactly with a marker's intensity gives Spearman r = 1. |
| `test_cell_outlier_flags_only_the_giant_cell` | Among normal cells, one cell 10× larger is the only one flagged. |
| `test_cell_level_summary_percent_outlier_cells` | 1 odd cell out of 10 gives 10% outlier cells in that sample and 0% in the other. The cell counts and percentiles are also checked. |

---

## What these tests do **not** check

- **Whether the plots look right.** Pictures are checked by eye. `test_app.py` checks that every plot file is created.
- **Whether your gating in FlowJo is correct.** The program trusts the gates it receives.
- **Real data.** The tests use tiny made-up tables. They prove that the calculations are correct, not that the biology is.
- **The scipy and numpy libraries themselves.** The tests assume that `scipy.stats.mannwhitneyu` and similar functions are correct (they are widely used and tested). The tests check that **our program calls them correctly and uses their results correctly**.

## When to run the tests

Run them **after every change** to `facs_vs_scrnaseq.py`, especially after changing SETTINGS or the statistics. If all tests pass, the calculations still give the right answers. If one fails, its name tells you which calculation broke, and pytest shows the expected and the actual value.
