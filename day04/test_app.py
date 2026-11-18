"""
test_app.py - automatic checks for the analysis script and the web app.

Run from the day04 folder:
    python -m pytest test_app.py -v

Every test works in a temporary folder, so nothing in day04 is changed.
"""

import io
import subprocess
import sys
import zipfile
from pathlib import Path

import pandas as pd
import pytest

import app as web
import facs_vs_scrnaseq as analysis

HERE = Path(__file__).resolve().parent
EXPECTED_FILES = [
    "A1_proportions_stacked.png", "A2_proportions_by_method.png", "A3_method_difference.png",
    "A4_composition_overview.png", "A_composition_overview.csv", "A_label_assignment.csv", "notes.txt",
    "A_method_comparison_stats.csv", "A_proportions.csv", "A_proportions_per_sample.csv",
    "B1_proportions_by_condition.png", "B_condition_tests_proportions.csv",
    "C1_protein_vs_rna.png", "C_protein_per_sample.csv", "C_rna_per_sample.csv",
    "C_protein_vs_rna.csv", "C_protein_vs_rna_stats.csv",
    "D1_imaging_features.png", "D_imaging_features.csv",
    "D_imaging_tests_between_conditions.csv", "D_imaging_tests_between_populations.csv",
    "E1_atlas_proportions.png", "E2_atlas_composition_difference.png", "E3_atlas_marker_profiles.png",
    "E_atlas_composition_difference.csv", "E_atlas_marker_profiles.csv",
    "E_atlas_proportion_tests.csv", "E_atlas_proportions.csv",
    "F1_imaging_feature_means.png", "F2_imaging_fingerprint.png", "F3_imaging_pca.png",
    "F4_imaging_outlier_samples.png", "F5_imaging_vs_markers.png",
    "F6_imaging_cell_distributions.png", "F7_imaging_cell_outliers.png", "F8_imaging_morphology_map.png",
    "F_imaging_cell_summary.csv", "F_imaging_marker_correlations.csv", "F_imaging_outlier_samples.csv",
    "F_imaging_pca.csv", "F_imaging_pca_loadings.csv", "F_imaging_population_summary.csv",
]
ATLAS_FILES = [f for f in EXPECTED_FILES if f.startswith("E")]
# Section F files that need the optional per-cell imaging export
CELL_LEVEL_FILES = ["F6_imaging_cell_distributions.png", "F7_imaging_cell_outliers.png",
                    "F8_imaging_morphology_map.png", "F_imaging_cell_summary.csv"]
SAMPLE_LEVEL_F_FILES = [f for f in EXPECTED_FILES if f.startswith("F") and f not in CELL_LEVEL_FILES]

# The demo data: scRNA-seq of 12 samples from 8 people (4 Healthy, 4 Disease) and
# FACS of OTHER people (one healthy pool + 6 Disease patients)
SEQ_CONDITIONS = ["Healthy", "Disease"]
FACS_CONDITIONS = ["Healthy pool", "Disease"]


@pytest.fixture(autouse=True)
def default_panel():
    """Every test starts and ends with the default FACS panel (some tests switch panels)."""
    analysis.use_panel(analysis.FACS_PANEL)
    yield
    analysis.use_panel(analysis.FACS_PANEL)


@pytest.fixture
def demo_data(tmp_path):
    folder = tmp_path / "data"
    analysis.make_demo_data(folder)
    return folder


@pytest.fixture
def client(tmp_path, monkeypatch):
    monkeypatch.setattr(web, "RUNS_DIR", tmp_path / "runs")
    web.app.config["TESTING"] = True
    return web.app.test_client()


# ---- the analysis script (business logic) ------------------------------------

def test_run_analysis_writes_all_files(demo_data, tmp_path):
    out = tmp_path / "results"
    result = analysis.run_analysis(demo_data, out, "condition", SEQ_CONDITIONS, FACS_CONDITIONS)
    assert sorted(p.name for p in out.iterdir()) == sorted(EXPECTED_FILES)
    assert (result["n_facs_samples"], result["n_seq_samples"], result["n_seq_people"]) == (7, 12, 8)
    # the statistics count PEOPLE: samples of the same person are combined
    assert set(result["comparison"]["n_seq"]) == {4}
    # every scRNA-seq condition against every FACS condition, for every population
    comparison = result["comparison"]
    assert len(comparison) == len(SEQ_CONDITIONS) * len(FACS_CONDITIONS) * len(analysis.POPULATIONS)
    pool = comparison[comparison["facs_group"] == "Healthy pool"]
    assert (pool["n_facs"] == 1).all()
    assert set(pool["test"]) == {"Wilcoxon: scRNA-seq vs single FACS value"}
    assert not result["by_condition"].empty


def test_run_analysis_one_condition_per_method_skips_condition_plot(demo_data, tmp_path):
    out = tmp_path / "results"
    result = analysis.run_analysis(demo_data, out, "condition", ["Disease"], ["Disease"])
    assert not (out / "B1_proportions_by_condition.png").exists()
    assert result["by_condition"].empty
    assert result["n_facs_samples"] == 6
    assert set(result["comparison"]["test"]) == {"Mann-Whitney U"}


def test_run_analysis_unknown_condition_stops_with_message(demo_data, tmp_path):
    with pytest.raises(SystemExit, match="not found"):
        analysis.run_analysis(demo_data, tmp_path / "results", "condition", ["Nope"], FACS_CONDITIONS)
    # "Healthy pool" is a FACS condition only
    with pytest.raises(SystemExit, match="scRNA-seq"):
        analysis.run_analysis(demo_data, tmp_path / "results", "condition", ["Healthy pool"], FACS_CONDITIONS)


def test_command_line_still_works(tmp_path):
    """Run the script from the terminal the old way, on a copy so day04 stays clean."""
    script = tmp_path / "facs_vs_scrnaseq.py"
    script.write_text((HERE / "facs_vs_scrnaseq.py").read_text())
    done = subprocess.run([sys.executable, str(script), "--demo", "--condition-column", "condition",
                           "--seq-conditions", *SEQ_CONDITIONS, "--facs-conditions", *FACS_CONDITIONS],
                          capture_output=True, text=True, cwd=tmp_path)
    assert done.returncode == 0, done.stderr
    assert "Done." in done.stdout
    assert len(list((tmp_path / "results_demo").iterdir())) == len(EXPECTED_FILES)


# ---- the web app ---------------------------------------------------------------

def test_start_page_shows_upload_form_and_settings(client):
    page = client.get("/").get_data(as_text=True)
    assert 'type="file" name="facs"' in page
    assert "Proximal tubule" in page and "Other CD45-" in page and "PTPRC" in page


def test_demo_run_end_to_end(client):
    response = client.post("/upload", data={"demo": "1"})
    assert response.status_code == 302
    setup_url = response.headers["Location"]
    run_id = setup_url.rstrip("/").split("/")[-1]

    setup_page = client.get(setup_url).get_data(as_text=True)
    assert "scRNA-seq samples" in setup_page and "FACS samples" in setup_page
    assert "Healthy pool" in setup_page and "(1 sample)" in setup_page
    # 6 Healthy Flex samples come from 4 people; the 6 FACS Disease samples are 6 people
    assert "(6 samples from 4 people)" in setup_page and "(6 samples)" in setup_page

    response = client.post(f"/run/{run_id}/analyse",
                           data={"condition_column": "condition", "seq_conditions": SEQ_CONDITIONS,
                                 "facs_conditions": FACS_CONDITIONS})
    assert response.status_code == 302

    page = client.get(f"/run/{run_id}/results").get_data(as_text=True)
    for name in EXPECTED_FILES:
        assert name in page
    # every section has its question, and every plot has a "how to read it" explanation
    assert page.count("<b>Question:</b>") == 6
    assert page.count("How to read it") == len([f for f in EXPECTED_FILES if f.endswith(".png")])
    assert all(plot in web.explanations.PLOT_HELP for plot in EXPECTED_FILES if plot.endswith(".png"))
    assert client.get(f"/run/{run_id}/file/A1_proportions_stacked.png").mimetype == "image/png"

    archive = zipfile.ZipFile(io.BytesIO(client.get(f"/run/{run_id}/download").data))
    assert sorted(archive.namelist()) == sorted(EXPECTED_FILES)


def test_upload_own_files(client, demo_data):
    files = {
        "metadata": (open(demo_data / "metadata.csv", "rb"), "my_samples.csv"),
        "facs": (open(demo_data / "facs_flowjo_export.xlsx", "rb"), "flowjo.xlsx"),
        "seq": (open(demo_data / "scrnaseq_cells.csv", "rb"), "cells.csv"),
    }
    response = client.post("/upload", data=files, content_type="multipart/form-data")
    assert response.status_code == 302
    run_id = response.headers["Location"].rstrip("/").split("/")[-1]

    # no condition ticked = all conditions of each method, like pressing Enter in the terminal
    client.post(f"/run/{run_id}/analyse", data={"condition_column": "condition"})
    page = client.get(f"/run/{run_id}/results").get_data(as_text=True)
    assert "Healthy, Disease" in page and "Healthy pool, Disease" in page
    # the atlas is optional: without it section E is skipped with a message
    assert "No reference atlas was given" in page
    assert not any(name in page for name in ATLAS_FILES)
    # without the per-cell imaging export: sample-level imaging plots only, plus a note
    assert all(name in page for name in SAMPLE_LEVEL_F_FILES)
    assert not any(name in page for name in CELL_LEVEL_FILES)
    assert "The cell-level plots (F6–F8) need the optional per-cell imaging export" in page


def test_setup_page_explains_the_condition_column(client, demo_data, tmp_path):
    # Add a column that only the scRNA-seq samples have: it is shown, but can't be chosen
    metadata = pd.read_csv(demo_data / "metadata.csv", dtype=str)
    metadata["sex"] = metadata["method"].map({"scRNA-seq": "F", "FACS": None})
    sheet = tmp_path / "with_sex.csv"
    metadata.to_csv(sheet, index=False)
    files = {"metadata": (open(sheet, "rb"), "samples.csv"),
             "facs": (open(demo_data / "facs_flowjo_export.xlsx", "rb"), "flowjo.xlsx"),
             "seq": (open(demo_data / "scrnaseq_cells.csv", "rb"), "cells.csv")}
    run_id = client.post("/upload", data=files, content_type="multipart/form-data").headers["Location"].split("/")[-1]
    page = client.get(f"/run/{run_id}").get_data(as_text=True)

    assert "Which column holds the groups?" in page and "Which groups to include?" in page
    assert "Which group each sample belongs to" in page                     # help text of "condition"
    assert "scRNA-seq: Healthy, Disease" in page and "FACS: Healthy pool, Disease" in page
    assert 'value="flex_probe_barcode"' not in page                         # a barcode is not a group
    assert "Can't be used: the FACS samples have no value here" in page
    response = client.post(f"/run/{run_id}/analyse", data={"condition_column": "sex"})
    assert response.status_code == 400
    assert "no values for the FACS samples" in response.get_data(as_text=True)


def test_upload_per_cell_imaging_export(client, demo_data):
    files = {
        "metadata": (open(demo_data / "metadata.csv", "rb"), "my_samples.csv"),
        "facs": (open(demo_data / "facs_flowjo_export.xlsx", "rb"), "flowjo.xlsx"),
        "seq": (open(demo_data / "scrnaseq_cells.csv", "rb"), "cells.csv"),
        "facs_cells": (open(demo_data / analysis.FACS_CELLS_FILE, "rb"), "single_cells.csv"),
    }
    response = client.post("/upload", data=files, content_type="multipart/form-data")
    run_id = response.headers["Location"].rstrip("/").split("/")[-1]
    client.post(f"/run/{run_id}/analyse", data={"condition_column": "condition"})
    page = client.get(f"/run/{run_id}/results").get_data(as_text=True)
    assert all(name in page for name in CELL_LEVEL_FILES)
    assert "need the optional per-cell imaging export" not in page


def test_robust_z_flags_an_outlier():
    values = [10, 11, 9, 10, 12, 10, 11, 40]          # 40 is far from the rest
    z = analysis.robust_z(values)
    assert abs(z[-1]) > analysis.OUTLIER_CUTOFF
    assert all(abs(v) < analysis.OUTLIER_CUTOFF for v in z[:-1])
    # no spread at all: outliers cannot be judged -> NaN instead of a crash
    assert all(pd.isna(analysis.robust_z([5, 5, 5])))


def test_no_imaging_features_skips_sections_d_and_f(demo_data, tmp_path, monkeypatch):
    monkeypatch.setattr(analysis, "IMAGING_FEATURES", [])
    out = tmp_path / "results"
    analysis.run_analysis(demo_data, out, "condition", SEQ_CONDITIONS, FACS_CONDITIONS)
    made = {p.name for p in out.iterdir()}
    assert not any(name.startswith(("D", "F")) for name in made)
    assert "A1_proportions_stacked.png" in made


def test_upload_own_files_with_atlas(client, demo_data):
    files = {
        "metadata": (open(demo_data / "metadata.csv", "rb"), "my_samples.csv"),
        "facs": (open(demo_data / "facs_flowjo_export.xlsx", "rb"), "flowjo.xlsx"),
        "seq": (open(demo_data / "scrnaseq_cells.csv", "rb"), "cells.csv"),
        "atlas": (open(demo_data / "reference_atlas.csv", "rb"), "urine_atlas.csv"),
    }
    response = client.post("/upload", data=files, content_type="multipart/form-data")
    run_id = response.headers["Location"].rstrip("/").split("/")[-1]
    client.post(f"/run/{run_id}/analyse", data={"condition_column": "condition"})
    page = client.get(f"/run/{run_id}/results").get_data(as_text=True)
    for name in ATLAS_FILES:
        assert name in page
    assert "No reference atlas was given" not in page


def test_upload_without_scrnaseq_file_uses_the_default_flex_data(client, demo_data, tmp_path, monkeypatch):
    # A small fake "default Flex" .h5ad made from the demo cells (the real one is on the X: drive)
    from test_business_logic import write_fake_h5ad
    cells = pd.read_csv(demo_data / analysis.SEQ_FILE, dtype={analysis.SAMPLE_COLUMN: str})
    metadata = pd.read_csv(demo_data / analysis.METADATA_FILE, dtype=str)
    seq_rows = metadata[metadata["method"] == "scRNA-seq"].set_index(analysis.SAMPLE_COLUMN)
    cells = cells.assign(leiden_atlas=cells[analysis.CELL_TYPE_COLUMN],
                         patient_id=cells[analysis.SAMPLE_COLUMN].map(seq_rows["donor_id"]),
                         role=cells[analysis.SAMPLE_COLUMN].map(seq_rows["condition"]), timepoint="T00")
    monkeypatch.setattr(analysis, "DEFAULT_FLEX_H5AD", write_fake_h5ad(tmp_path / "flex.h5ad", cells))
    monkeypatch.setattr(analysis, "DEFAULT_FLEX_CACHE_DIR", tmp_path / "cache")

    # Upload only the FACS part of the sample sheet and the FlowJo export, no scRNA-seq file
    facs_only = tmp_path / "facs_samples.csv"
    metadata[metadata["method"] == "FACS"].to_csv(facs_only, index=False)
    files = {"metadata": (open(facs_only, "rb"), "facs_samples.csv"),
             "facs": (open(demo_data / "facs_flowjo_export.xlsx", "rb"), "flowjo.xlsx")}
    response = client.post("/upload", data=files, content_type="multipart/form-data")
    assert response.status_code == 302
    run_id = response.headers["Location"].rstrip("/").split("/")[-1]
    setup_page = client.get(f"/run/{run_id}").get_data(as_text=True)
    assert "default Flex data (flex.h5ad)" in setup_page
    assert "12 scRNA-seq samples" in setup_page and "(6 samples from 4 people)" in setup_page

    client.post(f"/run/{run_id}/analyse", data={"condition_column": "condition"})
    page = client.get(f"/run/{run_id}/results").get_data(as_text=True)
    assert "12 samples of 8 people; default Flex data" in page
    assert "A4_composition_overview.png" in page


def test_missing_default_flex_data_is_reported(client, demo_data, tmp_path, monkeypatch):
    monkeypatch.setattr(analysis, "DEFAULT_FLEX_H5AD", tmp_path / "not_connected.h5ad")
    monkeypatch.setattr(analysis, "DEFAULT_FLEX_CACHE_DIR", tmp_path / "no_cache")
    files = {"metadata": (open(demo_data / "metadata.csv", "rb"), "metadata.csv"),
             "facs": (open(demo_data / "facs_flowjo_export.xlsx", "rb"), "flowjo.xlsx")}
    page = client.post("/upload", data=files, content_type="multipart/form-data",
                       follow_redirects=True).get_data(as_text=True)
    assert "default Flex data was not found" in page
    assert not any(web.RUNS_DIR.iterdir())      # the half-finished run was removed


def test_atlas_summary_table_gives_same_proportions(demo_data, tmp_path):
    # Turn the per-cell demo atlas into a summary table (cell type x donor -> n_cells)
    cells = pd.read_csv(demo_data / analysis.ATLAS_FILE)
    summary = (cells.groupby([analysis.ATLAS_DONOR_COLUMN, analysis.ATLAS_CELL_TYPE_COLUMN])
               .size().rename(analysis.ATLAS_COUNT_COLUMN).reset_index())
    summary.to_csv(tmp_path / "summary_atlas.csv", index=False)

    per_cell = analysis.build_atlas_proportions(analysis.load_atlas(demo_data / analysis.ATLAS_FILE))
    from_summary = analysis.build_atlas_proportions(analysis.load_atlas(tmp_path / "summary_atlas.csv"))
    pd.testing.assert_series_equal(per_cell["proportion"], from_summary["proportion"])

    # A summary table has no gene values, so the marker comparison (E3) is skipped
    (demo_data / analysis.ATLAS_FILE).unlink()
    (tmp_path / "summary_atlas.csv").rename(demo_data / analysis.ATLAS_FILE)
    out = tmp_path / "results"
    analysis.run_analysis(demo_data, out, "condition", SEQ_CONDITIONS, FACS_CONDITIONS)
    made = {p.name for p in out.iterdir()}
    assert "E1_atlas_proportions.png" in made and "E2_atlas_composition_difference.png" in made
    assert "E3_atlas_marker_profiles.png" not in made


def test_missing_file_is_reported(client, demo_data):
    files = {"metadata": (open(demo_data / "metadata.csv", "rb"), "metadata.csv")}
    response = client.post("/upload", data=files, content_type="multipart/form-data",
                           follow_redirects=True)
    assert "Please choose a file" in response.get_data(as_text=True)
    assert not any((web.RUNS_DIR).iterdir())     # the half-finished run was removed


def test_wrong_file_type_is_reported(client, demo_data):
    files = {
        "metadata": (open(demo_data / "metadata.csv", "rb"), "metadata.csv"),
        "facs": (open(demo_data / "metadata.csv", "rb"), "flowjo.csv"),
        "seq": (open(demo_data / "scrnaseq_cells.csv", "rb"), "cells.csv"),
    }
    response = client.post("/upload", data=files, content_type="multipart/form-data",
                           follow_redirects=True)
    assert "should be a .xlsx file" in response.get_data(as_text=True)


def test_bad_metadata_shows_script_error(client, tmp_path):
    bad = tmp_path / "bad.csv"
    pd.DataFrame({"name": ["a"], "condition": ["x"]}).to_csv(bad, index=False)
    files = {"metadata": (open(bad, "rb"), "bad.csv"),
             "facs": (io.BytesIO(b"x"), "f.xlsx"), "seq": (open(bad, "rb"), "s.csv")}
    page = client.post("/upload", data=files, content_type="multipart/form-data",
                       follow_redirects=True).get_data(as_text=True)
    assert "Column not found" in page and "sample_id" in page


def test_wrong_gate_name_shows_error_instead_of_crashing(client, monkeypatch):
    run_id = client.post("/upload", data={"demo": "1"}).headers["Location"].split("/")[-1]
    monkeypatch.setitem(analysis.POPULATIONS, "Immune",
                        {"gate": "Wrong/Gate", "seq_labels": ["T cell"]})
    response = client.post(f"/run/{run_id}/analyse", data={"condition_column": "condition"})
    assert response.status_code == 400
    assert "Wrong/Gate | Count" in response.get_data(as_text=True)


def test_two_marker_panel_end_to_end(demo_data, tmp_path):
    """With CD45 + EpCAM only: other populations and markers, same outputs."""
    out = tmp_path / "results"
    result = analysis.run_analysis(demo_data, out, "condition", SEQ_CONDITIONS, FACS_CONDITIONS,
                                   panel="CD45_EpCAM")
    assert sorted(p.name for p in out.iterdir()) == sorted(EXPECTED_FILES)
    assert result["panel"] == "CD45_EpCAM"
    assert list(dict.fromkeys(result["comparison"]["population"])) == [
        "Neutrophils", "Other immune", "EpCAM-low (proximal tubule)", "EpCAM-high epithelium"]
    assert list(result["rna_stats"]["marker"]) == ["CD45", "EpCAM", "All markers (scaled)"]
    # the pool, the Flex groups and the atlas side by side
    assert set(result["overview"]["label"]) == {"Healthy pool (FACS)", "Disease (FACS)", "Healthy (scRNA-seq)",
                                                "Disease (scRNA-seq)", "Reference atlas"}


def test_web_panel_choice(client):
    run_id = client.post("/upload", data={"demo": "1"}).headers["Location"].split("/")[-1]
    setup_page = client.get(f"/run/{run_id}").get_data(as_text=True)
    assert "CD45 + EpCAM (+ Fixable Yellow)" in setup_page and "CD45 + EpCAM + CD10 + CD13" in setup_page
    client.post(f"/run/{run_id}/analyse", data={"condition_column": "condition", "panel": "CD45_EpCAM"})
    page = client.get(f"/run/{run_id}/results").get_data(as_text=True)
    assert "FACS panel <b>CD45 + EpCAM (+ Fixable Yellow)</b>" in page
    assert "EpCAM-low (proximal tubule)" in page and "Composition overview" in page
    # re-running keeps the chosen panel selected
    setup_again = " ".join(client.get(f"/run/{run_id}").get_data(as_text=True).split())
    assert 'value="CD45_EpCAM" checked' in setup_again
    assert client.post(f"/run/{run_id}/analyse", data={"condition_column": "condition",
                                                       "panel": "nope"}).status_code == 400


def test_old_sample_sheet_without_method_column_still_runs(demo_data, tmp_path):
    """An older sheet (same samples measured by both methods, no method column) still works."""
    metadata = pd.read_csv(demo_data / analysis.METADATA_FILE, dtype=str)
    facs = metadata[metadata["method"] == "FACS"]
    seq_ids = metadata.loc[metadata["method"] == "scRNA-seq", analysis.SAMPLE_COLUMN].tolist()
    # Pretend the first 7 scRNA-seq samples were the FACS samples: rename the cells to match
    cells = pd.read_csv(demo_data / analysis.SEQ_FILE, dtype={analysis.SAMPLE_COLUMN: str})
    cells[analysis.SAMPLE_COLUMN] = cells[analysis.SAMPLE_COLUMN].map(
        dict(zip(seq_ids, facs[analysis.SAMPLE_COLUMN])))
    cells.dropna(subset=[analysis.SAMPLE_COLUMN]).to_csv(demo_data / analysis.SEQ_FILE, index=False)
    facs.drop(columns="method").to_csv(demo_data / analysis.METADATA_FILE, index=False)

    result = analysis.run_analysis(demo_data, tmp_path / "results", "condition", FACS_CONDITIONS, FACS_CONDITIONS)
    assert (result["n_facs_samples"], result["n_seq_samples"]) == (7, 7)


def run_id_of(response):
    """The run id from the redirect after an upload."""
    return response.headers["Location"].rstrip("/").split("/")[-1]


# ---- populations scRNA-seq cannot see, and what the results can show -------------

def test_population_without_scrnaseq_cells_is_left_out_of_the_method_comparison(demo_data, tmp_path):
    # Without the proximal tubule label, scRNA-seq has no proximal tubule cells at all
    out = tmp_path / "results"
    result = analysis.run_analysis(demo_data, out, "condition", SEQ_CONDITIONS, FACS_CONDITIONS,
                                   label_map={"Proximal tubule cell": "Other"})
    assert result["seq_missing_populations"] == ["Proximal tubule"]
    comparison = result["comparison"]
    left_out = comparison[comparison["population"] == "Proximal tubule"]
    assert len(left_out) == len(SEQ_CONDITIONS) * len(FACS_CONDITIONS)
    assert left_out["test"].str.startswith("not compared").all() and left_out["p_value"].isna().all()
    # Both methods' proportions are recalculated out of the populations both measure: they add up to 1
    compared = pd.read_csv(out / "A_proportions_compared_populations.csv")
    assert "Proximal tubule" not in set(compared["population"])
    sums = compared.groupby([analysis.SAMPLE_COLUMN, "method"])["proportion"].sum()
    assert sums.round(9).eq(1).all()
    # ... the atlas too (section E), while section B still has the FACS proximal tubule
    assert "Proximal tubule" not in set(result["atlas_tests"]["population"])
    assert "Proximal tubule" in set(result["by_condition"]["population"])
    assert "Proximal tubule" in " ".join(result["notes"])
    assert "Proximal tubule" in (out / "notes.txt").read_text()
    # The label table says which label was moved
    labels = result["labels"].set_index("label")
    assert labels.loc["Proximal tubule cell", "population"] == "Other" and labels.loc["Proximal tubule cell", "changed"]


def test_notes_say_when_a_comparison_can_never_be_significant(demo_data, tmp_path):
    result = analysis.run_analysis(demo_data, tmp_path / "results", "condition", SEQ_CONDITIONS, FACS_CONDITIONS)
    notes = " ".join(result["notes"])
    # 4 people against one pool: the smallest possible p is 2 / 2^4 = 0.125
    assert "'Healthy pool' is a single sample" in notes
    assert "smallest possible p-value is 0.125" in notes
    pool = result["comparison"][result["comparison"]["facs_group"] == "Healthy pool"]
    assert pool["smallest_possible_p"].eq(0.125).all()


def test_results_page_shows_the_notes_and_the_label_table(client):
    run_id = run_id_of(client.post("/upload", data={"demo": "1"}))
    client.post(f"/run/{run_id}/analyse", data={"condition_column": "condition"})
    page = client.get(f"/run/{run_id}/results").get_data(as_text=True)
    assert "Read this first" in page and "single sample" in page
    assert "Which population each scRNA-seq label was counted in" in page
    assert "smallest_possible_p" in page


# ---- checking and changing the label -> population assignment on the setup page ----

def test_setup_page_lists_every_scrnaseq_label_with_its_population(client, demo_data, tmp_path):
    cells = pd.read_csv(demo_data / analysis.SEQ_FILE, dtype={analysis.SAMPLE_COLUMN: str})
    cells.loc[cells.index[:30], analysis.CELL_TYPE_COLUMN] = "13: mystery cluster"     # a label not in SETTINGS
    cells.to_csv(tmp_path / "cells.csv", index=False)
    files = {"metadata": (open(demo_data / "metadata.csv", "rb"), "metadata.csv"),
             "facs": (open(demo_data / "facs_flowjo_export.xlsx", "rb"), "flowjo.xlsx"),
             "seq": (open(tmp_path / "cells.csv", "rb"), "cells.csv")}
    run_id = run_id_of(client.post("/upload", data=files, content_type="multipart/form-data"))
    page = client.get(f"/run/{run_id}").get_data(as_text=True)
    assert "Which FACS population is each scRNA-seq cell type?" in page
    for label in cells[analysis.CELL_TYPE_COLUMN].unique():
        assert label in page
    assert 'title="This label is not in the seq_labels of the settings">new' in page
    assert "<i>PTPRC</i>" in page and "<i>EPCAM</i>" in page
    # the default choice of the current panel comes from SETTINGS
    summary = analysis.seq_label_summary(tmp_path / "cells.csv")
    row = int(summary.index[summary["label"] == "Neutrophil"][0])
    compact = " ".join(page.split())
    assert f'name="label_{row}" data-panel="{analysis.FACS_PANEL}" data-default="Neutrophils"' in compact


def test_changed_label_assignment_is_used_and_remembered(client):
    run_id = run_id_of(client.post("/upload", data={"demo": "1"}))
    folder = web.RUNS_DIR / run_id
    summary = analysis.seq_label_summary(folder / "data" / analysis.SEQ_FILE)
    row = int(summary.index[summary["label"] == "Macrophage"][0])
    # Count macrophages as neutrophils (only the menus of the chosen panel are sent)
    response = client.post(f"/run/{run_id}/analyse", data={"condition_column": "condition",
                                                            f"label_{row}": "Neutrophils"})
    assert response.status_code == 302
    labels = pd.read_csv(folder / "results" / "A_label_assignment.csv").set_index("label")
    assert labels.loc["Macrophage", "population"] == "Neutrophils" and labels.loc["Macrophage", "changed"]
    assert "1 scRNA-seq label(s) were put in another population" in client.get(
        f"/run/{run_id}/results").get_data(as_text=True)
    # the setup page shows the choice again
    compact = " ".join(client.get(f"/run/{run_id}").get_data(as_text=True).split())
    menu = compact.split(f'name="label_{row}" data-panel="{analysis.FACS_PANEL}"')[1].split("</select>")[0]
    assert '<option value="Neutrophils" selected>' in menu
    # a population that is not in the panel is refused
    response = client.post(f"/run/{run_id}/analyse", data={"condition_column": "condition",
                                                            f"label_{row}": "Hepatocytes"})
    assert response.status_code == 400 and "not a population of this panel" in response.get_data(as_text=True)


# ---- scRNA-seq as an .h5ad file -------------------------------------------------------

def fake_h5ad_from_demo(demo_data, path):
    """The demo cells as a small .h5ad (like a Flex notebook's output), plus the FACS-only sample sheet."""
    from test_business_logic import write_fake_h5ad
    cells = pd.read_csv(demo_data / analysis.SEQ_FILE, dtype={analysis.SAMPLE_COLUMN: str})
    metadata = pd.read_csv(demo_data / analysis.METADATA_FILE, dtype=str)
    seq_rows = metadata[metadata["method"] == "scRNA-seq"].set_index(analysis.SAMPLE_COLUMN)
    cells = cells.assign(leiden_atlas=cells[analysis.CELL_TYPE_COLUMN],
                         patient_id=cells[analysis.SAMPLE_COLUMN].map(seq_rows["donor_id"]),
                         role=cells[analysis.SAMPLE_COLUMN].map(seq_rows["condition"]), timepoint="T00")
    write_fake_h5ad(path, cells)
    facs_only = path.parent / "facs_samples.csv"
    metadata[metadata["method"] == "FACS"].to_csv(facs_only, index=False)
    return path, facs_only


def test_upload_h5ad_choose_its_columns_and_run(client, demo_data, tmp_path):
    h5ad, facs_only = fake_h5ad_from_demo(demo_data, tmp_path / "my_flex.h5ad")
    files = {"metadata": (open(facs_only, "rb"), "facs_samples.csv"),
             "facs": (open(demo_data / "facs_flowjo_export.xlsx", "rb"), "flowjo.xlsx"),
             "seq": (open(h5ad, "rb"), "my_flex.h5ad")}
    response = client.post("/upload", data=files, content_type="multipart/form-data")
    assert response.headers["Location"].endswith("/h5ad")
    run_id = response.headers["Location"].split("/")[-2]
    # the setup page is not reachable before the columns are chosen
    assert client.get(f"/run/{run_id}").headers["Location"].endswith("/h5ad")

    page = client.get(f"/run/{run_id}/h5ad").get_data(as_text=True)
    assert "Which columns of my_flex.h5ad?" in page
    compact = " ".join(page.split())
    assert '<option value="sample_id" selected>' in compact          # guessed from the column names
    assert '<option value="patient_id" selected>' in compact and '<option value="role" selected>' in compact
    assert '<option value="log1p_cp10k" selected>' in compact

    response = client.post(f"/run/{run_id}/h5ad", data={
        "sample": "sample_id", "label": "leiden_atlas", "donor": "patient_id", "condition": "role",
        "extra": ["timepoint"], "layer": "log1p_cp10k"})
    assert response.status_code == 302
    folder = web.RUNS_DIR / run_id
    assert not (folder / "data" / analysis.SEQ_H5AD_FILE).exists()     # the big file is deleted again
    setup_page = client.get(f"/run/{run_id}").get_data(as_text=True)
    assert "my_flex.h5ad (cell types from column &#39;leiden_atlas&#39;" in setup_page
    assert "12 scRNA-seq samples" in setup_page and "(6 samples from 4 people)" in setup_page
    assert "<b>timepoint</b>" in setup_page

    client.post(f"/run/{run_id}/analyse", data={"condition_column": "condition"})
    page = client.get(f"/run/{run_id}/results").get_data(as_text=True)
    assert "12 samples of 8 people" in page and "A4_composition_overview.png" in page


def test_h5ad_column_that_describes_cells_is_refused(client, demo_data, tmp_path):
    h5ad, facs_only = fake_h5ad_from_demo(demo_data, tmp_path / "my_flex.h5ad")
    files = {"metadata": (open(facs_only, "rb"), "facs_samples.csv"),
             "facs": (open(demo_data / "facs_flowjo_export.xlsx", "rb"), "flowjo.xlsx")}
    # given as a path on this computer instead of an upload
    response = client.post("/upload", data={**files, "seq_path": f'"{h5ad}"'}, content_type="multipart/form-data")
    run_id = response.headers["Location"].split("/")[-2]
    # the cluster label differs between the cells of a sample: it cannot be the group of a sample
    response = client.post(f"/run/{run_id}/h5ad", data={"sample": "sample_id", "label": "leiden_atlas",
                                                         "condition": "leiden_atlas", "layer": "log1p_cp10k"})
    assert response.status_code == 400
    assert "has more than one value inside a sample" in response.get_data(as_text=True)
    assert h5ad.exists()                       # a file given as a path is never deleted


def test_scrnaseq_path_that_does_not_exist_is_reported(client, demo_data, tmp_path):
    files = {"metadata": (open(demo_data / "metadata.csv", "rb"), "metadata.csv"),
             "facs": (open(demo_data / "facs_flowjo_export.xlsx", "rb"), "flowjo.xlsx")}
    page = client.post("/upload", data={**files, "seq_path": str(tmp_path / "nothing.h5ad")},
                       content_type="multipart/form-data", follow_redirects=True).get_data(as_text=True)
    assert "was not found on this computer" in page
    assert not any(web.RUNS_DIR.iterdir())


def test_unknown_or_unsafe_run_ids_give_404(client):
    assert client.get("/run/123/results").status_code == 404
    assert client.get("/run/" + "0" * 32 + "/results").status_code == 404
    assert client.get("/run/../../etc/results").status_code == 404
