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
    "A1_proportions_stacked.png", "A2_proportions_FACS_vs_seq.png", "A3_bland_altman.png",
    "A_method_agreement_stats.csv", "A_proportions.csv",
    "B1_proportions_by_condition.png", "B_condition_tests_proportions.csv",
    "C1_protein_vs_rna.png", "C_protein_vs_rna.csv", "C_protein_vs_rna_stats.csv",
    "D1_imaging_features.png", "D_imaging_features.csv",
    "D_imaging_tests_between_conditions.csv", "D_imaging_tests_between_populations.csv",
]


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
    result = analysis.run_analysis(demo_data, out, "condition", ["Healthy", "Disease"])
    assert sorted(p.name for p in out.iterdir()) == sorted(EXPECTED_FILES)
    assert result["n_facs_samples"] == 10
    assert list(result["agreement"]["population"]) == list(analysis.POPULATIONS)
    assert not result["by_condition"].empty


def test_run_analysis_one_condition_skips_condition_plot(demo_data, tmp_path):
    out = tmp_path / "results"
    result = analysis.run_analysis(demo_data, out, "condition", ["Healthy"])
    assert not (out / "B1_proportions_by_condition.png").exists()
    assert result["by_condition"].empty
    assert result["n_facs_samples"] == 5


def test_run_analysis_unknown_condition_stops_with_message(demo_data, tmp_path):
    with pytest.raises(SystemExit, match="not found"):
        analysis.run_analysis(demo_data, tmp_path / "results", "condition", ["Nope"])


def test_command_line_still_works(tmp_path):
    """Run the script from the terminal the old way, on a copy so day04 stays clean."""
    script = tmp_path / "facs_vs_scrnaseq.py"
    script.write_text((HERE / "facs_vs_scrnaseq.py").read_text())
    done = subprocess.run([sys.executable, str(script), "--demo", "--condition-column", "condition",
                           "--conditions", "Healthy", "Disease"],
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
    assert "Healthy" in setup_page and "(5 samples)" in setup_page

    response = client.post(f"/run/{run_id}/analyse",
                           data={"condition_column": "condition", "conditions": ["Healthy", "Disease"]})
    assert response.status_code == 302

    page = client.get(f"/run/{run_id}/results").get_data(as_text=True)
    for name in EXPECTED_FILES:
        assert name in page
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

    # no condition ticked = all conditions, like pressing Enter in the terminal
    client.post(f"/run/{run_id}/analyse", data={"condition_column": "condition"})
    page = client.get(f"/run/{run_id}/results").get_data(as_text=True)
    assert "Healthy, Disease" in page


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


def test_unknown_or_unsafe_run_ids_give_404(client):
    assert client.get("/run/123/results").status_code == 404
    assert client.get("/run/" + "0" * 32 + "/results").status_code == 404
    assert client.get("/run/../../etc/results").status_code == 404
