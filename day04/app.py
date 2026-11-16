"""
app.py
======

Web interface for facs_vs_scrnaseq.py.

All the "business logic" (loading the files, building the tables, the
statistics and the plots) stays in facs_vs_scrnaseq.py. This file only:
  1. receives the three input files (+ the optional reference atlas and per-cell
     imaging export) through an
     upload form (or makes demo data)
  2. lets you pick the condition column and conditions with checkboxes
     (instead of typing numbers in the terminal)
  3. calls run_analysis() from facs_vs_scrnaseq.py
  4. shows the summary tables and plots in the browser and offers the CSVs
     and PNGs for download

-------------------------------------------------------------------------------
HOW TO RUN
-------------------------------------------------------------------------------
      python app.py
  then open http://127.0.0.1:5000 in your browser. Stop the server with Ctrl+C.

  Every analysis gets its own folder in runs/<run id>/ with:
      data/      the uploaded (or demo) input files
      results/   the tables and plots written by facs_vs_scrnaseq.py
      run.json   which conditions were chosen + the messages the analysis printed

  The SETTINGS (gate paths, channel names, cell-type labels, ...) are still
  edited in facs_vs_scrnaseq.py - the web page shows the current values.
"""

import contextlib
import io
import json
import os
import re
import shutil
import threading
import uuid
import zipfile
from pathlib import Path

import pandas as pd
from flask import (Flask, abort, flash, redirect, render_template, request,
                   send_file, send_from_directory, url_for)

import explanations   # help texts for the results page (how to read each table and plot)
import facs_vs_scrnaseq as analysis


# =============================================================================
# SETTINGS for the web app
# =============================================================================

APP_DIR = Path(__file__).resolve().parent
RUNS_DIR = APP_DIR / "runs"           # one sub-folder per analysis
MAX_UPLOAD_MB = 500                   # largest total upload accepted

# The upload fields (three required + two optional): form field name -> (file name the analysis expects, allowed extension)
UPLOADS = {
    "metadata": (analysis.METADATA_FILE, ".csv"),
    "facs": (analysis.FACS_FILE, ".xlsx"),
    "seq": (analysis.SEQ_FILE, ".csv"),
    "atlas": (analysis.ATLAS_FILE, ".csv"),
    "facs_cells": (analysis.FACS_CELLS_FILE, ".csv"),
}
# Upload fields that may be left empty: without an atlas, section E is skipped;
# without the per-cell imaging export, the cell-level plots of section F are skipped
OPTIONAL_UPLOADS = {"atlas", "facs_cells"}

# Headings for the plots and tables, by the letter their file name starts with
SECTIONS = {
    "A": "A. Proportions: FACS vs scRNA-seq",
    "B": "B. Proportions between conditions",
    "C": "C. Protein vs RNA",
    "D": "D. S8 imaging features",
    "E": "E. Comparison with the reference atlas",
    "F": "F. Imaging profiles per cell type",
}

app = Flask(__name__)
app.config["MAX_CONTENT_LENGTH"] = MAX_UPLOAD_MB * 1024 * 1024
app.secret_key = os.environ.get("SECRET_KEY") or os.urandom(16).hex()   # needed for flash messages

# matplotlib and the module-level OUTPUT_DIR in the analysis script are not
# safe to use from two requests at once, so analyses run one at a time.
ANALYSIS_LOCK = threading.Lock()
RUN_ID_PATTERN = re.compile(r"^[0-9a-f]{32}$")


# =============================================================================
# HELPERS
# =============================================================================

def run_folder(run_id):
    """Return the folder of an existing run; 404 for anything that is not one."""
    if not RUN_ID_PATTERN.match(run_id):
        abort(404)
    folder = RUNS_DIR / run_id
    if not folder.is_dir():
        abort(404)
    return folder


def call_analysis(function, *args, **kwargs):
    """
    Call a function from facs_vs_scrnaseq.py and capture what it prints.
    The script stops with SystemExit("explanation") when a file or column is
    wrong; here that becomes an error message for the page instead of a crash.
    Returns (result, printed_text, error_message).
    """
    log = io.StringIO()
    with ANALYSIS_LOCK, contextlib.redirect_stdout(log):
        try:
            return function(*args, **kwargs), log.getvalue(), None
        except SystemExit as stop:
            return None, log.getvalue(), str(stop.code).strip()


def condition_options(metadata):
    """For every metadata column that could hold conditions: its values and sample counts."""
    options = []
    for column in metadata.columns:
        if column in (analysis.SAMPLE_COLUMN, analysis.FACS_NAME_COLUMN):
            continue
        counts = metadata[column].value_counts(sort=False)
        options.append({"column": column, "values": [(v, int(n)) for v, n in counts.items()]})
    return options


def table_html(table, columns=None):
    """Turn a results table into an HTML table (numbers rounded to 3 significant digits)."""
    if table is None or table.empty:
        return None
    if columns:
        table = table[[c for c in columns if c in table.columns]]
    return table.to_html(index=False, classes="data", border=0, na_rep="–",
                         float_format=lambda x: f"{x:.3g}")


def read_result(results_dir, name):
    """Read one of the CSV tables the analysis wrote (None if missing or empty)."""
    path = results_dir / name
    if not path.exists() or path.stat().st_size == 0:
        return None
    try:
        return pd.read_csv(path)
    except pd.errors.EmptyDataError:
        return None


# =============================================================================
# PAGES
# =============================================================================

@app.get("/")
def index():
    """Start page: upload form, demo button and the current SETTINGS."""
    return render_template("index.html", uploads=UPLOADS, analysis=analysis)


@app.post("/upload")
def upload():
    """Save the uploaded files (or make demo data) into a new run folder."""
    run_id = uuid.uuid4().hex
    folder = RUNS_DIR / run_id
    data_dir = folder / "data"
    data_dir.mkdir(parents=True)

    if "demo" in request.form:
        call_analysis(analysis.make_demo_data, data_dir)
    else:
        for field, (file_name, extension) in UPLOADS.items():
            file = request.files.get(field)
            if not file or not file.filename:
                if field in OPTIONAL_UPLOADS:
                    continue          # optional file not given: its section is skipped
                shutil.rmtree(folder)
                flash(f"Please choose a file for '{file_name}'.")
                return redirect(url_for("index"))
            if not file.filename.lower().endswith(extension):
                shutil.rmtree(folder)
                flash(f"'{file.filename}' should be a {extension} file.")
                return redirect(url_for("index"))
            file.save(data_dir / file_name)    # saved under the name the analysis expects

    _, _, error = call_analysis(analysis.load_metadata, data_dir / analysis.METADATA_FILE)
    if error:
        shutil.rmtree(folder)
        flash(error)
        return redirect(url_for("index"))
    return redirect(url_for("setup", run_id=run_id))


@app.get("/run/<run_id>")
def setup(run_id, error=None):
    """Second page: choose the condition column and the conditions to include."""
    folder = run_folder(run_id)
    metadata, _, load_error = call_analysis(analysis.load_metadata,
                                            folder / "data" / analysis.METADATA_FILE)
    if load_error:
        abort(404)
    options = condition_options(metadata)
    return render_template("setup.html", run_id=run_id, options=options,
                           n_samples=len(metadata), error=error)


@app.post("/run/<run_id>/analyse")
def analyse(run_id):
    """Run the analysis with the chosen conditions, then show the results page."""
    folder = run_folder(run_id)
    column = request.form.get("condition_column")
    conditions = request.form.getlist("conditions")

    metadata = pd.read_csv(folder / "data" / analysis.METADATA_FILE, dtype=str)
    if column not in [option["column"] for option in condition_options(metadata)]:
        return setup(run_id, error="Please choose a condition column."), 400
    if not conditions:   # nothing ticked = all conditions (like pressing Enter in the terminal)
        conditions = list(metadata[column].value_counts(sort=False).index)

    results_dir = folder / "results"
    shutil.rmtree(results_dir, ignore_errors=True)    # re-running replaces old results
    result, log, error = call_analysis(analysis.run_analysis, folder / "data", results_dir,
                                       column, conditions)
    if error:
        return setup(run_id, error=error), 400

    (folder / "run.json").write_text(json.dumps({
        "condition_column": result["condition_column"],
        "conditions": result["conditions"],
        "n_facs_samples": int(result["n_facs_samples"]),
        "n_cells": int(result["n_cells"]),
        "n_seq_samples": int(result["n_seq_samples"]),
        "log": log,
    }, indent=2))
    return redirect(url_for("results", run_id=run_id))


@app.get("/run/<run_id>/results")
def results(run_id):
    """Results page: summary tables, all plots, and download links."""
    folder = run_folder(run_id)
    results_dir = folder / "results"
    if not (folder / "run.json").exists():
        return redirect(url_for("setup", run_id=run_id))
    info = json.loads((folder / "run.json").read_text())

    by_condition = read_result(results_dir, "B_condition_tests_proportions.csv")
    rna_stats = read_result(results_dir, "C_protein_vs_rna_stats.csv")
    # For each section: the summary table and which of its columns to show
    tables = {
        "A": (read_result(results_dir, "A_method_agreement_stats.csv"),
              ["population", "n", "spearman_r", "spearman_q_fdr", "pearson_r",
               "bland_altman_bias", "bland_altman_lower", "bland_altman_upper",
               "wilcoxon_q_fdr"]),
        "B": (by_condition,
              ["population", "method", "test", "comparison", "p_value",
               "q_value_fdr", "significant"]),
        "C": (rna_stats[rna_stats["population"] == "All populations"] if rna_stats is not None else None,
              ["marker", "n", "spearman_r", "spearman_q_fdr", "pearson_r", "pearson_q_fdr"]),
        "D": (read_result(results_dir, "D_imaging_tests_between_populations.csv"),
              ["feature", "test", "comparison", "p_value", "q_value_fdr", "significant"]),
        "E": (read_result(results_dir, "E_atlas_proportion_tests.csv"),
              ["population", "method", "group", "n_samples", "median_samples", "atlas_proportion",
               "difference", "n_atlas_donors", "test", "q_value_fdr", "significant"]),
        "F": (read_result(results_dir, "F_imaging_population_summary.csv"),
              ["population", "feature", "n_samples", "mean", "sd", "ci95_low", "ci95_high",
               "median", "cv_percent", "n_outlier_samples"]),
    }
    n_significant = int(by_condition["significant"].sum()) if by_condition is not None else 0

    files = sorted(p.name for p in results_dir.iterdir())
    sections = []
    for letter, title in SECTIONS.items():
        table, columns = tables[letter]
        shown = [] if table is None or table.empty else [c for c in columns if c in table.columns]
        sections.append({
            "letter": letter, "title": title,
            "table": table_html(table, columns),
            # explanations: the question of the section, notes, and the meaning of each column shown
            "help": explanations.SECTION_HELP.get(letter),
            "column_help": [(c, explanations.COLUMN_HELP[c]) for c in shown if c in explanations.COLUMN_HELP],
            # each plot comes with its own "how to read it" text (None if there is none)
            "plots": [(f, explanations.PLOT_HELP.get(f))
                      for f in files if f.startswith(letter) and f.endswith(".png")],
            "csvs": [f for f in files if f.startswith(letter) and f.endswith(".csv")],
        })
    return render_template("results.html", run_id=run_id, info=info, sections=sections,
                           n_significant=n_significant, alpha=analysis.ALPHA,
                           general_help=explanations.GENERAL_HELP)


@app.get("/run/<run_id>/file/<path:name>")
def result_file(run_id, name):
    """Serve one plot or table from the results folder."""
    return send_from_directory(run_folder(run_id) / "results", name)


@app.get("/run/<run_id>/download")
def download_all(run_id):
    """All tables and plots of a run in one zip file."""
    results_dir = run_folder(run_id) / "results"
    if not results_dir.is_dir():
        abort(404)
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(results_dir.iterdir()):
            archive.write(path, arcname=path.name)
    buffer.seek(0)
    return send_file(buffer, mimetype="application/zip", as_attachment=True,
                     download_name=f"facs_vs_scrnaseq_results_{run_id[:8]}.zip")


@app.errorhandler(413)
def too_large(_):
    flash(f"The upload is larger than {MAX_UPLOAD_MB} MB.")
    return redirect(url_for("index"))


if __name__ == "__main__":
    # debug=False: the debugger would let anyone who can reach the page run code.
    # host 127.0.0.1: the app is only reachable from this computer.
    app.run(host="127.0.0.1", port=5000, debug=False)
