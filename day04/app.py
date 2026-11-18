"""
app.py
======

Web interface for facs_vs_scrnaseq.py.

All the "business logic" (loading the files, building the tables, the
statistics and the plots) stays in facs_vs_scrnaseq.py. This file only:
  1. receives the input files through an upload form (or makes demo data):
     the sample sheet and FlowJo export, plus the optional scRNA-seq data
     (a cells .csv, an .h5ad file - uploaded or as a path on this computer -
     or left empty = the default Flex data from DEFAULT_FLEX_H5AD), reference
     atlas and per-cell imaging export
  1b. for an .h5ad file: asks which obs columns hold the sample, the cell type,
     the person and the groups, and which layer to use
  2. lets you pick the FACS panel (which antibodies were used), the condition
     column, and the conditions to include for the
     scRNA-seq samples and for the FACS samples, with checkboxes (instead of
     typing numbers in the terminal). The two methods can have different
     samples and conditions, e.g. scRNA-seq 'Disease' vs a FACS 'Healthy pool'.
     It also shows every scRNA-seq label (cells, marker genes) with the FACS
     population it is counted in, and lets you change that per label.
  3. calls run_analysis() from facs_vs_scrnaseq.py
  4. shows the warnings (notes.txt), the summary tables and plots in the
     browser and offers the CSVs and PNGs for download

-------------------------------------------------------------------------------
HOW TO RUN
-------------------------------------------------------------------------------
      python app.py
  then open http://127.0.0.1:5000 in your browser. Stop the server with Ctrl+C.

  Every analysis gets its own folder in runs/<run id>/ with:
      data/      the uploaded (or demo) input files (an uploaded .h5ad is deleted
                 once its cells are read: only the small cells table is kept)
      results/   the tables and plots written by facs_vs_scrnaseq.py
      run.json   which conditions and label assignments were chosen + the messages the analysis printed

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
MAX_UPLOAD_MB = 20000                 # largest total upload accepted (.h5ad files can be several GB;
                                      # a big .h5ad is quicker to give as a path, see below)

# The upload fields (two required + three optional):
# form field name -> (file name the analysis expects, allowed extensions)
UPLOADS = {
    "metadata": (analysis.METADATA_FILE, (".csv",)),
    "facs": (analysis.FACS_FILE, (".xlsx",)),
    "seq": (analysis.SEQ_FILE, (".csv", ".h5ad")),
    "atlas": (analysis.ATLAS_FILE, (".csv",)),
    "facs_cells": (analysis.FACS_CELLS_FILE, (".csv",)),
}
# Upload fields that may be left empty: without a scRNA-seq file, the default Flex
# data (DEFAULT_FLEX_H5AD in the analysis SETTINGS) is used; without an atlas,
# section E is skipped; without the per-cell imaging export, the cell-level plots
# of section F are skipped
OPTIONAL_UPLOADS = {"seq", "atlas", "facs_cells"}
SEQ_SOURCE_FILE = "seq_source.txt"    # in the run folder: where the scRNA-seq data came from
H5AD_PATH_FILE = "h5ad_path.txt"      # in the run folder: an .h5ad given as a path, waiting to be read
H5AD_NAME_FILE = "h5ad_name.txt"      # in the run folder: the original name of an uploaded .h5ad

# Headings for the plots and tables, by the letter their file name starts with
SECTIONS = {
    "A": "A. Proportions: scRNA-seq vs FACS",
    "B": "B. Proportions between conditions",
    "C": "C. Protein vs RNA per population",
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


# The two methods on the setup page: (method, form field name, heading)
METHOD_FIELDS = [("scRNA-seq", "seq_conditions", "scRNA-seq samples"),
                 ("FACS", "facs_conditions", "FACS samples")]


def condition_options(metadata):
    """
    For every metadata column that could hold conditions: for each method, its
    values and how many samples of that method have each value.
    metadata = the sample sheet as returned by analysis.load_metadata().
    """
    options = []
    for column in analysis.condition_candidates(metadata):
        methods = []
        for method, field, heading in METHOD_FIELDS:
            samples = analysis.samples_of(metadata, method)
            counts = samples[column].value_counts(sort=False)
            # How many different people those samples come from (donor_id; none = one person per sample)
            people = analysis.donor_of_each_row(samples, metadata).groupby(samples[column]).nunique()
            methods.append({"method": method, "field": field, "heading": heading,
                            "values": [(v, int(n), int(people[v])) for v, n in counts.items()]})
        # A column can only be used if BOTH methods have values in it (otherwise there is nothing to compare)
        empty = [entry["method"] for entry in methods if not entry["values"]]
        options.append({"column": column, "methods": methods, "usable": not empty,
                        "missing_for": empty,
                        "help": analysis.CONDITION_COLUMN_HELP.get(column, "A column of your sample sheet.")})
    # Show the usable columns first, and "condition" first among them
    options.sort(key=lambda o: (not o["usable"], o["column"] != "condition"))
    return options


def last_panel(folder):
    """The FACS panel this run used last time (when re-running), else the default from SETTINGS."""
    return previous_run(folder).get("panel", analysis.FACS_PANEL)


def seq_source(folder):
    """Where the scRNA-seq data of a run came from (uploaded file, default Flex data, or demo)."""
    path = folder / SEQ_SOURCE_FILE
    return path.read_text(encoding="utf-8") if path.exists() else "demo data"


def pending_h5ad(folder):
    """
    The .h5ad file of a run whose columns still have to be chosen (h5ad page),
    or None. It is either uploaded into data/ or given as a path on this computer.
    """
    if (folder / "data" / analysis.SEQ_FILE).exists():
        return None
    uploaded = folder / "data" / analysis.SEQ_H5AD_FILE
    if uploaded.exists():
        return uploaded
    path_file = folder / H5AD_PATH_FILE
    return Path(path_file.read_text(encoding="utf-8")) if path_file.exists() else None


def previous_run(folder):
    """The run.json of the last analysis of this run folder ({} before the first one)."""
    try:
        return json.loads((folder / "run.json").read_text())
    except (OSError, ValueError):
        return {}


def label_rows(summary, label_maps):
    """
    The label table of the setup page: one row per scRNA-seq label with its cells,
    its marker genes (value + shade 0-1 relative to the highest label) and, per
    panel, the populations to choose from and the one selected (the choice of the
    last run, else SETTINGS).
    """
    genes = [(name, info["gene"]) for name, info in analysis.ALL_MARKERS.items() if info["gene"] in summary.columns]
    highest = {gene: summary[gene].max() for _, gene in genes}
    rows = []
    for i, row in summary.iterrows():
        choices = {}
        for panel, info in analysis.PANELS.items():
            default = row[f"population_{panel}"]
            chosen = label_maps.get(panel, {}).get(row["label"], default)
            choices[panel] = {"options": [*info["populations"], "Other"], "default": default, "chosen": chosen}
        rows.append({
            "index": i, "label": row["label"], "n_cells": int(row["n_cells"]), "n_samples": int(row["n_samples"]),
            "in_settings": bool(row["in_settings"]), "choices": choices,
            "genes": [(row[gene], row[gene] / highest[gene] if highest[gene] > 0 else 0) for _, gene in genes],
        })
    return genes, rows


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

    def give_up(message):
        shutil.rmtree(folder)
        flash(message)
        return redirect(url_for("index"))

    if "demo" in request.form:
        call_analysis(analysis.make_demo_data, data_dir)
    else:
        for field, (file_name, extensions) in UPLOADS.items():
            file = request.files.get(field)
            if not file or not file.filename:
                if field in OPTIONAL_UPLOADS:
                    continue          # optional file not given: its section is skipped
                return give_up(f"Please choose a file for '{file_name}'.")
            if not file.filename.lower().endswith(extensions):
                return give_up(f"'{file.filename}' should be a {' or '.join(extensions)} file.")
            if field == "seq" and file.filename.lower().endswith(".h5ad"):
                file.save(data_dir / analysis.SEQ_H5AD_FILE)   # its columns are chosen on the next page
                (folder / H5AD_NAME_FILE).write_text(Path(file.filename).name, encoding="utf-8")
            else:
                file.save(data_dir / file_name)    # saved under the name the analysis expects

        # The scRNA-seq data can also be given as a path on this computer (quicker for a big .h5ad)
        seq_path = request.form.get("seq_path", "").strip().strip('"')
        if seq_path and not ((data_dir / analysis.SEQ_FILE).exists() or (data_dir / analysis.SEQ_H5AD_FILE).exists()):
            path = Path(seq_path)
            if not path.is_file():
                return give_up(f"The scRNA-seq path was not found on this computer:\n  {seq_path}")
            if path.suffix.lower() == ".h5ad":
                (folder / H5AD_PATH_FILE).write_text(str(path), encoding="utf-8")
            elif path.suffix.lower() == ".csv":
                shutil.copyfile(path, data_dir / analysis.SEQ_FILE)
                (folder / SEQ_SOURCE_FILE).write_text(f"{path.name} (from a path)", encoding="utf-8")
            else:
                return give_up(f"The scRNA-seq path should end in .h5ad or .csv:\n  {seq_path}")

        if pending_h5ad(folder) is None and not (folder / SEQ_SOURCE_FILE).exists():
            # No scRNA-seq data given: use the default Flex data (and add its samples to the sheet)
            used_default, _, error = call_analysis(analysis.add_default_flex, data_dir)
            if error:
                return give_up(error)
            source = (f"default Flex data ({Path(analysis.DEFAULT_FLEX_H5AD).name})" if used_default
                      else "uploaded file")
            (folder / SEQ_SOURCE_FILE).write_text(source, encoding="utf-8")

    _, _, error = call_analysis(analysis.load_metadata, data_dir / analysis.METADATA_FILE)
    if error:
        return give_up(error)
    if pending_h5ad(folder):
        return redirect(url_for("h5ad_columns", run_id=run_id))
    return redirect(url_for("setup", run_id=run_id))


@app.route("/run/<run_id>/h5ad", methods=["GET", "POST"])
def h5ad_columns(run_id):
    """For an .h5ad file: choose which obs columns hold the sample, cell type, person and groups."""
    folder = run_folder(run_id)
    path = pending_h5ad(folder)
    if path is None:
        return redirect(url_for("setup", run_id=run_id))
    name_file = folder / H5AD_NAME_FILE
    file_name = name_file.read_text(encoding="utf-8") if name_file.exists() else path.name
    try:
        overview, _, error = call_analysis(analysis.h5ad_overview, path)
    except (OSError, KeyError) as problem:          # not an .h5ad file at all (h5py cannot open it)
        overview, error = None, f"{file_name} could not be read as an .h5ad file ({problem})."
    if error or overview is None:
        return render_template("h5ad.html", run_id=run_id, file_name=file_name, overview=None,
                               error=error), 400
    if request.method == "GET":
        return render_template("h5ad.html", run_id=run_id, file_name=file_name, overview=overview,
                               chosen=overview["guesses"], error=None)

    # POST: read the chosen columns and convert the file into the cells table + sample-sheet rows
    form = request.form
    columns = {c["name"] for c in overview["columns"]}
    chosen = {"sample": form.get("sample"), "label": form.get("label"), "donor": form.get("donor") or None,
              "condition": form.get("condition") or None, "layer": form.get("layer", "X")}
    extra = [c for c in form.getlist("extra") if c in columns and c != chosen["condition"]]
    problem = None
    if chosen["sample"] not in columns or chosen["label"] not in columns:
        problem = "Please choose the sample column and the cell-type column."
    elif any(chosen[k] and chosen[k] not in columns for k in ("donor", "condition")):
        problem = "Please choose columns of the file."
    elif chosen["layer"] not in overview["layers"]:
        problem = "Please choose a layer of the file."
    if problem:
        return render_template("h5ad.html", run_id=run_id, file_name=file_name, overview=overview,
                               chosen=chosen, error=problem), 400
    groups = {**({"condition": chosen["condition"]} if chosen["condition"] else {}), **{c: c for c in extra}}
    n_cells, _, error = call_analysis(analysis.add_h5ad, folder / "data", path, chosen["sample"], chosen["label"],
                                      chosen["donor"], groups, chosen["layer"], overview["gene_column"])
    if error:
        return render_template("h5ad.html", run_id=run_id, file_name=file_name, overview=overview,
                               chosen=chosen, error=error), 400
    # Done: the small cells table is kept; an uploaded .h5ad (big, maybe patient data) is deleted
    (folder / "data" / analysis.SEQ_H5AD_FILE).unlink(missing_ok=True)
    (folder / H5AD_PATH_FILE).unlink(missing_ok=True)
    (folder / SEQ_SOURCE_FILE).write_text(
        f"{file_name} (cell types from column '{chosen['label']}', {n_cells} cells)", encoding="utf-8")
    return redirect(url_for("setup", run_id=run_id))


@app.get("/run/<run_id>")
def setup(run_id, error=None):
    """Second page: choose the condition column and the conditions to include."""
    folder = run_folder(run_id)
    if pending_h5ad(folder):
        return redirect(url_for("h5ad_columns", run_id=run_id))
    metadata, _, load_error = call_analysis(analysis.load_metadata,
                                            folder / "data" / analysis.METADATA_FILE)
    if load_error:
        abort(404)
    options = condition_options(metadata)
    n_samples = {method: len(analysis.samples_of(metadata, method)) for method in analysis.METHODS}
    # Every scRNA-seq label with its cells, marker genes and population (to check / change)
    summary, _, label_error = call_analysis(analysis.seq_label_summary, folder / "data" / analysis.SEQ_FILE)
    genes, labels = label_rows(summary, previous_run(folder).get("label_maps", {})) if summary is not None else ([], [])
    return render_template("setup.html", run_id=run_id, options=options,
                           n_samples=n_samples, error=error,
                           panels=analysis.PANELS, current_panel=last_panel(folder),
                           seq_source=seq_source(folder), value_help=analysis.CONDITION_VALUE_HELP,
                           genes=genes, labels=labels, label_error=label_error)


def chosen_label_map(folder, panel):
    """
    The population chosen on the setup page for every scRNA-seq label (form fields
    label_0, label_1, ... in the order of seq_label_summary). Returns
    (label_map, error); labels not in the form keep their SETTINGS population.
    """
    summary, _, error = call_analysis(analysis.seq_label_summary, folder / "data" / analysis.SEQ_FILE)
    if error:
        return None, error
    allowed = [*analysis.PANELS[panel]["populations"], "Other"]
    label_map = {}
    for i, label in enumerate(summary["label"]):
        population = request.form.get(f"label_{i}")
        if population is None:
            continue
        if population not in allowed:
            return None, f"'{population}' is not a population of this panel (label '{label}')."
        label_map[label] = population
    return label_map, None


@app.post("/run/<run_id>/analyse")
def analyse(run_id):
    """Run the analysis with the chosen conditions, then show the results page."""
    folder = run_folder(run_id)
    column = request.form.get("condition_column")
    panel = request.form.get("panel", analysis.FACS_PANEL)
    if panel not in analysis.PANELS:
        return setup(run_id, error="Please choose the FACS panel."), 400

    metadata, _, load_error = call_analysis(analysis.load_metadata, folder / "data" / analysis.METADATA_FILE)
    if load_error:
        abort(404)
    option = next((o for o in condition_options(metadata) if o["column"] == column), None)
    if option is None:
        return setup(run_id, error="Please choose a condition column."), 400
    if not option["usable"]:
        return setup(run_id, error=f"The column '{column}' has no values for the "
                                   f"{' and '.join(option['missing_for'])} samples, so there is nothing to compare. "
                                   "Choose another column."), 400
    # The ticked conditions of each method; nothing ticked = all conditions of that
    # method (like pressing Enter in the terminal)
    chosen = {}
    for entry in option["methods"]:
        ticked = request.form.getlist(entry["field"])
        chosen[entry["method"]] = ticked or [value for value, *_ in entry["values"]]

    label_map, error = chosen_label_map(folder, panel)
    if error:
        return setup(run_id, error=error), 400

    results_dir = folder / "results"
    shutil.rmtree(results_dir, ignore_errors=True)    # re-running replaces old results
    result, log, error = call_analysis(analysis.run_analysis, folder / "data", results_dir,
                                       column, chosen["scRNA-seq"], chosen["FACS"], panel, label_map)
    if error:
        return setup(run_id, error=error), 400

    # Remember the label choices per panel, so the setup page shows them again next time
    label_maps = {**previous_run(folder).get("label_maps", {}), panel: label_map}
    (folder / "run.json").write_text(json.dumps({
        "label_maps": label_maps,
        "n_labels_changed": int(result["labels"]["changed"].sum()),
        "seq_missing_populations": result["seq_missing_populations"],
        "notes": result["notes"],
        "panel": result["panel"],
        "panel_description": analysis.PANELS[result["panel"]]["description"],
        "condition_column": result["condition_column"],
        "conditions": result["conditions"],          # {"scRNA-seq": [...], "FACS": [...]}
        "n_facs_samples": int(result["n_facs_samples"]),
        "n_cells": int(result["n_cells"]),
        "n_seq_samples": int(result["n_seq_samples"]),
        "n_seq_people": int(result["n_seq_people"]),
        "seq_source": seq_source(folder),
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

    comparison = read_result(results_dir, "A_method_comparison_stats.csv")
    by_condition = read_result(results_dir, "B_condition_tests_proportions.csv")
    rna_stats = read_result(results_dir, "C_protein_vs_rna_stats.csv")
    # For each section: the summary table and which of its columns to show
    tables = {
        "A": (comparison,
              ["comparison", "population", "n_seq", "n_facs", "median_seq", "median_facs",
               "seq_minus_facs", "seq_minus_facs_ci95_low", "seq_minus_facs_ci95_high",
               "seq_to_facs_ratio", "test", "smallest_possible_p", "q_value_fdr", "significant"]),
        "B": (by_condition,
              ["population", "method", "test", "comparison", "p_value",
               "q_value_fdr", "significant"]),
        "C": (rna_stats, ["marker", "n", "spearman_r", "spearman_q_fdr", "pearson_r", "pearson_q_fdr"]),
        "D": (read_result(results_dir, "D_imaging_tests_between_populations.csv"),
              ["feature", "test", "comparison", "p_value", "q_value_fdr", "significant"]),
        "E": (read_result(results_dir, "E_atlas_proportion_tests.csv"),
              ["population", "method", "group", "n_samples", "median_samples", "atlas_proportion",
               "difference", "n_atlas_donors", "test", "smallest_possible_p", "q_value_fdr", "significant"]),
        "F": (read_result(results_dir, "F_imaging_population_summary.csv"),
              ["population", "feature", "n_samples", "mean", "sd", "ci95_low", "ci95_high",
               "median", "cv_percent", "n_outlier_samples"]),
    }
    n_significant = int(by_condition["significant"].sum()) if by_condition is not None else 0
    n_method_significant = int(comparison["significant"].fillna(False).astype(bool).sum()) if comparison is not None else 0

    # Section A also shows the composition overview: every group's average proportions side by side
    overview = read_result(results_dir, "A_composition_overview.csv")
    overview_html = table_html(analysis.composition_overview_wide(overview)) if overview is not None else None
    # ... and which population every scRNA-seq label was counted in
    labels_html = table_html(read_result(results_dir, "A_label_assignment.csv"))

    files = sorted(p.name for p in results_dir.iterdir())
    sections = []
    for letter, title in SECTIONS.items():
        table, columns = tables[letter]
        shown = [] if table is None or table.empty else [c for c in columns if c in table.columns]
        sections.append({
            "letter": letter, "title": title,
            "overview": overview_html if letter == "A" else None,
            "labels": labels_html if letter == "A" else None,
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
                           n_significant=n_significant, n_method_significant=n_method_significant,
                           alpha=analysis.ALPHA, min_samples=analysis.MIN_SAMPLES_PER_GROUP,
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
