# Day 04: a web interface for the FACS vs scRNA-seq program

In [day03](../day03) I built a command-line program, `facs_vs_scrnaseq.py`. It compares urine cell populations measured on the same samples in two ways:

- **BD FACSDiscover S8 sorting**, gated in FlowJo
- **10x Flex single-cell RNA-seq**

The populations are immune cells (CD45⁺), proximal tubule cells (CD45⁻ CD10⁺ CD13⁺), epithelial cells (CD45⁻ NOT CD10⁺CD13⁺, EpCAM⁺) and other CD45⁻ cells (CD45⁻ NOT CD10⁺CD13⁺, EpCAM⁻). The FACS gating is Cells → Single Cells → Live (Fixable Yellow⁻) → CD45⁺ / CD45⁻; CD45⁻ → CD10⁺CD13⁺ / NOT CD10⁺CD13⁺; NOT CD10⁺CD13⁺ → EpCAM⁺ / EpCAM⁻. CD10/CD13 is gated before EpCAM because proximal tubule cells can express EpCAM: with EpCAM⁺ and CD10⁺CD13⁺ side by side under CD45⁻, double-positive cells would be counted twice. With these gates every live cell is in exactly one population, so the proportions add up to 100%.

For day04 I copied that program into this folder and added a web application built with **Flask**. The web app reuses the program's "business logic" (loading the files, building the tables, statistics and plots) and puts a browser interface on top. I used **Claude Code** (Anthropic's AI coding assistant, running in VS Code) to do the work. This README explains what was done, how, and how I checked that it works.

---

## What's in this folder

| File | What it is |
|---|---|
| `facs_vs_scrnaseq.py` | The analysis program copied from day03, with one small refactor (see below). It still works from the terminal exactly as before. |
| `app.py` | The Flask web application. It contains no analysis code of its own and calls functions from `facs_vs_scrnaseq.py`. |
| `templates/` | The HTML pages: `base.html` (shared layout), `index.html` (upload), `setup.html` (choose conditions), `results.html` (tables and plots). |
| `static/style.css` | The page styling, with light and dark mode. |
| `test_app.py` | 12 automatic tests covering both the terminal program and the web app. |
| `.gitignore` | Keeps uploaded data and generated results out of git. |

---

## How it was done

### Step 1: Copy the program and look at its structure

Claude copied `day03/facs_vs_scrnaseq.py` into `day04/` and read through it to find what could be reused. The program was already split into functions:

- **Loading:** `load_metadata`, `load_facs`, `load_seq`
- **Building tables:** `build_proportions`, `build_protein_vs_rna`, `build_imaging`
- **Statistics:** `method_agreement_stats`, `condition_stats`, `protein_rna_stats`, `imaging_between_populations`
- **Plots:** the `plot_*` functions

This is the business logic. Only `main()` mixed two jobs: it read the terminal options and asked questions with `input()`, and it also ran all the analysis steps.

### Step 2: A small refactor so the logic can be called from somewhere else

Claude split `main()` in two:

- **`run_analysis(data_dir, output_dir, condition_column, conditions)`** runs steps 1–5 (load, build tables, statistics, save tables, plots) and returns the summary tables. It takes the input and output folders as arguments, so the web app can give every analysis its own folder.
- **`main()`** now only reads the command-line options, makes demo data if `--demo` is given, calls `run_analysis()` and prints the summary.

If you don't pass the condition column and conditions to `run_analysis()`, it still asks for them in the terminal, so the command-line behaviour is unchanged. Apart from a two-line note in the header comment about the web app, nothing else in the script was changed. The SETTINGS section, the statistics and the plots are identical to day03.

### Step 3: The Flask web app (`app.py`)

The web app is a thin layer over the analysis script. It does `import facs_vs_scrnaseq as analysis` and calls its functions. You go through three pages:

1. **Upload page (`/`):**
   - Upload the three input files: the sample sheet `.csv`, the FlowJo export `.xlsx` and the per-cell scRNA-seq `.csv`.
   - Or click **Use demo data** to make the made-up demo dataset (10 samples, 5 Healthy and 5 Disease) with the script's own `make_demo_data()`.
   - The page also shows the current SETTINGS from the script (populations, gate paths, markers and channels, imaging features), so you can see what the program expects.
2. **Setup page (`/run/<id>`):** this replaces typing numbers in the terminal.
   - Pick the condition column with radio buttons.
   - Tick the conditions you want; each one shows how many samples it has.
   - Ticking nothing includes everything, like pressing Enter in the terminal version.
3. **Results page (`/run/<id>/results`):** one section per analysis (A proportions, B conditions, C protein vs RNA, D imaging). Each section has:
   - a summary table
   - all the plots
   - links to the CSV tables
   - a **Download all results (.zip)** button

   The messages the script printed (for example warnings about samples missing from a file) are under "Messages from the analysis".

Some design decisions I want to point out:

- **Every analysis gets its own folder.** It's `runs/<random id>/`, holding `data/` (the inputs), `results/` (the outputs) and `run.json` (which conditions were chosen). Uploaded files are saved under the names the script expects, so `run_analysis()` can read them without any changes.
- **Errors show on the page instead of crashing.** When a column name is wrong, the original script stops with a message that lists the real column names. The web app catches that and shows the same message on the page. You can then fix SETTINGS and try again.
- **Analyses run one at a time.** A lock makes this happen, because matplotlib and the script's `OUTPUT_DIR` setting are shared by the whole program.
- **Safety:**
  - The server only listens on `127.0.0.1` (this computer) and has debug mode off.
  - Run IDs must be exactly 32 hex characters.
  - Files are served with Flask's `send_from_directory`, so nobody can read files outside a run's results folder.
  - The upload size is limited to 500 MB.

### Step 4: Checking that both systems work

I asked Claude to check both the old (day03) and the new (day04) versions. It did the following:

1. **Old system (day03), terminal.** Claude ran the original day03 script on the demo data, in a temporary copy so the results already committed in day03 wouldn't change. It tested both ways of choosing conditions:
   - with options: `--condition-column condition --conditions Healthy Disease`
   - interactively, typing `1` and then pressing Enter

   Both finished with `Done.` and wrote all 14 result files. They found 2 significant condition differences (immune cell proportions, Healthy vs Disease, in both FACS and scRNA-seq).
2. **New system (day04), terminal.** Claude ran the same commands with the refactored day04 script, then compared all 8 CSV tables with the day03 output. **All 8 were byte-for-byte identical**, so the refactor didn't change any results.
3. **New system (day04), web app.** Claude started the real server with `python app.py` and went through every page with HTTP requests:
   - the upload page loads (200 OK)
   - the demo upload leads to the setup page, which shows "5 samples" per condition
   - running the analysis leads to the results page, which shows 6 plots and "2 significant condition differences", the same as the terminal
   - a plot is served as `image/png`
   - the zip download contains all 14 files
   - a made-up run ID gives 404
4. **Automatic tests (`test_app.py`).** These can be re-run any time. All **12 tests pass** (about 40 seconds). They check:
   - `run_analysis()` writes all 14 files and finds 10 samples
   - with a single condition, the condition plot is skipped
   - an unknown condition stops with a clear message
   - the terminal command still works
   - the upload page shows the form and the settings
   - the full demo run works in the web app, including the zip download
   - uploading your own files works, and ticking no conditions means all of them
   - a missing file, a wrong file type, a sample sheet missing columns, and a wrong gate name each give an error message, not a crash
   - unknown or unsafe run IDs give 404

---

## How to use it

Requirements: Python with `pandas`, `numpy`, `scipy`, `matplotlib`, `openpyxl` and `flask`. Anaconda has all of them. On my Mac I use `/opt/anaconda3/bin/python`, because the system `python3` doesn't have pandas.

**Web interface:**
```bash
cd day04
python app.py
```
Then open **http://127.0.0.1:5000** in a browser. Stop the server with `Ctrl+C`.

**Terminal, as in day03:**
```bash
cd day04
python facs_vs_scrnaseq.py --demo                    # made-up demo data
python facs_vs_scrnaseq.py                           # real data from day04/data/
python facs_vs_scrnaseq.py --condition-column condition --conditions Healthy Disease
```

**Tests:**
```bash
cd day04
python -m pytest test_app.py -v
```

**Changing the settings:** gate paths, fluorochrome channel names, the cell-type labels for each population, imaging features and statistical thresholds are still edited in the SETTINGS section at the top of `facs_vs_scrnaseq.py`. Restart `app.py` after changing them. The start page shows the values it's using.

---

## Input files (same as day03)

1. **`metadata.csv`:** one row per sample, with `sample_id`, `facs_name` (the sample name as it appears in FlowJo) and any condition columns.
2. **FlowJo export `.xlsx`:** one row per sample, with the count of each population gate and the median of each marker and imaging feature in each population.
3. **`scrnaseq_cells.csv`:** one row per cell, with `cell_barcode, sample_id, cell_type, PTPRC, MME, ANPEP, EPCAM` (log-normalized). The script's header comment has a snippet to export this from Scanpy.

---

## Limitations and notes

- **Placeholder settings:** I don't have real data yet. The channel names (`FLUOR1`–`4`), gate paths and imaging feature names are placeholders. Everything has only been tested on the made-up demo data.
- **Local use only:** the web app is meant to run on my own computer. It has no login, and runs are kept in `runs/` until I delete them.
- **Sample size:** with fewer than 3 samples in a group, that test is skipped. With only 3 per group, a Mann–Whitney test can't reach p < 0.05, so I should aim for at least 5 samples per condition.
