# Day 04: a web interface for the FACS vs scRNA-seq program

In [day03](../day03) I built a command-line program, `facs_vs_scrnaseq.py`. It compares urine cell populations measured on the same samples in two ways:

- **BD FACSDiscover S8 sorting**, gated in FlowJo
- **10x Flex single-cell RNA-seq**

The populations are immune cells (CD45⁺), proximal tubule cells (CD45⁻ CD10⁺ CD13⁺), epithelial cells (CD45⁻ NOT CD10⁺CD13⁺, EpCAM⁺) and other CD45⁻ cells (CD45⁻ NOT CD10⁺CD13⁺, EpCAM⁻). The FACS gating is Cells → Single Cells → Live (Fixable Yellow⁻) → CD45⁺ / CD45⁻; CD45⁻ → CD10⁺CD13⁺ / NOT CD10⁺CD13⁺; NOT CD10⁺CD13⁺ → EpCAM⁺ / EpCAM⁻. CD10/CD13 is gated before EpCAM because proximal tubule cells can express EpCAM: with EpCAM⁺ and CD10⁺CD13⁺ side by side under CD45⁻, double-positive cells would be counted twice. With these gates every live cell is in exactly one population, so the proportions add up to 100%.
<img width="609" height="145" alt="image" src="https://github.com/user-attachments/assets/198757aa-7b32-4ba2-88fd-b27c083d29bf" />


For day04 I copied that program into this folder and added a web application built with **Flask**. The web app reuses the program's "business logic" (loading the files, building the tables, statistics and plots) and puts a browser interface on top. I used **Claude Code** (Anthropic's AI coding assistant, running in VS Code) to do the work. This README explains what was done, how, and how I checked that it works.

---

## Open the website (quick start)

1. Open a terminal: the Terminal app, or in VS Code **View → Terminal**.
2. Go to the day04 folder (inside the Assignments folder) and start the server:
   ```bash
   cd day04
   python -m flask --app app run --port 5050
   ```
   Use a Python that has the packages listed under [How to use it](#how-to-use-it). On my Mac that's Anaconda: `/opt/anaconda3/bin/python -m flask --app app run --port 5050`.
3. Wait until the terminal shows `Running on http://127.0.0.1:5050`.
4. Open **http://127.0.0.1:5050** in your browser (Chrome, Safari, ...).
5. Upload your three files (plus, optionally, a reference atlas and a single-cell imaging export), or click **Use demo data** to try it with made-up data.
6. When you are done, go back to the terminal and press **Ctrl+C** to stop the server.

Keep the terminal open while you use the website. Closing it, or pressing Ctrl+C, switches the website off.

**Why port 5050 and not `python app.py`?** `python app.py` starts the website on port 5000. On a Mac, port 5000 is often taken by the **AirPlay Receiver**, and the browser then shows an empty "403 Forbidden" page. Port 5050 avoids this. If you switch off AirPlay Receiver (System Settings → General → AirDrop & Handoff), `python app.py` and http://127.0.0.1:5000 work too.

### What does it mean that the website runs "on a server"?

A website always has two parts:

- **The server** is a program that waits for requests and answers them. Here the server is `app.py`, built with **Flask** (a Python library for making websites). When you start it, it keeps running in the terminal and listens for requests.
- **The client** is your web browser. When you open a page or click a button, the browser sends a *request* to the server ("give me the start page", "here are my 3 files", "run the analysis with Healthy and Disease"). The server runs the Python code for that request and sends back a *response*: an HTML page, a plot image, or a zip file. The browser then shows it.

What happens in this app:

1. You open the start page → the browser asks the server for `/` → Flask runs the `index()` function in `app.py` and sends back the upload page (made from `templates/index.html` and styled by `static/style.css`).
2. You upload the files → the browser sends them to `/upload` → the server saves them in a new folder `runs/<id>/data/`.
3. You choose the conditions and click **Run analysis** → the server calls `run_analysis()` from `facs_vs_scrnaseq.py`, the same code the terminal version uses. The tables and plots are saved in `runs/<id>/results/`.
4. The browser asks for the results page → the server builds it from the saved tables and plots, together with the explanations from `explanations.py`.

What the address **http://127.0.0.1:5050** means:

- **`127.0.0.1`** (also called *localhost*) is a special address that always means "this same computer". The server and the browser both run on my Mac, so nothing is sent over the internet.
- **`5050`** is the *port*: a number that tells the computer which program a request is for, because many programs can be listening on the same computer at once.
- Because the server only listens on `127.0.0.1`, **only this computer can open the website**. Other people on the network or the internet cannot reach it, and the data never leaves the computer. To put it online for others you would need to run it on a public server with a login, which this project does not do.

---

## What's in this folder

| File | What it is |
|---|---|
| `facs_vs_scrnaseq.py` | The analysis program copied from day03, with one small refactor (see below). It still works from the terminal exactly as before. |
| `app.py` | The Flask web application. It contains no analysis code of its own and calls functions from `facs_vs_scrnaseq.py`. |
| `HOW_THE_CODE_WORKS.md` | A plain-language explanation of the logic of the code, for readers who don't know much Python. |
| `explanations.py` | The help texts on the results page: the question each section answers, what each table column means, and how to read each plot. Edit the texts here without touching the code. |
| `templates/` | The HTML pages: `base.html` (shared layout), `index.html` (upload), `setup.html` (choose conditions), `results.html` (tables and plots). |
| `static/style.css` | The page styling, with light and dark mode. |
| `test_app.py` | 17 automatic tests covering both the terminal program and the web app. |
| `test_business_logic.py` | 46 unit tests of the calculations (proportions, statistics, atlas, imaging), each on a tiny hand-made table with a known answer. |
| `TESTS_README.md` | How the unit tests in `test_business_logic.py` were written and what each one checks. |
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
   - Upload the three input files: the sample sheet `.csv`, the FlowJo export `.xlsx` and the per-cell scRNA-seq `.csv`. Two more files are optional: a reference atlas `.csv` (see Step 5) and a single-cell imaging export `.csv` (see Step 6).
   - Or click **Use demo data** to make the made-up demo dataset (10 samples, 5 Healthy and 5 Disease, plus a reference atlas and single-cell imaging data) with the script's own `make_demo_data()`.
   - The page also shows the current SETTINGS from the script (populations, gate paths, markers and channels, imaging features), so you can see what the program expects.
2. **Setup page (`/run/<id>`):** this replaces typing numbers in the terminal.
   - Pick the condition column with radio buttons.
   - Tick the conditions you want; each one shows how many samples it has.
   - Ticking nothing includes everything, like pressing Enter in the terminal version.
3. **Results page (`/run/<id>/results`):** one section per analysis (A proportions, B conditions, C protein vs RNA, D imaging, E reference atlas when an atlas was given, F imaging profiles per cell type). Each section has:
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

### Step 5: Comparison with a reference single-cell atlas (section E)

Later I added an **optional** fourth input: a reference single-cell atlas of urine with known cell types, such as a published dataset. When it is given, a new section **E** compares my FACS and Flex results with the atlas. When it isn't, section E is skipped and everything else runs exactly as before.

**What was added:**
- **Settings:** `ATLAS_FILE` and the atlas column names. Each population also gets an `atlas_labels` list, because an atlas usually names cell types differently from my notebook (for example "Myeloid cell" → Immune).
- **Two atlas formats:**
  - one row per cell (`cell_type`, optional `donor_id`, optional gene columns `PTPRC, MME, ANPEP, EPCAM`);
  - a summary table with one row per cell type, or per cell type × donor, with an `n_cells` column (for example from a paper's supplementary table).
- **New analysis code:** a separate block in `facs_vs_scrnaseq.py` (section "4b"). It loads the atlas, turns it into proportions with the same rule as my samples, runs the statistics, and saves 4 tables and 3 plots:

  | Output | What it shows |
  |---|---|
  | `E1_atlas_proportions.png` | Each population per condition (FACS and scRNA-seq) next to the atlas donors, with the pooled atlas proportion as a dashed line. |
  | `E2_atlas_composition_difference.png` | For each sample, the % of cells that would have to change population to match the atlas, compared with how much the atlas donors differ among themselves. |
  | `E3_atlas_marker_profiles.png` | Marker patterns as heatmaps (atlas RNA vs Flex RNA vs FACS protein), each scaled 0–1, with the correlation to the atlas pattern. Needs a per-cell atlas with the gene columns. |
  | `E_atlas_proportion_tests.csv` | Per population × method × condition: median vs atlas, difference, test and FDR q-value. With ≥ 3 atlas donors: Mann-Whitney U of my samples vs the atlas donors. Otherwise: one-sample Wilcoxon vs the atlas proportion. |
  | `E_atlas_proportions.csv`, `E_atlas_composition_difference.csv`, `E_atlas_marker_profiles.csv` | The numbers behind the plots. |
- **Website:**
  - an optional "Reference atlas" upload field;
  - an "Atlas labels" column in the settings table;
  - a section E on the results page, with explanations of the question, the table columns and each plot.
- **Demo data:** now also includes a made-up atlas of 8 healthy donors.

**How it was checked:**
- **No existing outputs changed.** All 14 existing plots and tables (sections A–D) were compared byte for byte with the version before the atlas was added, and all were identical.
- **New tests:** 14 tests pass. The new ones check:
  - an upload without an atlas skips section E with a message;
  - an upload with an atlas shows all 7 E files;
  - a summary-table atlas gives exactly the same proportions as the per-cell atlas it was made from, and skips the marker heatmap.

### Step 6: Imaging profiles per cell type (section F)

Next I asked for everything the S8 imaging can tell about each cell type: mean morphological features, outliers, and anything else interesting. The S8 images give **morphological** features (size, shape, texture, where the signal sits in the cell), not physiological ones. Those need functional dyes. Section F describes each population with all its imaging features at once. It runs whenever `IMAGING_FEATURES` is set, and section D is unchanged.

**Two levels of data:**
- **Sample level (always):** the FlowJo table already has the **median** of each feature per sample and population. Plots F1–F5 use it.
- **Cell level (optional new file, `facs_imaging_cells.csv`):** one row per FACS cell, with its population and imaging features, exported from FlowJo. Plots F6–F8 need it to look at single cells.

**Outliers** are found with a **robust z-score**, `0.6745 × (value − median) / MAD`, where MAD is the median absolute deviation. A value above 3.5 counts as an outlier (`OUTLIER_CUTOFF` in SETTINGS). Unlike the usual mean/SD z-score, a few extreme values can't hide themselves by inflating the spread.

| Output | What it shows |
|---|---|
| `F1_imaging_feature_means.png` | Mean of each feature per cell type, with a 95% confidence interval and every sample as a dot. Outlier samples are ringed in red and named. |
| `F2_imaging_fingerprint.png` | A heatmap "fingerprint" of each cell type: the mean values, colored by how much higher or lower they are than in the other populations. |
| `F3_imaging_pca.png` | A PCA map of all features at once: do the cell types separate by morphology alone, and does a condition shift them? Arrows show what each feature contributes. |
| `F4_imaging_outlier_samples.png` | Robust z-score of every sample × population × feature; outliers are labelled. |
| `F5_imaging_vs_markers.png` | Spearman correlation between imaging features and FACS marker intensity across samples, within each cell type (* = significant after FDR). |
| `F6_imaging_cell_distributions.png` | *(per-cell file)* Violin plots of all single cells per cell type. |
| `F7_imaging_cell_outliers.png` | *(per-cell file)* The % of outlier cells per sample and cell type, for example doublets, clumps or debris. |
| `F8_imaging_morphology_map.png` | *(per-cell file)* Single cells on a map of the first two features, per condition, with outlier cells marked. |
| `F_imaging_population_summary.csv` | Per cell type × feature: n, mean, SD, 95% CI, median, CV%, min, max and the number of outlier samples (shown on the website). |
| `F_imaging_outlier_samples.csv`, `F_imaging_pca.csv`, `F_imaging_pca_loadings.csv`, `F_imaging_marker_correlations.csv`, `F_imaging_cell_summary.csv` | The numbers behind the plots. The cell summary has per sample × cell type × feature: n cells, mean, median, SD, CV%, 5th/95th percentile and % outlier cells. |

**Website:** a new optional upload field "Single-cell imaging export", and section F on the results page with explanations of every plot and table column. If the per-cell file is missing, a note says that F6–F8 need it. The demo data now also includes single-cell imaging data, with about 2% "odd" cells such as doublets.

**How it was checked:**
- **No existing outputs changed.** All 21 earlier plots and tables (sections A–E) are byte for byte identical to the version before this step.
- **Plots checked by eye:** the new plots were rendered and checked. One fix: in F3 the feature arrows could fall outside the plot, so they are now scaled to always fit.
- **New tests:** 17 tests pass. The new ones check:
  - without the per-cell file, F1–F5 appear and F6–F8 are skipped with a note;
  - with it, F6–F8 appear;
  - `robust_z` flags an injected outlier but not the normal values, and doesn't crash when all values are equal;
  - with no imaging features set, sections D and F are skipped and the rest still runs.

---

## How to use it

Requirements: Python with `pandas`, `numpy`, `scipy`, `matplotlib`, `openpyxl` and `flask`. Anaconda has all of them. On my Mac I use `/opt/anaconda3/bin/python`, because the system `python3` doesn't have pandas.

**Web interface:** see [Open the website (quick start)](#open-the-website-quick-start) at the top. In short:
```bash
cd day04
python -m flask --app app run --port 5050
```
Then open **http://127.0.0.1:5050** in a browser. Stop the server with `Ctrl+C`. (`python app.py` uses port 5000, which only works if AirPlay Receiver is switched off.)

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
python -m pytest test_app.py -v              # the whole program and the website
python -m pytest test_business_logic.py -v   # unit tests of the calculations (see TESTS_README.md)
```

**Changing the settings:** gate paths, fluorochrome channel names, the cell-type labels for each population, imaging features and statistical thresholds are still edited in the SETTINGS section at the top of `facs_vs_scrnaseq.py`. Restart `app.py` after changing them. The start page shows the values it's using.

---

## Input files (1–3 as in day03, 4 and 5 are new)

1. **`metadata.csv`:** one row per sample, with `sample_id`, `facs_name` (the sample name as it appears in FlowJo) and any condition columns.
2. **FlowJo export `.xlsx`:** one row per sample, with the count of each population gate and the median of each marker and imaging feature in each population.
3. **`scrnaseq_cells.csv`:** one row per cell, with `cell_barcode, sample_id, cell_type, PTPRC, MME, ANPEP, EPCAM` (log-normalized). The script's header comment has a snippet to export this from Scanpy.
4. **`reference_atlas.csv` (optional, new in day04):** a reference single-cell atlas of urine. Either:
   - one row per cell, with `cell_type`, an optional `donor_id` and optional `PTPRC, MME, ANPEP, EPCAM` columns;
   - or one row per cell type (× donor), with `cell_type`, an optional `donor_id` and `n_cells`.

   Put the atlas cell-type names in `atlas_labels` in SETTINGS. The script's header comment has a snippet to export an atlas from Scanpy.
5. **`facs_imaging_cells.csv` (optional, new in day04):** the single FACS cells with their imaging features. It has one row per cell, with `facs_name` (or `sample_id`), `population` (one of the names in `POPULATIONS`) and one column per imaging feature, named exactly as in `IMAGING_FEATURES`. In FlowJo, export the events of each population gate as CSV and combine them. The script's header comment has a snippet that does this.

---

## Limitations and notes

- **Placeholder settings:** I don't have real data yet. The channel names (`FLUOR1`–`4`), gate paths, imaging feature names and atlas labels are placeholders. Everything has only been tested on the made-up demo data.
- **Comparing with an atlas:** an atlas is made with other protocols and is usually from healthy people. Small differences from it are expected, and differences in a disease condition can be biology rather than an error. The marker heatmap (E3) only compares patterns (which population is high or low), not absolute levels.
- **Imaging and outliers:** the S8 imaging features describe morphology, not physiology. An outlier is a reason to look at the images of that sample or those cells, not an automatic reason to exclude them. With few samples, the imaging vs marker correlations (F5) are uncertain, so only strong, significant ones should be trusted.
- **Local use only:** the web app is meant to run on my own computer. It has no login, and runs are kept in `runs/` until I delete them.
- **Sample size:** with fewer than 3 samples in a group, that test is skipped. With only 3 per group, a Mann–Whitney test can't reach p < 0.05, so I should aim for at least 5 samples per condition.
