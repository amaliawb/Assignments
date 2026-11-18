# Day 04: a web interface for the FACS vs scRNA-seq program

In [day03](../day03) I built a command-line program, `facs_vs_scrnaseq.py`. It compares urine cell populations measured in two ways:

- **BD FACSDiscover S8 sorting**, gated in FlowJo
- **10x Flex single-cell RNA-seq**

The two methods **don't have to measure the same samples** (see [Step 7](#step-7-facs-and-scrna-seq-on-different-samples)). For example, the scRNA-seq results of a group of patients can be compared with FACS of a **pool of healthy urine** from several people, or with FACS of **other patients**.

You choose which antibodies were used on FACS (see [Step 8](#step-8-facs-panels-people-and-the-pool)). The **CD45 + EpCAM** panel and the **CD45 + EpCAM + CD10 + CD13** panel both use Fixable Yellow for live/dead. The gating is:

- **Both panels:** Cells → Single Cells → Live (Fixable Yellow⁻). CD45⁺ cells are immune cells: SSC-high = **neutrophils**, SSC-low = **other immune**.
- **CD45 + EpCAM:** CD45⁻ EpCAM-high = **distal tubule, collecting duct and urothelium**. CD45⁻ EpCAM-low = **proximal tubule** (the same cells as CD10⁺CD13⁺, plus other EpCAM-negative CD45⁻ cells).
- **CD45 + EpCAM + CD10 + CD13:** CD45⁻ CD10⁺CD13⁺ = **proximal tubule**; the rest splits into EpCAM-high (distal tubule, collecting duct, urothelium) and EpCAM-low (**other CD45⁻**).

CD10/CD13 is gated before EpCAM because proximal tubule cells can carry some EpCAM: double-positive cells would otherwise be counted twice. In both panels every live cell is in exactly one population, so the proportions add up to 100%. The diagram below shows the original day03 gating.
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
5. Upload the sample sheet and the FlowJo export (the sample sheet says which method measured each sample). The scRNA-seq file is optional: leave it empty to use the default Flex data (see [Step 9](#step-9-the-real-flex-data-as-default)). A reference atlas and a single-cell imaging export are optional too. Or click **Use demo data** to try it with made-up data.
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
| `test_app.py` | 23 automatic tests covering both the terminal program and the web app. |
| `test_business_logic.py` | 86 unit tests of the calculations (proportions, statistics, atlas, imaging), each on a tiny hand-made table with a known answer. |
| `TESTS_README.md` | How the unit tests in `test_business_logic.py` were written and what each one checks. |
| `.gitignore` | Keeps uploaded data and generated results out of git. |

---

## How it was done

### Step 1: Copy the program and look at its structure

Claude copied `day03/facs_vs_scrnaseq.py` into `day04/` and read through it to find what could be reused. The program was already split into functions:

- **Loading:** `load_metadata`, `load_facs`, `load_seq`
- **Building tables:** `build_proportions`, `build_protein_vs_rna`, `build_imaging`
- **Statistics:** `method_agreement_stats` (replaced by `method_comparison_stats` in Step 7), `condition_stats`, `protein_rna_stats`, `imaging_between_populations`
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
   - Tick the conditions you want; each one shows how many samples it has. (Since Step 7 there are two sets of checkboxes: one for the scRNA-seq samples and one for the FACS samples.)
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

### Step 7: FACS and scRNA-seq on different samples

Until now the program assumed that **every sample was measured by both methods**. Section A then paired each sample's FACS result with its own scRNA-seq result (correlation, Bland-Altman, paired Wilcoxon test). In practice the two methods are run on **different samples**. The scRNA-seq results of a group of people are compared with FACS of either:

- a **pool of healthy urine** from several people. This is ONE FACS sample.
- **several patient samples** measured by FACS, from other people than the scRNA-seq ones.

The program was changed so the two methods are compared as **groups of samples**, never sample by sample.

**What changed for the user:**
- **Sample sheet:** a new `method` column (`FACS` or `scRNA-seq`) says which method measured each sample. Only FACS samples need a `facs_name`. A pool is simply one FACS row. Ways of writing the method like `facs`, `Flex` or `seq` are understood. An older sheet without a `method` column still works: every sample is then taken to be measured by both methods.
- **Choosing conditions:** you choose the conditions separately for the scRNA-seq samples and for the FACS samples. In the terminal: `--seq-conditions Disease --facs-conditions "Healthy pool" Disease`. On the website there are two sets of checkboxes. **Every chosen scRNA-seq condition is compared with every chosen FACS condition**, for example "Disease (scRNA-seq) vs Healthy pool (FACS)" and "Disease (scRNA-seq) vs Disease (FACS)" in one run.

**What changed in the analysis:**

| Section | Before (same samples) | Now (different samples) |
|---|---|---|
| A | Paired: correlation, Bland-Altman, paired Wilcoxon | Unpaired, per population and per pair of groups: the medians, the difference (scRNA-seq − FACS) with a bootstrap 95% interval, and the ratio. Two groups of ≥ 3 samples → **Mann-Whitney U**. One side is a single sample (a pool) → **one-sample Wilcoxon** of the other side against that value. FDR-corrected. |
| B | Conditions within each method | Unchanged. A single pool can't be tested here (it's 1 sample). |
| C | Protein vs RNA per sample | **Per population**: the typical protein level (median over the FACS samples) vs the typical RNA level (median over the scRNA-seq samples). Spearman per marker across the populations, plus one "All markers (scaled)" test over all population × marker points. |
| D, F | All samples | The FACS samples and FACS conditions only (only FACS has images). |
| E | Per method and condition | The same, with each method's own conditions. A single FACS pool is tested by comparing the atlas donors with the pool's value. |

New or renamed outputs: `A2_proportions_by_method.png` (box per group), `A3_method_difference.png` (difference per population with its interval), `A_method_comparison_stats.csv`, `C_protein_per_sample.csv` and `C_rna_per_sample.csv`. These were removed: `A2_proportions_FACS_vs_seq.png`, `A3_bland_altman.png` and `A_method_agreement_stats.csv`.

**The demo data** now has different people per method: scRNA-seq of 10 people (5 Healthy, 5 Disease), and FACS of a healthy pool plus 6 other Disease patients. scRNA-seq is made to capture fewer epithelial cells, so section A has a method effect to find.

**Two statistical fixes found along the way:**
- **Exact Spearman p-values for few points.** scipy's p-value for a Spearman correlation is an approximation that is far too small for few points. For 4 points in perfect order it gives p = 0, while the smallest possible p is 2/24 = 0.083. For up to 8 points the program now calculates the exact p-value by trying every possible ordering. This matters for section C (4 populations) and for small groups everywhere.
- **No test when every value is identical.** If every value is identical (for example a population that is 0 in all samples), no test is reported, instead of a "test" with an empty p-value.

**Speed improvements** (the results stay the same):
- Only the needed columns are read from the scRNA-seq file, the atlas and the per-cell imaging export. Exports with many extra columns load faster and use less memory.
- The violin plots (F6) estimate their shape from at most 5,000 random cells per population (`MAX_CELLS_PER_VIOLIN`). The median line still uses all cells. With 437,000 cells this plot took 6.6 s; it now takes about 1 s, however many cells there are.
- Several tables are now built and summarised all at once instead of row by row, using pandas `melt` and `groupby`. This applies to the proportions, the marker tables, the per-cell summary and the imaging-marker correlations. All box plots of a panel are drawn in one call.
- Most of the remaining run time (~16 s for the demo) is matplotlib drawing the 17 figures. The statistics take under a second.

**Important for interpreting the results:**
- **A single pool has little statistical power.** 5 scRNA-seq samples tested against one pool can never reach p < 0.05: the smallest possible p is 2/32 = 0.0625. You need at least 6 samples, and more after FDR correction. With a pool, look mainly at the size of the difference and its interval.
- **Measuring the same pool again is not new samples.** Repeated FACS measurements of the same pool are not independent samples. Listing them as separate rows would make the tests too optimistic.
- **Method effect or biology?** With different people per method, a difference can be either. Comparing biologically similar groups (e.g. Healthy scRNA-seq vs a Healthy FACS pool) is the fairest test of the methods themselves. The reference atlas (section E) gives both methods a shared yardstick.

**How it was checked:**
- All 82 tests pass (64 unit tests, 18 app tests). The new tests check:
  - the `method` column and the old sheet format;
  - choosing conditions per method;
  - proportions for different samples per method;
  - the unpaired tests, with exact p-values worked out by hand (e.g. 2/32 for 5 samples against one pool);
  - the bootstrap interval;
  - the exact Spearman p-value;
  - the per-population protein vs RNA comparison;
  - the atlas test for a single pool;
  - an end-to-end run with an old-style sample sheet.
- The exact Spearman p-value was checked against scipy's own permutation test (identical results).
- The new plots were rendered and checked by eye.

### Step 8: FACS panels, people and the pool

Three more changes, from how the real experiment is set up.

**1. Choosing the FACS panel.** The antibodies used on FACS can be CD45 + EpCAM, or CD45 + EpCAM + CD10 + CD13. Both panels always use Fixable Yellow. The new setting `FACS_PANEL` in SETTINGS holds the default. It can be changed per run with `--panel CD45_EpCAM` in the terminal, or with a choice on the website's setup page. The panel decides:

- **the gating, and so which populations exist:**

  | Population | CD45 + EpCAM | CD45 + EpCAM + CD10 + CD13 |
  |---|---|---|
  | Neutrophils | CD45⁺ SSC-high | CD45⁺ SSC-high |
  | Other immune | CD45⁺ SSC-low | CD45⁺ SSC-low |
  | Proximal tubule | CD45⁻ EpCAM-low (+ other EpCAM⁻ cells) | CD45⁻ CD10⁺CD13⁺ |
  | EpCAM-high epithelium (distal tubule, collecting duct, urothelium) | CD45⁻ EpCAM-high | CD45⁻ NOT CD10⁺CD13⁺ EpCAM-high |
  | Other CD45⁻ | – (inside EpCAM-low) | CD45⁻ NOT CD10⁺CD13⁺ EpCAM-low |

- **which Flex labels each population is compared with.** Every Flex cell type is listed under the gate it would fall into on FACS. With CD45 + EpCAM, podocytes and fibroblasts (EpCAM⁻) land in the EpCAM-low gate together with the proximal tubule cells, so their labels are counted there too. Both methods then count the same cells.
- **which markers are compared with their genes** (section C): CD45 and EpCAM, plus CD10 and CD13 with the 4-marker panel.

If the FlowJo file doesn't match the chosen panel, the error message now also asks whether the right panel was chosen.

**2. Several samples per person.** The Flex data is 12 samples from 8 people, and samples of the same person are not independent. A new optional column `donor_id` in the sample sheet says which person each sample comes from. Before any statistics, the samples of one person are combined by adding up their cells, so every person counts once. Section C combines its marker levels the same way.

- **New table:** `A_proportions_per_sample.csv` keeps every single sample.
- **People, not samples:** `A_proportions.csv` and all tests now use one value per person. The terminal prints a note like "Combined 12 scRNA-seq samples into 8 people".

**3. The pool as one reference value.** The FACS pool is one sample. Section A now starts with a **composition overview**: the pool's proportions next to the Flex proportions (the average over the people of each condition) and the reference atlas. It comes as a table (the first table of section A on the website, and `A_composition_overview.csv`) and as a stacked-bar plot with the percentages written in (`A4_composition_overview.png`).

The tests of the Flex people against the pool's value are still there. With 8 people in a group, p can go down to 0.008. With 4 people per condition, the smallest possible p is 0.125.

**Demo data:**
- Built from individual cell types (neutrophils, T cells, macrophages, proximal tubule, distal tubule, collecting duct, urothelium, podocytes, fibroblasts), which are summed into the gates. The FlowJo export holds the gates of both panels, so the same demo works with either panel.
- scRNA-seq is 12 samples from 8 people.
- scRNA-seq captures fewer neutrophils, as is common because neutrophils have little RNA.

**How it was checked:**
- All 97 tests pass (77 unit tests, 20 app tests). The new tests check:
  - both panels: no double counting, every gate under the live gate, neutrophils = CD45⁺ SSC-high, every cell type in exactly one population;
  - switching panels;
  - podocytes counted with the EpCAM-low cells in the 2-marker panel;
  - combining samples per person (by hand: (3 + 7) / (10 + 30) = 0.25);
  - the composition overview (each group adds up to 100%);
  - a full run with the 2-marker panel;
  - the panel choice on the website.
- The new plots were checked by eye.

### Step 9: The real Flex data as default

The website now uses the real Flex data when the scRNA-seq upload is left empty:
`X:\Rotation_students\Amalia\1. data\results_v5_flex\data\v5\flex_v5_processed.h5ad` (60,296 cells, 12 samples from 8 people). Uploading a scRNA-seq file still overrides it.

**What happens:**
- **Reading the `.h5ad`:** the file is 4 GB, so it is never loaded whole. `read_h5ad_cells(...)` reads only the per-cell columns that are needed and the four marker genes (PTPRC, MME, ANPEP, EPCAM, layer `log1p_cp10k`) out of the sparse matrix. This takes about 20 seconds.
- **Cache:** the result is saved in `day04/default_flex/` and reused, so later runs start at once. The cache is rebuilt when the `.h5ad` changes, and is used alone when the X: drive is not connected. The folder is in `.gitignore`: the patient data never goes to GitHub.
- **Sample sheet:** if your sample sheet has no scRNA-seq rows, the 12 Flex samples are added to it, with `donor_id` (= `patient_id`) and two condition columns: `condition` (healthy_control / donor / recipient) and `timepoint` (control, T00, T07, T30). **Name the condition column of your FACS samples `condition` too**, so both methods can be chosen from one column. If your sheet already lists scRNA-seq samples, they are used as written.
- The setup and results pages say whether the scRNA-seq data was uploaded or the default.
- **Clearer setup page.** The "condition column" is now explained: it is the column of the sample sheet that says which group each sample belongs to. Each possible column is shown with a short explanation and its values for each method (e.g. `condition` → scRNA-seq: healthy_control, donor, recipient · FACS: Healthy pool, Disease). Barcode, file and ID columns such as `flex_probe_barcode` are no longer offered. A column that only one method has values in (e.g. `timepoint`) is greyed out. The explanations live in `CONDITION_COLUMN_HELP` and `CONDITION_VALUE_HELP` in SETTINGS.

**Flex labels → FACS gates.** The Flex annotation (`leiden_atlas`, e.g. `6: myeloid`) is coarser than the FACS gates: there is no "neutrophil" label. I assigned each cluster by its marker genes (mean log1p CP10k):

| Flex cluster | Marker genes | FACS population |
|---|---|---|
| 6, 7, 8, 10, 11 `myeloid` | CSF3R, FCGR3B, CXCR2 high | Neutrophils |
| 5, 9 `myeloid` | CD68, C1QA, CD163, LYZ high | Other immune |
| 4 `T_NK`, 12 `B_plasma` | CD3E / MS4A1 | Other immune |
| 3 `ascending_thick_LOH_kidney ?` | some EPCAM | EpCAM-high epithelium |
| 0, 1, 2 `epithelium_prostate_urethra ?` | KRT13, S100A8; **no** EPCAM, MME or ANPEP | EpCAM-low (CD45 + EpCAM) / Other CD45⁻ (4-marker panel) |

Clusters 0–2 look like squamous cells (KRT13), which are EpCAM-negative, so on FACS they would land in the CD45⁻ EpCAM-low gate. There is no proximal tubule cluster in the Flex data, so with the 4-marker panel the Flex "Proximal tubule" proportion is 0. Cluster numbers change if the notebook is re-clustered: then check the labels in `PANELS` again.

**How it was checked:**
- A full run with the real Flex data and the demo FACS files, with both panels: all 60,296 cells were assigned to a population, and the 12 samples were combined into 8 people.
- 10 new tests (on a tiny made-up `.h5ad`): the right genes are read from the sparse matrix; a missing gene, layer or file gives a clear message; the cache is used when the drive is missing; the Flex samples are added only to a sheet without scRNA-seq rows; an uploaded file is never replaced; every Flex cluster is in exactly one population; and a full website run without a scRNA-seq upload.

---

## How to use it

Requirements: Python with `pandas`, `numpy`, `scipy`, `matplotlib`, `openpyxl` and `flask`, plus `anndata` for the default Flex data (`requirements.txt` lists them). Anaconda has all of them. On my Mac I use `/opt/anaconda3/bin/python`, because the system `python3` doesn't have pandas.

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
python facs_vs_scrnaseq.py --condition-column condition --seq-conditions Healthy Disease --facs-conditions "Healthy pool" Disease
python facs_vs_scrnaseq.py --panel CD45_EpCAM       # FACS stained with CD45 + EpCAM only
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

1. **`metadata.csv`:** one row per sample (FACS samples **and** scRNA-seq samples), with:
   - `sample_id`;
   - `method` (`FACS` or `scRNA-seq`);
   - `facs_name`, for FACS samples only: the sample name as it appears in FlowJo;
   - `donor_id` (optional): the person the sample comes from. Samples of one person are combined;
   - any condition columns.

   For example:
   ```
   sample_id,method,facs_name,donor_id,condition
   S01,scRNA-seq,,P1,Disease
   S02,scRNA-seq,,P1,Disease
   S03,scRNA-seq,,P2,Disease
   POOL1,FACS,healthy_pool.fcs,,Healthy pool
   ```
2. **FlowJo export `.xlsx`:** one row per FACS sample, with the gates of the chosen panel, with the count of each population gate and the median of each marker and imaging feature in each population.
3. **`scrnaseq_cells.csv` (optional on the website):** without it, the default Flex data is used (see Step 9). Otherwise: one row per cell, with `cell_barcode, sample_id, cell_type, PTPRC, MME, ANPEP, EPCAM` (log-normalized). The script's header comment has a snippet to export this from Scanpy.
4. **`reference_atlas.csv` (optional, new in day04):** a reference single-cell atlas of urine. Either:
   - one row per cell, with `cell_type`, an optional `donor_id` and optional `PTPRC, MME, ANPEP, EPCAM` columns;
   - or one row per cell type (× donor), with `cell_type`, an optional `donor_id` and `n_cells`.

   Put the atlas cell-type names in `atlas_labels` in SETTINGS. The script's header comment has a snippet to export an atlas from Scanpy.
5. **`facs_imaging_cells.csv` (optional, new in day04):** the single FACS cells with their imaging features. It has one row per cell, with `facs_name` (or `sample_id`), `population` (one of the names in `POPULATIONS`) and one column per imaging feature, named exactly as in `IMAGING_FEATURES`. In FlowJo, export the events of each population gate as CSV and combine them. The script's header comment has a snippet that does this.

---

## Limitations and notes

- **Placeholder settings:** I don't have real FACS data yet. The channel names (`FLUOR1`–`4`), gate paths, imaging feature names and atlas labels are placeholders. The Flex side has been run on the real data (Step 9); the FACS side only on the made-up demo data.
- **Comparing with an atlas:** an atlas is made with other protocols and is usually from healthy people. Small differences from it are expected, and differences in a disease condition can be biology rather than an error. The marker heatmap (E3) only compares patterns (which population is high or low), not absolute levels.
- **Imaging and outliers:** the S8 imaging features describe morphology, not physiology. An outlier is a reason to look at the images of that sample or those cells, not an automatic reason to exclude them. With few samples, the imaging vs marker correlations (F5) are uncertain, so only strong, significant ones should be trusted.
- **Local use only:** the web app is meant to run on my own computer. It has no login, and runs are kept in `runs/` until I delete them.
- **Sample size:** with fewer than 3 samples in a group, that test is skipped. With only 3 per group, a Mann–Whitney test can't reach p < 0.05, so I should aim for at least 5 samples per condition. Against a single FACS pool, at least 6 scRNA-seq samples are needed before p < 0.05 is even possible.
- **Different samples per method:** a difference between FACS and scRNA-seq can come from the method or from the different people measured (see Step 7).
- **Placeholder gates and labels:** the SSC-high/low and EpCAM-high/low gate names are placeholders: copy the real names from FlowJo into `PANELS` in SETTINGS. The Flex cluster labels are real (Step 9), but they were assigned to gates from marker genes, not from sorted cells.
- **Default Flex data:** the terminal program does not use it; there, put `scrnaseq_cells.csv` in `day04/data/` (the cached `day04/default_flex/flex_cells.csv` can be copied).
