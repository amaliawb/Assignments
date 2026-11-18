# How the code works

This page explains, in plain language, what the code in this folder does and why it is built this way. You don't need to know Python to follow it. Short pieces of the real code are shown along the way, each with an explanation underneath.

For *how to run* the program, see [README.md](README.md). This page is about *how it works inside*.

---

## The big picture

The program answers one question: **do two different lab methods find the same cell populations in urine?**

- **Method 1: FACS** (BD FACSDiscover S8). Cells are stained with antibodies and sorted into populations with gates in FlowJo. FlowJo exports a table: how many cells are in each gate, and how bright each marker is.
- **Method 2: scRNA-seq** (10x Flex). Each cell's RNA is read and the cell gets a name (a "cell type") in a Jupyter notebook. The notebook exports a table with one row per cell.

The two methods **don't have to measure the same samples**. Usually they don't. For example, scRNA-seq was done on a group of patients, while FACS was done on a **pool of healthy urine** from several people, or on **other patients**. So the program compares the two methods as **groups of samples**, not sample by sample.

The program reads both tables, puts them into the same shape, compares them with statistics, and draws plots.

Optionally, you can add a third, independent source: a **reference single-cell atlas** of urine, a well-annotated dataset with known cell types (for example a published one). The program then also checks how both methods compare with that reference (section E).

```
 3 input files ──► 1. LOAD ──► 2. BUILD TABLES ──► 3. STATISTICS ──► 4. SAVE TABLES ──► 5. PLOTS
                  (read +       (count cells,       (compare methods     (.csv files)        (.png files)
                   check)        proportions)        groups and conditions)

 + optional atlas ──► E. load atlas ──► atlas proportions ──► compare with FACS and Flex ──► E tables + plots
 + imaging features ─► F. imaging profiles per cell type (+ optional single-cell export) ──► F tables + plots
```

There are two ways to use it, and both run **the same analysis code**:

- **In the terminal:** `python facs_vs_scrnaseq.py`
- **In the browser:** the website (`app.py`), which is a friendlier front door to the same steps

---

## The files and how they connect

| File | Its job | Comparable to |
|---|---|---|
| `facs_vs_scrnaseq.py` | The analysis: reading files, counting, statistics, plots | The lab protocol |
| `app.py` | The website: shows pages, receives uploads, calls the analysis | The receptionist who takes your samples and hands back results |
| `explanations.py` | The help texts shown on the results page | The figure legends |
| `templates/*.html` | The layout of each web page | Empty forms that get filled in |
| `static/style.css` | Colors, fonts, spacing of the website | The page design |
| `test_app.py` | Automatic checks that everything still works | Positive and negative controls |
| `test_business_logic.py` | Checks each calculation on tiny tables with a known answer | Standards of known concentration |

`app.py` contains **no analysis of its own**. It always calls `facs_vs_scrnaseq.py`, so the terminal and the website can never give different results.

---

## A few Python words you will meet

| Word | What it means |
|---|---|
| **variable** | A named box that holds a value. `ALPHA = 0.05` puts 0.05 in a box called `ALPHA`. |
| **list** `[ ... ]` | An ordered set of things: `["T cell", "Macrophage", "Neutrophil"]` |
| **dictionary** `{ key: value }` | A lookup table: `{"CD45": "PTPRC"}` means "CD45 → PTPRC". |
| **function** `def name(...):` | A named recipe. You give it ingredients (inputs) and it gives back a result. It only runs when it is *called*. |
| **`for` loop** | "Do this for every item": `for pop in POPULATIONS:` repeats the same steps for each population. |
| **`if`** | "Only do this when ...". |
| **DataFrame** | A table (like an Excel sheet) from the **pandas** library. Most of the work happens on these. |
| **`#` comment** | A note for humans. Python ignores everything after `#`. |
| **`"""docstring"""`** | A longer note at the start of a file or function explaining what it does. |
| **import** | Loads a library of ready-made tools: pandas (tables), numpy (maths), scipy (statistics), matplotlib (plots), Flask (websites). |

---

## Part 1: the analysis (`facs_vs_scrnaseq.py`)

The file is read from top to bottom, but most of it only *defines* functions (recipes). Nothing happens until the bottom of the file calls them. Its sections are numbered in the code itself.

### SETTINGS: everything you might change

At the top is one block of settings, so you never have to search through the code to adapt it to your experiment.

**The FACS panel.** The first choice is which antibodies were used on FACS:

```python
FACS_PANEL = "CD45_EpCAM_CD10_CD13"     # or "CD45_EpCAM"
```

Both panels use the Fixable Yellow live/dead dye. Each panel is described in a dictionary `PANELS`, which lists the panel's populations. For every population it gives:
- **`gate`**: where to find it in the FlowJo export (the gate path);
- **`seq_labels`**: which cell-type names from the scRNA-seq notebook belong to it;
- **`atlas_labels`**: the same for the reference atlas.

For example:

```python
"Neutrophils": {
    "gate": f"{LIVE_GATE}/CD45+/SSC high",
    "seq_labels": ["Neutrophil"],
    ...
},
```

This dictionary is the **bridge between the two methods**: it tells the program that the FlowJo gate `.../CD45+/SSC high` and the notebook label "Neutrophil" are the same population. The function **`use_panel(...)`** switches the whole program to one panel. It sets `POPULATIONS` (the gates and labels) and `MARKERS` (the antibodies and their genes), which every other function uses. `run_analysis(...)` calls it at the start, with the panel you chose.

The gates of each panel are designed so that **every live cell lands in exactly one population**:

```
Both panels:   Live (Fixable Yellow⁻)
               ├─ CD45+ ─┬─ SSC high                 → Neutrophils
               │         └─ SSC low                  → Other immune
               └─ CD45-  → (per panel, below)

CD45 + EpCAM:              CD45- ─┬─ EpCAM high      → EpCAM-high epithelium
                                  │                    (distal tubule, collecting duct, urothelium)
                                  └─ EpCAM low       → EpCAM-low (proximal tubule)

CD45 + EpCAM + CD10 + CD13: CD45- ─┬─ CD10+ CD13+    → Proximal tubule
                                   └─ NOT CD10+ CD13+ ─┬─ EpCAM high → EpCAM-high epithelium
                                                       └─ EpCAM low  → Other CD45-
```

This matters because scRNA-seq also gives each cell exactly one label. If a FACS cell could sit in two gates, it would be counted twice and the comparison would be unfair. For the same reason, every Flex label is listed under the gate that cell type would land in. For example, with only CD45 + EpCAM, podocytes and fibroblasts are EpCAM-low, so they are counted with the proximal tubule cells.

Other settings include:
- the marker ↔ gene pairs (`ALL_MARKERS`: CD45 ↔ PTPRC, CD10 ↔ MME, CD13 ↔ ANPEP, EpCAM ↔ EPCAM; a panel uses only its own markers);
- the imaging features;
- the statistics thresholds (`ALPHA = 0.05`, at least 3 samples per group).

### Step 1: Loading the data

**`get_column(...)`** is a small safety helper used everywhere. Before the program uses a column, it checks that the column exists. If not, it stops with a clear message that lists the real column names, so you can see the typo or the different gate name:

```
Column not found in the FACS file:
    'Cells/Single Cells/Live/CD45+ | Count'
Columns that ARE in the file:
    ...
-> Fix the matching name in the SETTINGS section of this script.
```

**`load_metadata(...)`** reads the sample sheet. It has one row per sample, and the `method` column says which method measured it (`FACS` or `scRNA-seq`). A pool of urine from several people is just one FACS row. An optional `donor_id` column says which person a sample comes from (one person can give several samples). Everything is read as text, so a sample ID like `01` doesn't turn into the number `1`. The function also:
- translates different ways of writing the method (`facs`, `Flex`, `seq`, ...) into the two standard names;
- stops with a clear message for an unknown method, for a FACS sample without a FlowJo name, or for a sample ID listed twice for one method;
- accepts an older sheet without a `method` column. Every sample is then listed once for each method, because it was measured by both.

**`choose_conditions(...)`** asks which column holds the conditions (e.g. "condition"). Then it asks **separately for each method** which conditions to include, for example scRNA-seq: Disease, and FACS: Healthy pool and Disease. In the terminal it asks you with numbered lists. The website instead passes your ticked checkboxes straight in, so nothing is asked. The result is a small dictionary: `{"scRNA-seq": ["Disease"], "FACS": ["Healthy pool", "Disease"]}`.

**`load_facs(...)`** reads the FlowJo Excel file, using the FACS rows of the sample sheet. FlowJo names samples by their `.fcs` file name, but the rest of the program uses sample IDs. The sample sheet has both, so the program builds a translation dictionary:

```python
name_to_id = dict(zip(metadata[FACS_NAME_COLUMN], metadata[SAMPLE_COLUMN]))
```

That means: "pair each FlowJo name with its sample ID", e.g. `POOL_urine.fcs → POOL`. Rows that aren't real samples (FlowJo adds rows like "Mean" and "SD") are thrown away, and you get a warning if a sample is missing.

**`load_seq(...)`** reads the per-cell scRNA-seq table and gives every cell its population, using the `seq_labels` from SETTINGS. It reads only the columns it needs (sample, cell type and the four marker genes), so a large export with many extra columns loads quickly:

```python
cells["population"] = cells[CELL_TYPE_COLUMN].map(label_to_population).fillna("Other")
```

In words: "look up each cell's label in the dictionary and write the population in a new column; if a label isn't in the dictionary, call it 'Other'". The program prints which labels ended up as "Other", so you notice if you forgot one in SETTINGS.

### Step 2: Building the comparison tables

The raw files have very different shapes: FlowJo has one row per **sample**, scRNA-seq has one row per **cell**. Step 2 turns both into the same "long" shape, with one measurement per row. That makes comparing and plotting simple. Every row belongs to one method and one of *its* samples. FACS and scRNA-seq rows are later matched by **population and condition**, never by sample.

**`build_proportions(...)`**: the main table.

- **FACS:** take each population's gate count from the FlowJo file.
- **scRNA-seq:** count the cells of each population in each sample (`pd.crosstab` makes a count table, like a pivot table in Excel).
- **Both:** divide each count by the sum of all populations in that sample. This gives a **proportion** (0 to 1), and the four proportions of a sample add up to 1.

The result looks like this:

| sample_id | method | population | count | proportion | condition |
|---|---|---|---|---|---|
| P5 | scRNA-seq | Neutrophils | 310 | 0.17 | Disease |
| POOL | FACS | Neutrophils | 5100 | 0.13 | Healthy pool |
| FACS01 | FACS | Neutrophils | 9400 | 0.24 | Disease |
| ... | ... | ... | ... | ... | ... |

Why proportions and not cell numbers? FACS counts tens of thousands of cells and scRNA-seq a few thousand, so the raw numbers can't be compared. Proportions can.

**`combine_proportions_per_donor(...)`**: samples of the **same person are combined** into one. The sample sheet's optional `donor_id` column says which person each sample comes from (e.g. 12 Flex samples from 8 people). Samples of one person are not independent, so the cells of all their samples are added up and divided again. From here on, every statistic counts each person once. The per-sample proportions are still saved, in `A_proportions_per_sample.csv`. Without a `donor_id`, every sample is its own person and nothing changes.

**`build_facs_markers(...)`** and **`build_rna_markers(...)`** make one table per method:
- **protein**: for each FACS sample, population and marker, the FACS median brightness of the antibody (e.g. CD45) in that gate;
- **RNA**: for each scRNA-seq sample, population and marker, the average expression of the matching gene (e.g. PTPRC) in that population's cells. If a population has fewer than 10 sequenced cells in a sample, its average is left empty, because an average of very few cells isn't reliable.

**`build_protein_vs_rna(...)`** puts the two side by side. The samples differ, so it uses each population's **typical** level: the median over the FACS samples (protein) and the median over the scRNA-seq samples (RNA), with the 25th–75th percentile and the number of samples behind each.

**`build_imaging(...)`**: the S8 imaging features (size, shape, ...) per FACS sample and population. Only FACS has images, so sections D and F use the FACS samples only.

### Step 3: Statistics

All tests are **non-parametric**: they work on ranks instead of assuming a bell curve. That suits small groups of samples and proportions that aren't normally distributed.

Small helper functions do the actual tests:

- **`correlate(x, y)`**: Spearman and Pearson correlation between two paired lists. It skips the test if there are fewer than 3 pairs or if all values are identical. For up to 8 pairs, the Spearman p-value is **exact**: `exact_spearman_p(...)` tries every possible way of pairing the values and counts how many give a correlation at least as strong. scipy's standard p-value is only an approximation, and it is far too small for few points (it gives p = 0 for 4 points in perfect order, while the true value is 0.083).
- **`compare_groups(...)`**:
  - 2 groups → **Mann-Whitney U** test
  - 3 or more → **Kruskal-Wallis** (is there any difference?) plus Mann-Whitney U for every pair
  - Groups that are too small are left out.
- **`compare_to_reference(values, reference)`**: the **unpaired** comparison of two groups of *different* samples:
  - both groups have ≥ 3 values → **Mann-Whitney U**;
  - one side is a single value (e.g. ONE pooled FACS sample) → **one-sample Wilcoxon**: are the other side's samples centred on that value?
  - every value identical → no test.
- **`bootstrap_median_difference(...)`**: a 95% interval for "median of group 1 − median of group 2". Both groups are resampled 2,000 times, and the middle 95% of the differences is kept. A single pool stays fixed.
- **`fdr_correct(...)`**: **Benjamini-Hochberg** correction. When you run many tests, some come out "significant" by luck. This turns each p-value into a q-value that accounts for the number of tests. The program uses q < 0.05 to call something significant.

These helpers are then used for the four questions of the analysis:

| Section | Function | Question |
|---|---|---|
| A | `composition_overview` | The average composition of every group side by side: each FACS group (a pool is its own value), each scRNA-seq group (average over people) and the reference atlas |
| A | `method_comparison_stats` | Do scRNA-seq and FACS find the same proportions? Every chosen scRNA-seq condition vs every chosen FACS condition, per population: medians, difference with bootstrap interval, ratio, and `compare_to_reference` |
| B | `condition_stats` | Does a population's proportion change between conditions, within each method? |
| C | `protein_rna_stats` | Do the populations come out in the same order for antibody brightness (FACS) and RNA of the same gene (scRNA-seq)? Per marker, plus all markers together after scaling each from 0 to 1 |
| D | `imaging_between_populations` (+ `condition_stats`) | Do the populations look different in the images, and do the image features change between the FACS conditions? |
| E | `atlas_comparison_stats` (+ `composition_difference`) | Do FACS and Flex look like the reference atlas? (only when an atlas is given, see below) |
| F | `imaging_population_summary`, `imaging_outlier_samples`, `imaging_pca`, `imaging_marker_correlations` | What does each cell type look like in the images, and which samples or cells are unusual? (see below) |

Writing the tests once as helpers and reusing them means every section is tested the same way. A fix in one place fixes it everywhere.

### Step 4: Saving the tables

Every table is saved as a `.csv` file (opens in Excel). The name starts with its section letter, e.g. `A_method_comparison_stats.csv`. The website uses that letter to put each file in the right section.

### Step 5: Plots

Each `plot_...` function draws one figure with matplotlib and saves it as a `.png` file. Some small helpers keep the plots consistent:

- `population_color(...)` gives each population the same color in every plot.
- `boxplot_with_points(...)` draws box plots with every sample as a dot on top (all boxes of a panel in one call, which is faster).
- `save(...)` writes the image to the results folder.

| File | What it shows |
|---|---|
| `A1_proportions_stacked.png` | One bar per sample: first the scRNA-seq samples, then the FACS samples |
| `A2_proportions_by_method.png` | One box per scRNA-seq group and per FACS group (a pool is one dot) |
| `A3_method_difference.png` | Per comparison and population: median scRNA-seq − median FACS, with its 95% interval and q-value |
| `A4_composition_overview.png` | One stacked bar per group (e.g. the FACS pool, each Flex condition, the reference atlas), with the % of each population written in |
| `B1_proportions_by_condition.png` | Proportions per condition, both methods |
| `C1_protein_vs_rna.png` | Per marker: each population's typical RNA (x) vs typical protein (y) |
| `D1_imaging_features.png` | Imaging features per population and condition |

### Optional step: comparing with a reference atlas (section E)

All the atlas code is in **one separate block** of the file (section "4b"). Sections A–D don't use anything from it, so adding the atlas didn't change any of their results. It only runs when a file called `reference_atlas.csv` is in the data folder, or was uploaded on the website.

**Why an atlas?** FACS and Flex are compared with each other in sections A–D. But if both had the same problem, for example a cell type that is always lost, they would still agree with each other. An atlas is an *independent* reference made by other people, so it can reveal problems that the two methods share. Because FACS and scRNA-seq usually measured different people, the atlas also gives both methods one shared yardstick.

The steps mirror the main analysis:

1. **`load_atlas(...)`** reads the atlas and gives every row a population, using `atlas_labels` from SETTINGS. This is the same idea as `seq_labels` for the notebook: a dictionary that translates names, e.g. `"Myeloid cell" → Other immune`. An atlas can come in two shapes, and the program accepts both:
   - one row per cell, so each row counts as 1 cell;
   - a summary table with one row per cell type and an `n_cells` column saying how many cells that row stands for.
2. **`build_atlas_proportions(...)`** counts cells per donor and population and divides them with **the same rule** as our samples. Proportions are only comparable if they are calculated the same way.
3. **`atlas_comparison_stats(...)`** compares, for each population, method and condition, our samples with the atlas. Which test it uses depends on what the atlas contains:
   - **3 or more atlas donors** → **Mann-Whitney U**, our samples vs the atlas donors. This takes into account that healthy people also differ from each other.
   - **an atlas without donors (or too few)** → **one-sample Wilcoxon**: are our samples spread around the single atlas proportion?
   - **a single sample, such as a FACS pool** → **one-sample Wilcoxon** the other way round: are the atlas donors spread around the pool's value?
   - This is the same `compare_to_reference(...)` helper that section A uses.
   - As everywhere, the p-values are FDR-corrected.
4. **`composition_difference(...)`** gives every sample one number: the **percentage of cells that would have to change population** for the sample to look exactly like the atlas. It is half the sum of the differences over all populations: 0% means identical, 100% means nothing in common. The atlas donors get the same score, which shows how big normal person-to-person variation is.
5. **`build_marker_profiles(...)`** works out the typical level of each marker in each population from three sources: atlas RNA, Flex RNA and FACS protein. RNA and antibody brightness have different units, so each marker is **rescaled from 0 (lowest population) to 1 (highest)**. Only the *pattern* is compared, for example "CD45 is high only in the immune populations".
6. **`run_atlas_comparison(...)`** runs these steps in order, saves the tables and draws the three plots:

| File | What it shows |
|---|---|
| `E1_atlas_proportions.png` | Our proportions per condition next to the atlas donors, with the atlas average as a dashed line |
| `E2_atlas_composition_difference.png` | How far each sample's whole composition is from the atlas, compared with the atlas donors themselves |
| `E3_atlas_marker_profiles.png` | Marker patterns as heatmaps: atlas RNA vs Flex RNA vs FACS protein (needs a per-cell atlas with gene columns) |

### Imaging profiles per cell type (section F)

All the section F code is in **its own block** of the file (section "4c"). It runs whenever imaging features are set in SETTINGS, and it uses the same "imaging" table that section D builds from the FlowJo export, so it covers the **FACS samples** only. Sections A–E don't use anything from it, so none of their results changed.

**What the S8 imaging measures.** Each cell is photographed while it is sorted, and the instrument turns each image into numbers, such as size (area), shape (eccentricity) and how the signal is spread (radial moment). These describe what cells *look like* (morphology), not what they are doing (physiology). Section F collects everything these numbers can tell about each cell type.

**Two levels of data:**
- **Sample level (always):** the FlowJo table gives one **median** per sample and population, so one number summarises all the cells of a population in a sample.
- **Cell level (optional):** if you also export the single cells (`facs_imaging_cells.csv`, one row per cell), the program can look at every cell and find unusual ones inside a sample.

**Finding outliers: the robust z-score.** A normal z-score says how many standard deviations a value is from the mean. But one extreme value pulls the mean towards itself and inflates the standard deviation, so it can "hide" itself. The **robust z-score** uses the median and the *median absolute deviation* (MAD) instead, which a few extreme values can't move:

```python
0.6745 * (value - median) / MAD
```

A value with a robust z above 3.5 (`OUTLIER_CUTOFF` in SETTINGS) counts as an outlier. The function `robust_z(...)` does this calculation, and every outlier check in section F uses it.

The functions, step by step:

1. **`imaging_outlier_samples(...)`** gives every sample a robust z-score for each population × feature, compared with all samples of that population, and marks the outliers.
2. **`imaging_population_summary(...)`** describes each population × feature over the samples: mean, SD, the **95% confidence interval** of the mean (the range the true average very likely lies in), median, **CV%** (spread as a % of the mean, so features with different units can be compared), min, max and the number of outlier samples.
3. **`imaging_pca(...)`** runs a **PCA** (principal component analysis). Each sample × population is described by all its features at once, and PCA finds the two directions that capture most of the differences, so they can be drawn on a flat map. The features are first standardized, so each one counts equally.
4. **`imaging_marker_correlations(...)`** asks, within each population and across samples, whether an imaging feature goes together with a marker's brightness (Spearman correlation, FDR-corrected).
5. **With the single-cell file:**
   - `load_facs_cells(...)` reads it (only the needed columns);
   - `cell_outlier_flags(...)` compares every cell with all cells of its population;
   - `cell_level_summary(...)` gives, per sample × population × feature: the number of cells, mean, median, SD, CV%, the 5th and 95th percentiles, and the **% of outlier cells**.
6. **`run_imaging_profiles(...)`** runs all of these, saves the tables and draws the plots:

| File | What it shows |
|---|---|
| `F1_imaging_feature_means.png` | Mean of each feature per cell type with its 95% confidence interval; each sample is a dot, and outlier samples are ringed in red and named |
| `F2_imaging_fingerprint.png` | A heatmap "fingerprint": the mean of every feature for every cell type, colored by how it compares with the other cell types |
| `F3_imaging_pca.png` | The PCA map: do cell types separate by morphology alone? Arrows show which feature pulls in which direction |
| `F4_imaging_outlier_samples.png` | The robust z-score of every sample × population × feature; outliers are labelled |
| `F5_imaging_vs_markers.png` | Correlations between imaging features and marker brightness, per cell type |
| `F6_imaging_cell_distributions.png` | *(single-cell file)* Violin plots: the shape of the distribution of the cells of each type (from at most 5,000 random cells per type, so it stays fast with big exports; the median uses all cells) |
| `F7_imaging_cell_outliers.png` | *(single-cell file)* The % of unusual cells in each sample, for example doublets, clumps or debris |
| `F8_imaging_morphology_map.png` | *(single-cell file)* Every cell as a dot on a map of the first two features, with outlier cells marked |

An outlier is a reason to **look** at that sample or those cells, not an automatic reason to remove them. It can be a technical problem, but it can also be real biology.

### Demo data

**`make_demo_data(...)`** writes made-up input files in exactly the right format. As in a real study, the methods measured different people: scRNA-seq of 12 samples from 8 people (4 Healthy, 4 Disease; some people gave two samples), and FACS of one pool of healthy urine plus 6 other Disease patients. The demo is built from individual cell types (`DEMO_CELL_TYPES`), which are added up into the gates of **both** panels, so it works with either panel. scRNA-seq is made to capture fewer neutrophils and urothelial cells, so section A has a method difference to find. The demo also has a reference atlas of 8 healthy donors and a single-cell imaging export (with about 2% odd cells, such as doublets, to find as outliers). It uses random numbers with a fixed "seed", so the demo data is the same every time. This lets you test the whole program before real data exists.

### Putting it together: `run_analysis(...)` and `main()`

**`run_analysis(...)`** is the conductor. It calls steps 1 to 5 in order, then section E if an atlas file is there, then section F if imaging features are set, and returns the summary tables. It is the **one function both the terminal and the website use**.

**`main()`** is only for the terminal. It:
1. reads the options you typed (`--demo`, `--panel`, `--condition-column`, `--seq-conditions`, `--facs-conditions`)
2. makes demo data if you asked for it
3. calls `run_analysis(...)`
4. prints a short summary

The last two lines of the file:

```python
if __name__ == "__main__":
    main()
```

mean "run `main()` only when this file is started directly". When the website *imports* the file to borrow its functions, `main()` does **not** run, so nothing is asked in a terminal.

---

## Part 2: the website (`app.py`)

The website is built with **Flask**. Each page is a Python function. A line like `@app.get("/")` above a function says "when the browser asks for this address, run this function".

### The path through the website

```
 Browser                          app.py                            facs_vs_scrnaseq.py
 ───────                          ──────                            ───────────────────
 open  /                 ──►  index()     shows upload page
 upload files / demo     ──►  upload()    saves files in runs/<id>/data ──► load_metadata (check)
 open  /run/<id>         ──►  setup()     shows condition checkboxes (scRNA-seq and FACS)
 click "Run analysis"    ──►  analyse()   ─────────────────────────────────► run_analysis(...)
 open  /run/<id>/results ──►  results()   reads the saved tables + plots,
                                          adds texts from explanations.py
 click a plot / zip      ──►  result_file() / download_all()
```

### Important ideas in `app.py`

- **Optional uploads.** The atlas and the single-cell imaging fields may be left empty. `upload()` then skips them, and `run_analysis()` simply doesn't find those files. The results page then shows a short note: "no reference atlas was given" in section E, or that the cell-level plots F6–F8 need the per-cell export in section F.
- **The FACS panel.** The setup page asks which antibodies were used (radio buttons). The choice is passed to `run_analysis(...)` and saved in `run.json`, so it stays selected when you re-run.
- **Conditions per method.** The setup page has two sets of checkboxes for the chosen column: one for the scRNA-seq samples and one for the FACS samples. Ticking none in a set means all conditions of that method.
- **A separate folder for every analysis.** Each upload gets a random ID (like `f7adccc8f906...`) and its own folder `runs/<id>/`, with `data/`, `results/` and `run.json` (the conditions you chose for each method). Two analyses can never mix up their files.
- **Errors become messages, not crashes.** When the analysis finds a problem, such as a missing column, it stops with an explanation. The helper `call_analysis(...)` catches that and shows the explanation on the page in red, so you can fix SETTINGS and try again.
- **One analysis at a time.** A "lock" makes sure two analyses don't run at the same moment, because the plotting library isn't built for that.
- **Templates.** The HTML files in `templates/` are page layouts with gaps. Flask fills the gaps with the real values: tables, plot names, your conditions. `base.html` is the frame shared by every page (header, fonts, style), and the other pages fill in the middle.
- **Safety.**
  - The server only accepts connections from this computer (`127.0.0.1`).
  - Run IDs must look exactly like a random ID.
  - Files are only served from inside a run's results folder.
  - Uploads are limited to 500 MB.

### The explanations on the results page

`explanations.py` is plain text in dictionaries: the question each section answers, what each table column means, and how to read each plot. `results()` picks the right texts for the tables and plots on the page and passes them to `templates/results.html`, which shows them under each plot. Because the texts are in their own file, you can improve the wording without touching any code.

---

## Part 3: the tests (`test_app.py`)

The tests are small programs that use the code and check the answers, like controls in an experiment. Run them with `python -m pytest test_app.py -v`. They check, among other things, that:

- a full demo analysis writes all 40 result files (19 for sections A–D, 7 for the atlas in section E, 14 for imaging in section F), with either FACS panel;
- the 12 Flex samples of the demo are combined into 8 people before the statistics;
- every scRNA-seq group is compared with every FACS group, and a single FACS pool is tested with the one-sample Wilcoxon test;
- an older sample sheet without a `method` column still works;
- without the single-cell imaging export, F6–F8 are skipped with a note; with it, they appear;
- `robust_z` flags an injected outlier but not the normal values;
- with no imaging features set, sections D and F are skipped and the rest still runs;
- without an atlas, section E is skipped with a message; with one, all its files appear;
- an atlas given as a summary table gives exactly the same proportions as the per-cell atlas it was made from;
- choosing only one condition per method skips the condition plot instead of failing;
- the terminal command still works;
- the website goes from upload → conditions → results → zip download;
- the FACS panel chosen on the website is used for the analysis, shown on the results page, and stays selected when you re-run;
- every section and plot on the results page has its explanation;
- wrong inputs (missing file, wrong file type, wrong gate name) give a clear message instead of a crash;
- made-up or unsafe addresses give "404 not found".

After any change, run the tests. If all of them pass, the change didn't break anything they cover.

A second test file, `test_business_logic.py`, checks the **calculations** one function at a time. Each test gives a function a tiny table whose correct answer was worked out by hand, like running a standard of known concentration. For example, counts of 10/20/30/40 must give proportions of 0.1/0.2/0.3/0.4. Run it with `python -m pytest test_business_logic.py -v`. `TESTS_README.md` explains every test.
