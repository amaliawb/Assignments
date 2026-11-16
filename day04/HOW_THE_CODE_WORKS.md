# How the code works

This page explains, in plain language, what the code in this folder does and why it is built this way. You don't need to know Python to follow it. Short pieces of the real code are shown along the way, each with an explanation underneath.

For *how to run* the program, see [README.md](README.md). This page is about *how it works inside*.

---

## The big picture

The program answers one question: **do two different lab methods count the same cells in the same urine samples?**

- **Method 1: FACS** (BD FACSDiscover S8). Cells are stained with antibodies and sorted into populations with gates in FlowJo. FlowJo exports a table: how many cells are in each gate, and how bright each marker is.
- **Method 2: scRNA-seq** (10x Flex). Each cell's RNA is read and the cell gets a name (a "cell type") in a Jupyter notebook. The notebook exports a table with one row per cell.

The program reads both tables, puts them into the same shape, compares them with statistics, and draws plots.

Optionally, you can add a third, independent source: a **reference single-cell atlas** of urine, a well-annotated dataset with known cell types (for example a published one). The program then also checks how both methods compare with that reference (section E).

```
 3 input files ──► 1. LOAD ──► 2. BUILD TABLES ──► 3. STATISTICS ──► 4. SAVE TABLES ──► 5. PLOTS
                  (read +       (count cells,       (compare methods     (.csv files)        (.png files)
                   check)        proportions)        and conditions)

 + optional atlas ──► E. load atlas ──► atlas proportions ──► compare with FACS and Flex ──► E tables + plots
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

At the top is one block of settings, so you never have to search through the code to adapt it to your experiment. For example:

```python
LIVE_GATE = "Cells/Single Cells/Live"     # live/dead gate on Fixable Yellow

POPULATIONS = {
    "Immune": {
        "gate": f"{LIVE_GATE}/CD45+",
        "seq_labels": ["T cell", "Macrophage", "Neutrophil"],
    },
    ...
}
```

This is a dictionary that says, for every population:
- **`gate`**: where to find it in the FlowJo export (the gate path).
- **`seq_labels`**: which cell-type names from the scRNA-seq notebook belong to it.

This dictionary is the **bridge between the two methods**: it tells the program that the FlowJo gate `.../CD45+` and the notebook labels "T cell", "Macrophage" and "Neutrophil" are the same population.

The gates are designed so that **every live cell lands in exactly one population**:

```
Live (Fixable Yellow⁻)
├─ CD45+                          → Immune
└─ CD45-
   ├─ CD10+ CD13+                 → Proximal tubule
   └─ NOT (CD10+ CD13+)
      ├─ EpCAM+                   → Epithelial
      └─ EpCAM-                   → Other CD45-
```

This matters because scRNA-seq also gives each cell exactly one label. If a FACS cell could sit in two gates, it would be counted twice and the comparison would be unfair.

Other settings include the marker ↔ gene pairs (`MARKERS`: CD45 ↔ PTPRC, CD10 ↔ MME, CD13 ↔ ANPEP, EpCAM ↔ EPCAM), the imaging features, and the statistics thresholds (`ALPHA = 0.05`, at least 3 samples per group).

### Step 1: Loading the data

**`get_column(...)`** is a small safety helper used everywhere. Before the program uses a column, it checks that the column exists. If not, it stops with a clear message that lists the real column names, so you can see the typo or the different gate name:

```
Column not found in the FACS file:
    'Cells/Single Cells/Live/CD45+ | Count'
Columns that ARE in the file:
    ...
-> Fix the matching name in the SETTINGS section of this script.
```

**`load_metadata(...)`** reads the sample sheet. Everything is read as text, so a sample ID like `01` doesn't turn into the number `1`.

**`choose_conditions(...)`** asks which column holds the conditions (e.g. "condition") and which ones to include (e.g. Healthy and Disease). In the terminal it asks you with numbered lists. The website instead passes your ticked checkboxes straight in, so nothing is asked.

**`load_facs(...)`** reads the FlowJo Excel file. FlowJo names samples by their `.fcs` file name, but the scRNA-seq file uses a sample ID. The sample sheet has both, so the program builds a translation dictionary:

```python
name_to_id = dict(zip(metadata[FACS_NAME_COLUMN], metadata[SAMPLE_COLUMN]))
```

That means: "pair each FlowJo name with its sample ID", e.g. `S01_urine.fcs → S01`. Rows that aren't real samples (FlowJo adds rows like "Mean" and "SD") are thrown away, and you get a warning if a sample is missing.

**`load_seq(...)`** reads the per-cell scRNA-seq table and gives every cell its population, using the `seq_labels` from SETTINGS:

```python
cells["population"] = cells[CELL_TYPE_COLUMN].map(label_to_population).fillna("Other")
```

In words: "look up each cell's label in the dictionary and write the population in a new column; if a label isn't in the dictionary, call it 'Other'". The program prints which labels ended up as "Other", so you notice if you forgot one in SETTINGS.

### Step 2: Building the comparison tables

The raw files have very different shapes: FlowJo has one row per **sample**, scRNA-seq has one row per **cell**. Step 2 turns both into the same "long" shape, with one measurement per row. That makes comparing and plotting simple.

**`build_proportions(...)`**: the main table.

- **FACS:** take each population's gate count from the FlowJo file.
- **scRNA-seq:** count the cells of each population in each sample (`pd.crosstab` makes a count table, like a pivot table in Excel).
- **Both:** divide each count by the sum of all populations in that sample. This gives a **proportion** (0 to 1), and the four proportions of a sample add up to 1.

The result looks like this:

| sample_id | population | method | count | proportion | condition |
|---|---|---|---|---|---|
| S01 | Immune | FACS | 5100 | 0.24 | Healthy |
| S01 | Immune | scRNA-seq | 310 | 0.21 | Healthy |
| ... | ... | ... | ... | ... | ... |

Why proportions and not cell numbers? FACS counts tens of thousands of cells and scRNA-seq a few thousand, so the raw numbers can't be compared. Proportions can.

**`build_protein_vs_rna(...)`**: for each sample, population and marker, it puts side by side:
- **protein**: the FACS median brightness of the antibody (e.g. CD45) in that gate
- **RNA**: the average expression of the matching gene (e.g. PTPRC) in that population's cells

If a population has fewer than 10 sequenced cells in a sample, its RNA average is left empty, because an average of very few cells isn't reliable.

**`build_imaging(...)`**: the same idea for the S8 imaging features (size, shape, ...), taken from the FACS file only.

### Step 3: Statistics

All tests are **non-parametric**: they work on ranks instead of assuming a bell curve. That suits small groups of samples and proportions that aren't normally distributed.

Small helper functions do the actual tests:

- **`correlate(x, y)`**: Spearman and Pearson correlation between two paired lists. It skips the test if there are fewer than 3 pairs or if all values are identical.
- **`compare_groups(...)`**:
  - 2 groups → **Mann-Whitney U** test
  - 3 or more → **Kruskal-Wallis** (is there any difference?) plus Mann-Whitney U for every pair
  - Groups that are too small are left out.
- **`fdr_correct(...)`**: **Benjamini-Hochberg** correction. When you run many tests, some come out "significant" by luck. This turns each p-value into a q-value that accounts for the number of tests. The program uses q < 0.05 to call something significant.

These helpers are then used for the four questions of the analysis:

| Section | Function | Question |
|---|---|---|
| A | `method_agreement_stats` | Do FACS and scRNA-seq give the same proportions? (correlation, Bland-Altman bias and limits of agreement, paired Wilcoxon test) |
| B | `condition_stats` | Does a population's proportion change between conditions, in each method? |
| C | `protein_rna_stats` | Does antibody brightness follow RNA expression of the same gene? |
| D | `imaging_between_populations` (+ `condition_stats`) | Do the populations look different in the images, and do the image features change between conditions? |
| E | `atlas_comparison_stats` (+ `composition_difference`) | Do FACS and Flex look like the reference atlas? (only when an atlas is given, see below) |

Writing the tests once as helpers and reusing them means every section is tested the same way. A fix in one place fixes it everywhere.

### Step 4: Saving the tables

Every table is saved as a `.csv` file (opens in Excel). The name starts with its section letter, e.g. `A_method_agreement_stats.csv`. The website uses that letter to put each file in the right section.

### Step 5: Plots

Each `plot_...` function draws one figure with matplotlib and saves it as a `.png` file. Some small helpers keep the plots consistent:

- `population_color(...)` gives each population the same color in every plot.
- `boxplot_with_points(...)` draws a box plot with every sample as a dot on top.
- `save(...)` writes the image to the results folder.

| File | What it shows |
|---|---|
| `A1_proportions_stacked.png` | Each sample's FACS bar next to its scRNA-seq bar |
| `A2_proportions_FACS_vs_seq.png` | FACS vs scRNA-seq proportion, one dot per sample, with the line of perfect agreement |
| `A3_bland_altman.png` | The difference between the methods for each sample |
| `B1_proportions_by_condition.png` | Proportions per condition, both methods |
| `C1_protein_vs_rna.png` | Protein vs RNA for each marker |
| `D1_imaging_features.png` | Imaging features per population and condition |

### Optional step: comparing with a reference atlas (section E)

All the atlas code is in **one separate block** of the file (section "4b"). Sections A–D don't use anything from it, so adding the atlas didn't change any of their results. It only runs when a file called `reference_atlas.csv` is in the data folder, or was uploaded on the website.

**Why an atlas?** FACS and Flex are compared with each other in sections A–D. But if both had the same problem, for example a cell type that is always lost, they would still agree with each other. An atlas is an *independent* reference made by other people, so it can reveal problems that the two methods share.

The steps mirror the main analysis:

1. **`load_atlas(...)`** reads the atlas and gives every row a population, using `atlas_labels` from SETTINGS. This is the same idea as `seq_labels` for the notebook: a dictionary that translates names, e.g. `"Myeloid cell" → Immune`. An atlas can come in two shapes, and the program accepts both:
   - one row per cell, so each row counts as 1 cell;
   - a summary table with one row per cell type and an `n_cells` column saying how many cells that row stands for.
2. **`build_atlas_proportions(...)`** counts cells per donor and population and divides them with **the same rule** as our samples. Proportions are only comparable if they are calculated the same way.
3. **`atlas_comparison_stats(...)`** compares, for each population, method and condition, our samples with the atlas. Which test it uses depends on what the atlas contains:
   - **3 or more atlas donors** → **Mann-Whitney U**, our samples vs the atlas donors. This takes into account that healthy people also differ from each other.
   - **otherwise** → **one-sample Wilcoxon**: are our samples spread around the single atlas proportion?
   - As everywhere, the p-values are FDR-corrected.
4. **`composition_difference(...)`** gives every sample one number: the **percentage of cells that would have to change population** for the sample to look exactly like the atlas. It is half the sum of the differences over the four populations: 0% means identical, 100% means nothing in common. The atlas donors get the same score, which shows how big normal person-to-person variation is.
5. **`build_marker_profiles(...)`** works out the typical level of each marker in each population from three sources: atlas RNA, Flex RNA and FACS protein. RNA and antibody brightness have different units, so each marker is **rescaled from 0 (lowest population) to 1 (highest)**. Only the *pattern* is compared, for example "CD45 is high only in immune cells".
6. **`run_atlas_comparison(...)`** runs these steps in order, saves the tables and draws the three plots:

| File | What it shows |
|---|---|
| `E1_atlas_proportions.png` | Our proportions per condition next to the atlas donors, with the atlas average as a dashed line |
| `E2_atlas_composition_difference.png` | How far each sample's whole composition is from the atlas, compared with the atlas donors themselves |
| `E3_atlas_marker_profiles.png` | Marker patterns as heatmaps: atlas RNA vs Flex RNA vs FACS protein (needs a per-cell atlas with gene columns) |

### Demo data

**`make_demo_data(...)`** writes made-up input files in exactly the right format: 10 samples (5 Healthy and 5 Disease), plus a reference atlas of 8 healthy donors. It uses random numbers with a fixed "seed", so the demo data is the same every time. This lets you test the whole program before real data exists.

### Putting it together: `run_analysis(...)` and `main()`

**`run_analysis(...)`** is the conductor. It calls steps 1 to 5 in order, then section E if an atlas file is there, and returns the summary tables. It is the **one function both the terminal and the website use**.

**`main()`** is only for the terminal. It:
1. reads the options you typed (`--demo`, `--condition-column`, `--conditions`)
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
 open  /run/<id>         ──►  setup()     shows condition checkboxes
 click "Run analysis"    ──►  analyse()   ─────────────────────────────────► run_analysis(...)
 open  /run/<id>/results ──►  results()   reads the saved tables + plots,
                                          adds texts from explanations.py
 click a plot / zip      ──►  result_file() / download_all()
```

### Important ideas in `app.py`

- **Optional uploads.** The atlas upload field may be left empty. `upload()` then skips it, `run_analysis()` finds no atlas file, and the results page shows a short "no reference atlas was given" note in section E.
- **A separate folder for every analysis.** Each upload gets a random ID (like `f7adccc8f906...`) and its own folder `runs/<id>/`, with `data/`, `results/` and `run.json` (the conditions you chose). Two analyses can never mix up their files.
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

- a full demo analysis writes all 21 result files (14 for sections A–D and 7 for the atlas, section E);
- without an atlas, section E is skipped with a message; with one, all its files appear;
- an atlas given as a summary table gives exactly the same proportions as the per-cell atlas it was made from;
- choosing only one condition skips the condition plot instead of failing;
- the terminal command still works;
- the website goes from upload → conditions → results → zip download;
- every section and plot on the results page has its explanation;
- wrong inputs (missing file, wrong file type, wrong gate name) give a clear message instead of a crash;
- made-up or unsafe addresses give "404 not found".

After any change, run the tests. If all of them pass, the change didn't break anything they cover.
