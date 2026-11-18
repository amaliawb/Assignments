# What the app does

This app compares urine cell populations measured in two ways:

- **FACS** on the BD FACSDiscover S8, gated in FlowJo: it counts cells by their antibody staining and takes an image of every cell;
- **single-cell RNA-seq** with 10x Flex: it reads the genes of every cell, and each cell gets a cell-type label in the analysis notebook.

It answers questions such as: *do both methods find the same proportion of neutrophils? Does that proportion change between patient groups, and do both methods see the change? Do our results look like a published reference atlas?*

This page describes what the app does with your files, what each choice in the app means and how to read the results. It does not describe the code; for that, see [HOW_THE_CODE_WORKS.md](HOW_THE_CODE_WORKS.md).

---

## In one picture

```
  Your files                     The app                                   You get
  ──────────                     ───────                                   ───────
  Sample sheet ─────┐
  FlowJo export ────┼──► every cell → one FACS population ──► proportions ──► A  scRNA-seq vs FACS
  scRNA-seq data ───┘    (FACS: by gate; scRNA-seq: by its     per person      B  between patient groups
  (or the default        label, which you can check)                           C  protein vs RNA
   Flex data)                                                                   D  imaging features
  Reference atlas (optional) ──────────────────────────────────────────────►  E  vs the reference atlas
  Single-cell imaging export (optional) ───────────────────────────────────►  F  imaging per cell type
                                                                                + "Read this first" warnings
```

The key idea: FACS and scRNA-seq describe cells differently, so the app first puts every cell of both methods into the **same set of populations**. Then it compares the **proportions** of those populations, the share of cells, not the cell numbers, because FACS counts tens of thousands of cells and scRNA-seq a few thousand.

---

## The files you give it

| File | Needed? | What it must contain | What the app does with it |
|---|---|---|---|
| **Sample sheet** (`.csv`) | Yes | One row per sample, with: the sample name; which method measured it (FACS or scRNA-seq); for FACS samples, the name FlowJo uses; optionally the person the sample comes from; and one or more columns with groups (e.g. `condition`: Healthy pool / Disease). | Links the files together, combines samples of the same person, and offers the group columns to compare. |
| **FlowJo export** (`.xlsx`) | Yes | One row per FACS sample: the cell count of each gate, the median brightness of each marker, and the imaging features in each gate. | Takes the population counts, marker levels and imaging features from it. |
| **scRNA-seq data** | No | One of three: a cells table (`.csv`: sample, cell-type label and the 4 marker genes per cell), **or** an `.h5ad` file from scanpy (upload it, or type its path on this computer), **or** nothing. | Without it, the app uses the **default Flex data** (12 samples from 8 people, on the X: drive) and adds those samples to your sample sheet. |
| **Reference atlas** (`.csv`) | No | A published single-cell dataset of urine with known cell types, either one row per cell or one row per cell type with a cell count. | Adds section E. Without it, section E is skipped. |
| **Single-cell imaging export** (`.csv`) | No | The single FACS cells with their population and imaging features, exported from FlowJo. | Adds the cell-level plots of section F. Without it, section F uses the per-sample values only. |

**FACS and scRNA-seq don't have to come from the same samples.** For example, Flex of patients can be compared with FACS of a *pool* of healthy urine. The sample sheet says which method measured each sample. A pool of urine from several people counts as **one** FACS sample.

**Nothing leaves your computer.** The app runs only on this computer (address 127.0.0.1), and every analysis gets its own folder. An uploaded `.h5ad` file is deleted once its cells have been read; a file given as a path is only read, never changed.

---

## Going through the app

### Page 1: Upload

- Choose the sample sheet and the FlowJo export (required).
- For the scRNA-seq data, choose one:
  - leave it empty → the default Flex data;
  - upload a cells `.csv` or an `.h5ad` file;
  - type the **path** of an `.h5ad` or `.csv` file. This is quicker for big files: the Flex `.h5ad` is 4 GB.
- Optionally add the reference atlas and the single-cell imaging export.
- **Use demo data** makes a made-up dataset that shows every feature. Use it to try the app.

### Page 1b: Which columns of the `.h5ad`? (only for an `.h5ad` file)

An `.h5ad` file has many per-cell columns. The app lists them, with how many different values each has and a few examples. Columns with common names are already selected. You choose:

| Choice | Meaning | Example |
|---|---|---|
| **Sample** | Which sample each cell comes from | `sample_id` |
| **Cell type** | The cluster or annotation label of each cell, from your notebook | `leiden_atlas` (e.g. "6: myeloid") |
| **Person** (optional) | Which person each sample comes from. Samples of one person are combined. | `patient_id` |
| **Groups to compare** (optional) | Becomes the `condition` column of the sample sheet | `role` (healthy_control / donor / recipient) |
| **More sample columns** (optional) | Copied too, so you can compare them later | `timepoint` |
| **Expression values** | Which version of the gene values to use: log-normalised, not raw counts | `log1p_cp10k` |

A person or group column must have **one value per sample**. The cluster column changes from cell to cell, so it can't be used as a group, and the page tells you if you try.

### Page 2: Choose what to compare

**a) The FACS panel: which antibodies were used.**
- *CD45 + EpCAM* gives 4 populations: neutrophils, other immune cells, EpCAM-low (proximal tubule and other EpCAM-negative cells), and EpCAM-high epithelium (distal tubule, collecting duct, urothelium).
- *CD45 + EpCAM + CD10 + CD13* gives 5 populations: neutrophils, other immune cells, proximal tubule (CD10⁺CD13⁺), EpCAM-high epithelium, and other CD45⁻ cells.
- Both panels use Fixable Yellow to remove dead cells, and every live cell lands in exactly one population.

**b) Which FACS population is each scRNA-seq cell type?**

This step links the two methods. **The app doesn't cluster the scRNA-seq data and doesn't decide what the cells are.** Each cell already has a label from your notebook. Here you say which FACS gate the cells of each label would land in. The table shows, for every label:
- the number of cells and samples;
- the average RNA level of the four marker genes: PTPRC (CD45), EPCAM, MME (CD10) and ANPEP (CD13). Darker means higher; compare down a column;
- a menu with the population. It starts from the settings, and labels the settings don't know are marked **new** and start as *Other*.

Use the marker genes to check each choice. A cluster high in PTPRC is immune; high EPCAM means EpCAM-high epithelium; high MME and ANPEP means proximal tubule. *Other* leaves a label out of all proportions. Changed rows are highlighted, and your choices are remembered for the next run.

> If no label is assigned to a FACS population, scRNA-seq can't measure that population (e.g. the Flex data has no proximal tubule cluster). The app then **leaves that population out of the method comparison**, instead of comparing a false 0%. See "Populations scRNA-seq can't see" below.

**c) Which column holds the groups?** (the "condition column")

This is the column of the sample sheet that says which group each sample belongs to. Every usable column is shown with its values for each method, for example `condition`: scRNA-seq *healthy_control, donor, recipient*; FACS *Healthy pool, Disease*. A column with values for only one method is greyed out, because there would be nothing to compare.

**d) Which groups to include?** Tick the scRNA-seq groups and the FACS groups to analyse. Unticked groups are left out completely. Every ticked scRNA-seq group is compared with every ticked FACS group (section A), and the groups within each method are compared with each other (section B).

### Page 3: Results

- **Read this first** at the top: warnings about what these results can and cannot show (see below).
- Then sections A–F. Each has the question it answers, notes on how to interpret it, the main table with an explanation of every column, and the plots, each with "how to read it" and "what it tells you".
- **Download all results (.zip)** gives every table (as `.csv`, opens in Excel) and every plot (`.png`).
- **Change conditions** goes back to page 2, for example to try the other panel or change a label.

---

## What happens inside, step by step

1. **Check the files.** Every sample in the sheet must be found in its file. A missing column or sample gives a message that names it, instead of a wrong result.
2. **Put every cell in a population.**
   - FACS: each population is one FlowJo gate, and its cell count is taken from the export.
   - scRNA-seq: each cell's label is translated into a population, as set on page 2.
3. **Count proportions.** For each sample, the cells of each population are divided by the cells in all populations of the panel, so the proportions add up to 100%. Both methods use the same rule.
4. **One value per person.** Samples of the same person (e.g. three timepoints) are not independent, so their cells are added together and the person counts **once**. The 12 Flex samples become 8 people.
5. **Populations scRNA-seq can't see.** If a population has no scRNA-seq cells at all, sections A and E leave it out. The proportions of **both** methods, and of the atlas, are then calculated again out of the remaining populations, so they stay comparable. Section B still shows it for FACS.
6. **Statistics.** All tests work on ranks (non-parametric), which suits small groups.
   - Two groups of people: **Mann-Whitney U test**.
   - A group against a single value, such as one pool: **one-sample Wilcoxon test**. It asks whether the people are centred on the pool's value.
   - Three or more groups: **Kruskal-Wallis test**, then every pair.
   - Correlations: **Spearman** and **Pearson**.
   - All p-values in a table are corrected for the number of tests (**FDR**, giving q-values). "Significant" means q < 0.05.
7. **Plots and tables** for each section, plus the warnings.

---

## What you get: sections A to F

| Section | Question | Main outputs |
|---|---|---|
| **A. scRNA-seq vs FACS** | Do both methods find the same proportions? | Average composition of every group side by side (incl. the pool and the atlas); per population and comparison: the medians, their difference with a 95% interval, the ratio and the test; the label table (which population each scRNA-seq label was counted in). |
| **B. Between groups** | Does a population change between patient groups, within each method? And do both methods see the same change? | Tests per population and method; box plots per group. |
| **C. Protein vs RNA** | Are the populations high or low in the same order for the antibody (FACS brightness) and its gene (RNA)? | Correlation per marker (CD45-PTPRC, EpCAM-EPCAM, CD10-MME, CD13-ANPEP). |
| **D. Imaging features** | Do the populations look different in the S8 images, and do the images change between FACS groups? | Tests per feature (size, shape, ...). |
| **E. Reference atlas** (optional) | Do our FACS and Flex results look like an independent atlas? | Proportions vs the atlas, how different each sample's whole composition is, marker patterns of atlas RNA vs Flex RNA vs FACS protein. |
| **F. Imaging per cell type** | What does each cell type look like, and which samples or cells are unusual? | Means with confidence intervals, a "fingerprint" heatmap, a PCA map, outlier samples, imaging vs marker links; with the single-cell export also cell distributions and % odd cells (e.g. doublets). |

---

## "Read this first": what the results can and cannot show

The app checks its own results and writes plain warnings at the top of the results page (also saved as `notes.txt`):

- **A population scRNA-seq can't see.** Example: *no scRNA-seq cells in Proximal tubule*. It was left out of sections A and E. If those cells are in your data under another label, assign that label on page 2.
- **A FACS pool is one value.** Comparing 8 people with one pool tests whether the people are centred on the pool's value. How much the pool itself would vary if measured again is unknown, so this **describes** the difference; it can't prove that the methods differ.
- **Comparisons that can never be significant.** With few samples, a test has a lowest possible p-value, whatever the data. The column `smallest_possible_p` shows it for every test. For example, 4 people against one pool can never go below p = 0.125. For such comparisons, look at the size of the difference and its interval instead of the p-value.
- **Too few samples for any test.** At least 3 values are needed on one side. Otherwise only the medians are shown.

Also keep in mind:
- **"Not significant" does not mean "no difference"**, especially with few samples.
- FACS and scRNA-seq measured **different people**, so a difference can come from the method *or* from the people. Comparing groups that should be biologically alike (healthy vs a healthy pool) is the fairest test of the methods.
- The comparison is only as good as the **labels**: the app can't split a cluster that mixes two cell types, or find a cell type the notebook didn't label.

---

## All options at a glance

| Where | Option | Default |
|---|---|---|
| Page 1 | scRNA-seq data: none / cells `.csv` / `.h5ad` (upload or path) | none = the default Flex data |
| Page 1 | Reference atlas, single-cell imaging export | not given = section E skipped / no cell-level plots |
| Page 1 | Use demo data | — |
| Page 1b | Sample, cell type, person, groups, extra columns, expression values (`.h5ad` only) | guessed from the column names |
| Page 2 | FACS panel | CD45 + EpCAM + CD10 + CD13 |
| Page 2 | Population of each scRNA-seq label | from the settings; unknown labels = Other |
| Page 2 | Condition column | `condition` |
| Page 2 | Groups to include, per method | all |
| Settings file (`facs_vs_scrnaseq.py`) | Gate names, FlowJo channel names, labels per population, imaging features, minimum samples per group (3), significance threshold (0.05), the default Flex file | — |

The settings file holds what stays the same for every run, such as your FlowJo gate names. The start page of the app shows the current values.
