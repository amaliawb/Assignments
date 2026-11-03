# FACS vs scRNA-seq: comparing urine cell populations

Im currently doing a rotation in which the project goal is to allow single cell analysis of urine samples, specifically looking for kidney shed cells in the urine, in order to understand the health status of the kidneys in a non invasive method. 
To do this we are trying to plan an experimental system that will allow to take urine samples, sort the cells in it to enrich the sample with only kidney cells, and then sequence the RNA expression profiles of these kidney cells. 
first, we need to know whether our RNA-seq data and FACS cell sorting data have similar cell type proportions. FACS sorts cells into types based on flourencence marked antibodies that bind to cell surface marker genes. RNA seq data allows for a computational annotation of cells in samples into clusters for each cell type, defined by the RNA expression of certain marker genes. 
I want to be able to compare the results of FACS and RNA seq.

I asked Claude to build a draft for a program that will do exactly this:
A Python script that compares urine cell populations measured on the same samples with two methods:

1. **Flow cytometry / imaging sorting** on the BD FACSDiscover S8, gated in FlowJo and exported as an Excel table
2. **10x Genomics Flex single-cell RNA-seq**, processed and annotated in a Jupyter notebook and exported as a CSV

The populations and their gate logic:

| Population | Gate |
|---|---|
| Immune cells | CD45⁺ |
| Epithelial cells | CD45⁻ EpCAM⁺ |
| Proximal tubule cells | CD45⁻ CD10⁺ CD13⁺ |

## What the program does

| | Analysis | Statistics |
|---|---|---|
| A | Proportions of each population per sample, FACS vs scRNA-seq | Spearman and Pearson correlation, paired Wilcoxon signed-rank test, Bland–Altman bias and limits of agreement |
| B | Proportions compared between conditions | 2 conditions: Mann–Whitney U; 3 or more: Kruskal–Wallis plus pairwise Mann–Whitney U |
| C | Protein vs RNA: FACS median intensity vs mean RNA expression (CD45/*PTPRC*, CD10/*MME*, CD13/*ANPEP*, EpCAM/*EPCAM*) | Correlation, all populations together and each population separately |
| D | S8 imaging features per population | Comparisons between populations and between conditions |

All p-values are corrected for multiple testing (Benjamini–Hochberg FDR). The results are saved as CSV tables and PNG plots in the `results/` folder.

When it runs, the program asks which metadata column holds the conditions and which conditions to include. This means it works with any number of samples and conditions.

## How to run

Requirements: `pandas`, `numpy`, `scipy`, `matplotlib` and `openpyxl` (all included in Anaconda).

```bash
python facs_vs_scrnaseq.py --demo     # test run on made-up demo data
python facs_vs_scrnaseq.py            # run on real data in data/
python facs_vs_scrnaseq.py --condition-column condition --conditions Healthy Disease   # choose conditions without the interactive questions
```

## Input files

1. **`metadata.csv`**: one row per sample, with the columns `sample_id` and `facs_name`. `facs_name` is the sample name exactly as it appears in the FlowJo export. You can add any other columns, such as condition, donor or Flex probe barcode.
2. **FlowJo Excel export (`.xlsx`)**: one row per sample, with these statistics:
   - the count of the total gate
   - the count of each population gate
   - the median of each marker channel and imaging feature in each population
3. **`scrnaseq_cells.csv`**: one row per cell, with the columns `cell_barcode, sample_id, cell_type, PTPRC, MME, ANPEP, EPCAM`. The cells are already demultiplexed and annotated, and gene values are log-normalized. A Scanpy snippet at the top of the script exports this file.

## Settings

Everything you can change is in the **SETTINGS** section at the top of the script:
- fluorochrome and channel names
- FlowJo gate paths
- cell-type labels for each population
- imaging feature names
- how proportions are calculated (`PROPORTION_DENOMINATOR`)

If a column name doesn't match, the program stops and lists the column names it found in your file.

## Notes

- **Sample size:** with 3 samples per group, a Mann–Whitney test can't reach p < 0.05. Use at least 5 samples per condition.
- **Proportions:** by default, each population is a fraction of the three populations combined, in both methods. This keeps debris and unlabelled cells from skewing the comparison. Set `PROPORTION_DENOMINATOR = "total"` to use a fraction of all cells instead.

## How this program was created

The script was written with Claude Code (Anthropic's AI coding assistant) using the following prompts.

**Prompt 1: the request**
> in the new folder day03, i want you to create a python file with code for a program that will recieve flow cytometry data from BD S8 FACS sorter which i would use to sort urine samples with gates with antibody markers for cd45, cd10, cd13, epcam intending to sort immune cells, epithelial cells, and kidney proximal tubule cells. it will also recieve a flex 10x single cell rna seq data of urine samples. i want to be able to compare cell type populations. what more information is needed to create this program?

Claude replied with questions about:
- data formats and gate logic
- imaging
- how the scRNA-seq data is processed
- the experimental design
- which statistics to use

**Prompt 2: answers to Claude's questions**
> 1. excel from flojo, fluorochrome names: , gate logic is correct, i want to use imaging
> 2. i will use a jupiter notbook to process and annotate all the seq data. whichever data file you think is best I can extract from my jupiter notebook just choose yourself and let me know, cell type labels will be already anotated and the samples are multiplexed with probe barcodes.
> 3. FACS and seq data are from the same samples , I still dont know how many samples and conditions so I need you to allow an option to select all needed conditions once u use the program. i will have a metadata file
> 4. include proportions, protein vs rna, and p value tests and any other statistical test u think is needed
> 5. i still dont have data so create a simple script and include many comments that will explain what each section of the code does and which parameters can be changed

Claude then:
- wrote the script
- generated demo data and ran the full analysis on it
- tested the interactive condition selection

## Status

The script has only been tested on made-up demo data so far. Before it runs on real data, the channel names, gate paths and imaging feature names in SETTINGS need to be replaced; they are placeholders for now.
