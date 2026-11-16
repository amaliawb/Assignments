"""
explanations.py
===============

The help texts shown on the results page of the web app (app.py):
  - GENERAL_HELP : how to read p-values, q-values and box plots (top of the page)
  - SECTION_HELP : for each section A-D, the question it answers and what to keep in mind
  - COLUMN_HELP  : what each column of the summary tables means
  - PLOT_HELP    : for each plot, how to read it and what it tells you

Kept in its own file so the texts are easy to edit without touching the code.
Plain text only (no HTML); the page template adds the layout.
"""

# -----------------------------------------------------------------------------
# General notes, shown once at the top of the results page
# -----------------------------------------------------------------------------
GENERAL_HELP = [
    ("p-value",
     "How likely a difference (or correlation) at least this large would be by chance if "
     "there were no real effect. Small p = unlikely to be chance."),
    ("q-value (FDR)",
     "The p-value corrected for running many tests at once (Benjamini-Hochberg false "
     "discovery rate). Use q, not p, to decide what is significant: here q < 0.05."),
    ("Not significant",
     "Does NOT prove there is no difference. With few samples the tests have little power: "
     "with 3 samples per group a Mann-Whitney test can never reach p < 0.05. Look at the "
     "dots in the plots as well."),
    ("Box plots",
     "The box holds the middle 50% of the samples, the black line is the median, the "
     "whiskers show the spread of the rest. Every dot is one sample - with small groups "
     "the dots are more informative than the box."),
    ("Spearman vs Pearson r",
     "Both go from -1 to +1. Spearman (r of ranks) asks 'do the samples come in the same "
     "order?'; Pearson asks 'do the values fall on a straight line?'. Spearman is less "
     "affected by a single extreme sample."),
]

# -----------------------------------------------------------------------------
# One block per section: the question it answers + notes on interpretation
# -----------------------------------------------------------------------------
SECTION_HELP = {
    "A": {
        "question": "Do FACS and scRNA-seq give the same cell-type proportions for the same sample?",
        "notes": [
            "Agreement needs two things: a high correlation (samples ranked the same way) AND "
            "a bias close to 0 with narrow limits of agreement (similar actual numbers).",
            "A high correlation with a large bias means the methods follow the same trends "
            "but one systematically finds more of a population - e.g. because some cells "
            "survive sorting or capture less well in one method.",
            "Proportions are compared, not absolute numbers: each population is divided by "
            "the sum of all populations in that sample, separately for each method.",
        ],
    },
    "B": {
        "question": "Does the proportion of each population change between the conditions - "
                    "and do both methods see the same change?",
        "notes": [
            "2 conditions: Mann-Whitney U test. 3 or more: Kruskal-Wallis test for any "
            "difference, then Mann-Whitney U for every pair of conditions.",
            "The strongest result is a change that is significant (or at least in the same "
            "direction) in BOTH FACS and scRNA-seq.",
            "If only one method shows a change, it may be a technical effect of that method.",
        ],
    },
    "C": {
        "question": "Does the protein measured with antibodies on FACS match the RNA of the "
                    "same gene measured by scRNA-seq?",
        "notes": [
            "Marker pairs: CD45 - PTPRC, CD10 - MME, CD13 - ANPEP, EpCAM - EPCAM.",
            "The table shows the correlation with all populations pooled. This is mostly "
            "high because the markers differ strongly BETWEEN populations (e.g. CD45 is high "
            "only on immune cells).",
            "The correlation INSIDE each population (in C_protein_vs_rna_stats.csv) asks a "
            "harder question: do samples with more RNA also have more protein? This is often "
            "weak, because RNA and protein levels are regulated separately.",
            "A population with fewer than 10 sequenced cells in a sample gets no RNA value "
            "for that sample (not enough cells for a reliable mean).",
        ],
    },
    "D": {
        "question": "Do the populations look different in the S8 images (size, shape, ...), "
                    "and do these image features change between conditions?",
        "notes": [
            "The table compares the populations with each other (Kruskal-Wallis + every pair). "
            "The tests between conditions are in D_imaging_tests_between_conditions.csv.",
            "Big differences between populations are expected (e.g. epithelial cells are "
            "larger than immune cells) - they mainly confirm that the gates pick the right cells.",
            "The feature names are set in SETTINGS (IMAGING_FEATURES) in facs_vs_scrnaseq.py.",
        ],
    },
    "E": {
        "question": "Do our FACS and Flex results look like an independent reference single-cell "
                    "atlas of urine with known cell types?",
        "notes": [
            "The atlas cell types are translated to the same four populations (atlas_labels in "
            "SETTINGS) and turned into proportions with the same rule as our samples.",
            "The table compares each population, method and condition with the atlas. With 3 or more "
            "atlas donors: Mann-Whitney U test of our samples vs the atlas donors. Otherwise: "
            "Wilcoxon test of our samples vs the single atlas proportion.",
            "A reference atlas is usually made from healthy people. Healthy samples that match the "
            "atlas support the gating and the annotation; a disease condition that differs from it "
            "may be biology.",
            "If only ONE method differs from the atlas while the other matches it, suspect a "
            "technical bias of that method (e.g. a cell type that sorts or sequences less well).",
            "Atlases are made with other protocols (sample handling, chemistry, annotation), so "
            "small differences are expected. Look for large, consistent ones.",
        ],
    },
    "F": {
        "question": "What does each cell type look like in the S8 images - its typical size, shape "
                    "and texture - and which samples or cells are unusual?",
        "notes": [
            "The S8 images give morphological features (size, shape, how the signal is spread). "
            "They describe what cells look like, not their physiology; for that, add functional dyes.",
            "F1-F5 use the FlowJo table: one MEDIAN per sample and population. F6-F8 need the optional "
            "single-cell imaging export and show the individual cells.",
            "Outliers are found with a robust z-score: how many 'typical spreads' (based on the median "
            "and MAD) a value is away from the median. Above 3.5 counts as an outlier. Unlike a normal "
            "z-score, extreme values cannot hide themselves by inflating the spread.",
            "An outlier is not automatically an error: check it. Common causes are a staining or "
            "acquisition problem, a gating problem, debris, doublets or clumps, or real biology.",
            "Sample outliers are found among all samples of a population together. A whole condition "
            "that is shifted (e.g. all Disease samples) is shown in section D, not as outliers here.",
        ],
    },
}

# -----------------------------------------------------------------------------
# What each table column means (only the columns shown in a table are listed)
# -----------------------------------------------------------------------------
COLUMN_HELP = {
    "population": "Cell population (from the FACS gates / scRNA-seq labels).",
    "method": "FACS (BD S8) or scRNA-seq (10x Flex).",
    "marker": "Antibody marker measured on FACS (paired with its gene).",
    "feature": "S8 imaging feature.",
    "n": "Number of samples used.",
    "spearman_r": "Rank correlation, -1 to +1. Near +1 = samples come in the same order in both measurements.",
    "spearman_q_fdr": "FDR-corrected p-value of the Spearman correlation. Below 0.05 = unlikely to be chance.",
    "pearson_r": "Linear correlation, -1 to +1. Near +1 = the values lie on a straight rising line.",
    "pearson_q_fdr": "FDR-corrected p-value of the Pearson correlation.",
    "bland_altman_bias": "Average difference scRNA-seq minus FACS. Positive = scRNA-seq finds MORE of this "
                         "population; negative = fewer. 0.05 = 5 percentage points.",
    "bland_altman_lower": "Lower limit of agreement: 95% of samples have a difference above this.",
    "bland_altman_upper": "Upper limit of agreement: 95% of samples have a difference below this. "
                          "A narrow range = the methods agree sample by sample.",
    "wilcoxon_q_fdr": "Paired test (FDR-corrected): is one method systematically higher than the other "
                      "in the same samples? Below 0.05 = yes.",
    "test": "Mann-Whitney U (two groups) or Kruskal-Wallis (three or more groups at once).",
    "comparison": "Which groups were compared.",
    "p_value": "Uncorrected p-value.",
    "q_value_fdr": "p-value corrected for all tests in this table. Use this one.",
    "significant": "True if q_value_fdr is below 0.05.",
    "group": "The condition of our samples that is compared with the atlas.",
    "n_samples": "Number of our samples used (in section E: in this condition).",
    "median_samples": "Median proportion of this population in our samples (this method and condition).",
    "atlas_proportion": "Proportion of this population in the atlas, with all atlas cells pooled.",
    "difference": "median_samples minus atlas_proportion. Positive = our samples have MORE of this "
                  "population than the atlas; 0.05 = 5 percentage points.",
    "n_atlas_donors": "Number of atlas donors. With 3 or more, our samples are tested against the "
                      "spread between donors; otherwise against the single pooled atlas proportion.",
    "mean": "Average of the sample medians (the bar height in F1).",
    "sd": "Standard deviation: how much the samples differ from each other.",
    "ci95_low": "Lower end of the 95% confidence interval of the mean: the true average is very likely "
                "above this.",
    "ci95_high": "Upper end of the 95% confidence interval of the mean. A narrow interval = a precise estimate.",
    "median": "The middle value of the samples.",
    "cv_percent": "Coefficient of variation: SD as a % of the mean. Lets you compare spread between "
                  "features with different units. Below ~10% = samples are very consistent.",
    "n_outlier_samples": "How many samples are outliers for this population and feature (robust z > 3.5). "
                         "They are named in plot F1 and shown in F4.",
}

# -----------------------------------------------------------------------------
# One entry per plot file: how to read it, and what it means
# -----------------------------------------------------------------------------
PLOT_HELP = {
    "A1_proportions_stacked.png": {
        "title": "Stacked bars: proportions per sample",
        "how": [
            "Each sample has two bars next to each other: FACS (left) and scRNA-seq ('seq', right).",
            "Each colored block is one population; its height is that population's proportion. "
            "Every bar adds up to 1 (100%).",
            "Samples are sorted by condition, so you can also see shifts between conditions.",
        ],
        "meaning": "If the two bars of a sample look alike, the methods agree for that sample. "
                   "If one color is always bigger in the seq bar, that method systematically finds "
                   "more of that population.",
    },
    "A2_proportions_FACS_vs_seq.png": {
        "title": "Scatter: FACS proportion vs scRNA-seq proportion",
        "how": [
            "One panel per population; each dot is one sample.",
            "x = proportion measured by FACS, y = proportion measured by scRNA-seq.",
            "The dashed diagonal is perfect agreement (same value in both methods).",
            "The title gives the Spearman r and its q-value.",
        ],
        "meaning": "Dots ON the diagonal = the methods agree. Dots ABOVE it = scRNA-seq finds more "
                   "of this population, BELOW = fewer. Dots on a tight line that is shifted off the "
                   "diagonal = the methods rank samples the same way but with a constant offset.",
    },
    "A3_bland_altman.png": {
        "title": "Bland-Altman: how big are the differences between the methods?",
        "how": [
            "One panel per population; each dot is one sample.",
            "x = the average of the two methods, y = scRNA-seq minus FACS.",
            "Grey line = 0 (no difference). Black line = the bias (average difference). "
            "Dashed lines = limits of agreement: 95% of samples should fall between them.",
        ],
        "meaning": "Black line far from 0 = one method is systematically higher. Dashed lines far "
                   "apart = the methods can disagree a lot for a single sample. Dots that slope or "
                   "spread out from left to right = the difference depends on how common the "
                   "population is.",
    },
    "B1_proportions_by_condition.png": {
        "title": "Proportions per condition, FACS and scRNA-seq",
        "how": [
            "One panel per population; the conditions are on the x-axis.",
            "For each condition: dark grey box = FACS, light grey box = scRNA-seq. Each dot is one sample.",
            "Compare the boxes of the same color between conditions.",
        ],
        "meaning": "If the median (black line) moves in the same direction in both methods, the change "
                   "between conditions is robust. Check the table above to see whether it is "
                   "significant. Only drawn when at least two conditions are chosen.",
    },
    "C1_protein_vs_rna.png": {
        "title": "Protein (FACS) vs RNA (scRNA-seq) per marker",
        "how": [
            "One panel per marker; each dot is one sample x one population, colored by population.",
            "x = mean RNA of the gene in that population (log-normalized), y = FACS median "
            "intensity of the antibody in that population's gate.",
            "The title gives the Spearman r with all populations pooled.",
        ],
        "meaning": "The population that should carry the marker (e.g. immune cells for CD45, "
                   "proximal tubule for CD10/CD13, epithelial for EpCAM) should sit top-right: high "
                   "in both methods. Separate color clouds = the marker separates the populations in "
                   "both methods. A rising trend inside one color = protein follows RNA between samples.",
    },
    "D1_imaging_features.png": {
        "title": "S8 imaging features per population and condition",
        "how": [
            "Rows = imaging features, columns = populations, x-axis = conditions.",
            "Each box shows the values of one condition; each dot is one sample "
            "(the median of that feature inside the population's gate).",
            "Every panel has its own y-axis scale - read the numbers before comparing heights "
            "between panels.",
        ],
        "meaning": "Within one panel: does the feature change between conditions in that population? "
                   "Along one row: how different are the populations for that feature (compare the "
                   "y-axis numbers)?",
    },
    "E1_atlas_proportions.png": {
        "title": "Proportions compared with the reference atlas",
        "how": [
            "One panel per population; our conditions on the x-axis, as in plot B1 "
            "(dark grey = FACS, light grey = scRNA-seq, each dot = one sample).",
            "The pink box on the right ('Atlas') shows the atlas donors, each pink dot = one donor "
            "(only drawn when the atlas has more than one donor).",
            "The pink dashed line is the atlas proportion with all atlas cells pooled.",
        ],
        "meaning": "Boxes around the dashed line, and overlapping the pink atlas box, mean our samples "
                   "have a similar amount of this population as the reference. A box clearly above or "
                   "below the line means more or fewer cells than in the atlas. Compare the two grey "
                   "boxes: if only one of them is far from the atlas, that method may be biased.",
    },
    "E2_atlas_composition_difference.png": {
        "title": "How different is each sample's whole composition from the atlas?",
        "how": [
            "Each dot is one sample (grey) or one atlas donor (pink).",
            "y = composition difference: the percentage of cells that would have to change "
            "population for that sample to look exactly like the pooled atlas. 0% = identical, "
            "100% = nothing in common.",
            "The pink 'Atlas donors' box shows how much reference individuals differ from their own "
            "pooled atlas: this is the normal person-to-person variation.",
        ],
        "meaning": "Samples at about the same height as the atlas donors are as close to the atlas as "
                   "healthy reference people are to each other. Samples clearly higher differ more "
                   "than normal variation. Comparing FACS and scRNA-seq shows which method is closer "
                   "to the reference.",
    },
    "E3_atlas_marker_profiles.png": {
        "title": "Marker patterns: atlas RNA, Flex RNA and FACS protein",
        "how": [
            "Three heatmaps: the reference atlas (RNA), our Flex scRNA-seq (RNA) and our FACS (protein).",
            "Rows = populations, columns = markers (CD45 = PTPRC, CD10 = MME, CD13 = ANPEP, "
            "EpCAM = EPCAM). Darker = higher.",
            "Each column is scaled from 0 (the population with the lowest level) to 1 (the highest), "
            "separately in each heatmap, because RNA and antibody brightness have different units.",
            "The titles give the Spearman correlation between that heatmap and the atlas heatmap. "
            "Only drawn when the atlas has one row per cell with the four gene columns.",
        ],
        "meaning": "If the three heatmaps show the same pattern (e.g. CD45 dark only in immune, CD10/CD13 "
                   "dark only in proximal tubule), our populations carry the same markers as the atlas "
                   "cell types, which supports the gating and the annotation. A correlation close to 1 "
                   "means the patterns match. A marker that is dark in a different population than in "
                   "the atlas points to a gating or annotation problem.",
    },
    "F1_imaging_feature_means.png": {
        "title": "Mean imaging features per cell type",
        "how": [
            "One panel per imaging feature; one bar per population (all chosen samples together).",
            "Bar height = mean of the sample medians. Error bar = 95% confidence interval of that mean.",
            "Each dot is one sample. A red ring and a name mark an outlier sample (robust z > 3.5).",
        ],
        "meaning": "Shows the typical size, shape and texture of each cell type, and how much samples vary. "
                   "Non-overlapping error bars mean the populations clearly differ for that feature. A "
                   "ringed sample is unusual for that population: check its gating and image quality.",
    },
    "F2_imaging_fingerprint.png": {
        "title": "Imaging fingerprint of each cell type",
        "how": [
            "Rows = populations, columns = imaging features. The number in each cell is the mean value.",
            "The color compares the populations within each column: red = higher than the other "
            "populations, blue = lower, white = average (z-score per feature).",
        ],
        "meaning": "Each row is a population's 'morphological fingerprint'. Read across a row to describe "
                   "a cell type (e.g. large, elongated, spread-out signal). Two rows with similar colors "
                   "look alike in the images; very different rows can be told apart by imaging alone.",
    },
    "F3_imaging_pca.png": {
        "title": "Imaging map (PCA)",
        "how": [
            "Each point is one sample x one population, placed by ALL its imaging features at once "
            "(principal component analysis). Color = population, shape = condition.",
            "Points close together look alike in the images. The axis labels say how much of the "
            "differences each axis captures.",
            "Arrows show the direction in which each feature increases. Arrows pointing the same way "
            "= features that carry similar information.",
        ],
        "meaning": "Separate color clouds mean the cell types can be told apart by morphology alone, a "
                   "label-free check of the gating. Shapes of one condition shifted within a color cloud "
                   "mean the condition changes that cell type's morphology. A lone point far from its "
                   "color is an unusual sample.",
    },
    "F4_imaging_outlier_samples.png": {
        "title": "Outlier samples (robust z-score)",
        "how": [
            "Rows = samples (with their condition), columns = population x feature.",
            "Color = robust z-score of the sample's median compared with all samples of that population: "
            "red = higher than usual, blue = lower, white = typical.",
            "Cells with |z| > 3.5 are outliers and show their z-value.",
        ],
        "meaning": "A row with many dark cells is a sample that looks unusual overall (e.g. an acquisition or "
                   "staining problem). One dark cell is a single unusual population x feature. Check the "
                   "images of that sample before excluding anything.",
    },
    "F5_imaging_vs_markers.png": {
        "title": "Imaging features vs marker intensity",
        "how": [
            "One heatmap per population. Rows = imaging features, columns = FACS markers.",
            "Each number is the Spearman correlation across samples between the feature and the marker's "
            "median intensity, within that population. Red = they rise together, blue = one rises as the "
            "other falls. * = significant after FDR correction.",
        ],
        "meaning": "Links shape to marker expression, e.g. whether samples with larger proximal tubule cells "
                   "also have brighter CD13. With few samples, correlations are uncertain; trust only strong, "
                   "starred ones. Larger cells can also look brighter simply because they are bigger.",
    },
    "F6_imaging_cell_distributions.png": {
        "title": "Distribution of single cells per cell type",
        "how": [
            "One panel per feature. Each violin shows ALL single cells of a population (all chosen samples).",
            "The width shows how many cells have that value; the black line is the median.",
            "Long thin tails are rare cells with extreme values.",
        ],
        "meaning": "Shows the full spread of cells, not just the sample medians. A violin with two bulges "
                   "means the population contains two kinds of cells (e.g. two subtypes, or cells and "
                   "doublets), which may be worth gating separately.",
    },
    "F7_imaging_cell_outliers.png": {
        "title": "% outlier cells per sample",
        "how": [
            "One panel per feature, plus 'Any feature' (outlier in at least one feature).",
            "Each dot is one sample: the % of its cells in that population with a robust z > 3.5 compared "
            "with all cells of that population.",
            "With normally spread data, only a few % of cells should be outliers.",
        ],
        "meaning": "A sample with a much higher % than the others has many unusual cells: possible doublets "
                   "or clumps (often large area), debris (small area), damaged cells, or a real biological "
                   "change. Look at that sample's images.",
    },
    "F8_imaging_morphology_map.png": {
        "title": "Morphology map of single cells",
        "how": [
            "Each dot is one cell, placed by the first two imaging features (x and y). One panel per "
            "condition; color = population; at most 4000 cells are drawn per panel.",
            "Black crosses are outlier cells (robust z > 3.5 in any feature).",
        ],
        "meaning": "Shows whether the populations occupy different regions of 'shape space' and how much they "
                   "overlap. Outlier cells far from every cloud are likely doublets, clumps or debris. "
                   "Compare panels to see a condition shifting cells to a new region.",
    },
}
