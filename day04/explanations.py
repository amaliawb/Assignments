"""
explanations.py
===============

The help texts shown on the results page of the web app (app.py):
  - GENERAL_HELP : how to read p-values, q-values and box plots (top of the page)
  - SECTION_HELP : for each section A-F, the question it answers and what to keep in mind
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
     "with 3 samples per group a Mann-Whitney test can never reach p < 0.05, and 5 samples "
     "compared with ONE pooled sample can never reach p < 0.05 either (the smallest possible "
     "p is 0.0625; 6 samples are needed). Look at the dots in the plots as well."),
    ("FACS panel",
     "Which antibodies were used on FACS: CD45 + EpCAM, or CD45 + EpCAM + CD10 + CD13 (always with "
     "the Fixable Yellow live/dead dye). The panel decides the populations: both split CD45+ cells "
     "into neutrophils (high side scatter) and other immune cells; CD45- cells are split by EpCAM "
     "(high = distal tubule, collecting duct, urothelium; low = proximal tubule) and, with the "
     "4-marker panel, first by CD10/CD13 (proximal tubule). Each population is compared with the "
     "Flex cell types that would fall into its gate."),
    ("People, not samples",
     "Several samples of the same person (donor_id column of the sample sheet) are combined into one "
     "before any statistics: their cells are added up. Samples of one person are not independent, so "
     "every person counts once (e.g. 12 Flex samples from 8 people = 8 values)."),
    ("Different samples in each method",
     "FACS and scRNA-seq were measured on different samples (e.g. a pool of healthy urine, or "
     "other patients). The methods are therefore compared as GROUPS (unpaired tests). A "
     "difference between them can come from the method OR from the different people: compare "
     "groups that should be biologically alike (e.g. Healthy vs Healthy pool) to judge the method."),
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
        "question": "Do scRNA-seq and FACS find the same cell-type proportions, when each method "
                    "measured its own samples?",
        "notes": [
            "The composition overview (first table, plot A4) simply shows the average proportions "
            "side by side: each FACS group (e.g. the healthy pool, whose numbers are its own "
            "proportions), each scRNA-seq group (average over the people) and the reference atlas.",
            "The test table below it compares every chosen scRNA-seq condition with every chosen FACS "
            "condition, per population (column 'comparison'). Untick conditions on the setup page to "
            "drop comparisons you do not need.",
            "Two groups (at least 3 each): Mann-Whitney U test. A FACS POOL is one sample: the "
            "scRNA-seq people are then tested against the pool's value with a one-sample Wilcoxon "
            "test. The smallest possible p is 2/2^n for n people (4 people: 0.125, 6: 0.031, 8: "
            "0.0078), so with few people look at the size of the difference and its interval.",
            "seq_minus_facs is the difference between the group medians: 0.05 = 5 percentage "
            "points more in scRNA-seq. The 95% interval comes from a bootstrap and is only rough "
            "with few samples.",
            "Because the samples differ, a difference can be biology (different people or "
            "conditions) or a method effect (e.g. some cells survive sorting or capture less "
            "well). Comparing biologically similar groups, e.g. Healthy scRNA-seq vs a Healthy "
            "pool on FACS, is the fairest test of the methods themselves.",
            "Proportions are compared, not absolute numbers: each population is divided by "
            "the sum of all populations in that sample, separately for each method.",
            "Neutrophils have little RNA and are fragile, so scRNA-seq often finds fewer of them "
            "than FACS does - a classic method effect to keep in mind.",
            "With the CD45 + EpCAM panel, the EpCAM-low gate holds the proximal tubule cells AND any "
            "other EpCAM-negative CD45- cells (e.g. podocytes, fibroblasts); the Flex labels of those "
            "cell types are counted in the same population so both methods count the same cells.",
        ],
    },
    "B": {
        "question": "Does the proportion of each population change between the conditions "
                    "within each method - and do both methods see the same change?",
        "notes": [
            "Conditions are only compared within one method (scRNA-seq with scRNA-seq, FACS with "
            "FACS). 2 conditions: Mann-Whitney U test. 3 or more: Kruskal-Wallis test for any "
            "difference, then Mann-Whitney U for every pair of conditions.",
            "A condition needs at least 3 samples of that method to be tested, so a single FACS "
            "pool never appears here (section A compares it with the scRNA-seq samples).",
            "The strongest result is a change that is significant (or at least in the same "
            "direction) in BOTH FACS and scRNA-seq.",
        ],
    },
    "C": {
        "question": "Do the populations that carry a marker protein on FACS also express its gene "
                    "in scRNA-seq?",
        "notes": [
            "Marker pairs: CD45 - PTPRC and EpCAM - EPCAM; with the 4-marker panel also "
            "CD10 - MME and CD13 - ANPEP.",
            "The samples differ between the methods, so each population's TYPICAL level is "
            "compared: the median over the FACS samples (protein) vs the median over the "
            "scRNA-seq samples (RNA).",
            "One row per marker: the correlation across the populations. With only 4 or 5 populations "
            "the p-value can never go below 0.083 (4) or 0.017 (5), so read r as 'same pattern or not' "
            "(r = 1: the populations come in exactly the same order).",
            "'All markers (scaled)' puts all markers together after scaling each one from 0 "
            "(lowest population) to 1 (highest) in each method. With all points together it is the best "
            "single test of whether the protein and RNA patterns match.",
            "A population with fewer than 10 sequenced cells in a sample gets no RNA value "
            "for that sample (not enough cells for a reliable mean).",
        ],
    },
    "D": {
        "question": "Do the populations look different in the S8 images (size, shape, ...), "
                    "and do these image features change between conditions?",
        "notes": [
            "Only the FACS samples have images, so this section (and F) uses the FACS samples "
            "and the FACS conditions only.",
            "The table compares the populations with each other (Kruskal-Wallis + every pair). "
            "The tests between conditions are in D_imaging_tests_between_conditions.csv.",
            "Big differences between populations are expected (e.g. urothelial cells are "
            "larger than immune cells, neutrophils are more granular) - they mainly confirm that the "
            "gates pick the right cells.",
            "The feature names are set in SETTINGS (IMAGING_FEATURES) in facs_vs_scrnaseq.py.",
        ],
    },
    "E": {
        "question": "Do our FACS and Flex results look like an independent reference single-cell "
                    "atlas of urine with known cell types?",
        "notes": [
            "The atlas cell types are translated to the same populations as our data (atlas_labels in "
            "SETTINGS) and turned into proportions with the same rule as our samples.",
            "The table compares each population, method and condition with the atlas. With 3 or more "
            "samples and 3 or more atlas donors: Mann-Whitney U test of our samples vs the atlas donors. "
            "With an atlas without donors: Wilcoxon test of our samples vs the single atlas proportion. "
            "For a single sample such as a FACS pool: Wilcoxon test of the atlas donors vs the pool's value.",
            "The atlas is a shared reference for both methods, which helps when FACS and scRNA-seq "
            "measured different people: if both match the atlas, they also match each other.",
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
            "Sample outliers are found among all FACS samples of a population together. A whole condition "
            "that is shifted (e.g. all Disease samples) is shown in section D, not as outliers here.",
            "Like section D, this section uses the FACS samples only (the scRNA-seq samples have no images).",
        ],
    },
}

# -----------------------------------------------------------------------------
# What each table column means (only the columns shown in a table are listed)
# -----------------------------------------------------------------------------
COLUMN_HELP = {
    "population": "Cell population (from the FACS gates / scRNA-seq labels).",
    "method": "scRNA-seq (10x Flex) or FACS (BD S8).",
    "marker": "Antibody marker measured on FACS (paired with its gene); 'All markers (scaled)' = all "
              "markers together after scaling each from 0 to 1.",
    "feature": "S8 imaging feature.",
    "n": "Number of values used (in section C: the number of population points).",
    "spearman_r": "Rank correlation, -1 to +1. Near +1 = the values come in the same order in both measurements.",
    "spearman_q_fdr": "FDR-corrected p-value of the Spearman correlation (exact for up to 8 points). "
                      "Below 0.05 = unlikely to be chance.",
    "pearson_r": "Linear correlation, -1 to +1. Near +1 = the values lie on a straight rising line.",
    "pearson_q_fdr": "FDR-corrected p-value of the Pearson correlation.",
    "comparison": "Which groups were compared, e.g. 'Disease (scRNA-seq) vs Healthy pool (FACS)'.",
    "n_seq": "Number of scRNA-seq PEOPLE in this comparison (samples of the same person are combined).",
    "n_facs": "Number of FACS samples in this comparison (1 = a single sample, e.g. a pool).",
    "median_seq": "Median proportion of this population in the scRNA-seq samples.",
    "median_facs": "Median proportion of this population in the FACS samples (for one pool: its value).",
    "seq_minus_facs": "median_seq minus median_facs. Positive = scRNA-seq finds MORE of this population; "
                      "0.05 = 5 percentage points.",
    "seq_minus_facs_ci95_low": "Lower end of the bootstrap 95% interval of the difference.",
    "seq_minus_facs_ci95_high": "Upper end of the bootstrap 95% interval of the difference. An interval that "
                                "does not include 0 points to a real difference (rough with few samples).",
    "seq_to_facs_ratio": "median_seq divided by median_facs. 2 = twice as much in scRNA-seq, 0.5 = half.",
    "test": "Mann-Whitney U (two groups of samples), Kruskal-Wallis (three or more groups at once) or "
            "one-sample Wilcoxon (a group of samples against a single value, e.g. a pool). "
            "Empty = too few samples for a test.",
    "p_value": "Uncorrected p-value.",
    "q_value_fdr": "p-value corrected for all tests in this table. Use this one.",
    "significant": "True if q_value_fdr is below 0.05.",
    "group": "The condition of our samples that is compared with the atlas (for this method).",
    "n_samples": "Number of our samples used (in section E: of this method and condition).",
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
            "One bar per scRNA-seq PERSON (their samples combined) and per FACS sample (see the labels "
            "above the bars). Each axis label gives the person / sample and its condition. The proportions "
            "of every single sample are in A_proportions_per_sample.csv.",
            "Each colored block is one population; its height is that population's proportion. "
            "Every bar adds up to 1 (100%).",
        ],
        "meaning": "Compare the scRNA-seq bars with the FACS bars: if one color is clearly bigger in "
                   "all bars of one method, that method finds more of that population. Within one "
                   "method, you can see how much the samples differ from each other.",
    },
    "A2_proportions_by_method.png": {
        "title": "Proportions per group: scRNA-seq groups and FACS groups",
        "how": [
            "One panel per population. Each box is one group: the scRNA-seq conditions (light grey) "
            "and the FACS conditions (dark grey). Each dot is one sample.",
            "A FACS pool is a single sample, so it is drawn as one dot with a flat line.",
        ],
        "meaning": "If the scRNA-seq and FACS boxes of biologically similar groups sit at the same "
                   "height, the methods agree for that population. A gap between them is a "
                   "difference between the methods (or between the people measured).",
    },
    "A4_composition_overview.png": {
        "title": "Composition overview: FACS, scRNA-seq and the reference atlas",
        "how": [
            "One stacked bar per group: the FACS groups (e.g. the healthy pool), the scRNA-seq groups and "
            "the reference atlas. Each block is a population, with its percentage written in it.",
            "A pool is a single sample, so its bar is its own composition. For scRNA-seq and the atlas the bar "
            "is the average over the people / donors (n is written under each bar).",
        ],
        "meaning": "The quickest way to see whether the pool, the Flex data and the atlas contain the cell "
                   "types in similar amounts. A block that is clearly bigger or smaller in one bar points to a "
                   "biological difference or a method effect (e.g. neutrophils are often lower in scRNA-seq). "
                   "The spread between people is shown in A2, B1 and E1.",
    },
    "A3_method_difference.png": {
        "title": "How much do scRNA-seq and FACS differ?",
        "how": [
            "One panel per comparison (a scRNA-seq group vs a FACS group); one row per population.",
            "Dot = median scRNA-seq minus median FACS. Line = bootstrap 95% interval. "
            "The dashed line at 0 = no difference. The q-value is written on the right.",
        ],
        "meaning": "Dots to the right of 0: scRNA-seq finds more of that population; to the left: fewer. "
                   "A line that does not cross 0 points to a real difference. The same population "
                   "shifted the same way in every comparison suggests a method effect rather than biology.",
    },
    "B1_proportions_by_condition.png": {
        "title": "Proportions per condition, scRNA-seq and FACS",
        "how": [
            "One panel per population; the conditions of both methods are on the x-axis.",
            "For each condition: light grey box = scRNA-seq, dark grey box = FACS. Each dot is one "
            "sample. A method without samples in a condition has no box there.",
            "Compare the boxes of the same color between conditions.",
        ],
        "meaning": "If the median (black line) moves in the same direction in both methods, the change "
                   "between conditions is robust. Check the table above to see whether it is "
                   "significant. Only drawn when a method has at least two chosen conditions.",
    },
    "C1_protein_vs_rna.png": {
        "title": "Protein (FACS) vs RNA (scRNA-seq) per marker",
        "how": [
            "One panel per marker; each dot is one population.",
            "x = typical RNA of the gene in that population (median over the scRNA-seq samples), "
            "y = typical FACS median intensity of the antibody (median over the FACS samples). "
            "The thin lines show the 25th-75th percentile over the samples of each method.",
            "The title gives the Spearman r across the populations.",
        ],
        "meaning": "The population that should carry the marker (e.g. neutrophils and other immune cells "
                   "for CD45, proximal tubule for CD10/CD13, EpCAM-high epithelium for EpCAM) should sit top-right: high "
                   "in both methods. Dots rising from bottom-left to top-right = protein and RNA put "
                   "the populations in the same order.",
    },
    "D1_imaging_features.png": {
        "title": "S8 imaging features per population and condition",
        "how": [
            "Rows = imaging features, columns = populations, x-axis = conditions.",
            "Each box shows the FACS samples of one condition; each dot is one sample "
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
            "(light grey = scRNA-seq, dark grey = FACS, each dot = one sample).",
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
            "Rows = populations, columns = the markers of the panel (CD45 = PTPRC, CD10 = MME, CD13 = ANPEP, "
            "EpCAM = EPCAM). Darker = higher.",
            "Each column is scaled from 0 (the population with the lowest level) to 1 (the highest), "
            "separately in each heatmap, because RNA and antibody brightness have different units.",
            "The titles give the Spearman correlation between that heatmap and the atlas heatmap. "
            "Only drawn when the atlas has one row per cell with the four gene columns.",
        ],
        "meaning": "If the three heatmaps show the same pattern (e.g. CD45 dark only in the immune populations, CD10/CD13 "
                   "dark only in proximal tubule), our populations carry the same markers as the atlas "
                   "cell types, which supports the gating and the annotation. A correlation close to 1 "
                   "means the patterns match. A marker that is dark in a different population than in "
                   "the atlas points to a gating or annotation problem.",
    },
    "F1_imaging_feature_means.png": {
        "title": "Mean imaging features per cell type",
        "how": [
            "One panel per imaging feature; one bar per population (all chosen FACS samples together).",
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
            "One panel per feature. Each violin shows the single cells of a population (all chosen FACS samples; "
            "with very many cells the shape is drawn from 5000 random cells, the median from all).",
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
