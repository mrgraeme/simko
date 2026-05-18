import streamlit as st

st.markdown("""
## SimKO Documentation 👓

## Overview

This tool is an interactive Streamlit application for exploring the downstream effects of simulated protein knockouts (KO) using proteomics data from the **Gygi Lab**. 
By selecting one or more proteins of interest, the tool identifies cell lines with low versus median protein abundance and tests which other proteins change significantly between those groups — providing a proxy for the consequences of depleting a given protein.

---

## Data

### Source and Citation

Proteomics abundance data was sourced via:

> Nusinow DP, Szpyt J, Ghandi M, Rose CM, McDonald ER 3rd, Kalocsay M, Jané-Valbuena J, Gelfand E, Schweppe DK, Jedrychowski M, Golji J, Porter DA, Rejtar T, Wang YK, Kryukov GV, Stegmeier F, Erickson BK, Garraway LA, Sellers WR, Gygi SP. **Quantitative Proteomics of the Cancer Cell Line Encyclopedia.** *Cell.* 2020 Jan 23;180(2):387-402.e16. doi: [10.1016/j.cell.2019.12.023](https://doi.org/10.1016/j.cell.2019.12.023)

The dataset was produced by the **Gygi Lab at Harvard Medical School** in collaboration with the Broad Institute. Normalised data files are freely available via the [Gygi Lab website](https://gygi.hms.harvard.edu/publications/ccle.html) and deposited in the [MassIVE repository](https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?task=02cd1b6a7c674f3ebdbed300b5d9aa57).

### Experimental Design

#### Mass Spectrometry Approach

Protein abundance was measured using **Tandem Mass Tag (TMT) multiplex mass spectrometry**, to allow relative quantification of proteins across multiple samples simultaneously:

- Cell lines were processed in **10-plex format**: 9 biological samples per run, plus one common reference sample used to normalise between runs.
- A total of **42 multiplex experiments** were performed, comprising **504 individual mass spectrometer runs** and over 1,500 hours of instrument time.
- On average, **over 9,000 proteins** were quantified per experiment.
- The first two multiplex experiments included **biological triplicates** to assess reproducibility.

Protein identification used a **UniProt database search**, and protein identifiers in the dataset follow UniProt annotations.

#### Normalisation

1. Within each 10-plex run, samples are normalised to the common reference channel to generate relative abundance values.
2. Cross-run normalisation is applied to make values comparable across all 42 experiments.
3. The resulting values are provided as **log2-transformed normalised protein quantities**, where positive values indicate higher-than-reference abundance and negative values indicate lower-than-reference abundance.

A detailed guide to the normalisation procedure is available in the companion preprint: [Nusinow & Gygi, 2020 (bioRxiv)](https://www.biorxiv.org/content/10.1101/2020.02.03.932384v1).
            
### Pre-processing of the Abundance Matrix

The following steps were applied to the raw Nusinow et al. normalised protein quantitation table before loading into the tool:

1. **Sparse protein filtering** — any protein row with valid (non-missing) values in fewer than 50% of cell lines was removed. This avoids unreliable group statistics for proteins detected in only a minority of the panel.

2. **Mean imputation** — for proteins passing the coverage filter, remaining missing values were imputed with the **mean abundance of that protein across all cell lines**. This is a conservative imputation strategy that pulls missing values toward the centre of the distribution, minimising their influence on differential analysis while avoiding the data loss of complete-case analysis.

> As such all proteins in the loaded abundance matrix have full coverage across cell lines, and no missing value handling is performed at runtime by the app itself.


#### Key Characteristics of the Dataset

- **375 cancer cell lines** across 22 lineages including lung, breast, colorectal, haematopoietic, and others.
- **~12,000 unique proteins** quantified across the full dataset (coverage varies per cell line).
- **Missing values** are present — not all proteins are detected in every cell line; the published normalised matrix retains these as blanks/NaN.
- **RNA–protein correlation** averages ~0.5 across proteins, meaning the proteome captures substantial post-transcriptional regulation not visible in RNA data.

---

### Data Assumptions

- Abundance and expression values are treated as **continuous, approximately normally distributed** measurements, which is appropriate given the log2-transformed normalised values provided by the dataset.
- Mean imputation of missing abundance values is appropriate for the t-test approach used here, though imputed values will slightly compress variance estimates for proteins that required substantial imputation. Proteins near the 50% missingness threshold should be interpreted with additional caution.
- Mutation values are treated as **binary counts** (mutated / not mutated) and are summed rather than averaged when comparing groups.

---

## Suitability of the Data for This Tool

The Nusinow et al. dataset is well-suited to the type of analysis this tool performs, with some important caveats worth noting.

**Breadth of coverage.** With 375 cell lines across 22 tissue types and ~9,000+ proteins per experiment, there is sufficient statistical power to meaningfully compare groups of 20 cell lines and detect differential protein abundance.

**Quantitative, normalised values.** The log2-normalised TMT values are designed for exactly this kind of between-sample comparison. The normalisation strategy corrects for run-to-run technical variation, making abundance differences interpretable as biological signal rather than measurement artefact.

**Matched multi-omic data.** Because the cell lines belong to the CCLE, matched RNA expression and mutation data are available for the same samples. This allows the tool to go beyond protein abundance and ask whether genes that change at the protein level also show concordant changes in RNA or mutation burden — a meaningful triangulation.

**Natural variation as a proxy for depletion.** The tool leverages the fact that protein abundance varies naturally across cell lines due to differences in gene expression, copy number, mutation, and post-transcriptional regulation. Cell lines at the low end of abundance for a given protein serve as a natural model for its partial depletion, which is a widely used approach in correlative proteomics analysis.

## Important Caveats

**This is not a true knockout experiment.** The low-abundance cell lines were not experimentally manipulated. Their low protein levels reflect the full complexity of their cancer biology — including co-occurring mutations, lineage effects, and other confounders. Any protein that co-varies with your protein of interest may do so for reasons unrelated to a functional dependency.

**Tissue type is a dominant source of variation.** The original paper noted a striking separation between haematopoietic/lymphoid lineages and solid organ lineages at the proteome level. Running the analysis across all tissue types simultaneously risks confounding functional signal with lineage effects. The tissue filter in the tool is strongly recommended when working with proteins of known tissue-specific biology.

**RNA–protein correlation is imperfect (~0.5 on average).** A key finding of Nusinow et al. is that protein and RNA abundances are only moderately correlated, particularly for members of protein complexes. This is actually a strength of using proteomics data — it captures regulation invisible to RNA — but it means that RNA expression changes seen in the tool should be interpreted as complementary evidence rather than a confirmation of the protein-level findings.

**Mean imputation introduces a subtle bias.** Missing values in the abundance matrix were imputed with each protein's cross-cell-line mean prior to loading. This means imputed values contribute no variance of their own — they pull group means toward the global average and slightly deflate standard deviations, which can modestly inflate t-statistics for proteins that required heavy imputation. Proteins that were close to the 50% missingness threshold (and therefore retained with a meaningful proportion of imputed values) deserve additional scrutiny in the results.

**No multiple testing correction.** See the Limitations section below.

---

## Methodology

### Step 1 — Classifying Cell Lines by Protein Abundance

For the selected protein(s), the tool:

1. Extracts the abundance values across all cell lines (optionally filtered by tissue type).
2. Calculates the **mean abundance** of the selected protein(s) across each cell line.
3. Sorts cell lines by this mean value and assigns them to two classes:

| Class | Definition |
|---|---|
| **Median** | The `n` cell lines closest to the middle of the distribution |
| **Low** | The `n` cell lines with the lowest abundance |

The number of cell lines per class (`n`) is set automatically:
- If fewer than 60 cell lines are available: `n = (total_cell_lines - 1) // 3`
- Otherwise: `n = 20`

This design is intentional — the **median group acts as a baseline**, avoiding the confounding effects that a "high vs low" comparison can introduce (e.g. comparing extremes may capture unrelated biology). Using median expressors as the reference group provides a more conservative and interpretable comparison.

### Step 2 — Differential Analysis

For each of the three data types (abundance, expression, mutation), the tool computes differences between the **low** and **median** groups:

**Abundance and Expression (continuous data):**

- Group means and standard deviations are calculated for each protein.
- A **two-sample t-test** (pooled variance, two-tailed) is performed using the formula:

$$t = \frac{\bar{x}_{low} - \bar{x}_{median}}{s_p \sqrt{\frac{1}{n} + \frac{1}{n}}}$$

where $s_p$ is the pooled standard deviation and degrees of freedom = $2n - 2$.

- The resulting **p-value** indicates whether the abundance/expression of a given protein differs significantly between the low and median cell line groups.
- Results are sorted by **fold change** (low − median), with the most depleted proteins at the top.

**Mutation (binary data):**

- Mutation counts are **summed** within each group rather than averaged.
- The **difference in mutation counts** (low − median) is reported. No significance test is applied to mutation data.

### Step 3 — Reporting

The tool reports:

- A **heatmap** of the mean abundance of selected proteins across the classified cell lines.
- A **boxplot** of abundance/expression values split by class (low vs median) for proteins of interest.
- A **summary table** combining abundance, expression, and mutation differentials for selected proteins.
- A ranked list of **top differentially abundant proteins** (filtered to p < 0.01) — those most affected in cell lines where your protein of interest is depleted.

---

## Using the Tool

### 1. Select Proteins for Simulated KO

Use the **"Proteins for KO"** multiselect to choose one or more proteins. The tool will classify cell lines based on the mean abundance of all selected proteins combined.

### 2. Filter by Tissue (Optional)

Use **"Filter for Tissue"** to restrict the analysis to specific tissue types. This is useful when you want to avoid mixing tissue-specific biology and focus on a more homogeneous cell line panel.

### 3. Explore the Cell Line Classification

The heatmap tab shows how the classified cell lines (median vs low) compare in mean abundance. Use the **"Figure Values"** tab to inspect the raw values.

### 4. Investigate Additional Proteins

Use **"Select additional proteins differences to view"** to add proteins of interest to the summary boxplot and table — for example, known interactors, pathway members, or candidates you want to check manually.

### 5. Explore Top Differentially Abundant Proteins

Adjust the **"Number of proteins by median differential"** slider to control how many top-ranked proteins (by fold change, filtered to p < 0.01) are shown in the results table.

### 6. Download Results

Two download buttons are available:

- **Download All Abundance Foldchange** — full differential abundance table for all proteins, with mean values, fold change, and p-values.
- **Download All Abundance** — the raw abundance matrix restricted to the classified cell lines, sorted by fold change.

File names are automatically tagged with the selected protein names.

---

## Interpreting Results

| Column | Meaning |
|---|---|
| `median` | Mean abundance/expression in the median cell line group |
| `low` | Mean abundance/expression in the low cell line group |
| `diff` | Difference (low − median); negative = depleted in low group |
| `p` | Two-tailed t-test p-value; values < 0.01 are highlighted in downstream tables |
| `Mutation - diff` | Difference in mutation count between groups (not a statistical test) |

A **negative `diff`** for a protein means it tends to be lower in abundance in the cell lines where your KO protein is also low — suggestive of co-dependency or co-regulation. A **positive `diff`** suggests the protein is relatively higher in those cell lines, which could indicate compensatory upregulation.

---

## Limitations and Caveats

- This analysis is **correlational**, not causal. Classifying cell lines by natural abundance variation is not equivalent to a controlled genetic knockout experiment.
- The **median group as baseline** is a pragmatic choice; it reduces but does not eliminate the risk of comparing against a biologically extreme group.
- **No multiple testing correction** is currently applied to p-values. With many proteins tested simultaneously, a threshold of p < 0.01 will still produce false positives. Treat results as hypothesis-generating rather than definitive.
- Results may vary depending on the **tissue filter** applied, as tissue type is a major source of biological variation in proteomics data.
- Mutation data is reported as raw count differences and **no statistical test** is applied; interpret with caution.

---

## Dependencies

```
streamlit
pandas
numpy
scipy
matplotlib
seaborn
```

---
            
""")