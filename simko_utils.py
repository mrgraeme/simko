"""
Shared data loading, statistics, and plotting utilities for the SimKO
Streamlit app. Imported by each of the three per-layer pages (Abundance,
Expression, Mutation) so the statistical framework can't drift between them.
"""
import streamlit as st
import pandas as pd
import numpy as np
from scipy.stats import t, fisher_exact
from scipy.interpolate import UnivariateSpline
import matplotlib.pyplot as plt
import seaborn as sns


# ---------------------------------------------------------------------------
# Data loading (cached)
# ---------------------------------------------------------------------------

@st.cache_data
def get_abundance_data():
    abundance = pd.read_csv('./data/abundance.csv')
    return abundance.set_index('protein')

@st.cache_data
def get_expression_data():
    expression = pd.read_csv('./data/expression.csv')
    return expression.set_index('protein')

@st.cache_data
def get_mutation_data():
    mutation = pd.read_csv('./data/mutation.csv')
    return mutation.set_index('protein')

@st.cache_data
def get_dependency_data():
    dependency = pd.read_csv('./data/dependency.csv')
    return dependency.set_index('protein')


# ---------------------------------------------------------------------------
# Multiple-testing correction
# ---------------------------------------------------------------------------

def benjamini_hochberg(pvals):
    """Benjamini-Hochberg FDR correction, implemented directly on numpy so we
    don't need statsmodels as a dependency. Returns q-values in the same
    order as the input array. NaN p-values (e.g. degenerate cohorts) are
    treated as 1.0 so they never rank as significant.
    """
    pvals = np.asarray(pvals, dtype=float)
    pvals = np.nan_to_num(pvals, nan=1.0)
    n = len(pvals)
    if n == 0:
        return pvals
    order = np.argsort(pvals)
    ranked = pvals[order]
    ranks = np.arange(1, n + 1)
    adjusted = ranked * n / ranks
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    out = np.empty(n)
    out[order] = adjusted
    return out


def estimate_pi0(pvals, lambdas=None):
    """Estimate pi0, the proportion of tested features that are truly null,
    following Storey & Tibshirani (2003). For each candidate lambda, the
    p-values exceeding lambda should be almost entirely null features (real
    effects concentrate at small p), so counting them and scaling up
    estimates pi0; pi0(lambda) is then smoothed with a cubic spline and read
    off as lambda -> 1. Falls back to the median of the raw pi0(lambda)
    estimates if the spline fit fails (e.g. too few distinct p-values).
    """
    pvals = np.asarray(pvals, dtype=float)
    # Below ~100 features, pi0(lambda) is too noisy to trust (verified this
    # empirically -- pure-null pi0_hat bounced between ~0.5 and 1.0 at
    # n=2..50 across random seeds). Fall back to pi0=1, i.e. the same
    # (maximally conservative) assumption plain BH makes, rather than risk
    # an unstable estimate making a small scan look more permissive than
    # it should.
    if len(pvals) < 100:
        return 1.0
    if lambdas is None:
        lambdas = np.arange(0.0, 0.91, 0.05)
    pi0_hats = np.array([np.mean(pvals > lam) / (1 - lam) for lam in lambdas])
    try:
        spline = UnivariateSpline(lambdas, pi0_hats, k=3, s=len(lambdas))
        pi0 = float(spline(lambdas[-1]))
    except Exception:
        pi0 = float(np.median(pi0_hats))
    return float(np.clip(pi0, 0.0, 1.0))


def storey_qvalue(pvals, pi0=None):
    """Storey q-value: a more permissive alternative to Benjamini-Hochberg.

    Both are rank/order-statistic-based FDR procedures -- the underlying
    formula is the same shape as BH -- but standard BH implicitly assumes
    every tested feature COULD be truly null (pi0 = 1), which is the most
    conservative possible assumption. Storey's method instead estimates the
    actual proportion of null features (pi0) from the data via
    estimate_pi0(): when there is real signal (an excess of small
    p-values), pi0_hat < 1 and the correction relaxes accordingly, so more
    features pass a given threshold for the same underlying evidence. When
    there is little or no real signal, pi0_hat is close to 1 and Storey's
    q-values converge to something very close to BH's adjusted p-values.

    Pass pi0 explicitly to reuse an estimate across many similar scans
    (e.g. many KO targets) rather than re-estimating it every time.
    """
    pvals = np.asarray(pvals, dtype=float)
    pvals = np.nan_to_num(pvals, nan=1.0)
    m = len(pvals)
    if m == 0:
        return pvals
    if pi0 is None:
        pi0 = estimate_pi0(pvals)

    order = np.argsort(pvals)
    ranked = pvals[order]
    ranks = np.arange(1, m + 1)
    q = pi0 * m * ranked / ranks
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    out = np.empty(m)
    out[order] = q
    return out


# ---------------------------------------------------------------------------
# Cohort classification
# ---------------------------------------------------------------------------

def get_classes_by_joint_value(protein_list, value_df, n, mode='max'):
    """
    Partition cell lines into 'low' and 'median' cohorts from one or more KO
    targets' values in a continuous matrix (abundance or expression).

    mode='max' (default, "strict"): each target is rank-normalised to a
        percentile across cell lines, and the cohort ranking uses the MAX
        percentile across targets -- a cell line only scores as jointly low
        if it is low in EVERY selected target (AND). Use this for genuine
        double-knockout / genetic-interaction-style stratification.
    mode='mean' ("compensatory"): averages raw values across targets, so a
        cell line can qualify by being very low in one target and only
        average in another (OR-ish, partial redundancy model).

    With a single target, mode makes no difference to cohort membership,
    since percentile rank is monotonic in the raw value.
    """
    target_vals = value_df.loc[value_df.index.isin(protein_list)].T

    if mode == 'max':
        pct_ranks = target_vals[protein_list].rank(pct=True, axis=0)
        target_vals['joint_low_score'] = pct_ranks.max(axis=1)
    else:
        target_vals['joint_low_score'] = target_vals[protein_list].mean(axis=1)

    target_vals = target_vals.sort_values('joint_low_score', ascending=True)
    total = target_vals.shape[0]

    low = target_vals.head(n).copy()
    low['class'] = 'low'

    median_start = max(round(total / 2 - n / 2), 0)
    median = target_vals.iloc[median_start: median_start + n].copy()
    median['class'] = 'median'

    return pd.concat([median, low])


def get_classes_by_mutation_status(protein_list, mutation, mode='any'):
    """
    Partition cell lines into a target-mutated cohort and a target
    wild-type cohort, from non-synonymous mutation status of one or more KO
    target genes.

    Internal class labels are kept as 'low'/'median' so this plugs into the
    same downstream helpers as the abundance/expression classifiers: 'low'
    means "target(s) mutated", 'median' means "target(s) wild-type" -- it is
    not an abundance ranking. Pages using this should relabel the legend for
    the user (see class_labels on get_differentials_boxplot).

    mode='any' (default): a cell line counts as mutated if AT LEAST ONE
        selected target carries a non-synonymous mutation (OR).
    mode='all': a cell line counts as mutated only if EVERY selected target
        is mutated simultaneously (AND) -- for double-mutant /
        genetic-interaction-style questions.
    """
    target_mat = mutation.loc[mutation.index.isin(protein_list)]
    if mode == 'all':
        is_mutated = target_mat.sum(axis=0) == len(protein_list)
    else:
        is_mutated = target_mat.sum(axis=0) >= 1

    return pd.DataFrame(
        {'class': np.where(is_mutated, 'low', 'median')},
        index=target_mat.columns,
    )


# ---------------------------------------------------------------------------
# Tissue-variance annotation
# ---------------------------------------------------------------------------

def get_tissue_of(columns):
    """
    Map each cell line column name to its tissue label, using the same
    parsing convention as the rest of the app (cell line name =
    <LINE>_<TISSUE>, tissue is everything after the FIRST underscore, so
    multi-word tissues like UPPER_AERODIGESTIVE_TRACT are handled
    correctly). Returns a dict; cell lines with no underscore are omitted.
    """
    return {c: c.split('_', 1)[1] for c in columns if '_' in c}


def compute_tissue_variance_explained(data_df, tissue_of=None):
    """
    For every feature (row) in data_df, compute eta-squared from a
    one-way ANOVA against tissue-of-origin: the proportion of
    across-cell-line variance explained by which tissue a cell line comes
    from. A closed-form effect size (no permutation needed), so it can be
    computed for every protein/gene in a matrix in a couple of seconds --
    unlike the permutation-based, per-target tissue-bias validation this
    approximates (see the validation script), which is far too slow to
    run live for more than a handful of targets at a time.

    Interpretation: eta-squared near 0 means this feature's abundance is
    essentially unrelated to tissue-of-origin; near 1 means tissue
    membership explains nearly all of its variance. There's no universally
    agreed "high" cutoff -- treat it as a continuous warning signal (e.g.
    display it, let users judge, rather than hard-coding a threshold) --
    but as a loose reference, conventional ANOVA effect-size guidance
    treats ~0.01/0.06/0.14 as small/medium/large.

    tissue_of: dict/Series mapping cell line name -> tissue label; if
    None, computed automatically from data_df's own columns via
    get_tissue_of(). Cell lines not present in tissue_of are excluded.

    Returns a pandas Series indexed like data_df.index. NaN for features
    with zero total variance (e.g. entirely missing or constant rows).
    """
    if tissue_of is None:
        tissue_of = get_tissue_of(data_df.columns)

    cols_with_tissue = [c for c in data_df.columns if c in tissue_of]
    sub = data_df[cols_with_tissue]
    tissue_labels = pd.Series({c: tissue_of[c] for c in cols_with_tissue})

    overall_mean = sub.mean(axis=1)
    sst = sub.sub(overall_mean, axis=0).pow(2).sum(axis=1)

    ssb = pd.Series(0.0, index=sub.index)
    for tissue, group_cols in tissue_labels.groupby(tissue_labels).groups.items():
        group_cols = list(group_cols)
        n_t = len(group_cols)
        group_mean = sub[group_cols].mean(axis=1)
        ssb = ssb + n_t * (group_mean - overall_mean) ** 2

    with np.errstate(divide='ignore', invalid='ignore'):
        eta_sq = ssb / sst
    return eta_sq.where(sst > 0)


# ---------------------------------------------------------------------------
# Differential statistics
# ---------------------------------------------------------------------------

def ttest_from_sample_stats(row, n_low, n_median):
    """Pooled-variance, two-tailed, two-sample t-test from precomputed
    per-cohort mean/SD.

    n_low and n_median are passed separately, not assumed equal, because
    cohorts are equal-sized by construction on the Abundance/Expression
    pages but are naturally UNEQUAL on the Mutation page -- however many
    cell lines happen to carry the mutation, which can be a small minority.
    """
    if n_low < 2 or n_median < 2:
        return np.nan
    pooled_sd = np.sqrt(
        (((n_low - 1) * (row['low_std'] ** 2)) + ((n_median - 1) * (row['median_std'] ** 2)))
        / (n_low + n_median - 2)
    )
    if pooled_sd == 0:
        return 1.0
    t_stat = (row['low'] - row['median']) / (pooled_sd * np.sqrt(1 / n_low + 1 / n_median))
    return 2 * (1 - t.cdf(abs(t_stat), (n_low + n_median - 2)))


def get_correlation_scores(target_gene, data_df, active_cell_lines=None, method='spearman'):
    """
    A complementary alternative to the cohort-based differential above:
    correlates target_gene's raw value against every other gene/protein's
    raw value, across all cell lines, with no low/median split at all.

    Added because validation work found this genuinely complements the
    cohort-based approach rather than duplicating it -- it is not simply
    a simpler fallback. IMPORTANT SCOPE NOTE: this evidence is
    transcriptomic only -- ARID1A's Expression layer (CCLE RNA-seq)
    benchmarked against real HCT116 ARID1A-knockout RNA-seq (GSE101966).
    The cohort-based method recovered real, permutation-verified signal
    among genes that go UP in low-ARID1A cell lines (17 real hits vs ~7
    expected by chance, empirical p=0.02) while completely missing genes
    that go DOWN (1 hit, indistinguishable from chance); this correlation
    method showed exactly the reverse pattern (34 real DOWN hits,
    empirical p<0.0001, but only 1 real UP hit). Taking each method's
    stronger direction together captured 51 real hits against a
    chance-level null averaging ~13 (empirical p=0.0000) -- each method
    catches real signal the other structurally cannot. See
    validation/benchmark_vs_correlation.py for the full analysis.

    This function is nonetheless wired into Abundance and Dependency as
    well as Expression in the app (see get_diff_summary), on the
    assumption that the same logic -- a correlation is structurally
    different information from a cohort split, regardless of data type --
    likely generalises. That assumption is UNTESTED for Abundance/
    Dependency specifically: no proteomic or dependency ground-truth
    comparison has been run. Re-run the same benchmark against real
    proteomic knockout data (see validation/benchmark_vs_correlation.py)
    once available, to confirm or correct this before treating the
    Abundance/Dependency correlation columns as equally well-evidenced.

    Not applicable to Mutation at all -- binary status has no natural
    correlation analog.

    Returns a DataFrame indexed by gene/protein with 'corr_rho' (the
    correlation coefficient -- the effect-size column to read) and
    'corr_p' (its own nominal p-value). Like the cohort method's own p,
    this nominal p-value can be inflated by shared tissue/lineage
    structure across cell lines -- treat corr_rho as a ranking signal,
    not corr_p as a trustworthy standalone significance value.
    target_gene's own row is dropped.
    """
    if active_cell_lines is not None:
        data_df = data_df.filter(active_cell_lines)

    target_vals = data_df.loc[target_gene].values
    other = data_df.drop(index=target_gene)

    if method == 'spearman':
        target_ranks = pd.Series(target_vals).rank().values
        other_ranks = other.rank(axis=1).to_numpy()
    else:
        target_ranks = target_vals
        other_ranks = other.to_numpy()

    n = len(target_vals)
    x = target_ranks - target_ranks.mean()
    y = other_ranks - other_ranks.mean(axis=1, keepdims=True)
    cov = (y * x).sum(axis=1)
    denom = np.sqrt((y ** 2).sum(axis=1)) * np.sqrt((x ** 2).sum())
    with np.errstate(divide='ignore', invalid='ignore'):
        corr = cov / denom
        t_stat = corr * np.sqrt((n - 2) / (1 - corr ** 2))
    pvals = 2 * t.sf(np.abs(t_stat), df=n - 2)

    return pd.DataFrame({'corr_rho': corr, 'corr_p': pvals}, index=other.index)


def get_differentials(class_df, data_df, n_low=None, n_median=None,
                       prefilter_mode=None, prefilter_min=0.0, prefilter_protein_list=None):
    """Per-feature mean/SD/diff/t-test/q between the 'low' and 'median'
    cohorts in class_df, for a continuous data matrix (abundance,
    expression, or dependency). n_low/n_median default to the actual cohort
    sizes recorded in class_df if not passed explicitly -- always pass them
    explicitly if class_df's cohorts might not match data_df's columns
    exactly (e.g. after further filtering).

    prefilter_mode restricts the SET OF FEATURES the multiple-testing
    correction is computed over (i.e. shrinks m), rather than filtering
    AFTER q-values are already computed genome-wide (that's what
    select_outliers does). Three modes, with different statistical standing:

    prefilter_mode=None (default): no restriction -- correction computed
        across every feature in data_df.
    prefilter_mode='protein_list': restrict to exactly the genes/proteins in
        prefilter_protein_list -- e.g. a user's pre-specified tracked
        protein list, chosen independently of this run's p-values. This is
        the cleanest of the three: a genuinely a priori restriction, so
        m = (number of that list actually present in data_df) and the
        resulting q-values carry the full formal FDR guarantee, the same as
        restricting to a pre-registered candidate gene list in any other
        study.
    prefilter_mode='overall_range': restrict to features whose standard
        deviation ACROSS ALL CELL LINES (not just the two cohorts) is >=
        prefilter_min. Independent of cohort assignment under the null --
        the formally-justified "independent filtering" approach (Bourgon,
        Gentleman & Huber, PNAS 2010), the same logic DESeq2 uses when
        filtering on mean expression rather than fold-change.
    prefilter_mode='diff': restrict to features with |diff| >=
        prefilter_min. Less formally justified: diff is the same quantity
        driving the t-statistic's numerator, so it's correlated with p even
        under pure noise. Treat correction under this mode as indicative
        rather than a strict FDR guarantee if reporting it -- this is the
        tool's current default for the genome-wide "most changed" scan,
        given how conservative the alternatives are at typical cohort
        sizes; see Design and Implementation for the caveat as written up.
    """
    median_class = list(class_df.loc[class_df['class'] == 'median'].index)
    low_class = list(class_df.loc[class_df['class'] == 'low'].index)
    if n_low is None:
        n_low = len(low_class)
    if n_median is None:
        n_median = len(median_class)

    working_df = data_df
    if prefilter_mode == 'protein_list' and prefilter_protein_list:
        working_df = data_df.loc[data_df.index.isin(prefilter_protein_list)]
    elif prefilter_mode == 'overall_range' and prefilter_min > 0:
        overall_std = data_df.filter(median_class + low_class).std(axis=1)
        working_df = data_df.loc[overall_std >= prefilter_min]

    diff_df = pd.DataFrame(index=working_df.index)
    diff_df['median'] = working_df.filter(median_class).mean(axis=1)
    diff_df['median_std'] = working_df.filter(median_class).std(axis=1)
    diff_df['low'] = working_df.filter(low_class).mean(axis=1)
    diff_df['low_std'] = working_df.filter(low_class).std(axis=1)

    diff_df['diff'] = diff_df['low'] - diff_df['median']

    if prefilter_mode == 'diff' and prefilter_min > 0:
        diff_df = diff_df[diff_df['diff'].abs() >= prefilter_min]

    diff_df['p'] = diff_df.apply(lambda row: ttest_from_sample_stats(row, n_low, n_median), axis=1)
    # Storey's q-value by default (see storey_qvalue docstring) -- less
    # conservative than plain BH when there's real signal in the scan.
    # benjamini_hochberg() is still available if you want to report both
    # side by side for the paper's statistical-robustness comparison.
    diff_df['p_adj'] = storey_qvalue(diff_df['p'].values)
    return diff_df.drop(columns=['low_std', 'median_std']).sort_values('diff', ascending=True)


def get_mutation_differentials(class_df, mutation_df, prefilter_protein_list=None):
    """Per-gene mutation-count differential between the 'low' and 'median'
    cohorts, with a Fisher's exact test on the 2x2 (mutated/wild-type x
    cohort) table and BH correction.

    This replaces the tool's original count-only behaviour: previously the
    Mutation layer had a diff but no p-value at all, which didn't match the
    tool's own documentation (which claimed a p<0.05 threshold applied
    here). Used identically across all three pages so "Mutation" always
    means the same thing regardless of which page's cohort split produced
    it.

    prefilter_protein_list: if given, restricts the correction to exactly
    this set of genes (e.g. a pre-specified tracked protein list), so
    m = (count of that list actually present in mutation_df) rather than
    every gene in mutation_df -- the same "protein_list" idea as
    get_differentials, applied here to the Fisher's exact test instead of
    the t-test.
    """
    median_class = list(class_df.loc[class_df['class'] == 'median'].index)
    low_class = list(class_df.loc[class_df['class'] == 'low'].index)

    working_df = mutation_df
    if prefilter_protein_list:
        working_df = mutation_df.loc[mutation_df.index.isin(prefilter_protein_list)]

    low_mat = working_df.filter(low_class)
    median_mat = working_df.filter(median_class)
    n_low = low_mat.shape[1]
    n_median = median_mat.shape[1]

    diff_df = pd.DataFrame(index=working_df.index)
    diff_df['median'] = median_mat.sum(axis=1)
    diff_df['low'] = low_mat.sum(axis=1)
    diff_df['diff'] = diff_df['low'] - diff_df['median']

    pvals = np.empty(len(diff_df))
    low_counts = diff_df['low'].to_numpy()
    median_counts = diff_df['median'].to_numpy()
    for i in range(len(diff_df)):
        a = int(low_counts[i])          # low/mutated-target cohort, this gene mutated
        b = n_low - a                    # low/mutated-target cohort, this gene WT
        c = int(median_counts[i])       # median/WT-target cohort, this gene mutated
        d = n_median - c                 # median/WT-target cohort, this gene WT
        try:
            _, p = fisher_exact([[a, b], [c, d]])
        except ValueError:
            p = 1.0
        pvals[i] = p

    diff_df['p'] = pvals
    diff_df['p_adj'] = storey_qvalue(diff_df['p'].values)
    return diff_df.sort_values('diff', ascending=True)


def select_outliers(diff_df, filter_col='p', threshold=0.01, n_outliers=5, min_abs_diff=0.0, sort_by='diff'):
    """
    Filter diff_df to features passing BOTH a significance threshold
    (filter_col < threshold) AND an effect-size floor (|diff| >=
    min_abs_diff), then return up to 2*n_outliers of them.

    sort_by='diff' (default, backward-compatible): the classic head/tail
        split -- the n_outliers most negative and n_outliers most positive
        differentials among those that pass both filters.
    sort_by='significance': ignores the +/- split and returns the
        2*n_outliers features with the smallest filter_col values (i.e.
        genuinely ranked by p or q, mixing positive and negative
        differentials according to which is more significant -- this is
        what "order these by p-value" actually means, as opposed to the
        diff-magnitude ordering the tool used previously).

    min_abs_diff's sensible range depends entirely on which layer/factor is
    selected (typical |diff| is well under 1 for Abundance/Expression/
    Dependency, but can be several counts for the raw Mutation differential)
    -- there's no one sensible default across layers, so callers should
    expose this as an adjustable control rather than hardcoding a value.

    Returns (protein_list, n_pass, n_total).
    """
    total = len(diff_df)
    filtered = diff_df
    if filter_col in filtered.columns:
        filtered = filtered[filtered[filter_col] < threshold]
    if min_abs_diff > 0 and 'diff' in filtered.columns:
        filtered = filtered[filtered['diff'].abs() >= min_abs_diff]
    n_pass = len(filtered)

    if sort_by == 'significance' and filter_col in filtered.columns:
        ranked = filtered.sort_values(filter_col, ascending=True)
        proteins = list(ranked.head(2 * n_outliers).index)
    else:
        ranked = filtered.sort_values('diff', ascending=True)
        proteins = list(ranked.head(n_outliers).index) + list(ranked.tail(n_outliers).index)

    return proteins, n_pass, total


# ---------------------------------------------------------------------------
# Summary / raw table helpers
# ---------------------------------------------------------------------------

TISSUE_VARIANCE_COL = 'Tissue variance (\u03b7\u00b2)'


def get_tissue_variance_column_config(medium_cutoff=0.20, large_cutoff=0.30):
    """
    st.column_config entry for TISSUE_VARIANCE_COL, giving it a "?" help
    tooltip on the column header explaining what the value and colours
    mean. Pass this in the column_config= argument of any st.dataframe()
    call that displays a table built by get_diff_summary()/
    render_raw_table() with eta_sq set. Thresholds should match whatever
    was passed to style_results_table() for the same table, so the
    tooltip text stays accurate.
    """
    return {
        TISSUE_VARIANCE_COL: st.column_config.NumberColumn(
            TISSUE_VARIANCE_COL,
            help=(
                "Proportion of this protein's cell-line-to-cell-line variance explained by "
                "tissue-of-origin (eta-squared, one-way ANOVA) -- not a p-value, just a 0-1 "
                "effect size. High values mean pan-cancer results for this protein may be "
                "confounded by cancer type; consider the tissue filter (Global Settings) or "
                "interpret with extra caution.\n\n"
                f"Green: < {medium_cutoff:g}  \u2022  Yellow: {medium_cutoff:g}\u2013{large_cutoff:g}  "
                f"\u2022  Red: \u2265 {large_cutoff:g}  \u2022  Grey: not covered by this matrix."
            ),
            format="%.3f",
        )
    }


def _tissue_variance_style(value, medium_cutoff=0.20, large_cutoff=0.30):
    """Excel-style traffic-light conditional formatting for a tissue-variance
    (eta-squared) cell: green = low bias (<medium_cutoff), yellow = some
    bias (medium_cutoff to large_cutoff), red = heavy bias (>=large_cutoff),
    grey = not covered by this matrix (NaN).

    Defaults (0.20 / 0.30) follow Cohen's (1988) small/medium/large
    guidance as applied in biological effect-size literature -- e.g.
    Jaehne et al., Prenatal and neonatal housing conditions affect
    anxiety-like behavior in adulthood in rats..., cites Cohen 1988/1992
    as: small >0.20, medium 0.20-0.30, large >0.30. This is looser than
    Cohen's original PSYCHOLOGY-specific eta-squared convention
    (0.01/0.06/0.14, used in the first version of this function), which we
    found flagged far more of a real proteomics panel than was useful in
    practice. Neither is a threshold validated specifically for tissue
    variance in cancer cell line abundance data -- adjust freely if these
    still don't feel right for your panel.
    """
    if pd.isna(value):
        return 'background-color: #d9d9d9; color: #595959'
    elif value < medium_cutoff:
        return 'background-color: #c6efce; color: #006100'
    elif value < large_cutoff:
        return 'background-color: #ffeb9c; color: #9c6500'
    else:
        return 'background-color: #ffc7ce; color: #9c0006'


def _style_map_compat(styler, func, subset):
    """Styler.map() (pandas >=2.1) with an .applymap() fallback for older
    pandas, so this doesn't break depending on which pandas version is
    installed."""
    if hasattr(styler, 'map'):
        return styler.map(func, subset=subset)
    return styler.applymap(func, subset=subset)


def style_results_table(df, cmap, vmin=-6, vmax=6, tissue_medium_cutoff=0.20, tissue_large_cutoff=0.30):
    """
    Style a results table for display: continuous background_gradient on
    every column EXCEPT TISSUE_VARIANCE_COL (if present), and the discrete
    green/yellow/red _tissue_variance_style on that column instead. Use
    this in place of calling .style.background_gradient(...) directly
    wherever a table might include the tissue-variance column added by
    get_diff_summary() or render_raw_table() below.

    tissue_medium_cutoff/tissue_large_cutoff: pass through to
    _tissue_variance_style if you want a different green/yellow/red split
    than the default -- see that function's docstring for the reasoning
    behind the current defaults.
    """
    if TISSUE_VARIANCE_COL in df.columns:
        data_cols = [c for c in df.columns if c != TISSUE_VARIANCE_COL]
        styler = df.style.background_gradient(cmap=cmap, vmin=vmin, vmax=vmax, subset=data_cols)
        color_func = lambda v: _tissue_variance_style(v, tissue_medium_cutoff, tissue_large_cutoff)
        styler = _style_map_compat(styler, color_func, subset=[TISSUE_VARIANCE_COL])
        return styler.format("{:.3f}")
    return df.style.background_gradient(cmap=cmap, vmin=vmin, vmax=vmax, axis=None).format("{:.3f}")


def get_diff_summary(diff_abund_df, diff_exp_df, diff_mut_df, protein_list, diff_dep_df=None, eta_sq=None,
                      corr_abund_df=None, corr_exp_df=None, corr_dep_df=None):
    """
    eta_sq: optional Series (protein -> eta-squared, from
    compute_tissue_variance_explained) -- if given, inserted as the first
    column (TISSUE_VARIANCE_COL). Conventionally the ABUNDANCE-layer
    eta-squared, since that's what this tool's tissue-confounding
    validation was run on -- pass whichever layer's eta_sq is most
    relevant to how you're using this summary if you want something else.

    corr_abund_df/corr_exp_df/corr_dep_df: optional DataFrames from
    get_correlation_scores() (columns 'corr_rho'/'corr_p'), one per
    continuous layer -- if given, appended alongside that layer's
    cohort-based diff/p columns. See get_correlation_scores()'s docstring
    for why both views are shown together rather than one replacing the
    other. No Mutation equivalent -- binary status has no natural
    correlation analog.
    """
    ab_sub = diff_abund_df.loc[diff_abund_df.index.isin(protein_list)].copy()
    if corr_abund_df is not None:
        ab_sub = ab_sub.join(corr_abund_df.loc[corr_abund_df.index.isin(protein_list)])
    ab_sub.columns = ['Abundance - ' + col for col in ab_sub.columns]

    ex_sub = diff_exp_df.loc[diff_exp_df.index.isin(protein_list)].copy()
    if corr_exp_df is not None:
        ex_sub = ex_sub.join(corr_exp_df.loc[corr_exp_df.index.isin(protein_list)])
    ex_sub.columns = ['Expression - ' + col for col in ex_sub.columns]

    mu_sub = diff_mut_df.loc[diff_mut_df.index.isin(protein_list)].copy()
    mu_sub.columns = ['Mutation - ' + col for col in mu_sub.columns]

    frames = []
    if diff_dep_df is not None:
        dep_sub = diff_dep_df.loc[diff_dep_df.index.isin(protein_list)].copy()
        if corr_dep_df is not None:
            dep_sub = dep_sub.join(corr_dep_df.loc[corr_dep_df.index.isin(protein_list)])
        dep_sub.columns = ['Dependency - ' + col for col in dep_sub.columns]
        frames.append(dep_sub)  # Dependency goes first / leftmost

    frames += [ab_sub, ex_sub, mu_sub]
    result = pd.concat(frames, axis=1)

    if eta_sq is not None:
        result.insert(0, TISSUE_VARIANCE_COL, [eta_sq.get(p, np.nan) for p in result.index])

    return result


def render_raw_table(matrix, protein_subset, class_df, cmap, label, filename_prefix, eta_sq=None):
    """Render a raw-matrix dataframe (filtered to protein_subset x classified
    cell lines) plus a download button.

    eta_sq: optional Series (protein -> eta-squared, from
    compute_tissue_variance_explained on THIS SAME matrix -- i.e. abundance
    eta_sq for an abundance table, expression eta_sq for an expression
    table, etc., so the flagged bias reflects the layer actually being
    shown) -- if given, inserted as the first column with traffic-light
    styling.
    """
    subset = matrix.loc[matrix.index.isin(protein_subset)]
    subset = subset.filter(class_df.index).head(100)

    if eta_sq is not None:
        subset = subset.copy()
        subset.insert(0, TISSUE_VARIANCE_COL, [eta_sq.get(p, np.nan) for p in subset.index])

    col_config = get_tissue_variance_column_config() if eta_sq is not None else None
    st.dataframe(style_results_table(subset, cmap), column_config=col_config)
    st.download_button(
        f'Download Filtered {label}',
        subset.to_csv(index=True).encode('utf-8'),
        f"{filename_prefix}.csv",
        "text/csv",
        key=f"download_{filename_prefix}"
    )


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def get_differentials_boxplot(class_df, data_df, protein_list, diff_df=None, title=None,
                               ylabel='value', use_global_order=True, legend_labels=None):
    """
    legend_labels: optional dict, e.g. {'low': 'Mutated', 'median': 'Wild-type'},
    to relabel the legend for pages where 'low'/'median' don't mean an
    abundance ranking (i.e. the Mutation page). The underlying data values
    stay 'low'/'median' either way -- only the displayed legend text changes.

    Note: the original version of this function took an unused `n`
    parameter (cohort size), which has been dropped here since it was never
    referenced in the plotting logic.
    """
    if use_global_order:
        global_order = st.session_state.get('saved_protein_list', [])
        plot_order = [p for p in global_order if p in protein_list]
        if not plot_order:
            plot_order = protein_list
    else:
        plot_order = protein_list

    data_df = data_df.reset_index()
    data_df = data_df.loc[data_df['protein'].isin(plot_order)]
    data_df = data_df.melt(id_vars='protein')

    class_df = class_df[['class']].reset_index()
    box_data = data_df.merge(class_df, how='left', left_on='variable', right_on='index').dropna()

    ggplot2_palette = {'median': '#619CFF', 'low': '#F8766D'}

    fig = plt.figure(figsize=(10, 4))
    ax = sns.boxplot(
        data=box_data,
        x="protein",
        y="value",
        hue='class',
        order=plot_order,
        hue_order=['median', 'low'],
        palette=ggplot2_palette
    )

    plt.xticks(rotation=45)
    plt.ylabel(ylabel)
    if title:
        plt.title(title)

    ax.legend(
        title='class',
        loc='upper left',
        bbox_to_anchor=(1.01, 1),
        borderaxespad=0,
        fontsize=8,
        title_fontsize=8,
        frameon=False
    )
    if legend_labels:
        for text_obj, key in zip(ax.legend_.texts, ['median', 'low']):
            if key in legend_labels:
                text_obj.set_text(legend_labels[key])

    if diff_df is not None and 'p' in diff_df.columns and len(box_data) > 0:
        y_max = box_data['value'].max()
        y_min = box_data['value'].min()
        pad = (y_max - y_min) * 0.08 if y_max > y_min else abs(y_max) * 0.08 or 1
        ax.set_ylim(top=y_max + pad * 3)
        for i, protein in enumerate(plot_order):
            if protein in diff_df.index:
                if 'p_adj' in diff_df.columns:
                    q_val = diff_df.loc[protein, 'p_adj']
                    label = 'q<0.001' if q_val < 0.001 else f'q={q_val:.3f}'
                else:
                    p_val = diff_df.loc[protein, 'p']
                    label = 'p<0.001' if p_val < 0.001 else f'p={p_val:.3f}'
                ax.text(i, y_max + pad, label, ha='center', va='bottom', fontsize=8)

    plt.tight_layout()
    fig.subplots_adjust(right=0.82)
    return fig