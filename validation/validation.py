"""
SimKO validation script.

Standalone, importable functions for validating SimKO's outputs -- run from a
script, a notebook, or the __main__ block at the bottom of this file. Not a
Streamlit page: nothing here depends on session state, so it can be run
completely outside the app.

Organised by phase, matching the validation plan. Currently implemented:

    Phase 0 -- Statistical foundation
        Permutation-based empirical FDR: shuffle which cell lines are
        'low' vs 'median' many times, recompute every feature's p-value
        under each shuffled (null) labelling, and compare the real hit
        counts against what pure chance produces. Answers directly: "at
        nominal p < X, how many of my hits would I expect to see even if
        there were no real relationship at all?"

Later phases (essential-gene positive/negative controls, Perturb-seq
concordance, paralog synthetic-lethality recovery, etc.) will be added here
as they're built -- see the validation plan document for the full roadmap.

No additional data is required for Phase 0: it only resamples your existing
abundance/expression/mutation/dependency matrices. Later phases will need
external reference data (documented when we get there).

This script lives inside validation/, one level below the app root where
simko_utils.py and data/ actually sit. All default paths below are resolved
relative to THIS FILE's own location, not the current working directory --
so it can be run as `python validate_simko.py` from inside validation/,
`python validation/validate_simko.py` from the app root, or imported from a
notebook anywhere, and will find simko_utils.py and the data either way.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import t as t_dist
import matplotlib
matplotlib.use('Agg')  # headless -- this script may run outside a GUI/notebook context
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.patches import Patch

# simko_utils.py sits one directory up from this file (the app root) --
# add it to the import path using this file's own location, not the
# current working directory, so this works regardless of where the script
# is invoked from.
SCRIPT_DIR = Path(__file__).resolve().parent
APP_ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(APP_ROOT))

from simko_utils import (
    get_classes_by_joint_value, get_classes_by_mutation_status,
    get_differentials, get_mutation_differentials,
    storey_qvalue, benjamini_hochberg,
)

# ---------------------------------------------------------------------------
# Data loading (plain -- no Streamlit caching needed for a one-off script)
# ---------------------------------------------------------------------------

def load_matrices(data_dir=None):
    """Load the four SimKO data matrices directly. Expects the same CSVs the
    live app uses (abundance.csv, expression.csv, mutation.csv,
    dependency.csv), each with a 'protein' column used as the index.

    data_dir defaults to <app root>/data (i.e. APP_ROOT / 'data'), resolved
    from this script's own location -- pass an explicit path to override.
    """
    if data_dir is None:
        data_dir = APP_ROOT / 'data'
    data_dir = Path(data_dir)
    abundance = pd.read_csv(data_dir / 'abundance.csv').set_index('protein')
    expression = pd.read_csv(data_dir / 'expression.csv').set_index('protein')
    mutation = pd.read_csv(data_dir / 'mutation.csv').set_index('protein')
    dependency = pd.read_csv(data_dir / 'dependency.csv').set_index('protein')
    return abundance, expression, mutation, dependency


def load_essential_gene_lists(validation_dir=None):
    """Load the Hart lab core-essential (CEGv2) and non-essential (NEGv1)
    reference gene lists from essential_genes/ -- see that folder's README
    for provenance. Returns two plain lists of gene symbols
    (core_essential, non_essential), ready to use as Phase 1 positive/
    negative controls.

    validation_dir defaults to this script's own directory (i.e. wherever
    validate_simko.py itself lives, since essential_genes/ is a sibling
    folder) -- pass an explicit path to override.
    """
    if validation_dir is None:
        validation_dir = SCRIPT_DIR
    validation_dir = Path(validation_dir)
    ceg = pd.read_csv(validation_dir / 'essential_genes' / 'CEGv2.txt', sep='\t')
    neg = pd.read_csv(validation_dir / 'essential_genes' / 'NEGv1.txt', sep='\t')
    # NEGv1.txt has at least one internal duplicate gene symbol (SLC34A1) --
    # dedupe here (preserving order) so downstream sampling/selection never
    # has to deal with the same gene name appearing twice in one list.
    return list(dict.fromkeys(ceg['GENE'])), list(dict.fromkeys(neg['GENE']))


def find_covered_genes(candidate_genes, available_index, n=None):
    """
    Filter a candidate gene list (e.g. the Hart CEG2/NEG reference lists)
    down to just the ones actually present in a given matrix's index (e.g.
    abundance.index) -- preserves the candidate list's original order.
    Use this instead of guessing individual gene names against your real
    data (tissue-restricted or lowly-expressed genes are often simply
    absent from the CCLE proteomics panel, e.g. G6PC). Pass n to cap the
    number of covered genes returned.

    Also useful directly for Phase 1 (essential-gene positive/negative
    controls), which needs exactly this same coverage check.
    """
    covered = [g for g in candidate_genes if g in available_index]
    return covered[:n] if n else covered


def find_covered_genes_with_variance(candidate_genes, data_df, min_unique_frac=0.5, n=None):
    """
    Like find_covered_genes, but also excludes genes that are effectively
    constant across the panel -- a real, observed failure mode: a gene
    that's heavily imputed (e.g. tissue-restricted genes like DDX4/VASA,
    essentially germline-specific, largely absent from a cancer cell line
    panel) can have a large block of cell lines sharing the exact same
    imputed value. Ranking on such a column doesn't produce a real
    low/median split -- it produces a split dominated by pandas' tie-break
    order (stable, i.e. deterministic by original row order), which is a
    QUALITATIVELY different (non-random) partition than any of the
    genuinely-randomised permutation shuffles, and can distort the
    permutation-based empirical FDR in either direction as a result.

    min_unique_frac: minimum fraction of cell lines that must have a
    distinct value (default 0.5, i.e. at least half the values are unique)
    for a gene to be considered usable as a clean target -- not a
    validated threshold, just a reasonable screen; inspect a candidate's
    value_counts() directly if you want to be sure.

    Also skips (rather than crashes on) any gene symbol that appears MORE
    THAN ONCE in data_df's index -- a real, known issue in proteomics
    matrices, where multiple detected protein groups/isoforms can map to
    the same gene symbol. data_df.loc[gene] returns a DataFrame rather than
    a Series in that case, which there's no principled way to reduce to a
    single set of values here -- so such genes are excluded rather than
    guessed at.
    """
    covered = find_covered_genes(candidate_genes, data_df.index)
    duplicated_in_matrix = set(data_df.index[data_df.index.duplicated(keep=False)])
    clean = []
    skipped_duplicates = []
    for g in covered:
        if g in duplicated_in_matrix:
            skipped_duplicates.append(g)
            continue
        vals = data_df.loc[g]
        if vals.nunique() / len(vals) >= min_unique_frac:
            clean.append(g)
    if skipped_duplicates:
        print(f"find_covered_genes_with_variance: skipped {len(skipped_duplicates)} gene(s) with "
              f"duplicate rows in the matrix index: {skipped_duplicates}")
    return clean[:n] if n else clean


def get_tissue_cell_lines(columns, tissue):
    """
    Return cell lines belonging to a given tissue, using the exact same
    parsing convention as the live app's Global Settings tissue filter
    (cell line name = <LINE>_<TISSUE>, tissue is everything after the
    FIRST underscore -- so multi-word tissues like
    UPPER_AERODIGESTIVE_TRACT are handled correctly). Case-sensitive,
    exact match (not substring).
    """
    return [c for c in columns if '_' in c and c.split('_', 1)[1] == tissue]


def compute_cohort_n(n_cell_lines):
    """Matches the live tool's dynamic cohort-size formula exactly (see
    Design and Implementation), so validation results are directly
    comparable to what the app would show for the same target/filter
    combination: n=20 once 60+ cell lines are available, otherwise a
    smaller cohort scaled to the available panel.
    """
    return ((n_cell_lines - 1) // 3) if n_cell_lines < 60 else 20


# ---------------------------------------------------------------------------
# Phase 0 -- Statistical foundation: permutation-based empirical FDR
# ---------------------------------------------------------------------------

def _vectorized_ttest_p(low_mean, low_std, median_mean, median_std, n_low, n_median):
    """Vectorized equivalent of simko_utils.ttest_from_sample_stats, for use
    inside permutation loops where calling the row-wise version thousands of
    times (once per permutation) would be far too slow. Operates on whole
    numpy arrays (one value per feature) rather than one row at a time.
    """
    pooled_sd = np.sqrt((((n_low - 1) * (low_std ** 2)) + ((n_median - 1) * (median_std ** 2))) / (n_low + n_median - 2))
    with np.errstate(divide='ignore', invalid='ignore'):
        t_stat = (low_mean - median_mean) / (pooled_sd * np.sqrt(1 / n_low + 1 / n_median))
        p = 2 * (1 - t_dist.cdf(np.abs(t_stat), n_low + n_median - 2))
    p = np.where(pooled_sd == 0, 1.0, p)
    p = np.nan_to_num(p, nan=1.0)
    return p


def permutation_null_pvalues(data_df, low_class, median_class, n_permutations=200, seed=None):
    """
    Build an empirical null distribution of per-feature p-values by
    repeatedly shuffling which of the (low_class + median_class) cell lines
    are labelled 'low' vs 'median' -- preserving cohort sizes exactly, but
    breaking any real relationship between cohort assignment and the data --
    and recomputing every feature's p-value under each shuffled labelling.

    Returns an array of shape (n_permutations, n_features). Pooling across
    both features and permutations to build the empirical null (rather than
    using each feature's own handful of permuted values in isolation) is
    standard practice in permutation-based genomics methods, e.g. SAM
    [Tusher VG, Tibshirani R, Chu G. Significance analysis of microarrays
    applied to the ionizing radiation response. PNAS. 2001;98(9):5116-5121],
    and is valid here because every feature's statistic is already
    standardised to a t-statistic before conversion to a p-value, so under
    the null all features are expected to share the same reference
    distribution regardless of their individual variance.
    """
    rng = np.random.default_rng(seed)
    all_cls = list(low_class) + list(median_class)
    n_low = len(low_class)
    n_median = len(all_cls) - n_low

    working = data_df.filter(all_cls).to_numpy()  # features x cell lines, columns in all_cls order
    n_features = working.shape[0]
    n_cls = len(all_cls)

    null_p = np.empty((n_permutations, n_features))
    for i in range(n_permutations):
        shuffled = rng.permutation(n_cls)
        low_idx = shuffled[:n_low]
        median_idx = shuffled[n_low:]

        low_block = working[:, low_idx]
        median_block = working[:, median_idx]
        low_mean, low_std = low_block.mean(axis=1), low_block.std(axis=1, ddof=1)
        median_mean, median_std = median_block.mean(axis=1), median_block.std(axis=1, ddof=1)

        null_p[i, :] = _vectorized_ttest_p(low_mean, low_std, median_mean, median_std, n_low, n_median)

    return null_p


def stratified_permutation_null_pvalues(data_df, low_class, median_class, tissue_of, n_permutations=200, seed=None):
    """
    Like permutation_null_pvalues, but shuffles labels ONLY WITHIN each
    tissue stratum -- preserving the exact number of each tissue's cell
    lines in the 'low' cohort as the real split has, in every single
    permutation, rather than reshuffling freely across the whole cohort.

    Why: holding tissue composition fixed across all permutations means
    any differential that is purely a consequence of low/median having
    different tissue mixes gets baked into the null itself -- it's present
    in the real split AND in every permuted one equally -- so it stops
    contributing to apparent significance. A gene whose apparent signal is
    genuinely target-specific (tissue-independent) is unaffected by this,
    since its real effect operates on top of whatever tissue structure
    exists regardless of how permutations are constrained; a gene whose
    apparent signal is actually tissue composition riding along with the
    ranking will look dramatically less significant once tissue is
    controlled for this way. Comparing this against the standard (global)
    permutation null on the SAME real cohort split isolates the tissue
    contribution using the full pan-cancer dataset -- no subsetting, so no
    confound with reduced statistical power the way comparing against a
    single restricted tissue would have.

    tissue_of: dict or Series mapping each cell line name to its tissue
    label (e.g. {cl: cl.split('_', 1)[1] for cl in columns if '_' in cl},
    matching get_tissue_cell_lines' parsing convention). Cell lines not
    present in tissue_of are treated as their own single-member stratum
    (never reassigned to a different label, since there's nothing to pool
    them with).
    """
    rng = np.random.default_rng(seed)
    all_cls = list(low_class) + list(median_class)
    n_low = len(low_class)

    working = data_df.filter(all_cls).to_numpy()
    n_features = working.shape[0]
    n_cls = len(all_cls)

    tissue_labels = [tissue_of.get(c, f'__singleton_{i}__') for i, c in enumerate(all_cls)]
    strata = {}
    for i, t in enumerate(tissue_labels):
        strata.setdefault(t, []).append(i)

    real_low_positions = set(range(n_low))
    stratum_low_counts = {t: sum(1 for i in idxs if i in real_low_positions) for t, idxs in strata.items()}

    null_p = np.empty((n_permutations, n_features))
    for p in range(n_permutations):
        low_idx = []
        for t, idxs in strata.items():
            k = stratum_low_counts[t]
            if k > 0:
                low_idx.extend(rng.choice(idxs, size=k, replace=False))
        low_idx = np.array(sorted(low_idx), dtype=int)
        low_set = set(low_idx.tolist())
        median_idx = np.array([j for j in range(n_cls) if j not in low_set], dtype=int)

        low_block = working[:, low_idx]
        median_block = working[:, median_idx]
        low_mean, low_std = low_block.mean(axis=1), low_block.std(axis=1, ddof=1)
        median_mean, median_std = median_block.mean(axis=1), median_block.std(axis=1, ddof=1)

        null_p[p, :] = _vectorized_ttest_p(low_mean, low_std, median_mean, median_std, len(low_idx), len(median_idx))

    return null_p


def quantify_tissue_bias(target_genes, layer_matrix, n=None, mode='max', n_permutations=200, seed=42,
                          thresholds=None, exclude_targets_from_stats=True):
    """
    For a chosen target (or list of targets, run one at a time), compute
    BOTH the standard (global) permutation empirical FDR and the
    tissue-STRATIFIED permutation empirical FDR on the exact same
    pan-cancer cohort split, and report both side by side with a
    tissue_bias_ratio column (stratified / standard).

    Interpretation: tissue_bias_ratio near 1 (or below) means tissue
    composition explains little of this target's apparent signal --
    controlling for it barely changes anything, consistent with genuine,
    tissue-independent biology. tissue_bias_ratio well above 1 means a
    substantial part of the apparent pan-cancer signal disappears once
    tissue is controlled for -- i.e. tissue composition, not target
    biology, was driving much of what looked like a hit.

    Uses ALL cell lines (pan-cancer) for every target -- this is
    specifically designed to isolate tissue's contribution WITHOUT the
    reduced-power confound that comparing against a single restricted
    tissue subset would introduce.
    """
    tissue_of = {c: c.split('_', 1)[1] for c in layer_matrix.columns if '_' in c}
    rows = []

    for i, target in enumerate(target_genes):
        target_list = [target] if isinstance(target, str) else list(target)
        target_label = target if isinstance(target, str) else '+'.join(target)

        if n is None:
            n_use = compute_cohort_n(layer_matrix.shape[1])
        else:
            n_use = n
        class_df = get_classes_by_joint_value(target_list, layer_matrix, n_use, mode=mode)
        low_class = list(class_df.loc[class_df['class'] == 'low'].index)
        median_class = list(class_df.loc[class_df['class'] == 'median'].index)

        observed = get_differentials(class_df, layer_matrix, len(low_class), len(median_class))
        if exclude_targets_from_stats:
            observed_for_stats = observed.loc[~observed.index.isin(target_list)]
        else:
            observed_for_stats = observed

        global_null = permutation_null_pvalues(layer_matrix, low_class, median_class,
                                                 n_permutations=n_permutations, seed=seed + i)
        strat_null = stratified_permutation_null_pvalues(layer_matrix, low_class, median_class, tissue_of,
                                                           n_permutations=n_permutations, seed=seed + i)

        global_curve = empirical_fdr_curve(observed_for_stats['p'].values, global_null, thresholds=thresholds)
        strat_curve = empirical_fdr_curve(observed_for_stats['p'].values, strat_null, thresholds=thresholds)
        strat_curve = strat_curve.drop(columns=['observed_hits'])  # identical to global_curve's -- same real data either way

        merged = global_curve.merge(strat_curve, on='nominal_p_threshold', suffixes=('_global', '_stratified'))
        merged.insert(0, 'target', target_label)
        merged['tissue_bias_ratio'] = merged['empirical_fdr_stratified'] / merged['empirical_fdr_global']
        rows.append(merged)

    return pd.concat(rows, ignore_index=True)


def empirical_pvalue_from_permutations(observed_p, null_p_pool):
    """
    For each feature's observed (real-label) nominal p-value, compute an
    empirical p-value: the fraction of the POOLED permuted (null) p-values
    -- across every feature and every permutation -- that are as small or
    smaller. Borrows strength across the whole matrix rather than relying on
    each feature's own small number of permuted values alone.
    """
    pooled = np.sort(null_p_pool.ravel())
    ranks = np.searchsorted(pooled, observed_p, side='right')
    return ranks / len(pooled)


def empirical_fdr_curve(observed_p, null_p_pool, thresholds=None):
    """
    For a range of nominal-p thresholds, compare the number of REAL features
    passing that threshold against the AVERAGE number of features that pass
    it per permutation (the expected false positives under a true null),
    giving an empirical FDR = expected_false_positives / observed_hits at
    each threshold.

    This is the direct, data-driven answer to "how much can the nominal
    p<0.01 threshold actually be trusted, given this specific scan's real
    sample size and correlation structure" -- rather than relying on the
    theoretical BH/Storey assumptions alone.
    """
    if thresholds is None:
        thresholds = [0.1, 0.05, 0.01, 0.005, 0.001, 0.0001]
    rows = []
    for thresh in thresholds:
        observed_hits = int(np.sum(observed_p < thresh))
        permuted_hits_per_run = (null_p_pool < thresh).sum(axis=1)
        expected_false_positives = float(permuted_hits_per_run.mean())
        empirical_fdr = expected_false_positives / observed_hits if observed_hits > 0 else np.nan
        rows.append({
            'nominal_p_threshold': thresh,
            'observed_hits': observed_hits,
            'expected_false_positives_per_scan': round(expected_false_positives, 2),
            'empirical_fdr': round(empirical_fdr, 4) if not np.isnan(empirical_fdr) else np.nan,
        })
    return pd.DataFrame(rows)


def run_permutation_validation(target_genes, layer_matrix, n=None, mode='max',
                                n_permutations=200, seed=42, thresholds=None,
                                active_cell_lines=None):
    """
    End-to-end Phase 0 permutation validation for one KO target (or joint
    target list) against one continuous data layer (pass whichever matrix
    you want to test -- abundance, expression, or dependency).

    n: cohort size; if None, computed automatically via compute_cohort_n()
       on the (optionally tissue/cell-line-filtered) matrix, matching what
       the live app would use for the same filter settings.
    active_cell_lines: optional list to restrict the matrix to before
       classifying cohorts (mirrors the app's tissue/cell-line filter).

    Returns a dict:
      'class_df'   -- the low/median cohort assignment used
      'observed'   -- get_differentials() result, genome-wide, no prefilter
                       (this validation wants the true, unrestricted p-value
                       distribution, not a diff-prefiltered one), with an
                       added 'empirical_p' column
      'null_p_pool'-- the raw (n_permutations x n_features) permuted p-values,
                       kept in case you want to inspect it directly
      'fdr_curve'  -- output of empirical_fdr_curve()
    """
    if active_cell_lines is not None:
        layer_matrix = layer_matrix.filter(active_cell_lines)

    if n is None:
        n = compute_cohort_n(layer_matrix.shape[1])
    else:
        # get_classes_by_joint_value's low/median bands can overlap if n is
        # too large relative to the number of cell lines available (this is
        # a property of the underlying band-based cohort split, not
        # specific to this validation script) -- catch it here with a clear
        # error rather than silently reusing cell lines in both cohorts.
        total = layer_matrix.shape[1]
        median_start = max(round(total / 2 - n / 2), 0)
        if median_start < n:
            raise ValueError(
                f"n={n} is too large for {total} available cell lines -- the low and median "
                f"cohorts would overlap (median band starts at row {median_start}, low band "
                f"covers rows 0-{n-1}). Use compute_cohort_n({total})={compute_cohort_n(total)} "
                f"or smaller, or pass n=None to compute it automatically."
            )

    class_df = get_classes_by_joint_value(target_genes, layer_matrix, n, mode=mode)
    low_class = list(class_df.loc[class_df['class'] == 'low'].index)
    median_class = list(class_df.loc[class_df['class'] == 'median'].index)

    observed = get_differentials(class_df, layer_matrix, len(low_class), len(median_class))
    null_p = permutation_null_pvalues(layer_matrix, low_class, median_class,
                                       n_permutations=n_permutations, seed=seed)

    observed = observed.copy()
    observed['empirical_p'] = empirical_pvalue_from_permutations(observed['p'].values, null_p)

    fdr_curve = empirical_fdr_curve(observed['p'].values, null_p, thresholds=thresholds)

    return {
        'class_df': class_df,
        'observed': observed,
        'null_p_pool': null_p,
        'fdr_curve': fdr_curve,
    }


def compare_correction_methods(observed_diff_df, target_label="", thresholds=None):
    """
    Side-by-side hit counts at a few standard thresholds for: the nominal
    (uncorrected) p, Storey's q (already in observed_diff_df['p_adj'] --
    theoretical), and the permutation-based empirical p (if present as
    'empirical_p' -- data-driven). This is the direct "nominal vs
    Storey-adjusted vs permutation-based hit lists" comparison called for in
    the validation plan's Phase 0.

    Run this across a handful of representative targets and include the
    combined table as a supplementary comparison in the manuscript --
    the point isn't that any one number is "correct", it's showing whether
    the theoretical correction (Storey) roughly agrees with the empirical,
    assumption-free one (permutation), or diverges from it.
    """
    if thresholds is None:
        thresholds = [0.05, 0.01, 0.001]
    rows = []
    for thresh in thresholds:
        row = {
            'target': target_label,
            'threshold': thresh,
            'nominal_p_hits': int((observed_diff_df['p'] < thresh).sum()),
            'storey_q_hits': int((observed_diff_df['p_adj'] < thresh).sum()),
        }
        if 'empirical_p' in observed_diff_df.columns:
            row['permutation_empirical_p_hits'] = int((observed_diff_df['empirical_p'] < thresh).sum())
        rows.append(row)
    return pd.DataFrame(rows)


def run_multi_target_validation(targets, layer_matrix, n=None, mode='max',
                                 n_permutations=200, seed=42, thresholds=None,
                                 active_cell_lines=None, exclude_targets_from_stats=True):
    """
    Run run_permutation_validation() once per target in `targets` (each gets
    its own cohort split), and return combined FDR-curve and
    correction-method-comparison tables (one row per target x threshold),
    ready to drop into a manuscript table or figure.

    exclude_targets_from_stats (default True): a KO target's own
    abundance/expression is trivially different between its 'low' and
    'median' cohorts BY CONSTRUCTION -- that's how the cohort was built --
    so it will always show as a "hit" regardless of any real downstream
    biology. Counting it toward observed_hits inflates detection-power
    statistics with a non-finding. This excludes each target's own gene
    from its row's hit counts before computing the FDR curve and
    correction-method comparison. The permutation null itself is
    unaffected either way (under shuffled labels the target's self-effect
    is destroyed, so it isn't a source of bias there).

    Each target gets a different permutation seed (seed + index) so the
    four runs aren't just repeating the same shuffle pattern.
    """
    fdr_rows = []
    comparison_rows = []

    for i, target in enumerate(targets):
        target_list = [target] if isinstance(target, str) else list(target)
        target_label = target if isinstance(target, str) else '+'.join(target)

        result = run_permutation_validation(
            target_list, layer_matrix, n=n, mode=mode,
            n_permutations=n_permutations, seed=seed + i,
            thresholds=thresholds, active_cell_lines=active_cell_lines,
        )
        observed = result['observed']
        observed_for_stats = (
            observed.loc[~observed.index.isin(target_list)] if exclude_targets_from_stats else observed
        )

        fdr = empirical_fdr_curve(observed_for_stats['p'].values, result['null_p_pool'], thresholds=thresholds)
        fdr.insert(0, 'target', target_label)
        fdr_rows.append(fdr)

        comp = compare_correction_methods(observed_for_stats, target_label=target_label, thresholds=thresholds)
        comparison_rows.append(comp)

    return pd.concat(fdr_rows, ignore_index=True), pd.concat(comparison_rows, ignore_index=True)


def find_negative_controls_from_dependency(target_matrix, dependency_df, max_abs_mean_dependency=0.1,
                                            min_unique_frac=0.5, n=None):
    """
    Alternative to find_covered_genes_with_variance(neg_list, ...): define
    candidate negative-control genes directly from this dataset's OWN
    dependency layer, rather than the external NEGv1 reference list.

    Motivation: NEGv1 is curated purely on CRISPR fitness phenotype, with
    no requirement that the gene's protein product be detectable by mass
    spec. Genes with no fitness phenotype are disproportionately genes
    with narrow or absent tissue expression -- the same property that
    makes a gene undetectable in TMT proteomics. In practice this can
    produce a severe coverage mismatch (observed: as few as 6 of 927
    NEGv1 genes present in a real CCLE abundance panel). Defining
    "non-essential" from this dataset's own dependency scores instead
    sidesteps the mismatch entirely, since the candidate pool is built by
    starting from target_matrix's own genes (guaranteed covered) and
    filtering by dependency, rather than starting from an external list
    and hoping enough of it happens to be covered.

    max_abs_mean_dependency: a gene's mean dependency score across all
    cell lines must fall within +/- this value of zero to count as
    "broadly non-essential in this panel" -- not a validated threshold
    (DepMap's own common-essential/non-essential classification, built
    from a two-component distribution fit, would be more rigorous; this
    is a quick, transparent screen, not a substitute for that if you want
    a fully rigorous version later).
    """
    common_genes = target_matrix.index.intersection(dependency_df.index)
    mean_dep = dependency_df.loc[common_genes].mean(axis=1)
    candidates = list(mean_dep[mean_dep.abs() <= max_abs_mean_dependency].index)
    return find_covered_genes_with_variance(candidates, target_matrix, min_unique_frac=min_unique_frac, n=n)


def run_random_target_panel(layer_matrix, n_targets=10, min_unique_frac=0.5, n=None, mode='max',
                             n_permutations=200, seed=42, thresholds=None, active_cell_lines=None,
                             random_state=42, exclude_genes=None):
    """
    Run the permutation validation across a panel of genes sampled
    completely at random from layer_matrix -- subject only to coverage
    and a basic variance screen (find_covered_genes_with_variance), no
    essentiality claim of any kind.

    Why this exists alongside (not instead of) run_negative_control_panel():
    essentiality (whether a gene has a CRISPR fitness phenotype) has no
    logical bearing on whether ranking cell lines by that gene's abundance
    produces a real downstream proteomic signature. A non-essential gene
    can still drive large, real, reproducible changes elsewhere --
    redundant pathway components and compensable regulators are routinely
    non-essential. And by the textbook definition of synthetic lethality,
    BOTH partners in a real SL pair are individually non-essential --
    meaning a gene set curated purely on "no growth phenotype" (NEGv1) is,
    if anything, enriched for exactly the kind of genuine buffering
    relationship that would produce a real signal, not free of one. Using
    such a set to characterise whether the METHOD produces spurious hits
    imports that flawed assumption. A panel that makes no essentiality
    claim at all sidesteps it, and isolates what we actually want to
    characterise here: whether cohort construction itself (ranking by an
    arbitrary gene, independent of what that gene does biologically)
    tends to produce lineage-driven apparent hits.

    exclude_genes: optional list to exclude from the random draw (e.g.
    your named biologically-motivated targets, so this panel can't
    accidentally re-select one of them).

    Returns a dict: 'targets', 'fdr_curves', 'comparisons', 'summary' --
    same shape as run_negative_control_panel(), minus the
    passed_variance_filter/source/summary_clean_only fields, since there's
    no backfill tier here (the whole-panel candidate pool is large enough
    that running short is unlikely).
    """
    candidates = find_covered_genes_with_variance(
        list(layer_matrix.index), layer_matrix, min_unique_frac=min_unique_frac
    )
    if exclude_genes:
        candidates = [g for g in candidates if g not in exclude_genes]
    if len(candidates) < n_targets:
        raise ValueError(
            f"Only {len(candidates)} genes pass the variance filter in this panel "
            f"(after excluding {exclude_genes}); requested {n_targets}."
        )
    rng = np.random.default_rng(random_state)
    chosen = list(rng.choice(candidates, size=n_targets, replace=False))

    fdr_curves, comparisons = run_multi_target_validation(
        chosen, layer_matrix, n=n, mode=mode, n_permutations=n_permutations,
        seed=seed, thresholds=thresholds, active_cell_lines=active_cell_lines,
    )

    summary = (
        fdr_curves.groupby('nominal_p_threshold')['empirical_fdr']
        .agg(
            median='median', mean='mean',
            q25=lambda s: s.quantile(0.25), q75=lambda s: s.quantile(0.75),
            min='min', max='max', n_genes='count',
        )
        .reset_index()
    )

    return {'targets': chosen, 'fdr_curves': fdr_curves, 'comparisons': comparisons, 'summary': summary}


def run_negative_control_panel(layer_matrix, n_targets=10, min_unique_frac=0.5, n=None, mode='max',
                                n_permutations=200, seed=42, thresholds=None, active_cell_lines=None,
                                validation_dir=None, random_state=42, backfill=True, dependency_df=None):
    """
    Run the permutation validation across a PANEL of randomly-selected
    non-essential candidate genes, rather than trusting a single gene, and
    return both the per-gene results and aggregate (median/mean/quartile/
    min/max) empirical FDR statistics per threshold.

    Why a panel: a single-gene negative control can't distinguish an
    idiosyncratic property of that one gene from a systematic property of
    the validation method -- e.g. DDX4's residual within-lineage signal
    might be specific to that gene, or might reflect something broader
    about how NEGv1 genes behave in this kind of comparison. n=1 can't
    tell you which. This mirrors treating the positive-control side
    (multiple independent targets) as several checks rather than one, and
    directly tests whether a fix (e.g. tissue restriction) generalises
    across genes rather than resolving one specific case.

    Candidate selection, in order of preference:
      1. NEGv1 genes present in layer_matrix and passing the variance
         filter (find_covered_genes_with_variance).
      2. If dependency_df is given and (1) doesn't reach n_targets: extra
         candidates from find_negative_controls_from_dependency(), which
         defines "non-essential" from THIS dataset's own dependency scores
         rather than the external NEGv1 list -- avoids the coverage
         mismatch between NEGv1 (curated on CRISPR fitness, no
         detectability requirement) and a proteomics panel (observed: as
         few as 6 of 927 NEGv1 genes covered at all in a real abundance
         matrix). Pass your loaded dependency matrix here if the NEGv1
         pool alone is too small -- this is usually the fix, not backfill.
      3. If backfill=True and still short: top up from the full covered-
         but-not-variance-filtered NEGv1 pool.
    All non-strict-NEGv1 genes (from source 2 or 3) are flagged False in
    'passed_variance_filter' EXCEPT genes from source 2, which get their
    own True/False based on whether THEY passed the variance filter --
    see 'source' in the per-gene columns for which pool each gene came
    from.

    Returns a dict:
      'targets'              -- the full gene list actually used
      'passed_variance_filter' -- dict mapping each target to True/False
      'source'                -- dict mapping each target to 'negv1_strict',
                                  'dependency_derived', or 'negv1_backfill'
      'fdr_curves'            -- per-gene FDR curve rows, with
                                  'passed_variance_filter' and 'source' columns
      'comparisons'           -- per-gene correction-method comparison rows, same columns
      'summary'               -- aggregate stats across the FULL panel
      'summary_clean_only'    -- the same aggregate stats using ONLY genes
                                  that passed the variance filter (None if
                                  none did) -- compare before deciding
                                  which to report
    """
    _, neg_list = load_essential_gene_lists(validation_dir)
    strict_candidates = find_covered_genes_with_variance(neg_list, layer_matrix, min_unique_frac=min_unique_frac)
    rng = np.random.default_rng(random_state)

    chosen = []
    passed_filter = {}
    source = {}

    if len(strict_candidates) >= n_targets:
        chosen = list(rng.choice(strict_candidates, size=n_targets, replace=False))
        passed_filter = {g: True for g in chosen}
        source = {g: 'negv1_strict' for g in chosen}
    else:
        chosen = list(strict_candidates)
        passed_filter.update({g: True for g in chosen})
        source.update({g: 'negv1_strict' for g in chosen})
        n_needed = n_targets - len(chosen)

        if dependency_df is not None and n_needed > 0:
            dep_candidates = find_negative_controls_from_dependency(
                layer_matrix, dependency_df, min_unique_frac=min_unique_frac
            )
            dep_pool = [g for g in dep_candidates if g not in chosen]
            n_from_dep = min(n_needed, len(dep_pool))
            if n_from_dep > 0:
                dep_chosen = list(rng.choice(dep_pool, size=n_from_dep, replace=False))
                print(f"Only {len(strict_candidates)} NEGv1 genes passed the variance filter; "
                      f"adding {n_from_dep} dependency-derived non-essential candidates instead "
                      f"of falling back to unfiltered NEGv1 backfill: {dep_chosen}")
                chosen += dep_chosen
                passed_filter.update({g: True for g in dep_chosen})  # these DID pass the variance filter, just via a different candidate pool
                source.update({g: 'dependency_derived' for g in dep_chosen})
                n_needed = n_targets - len(chosen)

        if n_needed > 0:
            if not backfill:
                raise ValueError(
                    f"Only {len(chosen)} candidates found (NEGv1 strict"
                    + (" + dependency-derived" if dependency_df is not None else "")
                    + f"); requested {n_targets} and backfill=False. Lower n_targets, pass "
                    "dependency_df for more candidates, or set backfill=True."
                )
            all_covered = find_covered_genes(neg_list, layer_matrix.index)
            remaining_pool = [g for g in all_covered if g not in chosen]
            if len(remaining_pool) < n_needed:
                raise ValueError(
                    f"Only {len(chosen)} candidates found in total (of {len(all_covered)} NEGv1 "
                    f"genes covered in this panel at all) -- cannot reach n_targets={n_targets} "
                    "even with backfill. Lower n_targets, or pass dependency_df for a much larger "
                    "candidate pool that doesn't depend on NEGv1 coverage."
                )
            backfilled = list(rng.choice(remaining_pool, size=n_needed, replace=False))
            print(f"Backfilling {n_needed} more from the full NEGv1-covered pool WITHOUT the "
                  f"variance screen: {backfilled}. These may include near-constant/heavily-imputed "
                  f"genes -- check summary_clean_only before treating the full panel as uniformly clean.")
            chosen += backfilled
            passed_filter.update({g: False for g in backfilled})
            source.update({g: 'negv1_backfill' for g in backfilled})

    fdr_curves, comparisons = run_multi_target_validation(
        chosen, layer_matrix, n=n, mode=mode, n_permutations=n_permutations,
        seed=seed, thresholds=thresholds, active_cell_lines=active_cell_lines,
    )
    fdr_curves['passed_variance_filter'] = fdr_curves['target'].map(passed_filter)
    fdr_curves['source'] = fdr_curves['target'].map(source)
    comparisons['passed_variance_filter'] = comparisons['target'].map(passed_filter)
    comparisons['source'] = comparisons['target'].map(source)

    def _summarize(df):
        return (
            df.groupby('nominal_p_threshold')['empirical_fdr']
            .agg(
                median='median', mean='mean',
                q25=lambda s: s.quantile(0.25), q75=lambda s: s.quantile(0.75),
                min='min', max='max', n_genes='count',
            )
            .reset_index()
        )

    summary = _summarize(fdr_curves)
    clean_only = fdr_curves[fdr_curves['passed_variance_filter']]
    summary_clean_only = _summarize(clean_only) if len(clean_only) else None

    return {
        'targets': chosen,
        'passed_variance_filter': passed_filter,
        'source': source,
        'fdr_curves': fdr_curves,
        'comparisons': comparisons,
        'summary': summary,
        'summary_clean_only': summary_clean_only,
    }


def run_phase0_tissue_comparison(layer_matrix, dependency_df, real_targets, n_permutations=200, seed=42,
                                  n_panel_targets=10, n_tissues=2, min_tissue_size=20):
    """
    Run the full Phase 0 comparison -- real (biologically-motivated)
    targets, the essentiality-based negative-control panel, and the
    random-gene panel -- pan-cancer AND restricted to each of the
    n_tissues largest tissue lineages with at least min_tissue_size cell
    lines. Returns one combined comparison table (rows = nominal p
    threshold, one column per panel-type x scope) plus the full underlying
    results per scope for further inspection.

    Why more than one tissue: a tissue-restricted improvement in the
    non-specific panels (empirical FDR moving toward the ~1.0 pure-noise
    baseline) could be a genuine, general effect of removing lineage
    structure, or it could be a quirk of that one tissue's particular
    composition/subtype makeup. Testing a second, biologically distinct
    tissue and seeing whether the same qualitative pattern holds is the
    direct way to tell those apart.

    Returns (comparison_table, results) where results is a dict keyed by
    scope ('pancancer', then each tissue name) with keys 'real_median',
    'neg_control_median', 'random_median' (per-threshold Series) and 'raw'
    (the full underlying run outputs, in case you want to inspect a
    specific scope/panel further).
    """
    tissue_counts = pd.Series([c.split('_', 1)[1] for c in layer_matrix.columns if '_' in c]).value_counts()
    tissue_counts = tissue_counts[tissue_counts >= min_tissue_size]
    tissues = list(tissue_counts.index[:n_tissues])
    print(f"Tissue composition (top {n_tissues} used, min {min_tissue_size} cell lines):")
    print(tissue_counts.head(10))
    print(f"\nUsing tissues: {tissues}")

    def _run_scope(active_cell_lines, label):
        print(f"\n{'=' * 70}\n{label}\n{'=' * 70}")
        fdr_curves, _ = run_multi_target_validation(
            real_targets, layer_matrix, n_permutations=n_permutations, seed=seed, active_cell_lines=active_cell_lines
        )
        print("Real targets:")
        print(fdr_curves.to_string(index=False))

        neg_panel = run_negative_control_panel(
            layer_matrix, n_targets=n_panel_targets, n_permutations=n_permutations, seed=seed,
            dependency_df=dependency_df, active_cell_lines=active_cell_lines,
        )
        print(f"\nNegative-control panel: {neg_panel['targets']}")
        print(neg_panel['summary'].to_string(index=False))

        rand_panel = run_random_target_panel(
            layer_matrix, n_targets=n_panel_targets, n_permutations=n_permutations, seed=seed,
            exclude_genes=real_targets, active_cell_lines=active_cell_lines,
        )
        print(f"\nRandom panel: {rand_panel['targets']}")
        print(rand_panel['summary'].to_string(index=False))

        return {
            'real_median': fdr_curves.groupby('nominal_p_threshold')['empirical_fdr'].median(),
            'neg_control_median': neg_panel['summary'].set_index('nominal_p_threshold')['median'],
            'random_median': rand_panel['summary'].set_index('nominal_p_threshold')['median'],
            'raw': {'real': fdr_curves, 'neg_control': neg_panel, 'random': rand_panel},
        }

    results = {'pancancer': _run_scope(None, "PAN-CANCER")}
    for tissue in tissues:
        tissue_cls = get_tissue_cell_lines(layer_matrix.columns, tissue)
        results[tissue] = _run_scope(tissue_cls, f"TISSUE-RESTRICTED: {tissue} ({len(tissue_cls)} cell lines)")

    cols = {}
    for scope, res in results.items():
        key = scope.lower()
        cols[f'real_{key}'] = res['real_median']
        cols[f'negctrl_{key}'] = res['neg_control_median']
        cols[f'random_{key}'] = res['random_median']
    comparison_table = pd.DataFrame(cols)

    return comparison_table, results


# ---------------------------------------------------------------------------
# Phase 0 figures -- built directly from the DataFrames the functions above
# return, so figures always reflect whatever data was actually run rather
# than needing separate hardcoded numbers kept in sync by hand.
# ---------------------------------------------------------------------------

_PANEL_COLORS = {
    'Real targets': '#1b9e77',
    'Essentiality-curated panel': '#7570b3',
    'Random panel': '#d95f02',
}


def plot_real_target_fdr_curves(fdr_curves, output_path, noise_baseline=(0.85, 0.96), title=None):
    """
    Line plot of empirical FDR vs. nominal p threshold, one line per
    target, built directly from a run_multi_target_validation() (or
    run_permutation_validation()'s fdr_curve, renamed to 'target') output.

    noise_baseline: shaded reference band showing the empirical FDR range
    a pure-noise target should occupy -- default (0.85, 0.96) is what
    synthetic pure-noise testing gave during development; rerun that
    yourself and pass your own range if you want it re-derived.

    How to read this figure: each line traces one target's empirical FDR
    as the significance threshold tightens, left to right. A line sitting
    well below the shaded band means that target's nominal hits are
    largely real rather than chance -- this is the direct evidence that
    SimKO's p-values reflect genuine downstream biology for that target,
    not an artefact of the testing framework.
    """
    fig, ax = plt.subplots(figsize=(6.5, 5))
    ax.axhspan(*noise_baseline, color='#cccccc', alpha=0.5,
               label=f'Pure-noise baseline\n(synthetic testing, {noise_baseline[0]:g}\u2013{noise_baseline[1]:g})')

    palette = plt.cm.Dark2(np.linspace(0, 1, fdr_curves['target'].nunique()))
    for color, (target, sub) in zip(palette, fdr_curves.groupby('target', sort=False)):
        sub = sub.sort_values('nominal_p_threshold', ascending=False)
        ax.plot(sub['nominal_p_threshold'], sub['empirical_fdr'], marker='o',
                label=target, color=color, linewidth=2, markersize=5)

    ax.set_xscale('log')
    ax.invert_xaxis()
    ax.set_xlabel('Nominal p-value threshold')
    ax.set_ylabel('Empirical FDR (permutation-based)')
    ax.set_title(title or 'Empirical FDR by target')
    ax.set_ylim(0, max(1.0, fdr_curves['empirical_fdr'].max() * 1.05))
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f'{x:g}'))
    ax.legend(loc='upper right', fontsize=8.5, frameon=False)
    ax.grid(alpha=0.25, linewidth=0.5)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    return output_path


def plot_multipanel_multitissue_comparison(comparison_table, output_path, scope_labels=None, title=None):
    """
    Faceted line plot (one panel per tissue scope) of median empirical FDR
    vs. threshold, one line per panel type, built directly from
    run_phase0_tissue_comparison()'s comparison_table (columns named
    '<panel>_<scope>', e.g. 'real_pancancer', 'negctrl_lung' -- exactly
    what that function returns).

    scope_labels: optional dict mapping scope key (e.g. 'pancancer',
    'lung') to a display title (e.g. 'Lung (77 lines)') -- if omitted,
    scope keys are shown as-is.

    How to read this figure: within each panel (tissue scope), the 'Real
    targets' line sitting below the non-specific panels is the core
    validation claim -- SimKO recovers genuine target-specific signal
    rather than an artefact of how cohorts are constructed. Comparing
    panels left-to-right shows what happens as tissue diversity shrinks:
    non-specific panels moving toward the horizontal dotted line (1.0,
    chance level) while Real targets stays separated is evidence that
    tissue restriction is doing its job without destroying real signal.
    Values noticeably above 1.0, especially at the strictest thresholds,
    usually reflect small observed-hit counts in that scope rather than a
    new artefact -- check the underlying counts before over-interpreting
    a single spike.
    """
    panel_keys = {'real': 'Real targets', 'negctrl': 'Essentiality-curated panel', 'random': 'Random panel'}
    scopes = []
    for col in comparison_table.columns:
        _, scope_key = col.split('_', 1)
        if scope_key not in scopes:
            scopes.append(scope_key)

    fig, axes = plt.subplots(1, len(scopes), figsize=(4.3 * len(scopes), 4.5), sharey=True)
    if len(scopes) == 1:
        axes = [axes]

    for ax, scope_key in zip(axes, scopes):
        ax.axhline(1.0, color='grey', linestyle=':', linewidth=1, zorder=0)
        for panel_key, panel_label in panel_keys.items():
            col = f'{panel_key}_{scope_key}'
            if col not in comparison_table.columns:
                continue
            series = comparison_table[col].sort_index(ascending=False)
            ax.plot(series.index, series.values, marker='o', label=panel_label,
                     color=_PANEL_COLORS[panel_label], linewidth=2, markersize=4)
        ax.set_xscale('log')
        ax.invert_xaxis()
        label = (scope_labels or {}).get(scope_key, scope_key)
        ax.set_title(label, fontsize=10.5)
        ax.set_xlabel('Nominal p threshold')
        ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f'{x:g}'))
        ax.grid(alpha=0.25, linewidth=0.5)

    axes[0].set_ylabel('Median empirical FDR')
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 1.06), ncol=3, frameon=False, fontsize=9.5)
    fig.suptitle(title or 'Panel comparison across tissue scopes', y=1.14, fontsize=12)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    return output_path


def plot_tissue_bias_ratio(tissue_bias_df, target_categories, output_path, threshold=0.01, title=None):
    """
    Horizontal bar chart of tissue-bias ratio (stratified / global
    empirical FDR) per gene at a chosen threshold, built directly from
    quantify_tissue_bias() output.

    target_categories: dict mapping gene name -> category label (e.g.
    'Real targets' / 'Essentiality-curated panel' / 'Random panel') used
    only to colour the bars -- build this from whatever target lists you
    actually ran (see the __main__ block below for the pattern).

    How to read this figure: this is a PER-TARGET diagnostic, not a
    real-vs-fake test (that question is answered by
    plot_multipanel_multitissue_comparison above). Bars extending past
    the dotted reference line at 1.0 indicate genes whose own abundance is
    itself tissue-structured -- not a validity concern for SimKO, but a
    signal that pan-cancer results for that specific target should be
    interpreted with the tissue filter in mind. Genuinely real,
    literature-motivated targets CAN and did show high ratios here (e.g.
    proliferation-linked genes), so a high bar for a real target is not
    itself evidence against the tool -- it's a reason to re-run that
    particular target with the tissue filter before trusting its
    pan-cancer output.
    """
    sub = tissue_bias_df[tissue_bias_df['nominal_p_threshold'] == threshold].copy()
    sub = sub.sort_values('tissue_bias_ratio', ascending=True)

    fig, ax = plt.subplots(figsize=(7, max(3, 0.5 * len(sub) + 1)))
    y_pos = np.arange(len(sub))
    colors = [_PANEL_COLORS.get(target_categories.get(t, 'Essentiality-curated panel'), '#7570b3') for t in sub['target']]
    ax.barh(y_pos, sub['tissue_bias_ratio'], color=colors, edgecolor='white', height=0.7)
    ax.axvline(1.0, color='grey', linestyle=':', linewidth=1.2)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(sub['target'])
    ax.set_xlabel(f'Tissue-bias ratio (stratified / global empirical FDR), p<{threshold:g}')
    ax.set_title(title or 'Per-gene tissue-bias ratio')
    for y, r in zip(y_pos, sub['tissue_bias_ratio']):
        if pd.notna(r):
            ax.text(r + sub['tissue_bias_ratio'].max() * 0.02, y, f'{r:.2f}', va='center', fontsize=8.5)

    used_categories = sorted(set(target_categories.get(t, 'Essentiality-curated panel') for t in sub['target']),
                              key=lambda c: list(_PANEL_COLORS).index(c) if c in _PANEL_COLORS else 99)
    legend_elems = [Patch(facecolor=_PANEL_COLORS.get(c, '#7570b3'), label=c) for c in used_categories]
    ax.legend(handles=legend_elems, loc='lower right', fontsize=8.5, frameon=False)
    finite_max = sub['tissue_bias_ratio'].replace([np.inf, -np.inf], np.nan).max()
    ax.set_xlim(0, (finite_max if pd.notna(finite_max) else 1) * 1.25)
    ax.grid(axis='x', alpha=0.25, linewidth=0.5)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    return output_path


# ---------------------------------------------------------------------------
# Demo / example usage
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    abundance, expression, mutation, dependency = load_matrices()
    LAYER_NAME = 'abundance'
    layer = {'abundance': abundance, 'expression': expression, 'dependency': dependency}[LAYER_NAME]

    # Four biologically-motivated targets: PCNA (CEGv2 positive control,
    # expect strong real signal), ARID1B and SMARCA4 (SWI/SNF chromatin
    # remodelling complex, expect real but more modest signal), WRN
    # (Discussion case-study candidate).
    TARGETS = ['PCNA', 'ARID1B', 'WRN', 'SMARCA4']
    missing = [t for t in TARGETS if t not in layer.index]
    if missing:
        print(f"Not covered in this panel, skipping: {missing}")
    targets_to_run = [t for t in TARGETS if t not in missing]

    # Full Phase 0 comparison: real targets, essentiality-based negative-
    # control panel, and the essentiality-free random-gene panel, each run
    # pan-cancer AND restricted to the top n_tissues largest tissue
    # lineages. Testing more than one tissue is the key check: if a
    # tissue-restricted improvement in the non-specific panels holds
    # across genuinely different tissues, that's evidence of a general
    # lineage-structure effect rather than a quirk of one tissue's
    # particular composition.
    comparison_table, results = run_phase0_tissue_comparison(
        layer, dependency, targets_to_run, n_permutations=200, seed=42,
        n_panel_targets=10, n_tissues=2, min_tissue_size=20,
    )

    print(f"\n{'=' * 70}\nCombined comparison (median empirical FDR, all scopes)\n{'=' * 70}")
    print(comparison_table.to_string())
    comparison_table.to_csv('phase0_full_comparison.csv')
    print("\nSaved phase0_full_comparison.csv -- this table is the headline result for Phase 0.")

    # Save each scope's underlying data too, in case you want to inspect a
    # specific panel/scope beyond the summary table above.
    for scope, res in results.items():
        suffix = scope.lower()
        res['raw']['real'].to_csv(f'phase0_real_targets_{suffix}.csv', index=False)
        res['raw']['neg_control']['fdr_curves'].to_csv(f'phase0_negative_control_panel_{suffix}.csv', index=False)
        res['raw']['random']['fdr_curves'].to_csv(f'phase0_random_panel_{suffix}.csv', index=False)
    print(f"Saved per-scope detail files for: {list(results.keys())}")

    # --- Figures (all built live from the results just computed above) ---
    FIG_DIR = SCRIPT_DIR / 'figures'
    FIG_DIR.mkdir(exist_ok=True)

    scope_labels = {'pancancer': f'Pan-cancer ({layer.shape[1]} lines)'}
    for scope, res in results.items():
        if scope == 'pancancer':
            continue
        tissue_cls = get_tissue_cell_lines(layer.columns, scope)
        scope_labels[scope.lower()] = f'{scope.title()} ({len(tissue_cls)} lines)'

    fig1_path = plot_real_target_fdr_curves(
        results['pancancer']['raw']['real'], FIG_DIR / 'fig1_real_target_fdr_curves.png',
        title='Empirical FDR for real targets (Abundance layer, pan-cancer)',
    )
    print(f"\nSaved {fig1_path}")

    fig2_path = plot_multipanel_multitissue_comparison(
        comparison_table, FIG_DIR / 'fig2_multipanel_multitissue_comparison.png',
        scope_labels=scope_labels, title='Real targets vs. non-specific panels, across tissue scopes',
    )
    print(f"Saved {fig2_path}")

    # Figure 3 needs quantify_tissue_bias(), which isn't otherwise part of
    # this pipeline (it's pan-cancer-only by construction -- stratified
    # permutation is degenerate once you've already restricted to one
    # tissue, since there's then only one stratum). Runs BOTH a global and
    # a stratified permutation null per target, so this step roughly
    # doubles this section's runtime.
    tissue_bias_targets = targets_to_run + results['pancancer']['raw']['neg_control']['targets']
    print(f"\n{'=' * 70}\nQuantifying tissue bias for {tissue_bias_targets}\n{'=' * 70}")
    tissue_bias_df = quantify_tissue_bias(tissue_bias_targets, layer, n_permutations=200, seed=42)
    tissue_bias_df.to_csv('phase0_tissue_bias_quantification_full.csv', index=False)
    print(tissue_bias_df.to_string(index=False))

    target_categories = {t: 'Real targets' for t in targets_to_run}
    target_categories.update({t: 'Essentiality-curated panel' for t in results['pancancer']['raw']['neg_control']['targets']})

    fig3_path = plot_tissue_bias_ratio(
        tissue_bias_df, target_categories, FIG_DIR / 'fig3_tissue_bias_ratio_per_gene.png',
        threshold=0.01, title='Per-gene tissue-bias ratio, pan-cancer (Abundance layer)',
    )
    print(f"\nSaved {fig3_path}")
    print(f"\nAll three Phase 0 figures saved to {FIG_DIR}")