"""
SimKO consolidated validation suite.

Everything from the validation programme in one place: shared utility and
statistical functions first, then a RUN PIPELINE section at the bottom
made of independent `# %%` cells. Open this file in VS Code (or similar)
and run cells with Shift+Enter, or run the whole file as a script -- no
CLI, no argparse. Each cell is self-contained; run only the ones you
need, in order where one cell's variables feed the next (noted in each
cell's comment).

CELL MAP -- which cell makes which figure
-----------------------------------------
Figures are only WRITTEN by the last cell (Cell 17); the earlier cells compute
the results each figure is drawn from. Run the feeder cell(s), then Cell 17.

  Figure 1, 2, 3   <- Cell 2   (Phase 0)
  Figure 4         <- Cell 5 (panel A, RPE1 protein) and Cell 7 (panel B, Abundance vs HCT116 RNA; needs Cells 3 and 6)
  Figure 5         <- Cells 3 and 8   (Abundance vs Expression layer, on real RPE1 protein data)
  Figure 6         <- Cell 10  (BRCA1/BRCA2 -> PARP1)
  Figure 7, 8, 8b  <- Cell 11  (BioGRID)
  Figure 9         <- Cell 11b (optional; BioGRID, protein vs RNA vs plain correlation; needs Cell 11)
  No figure (text results only): Cells 4, 7b, 7c and 7d (optional quick looks), 9, 11c (optional; SynLethDB, needs Cell 11b), 11d (optional; reference width + tissue-stratified null, needs Cell 11b), 11e (optional; readout-matched controls + SynLethDB high-confidence, needs Cells 11b and 11c), 12-16

Every cell starts with a banner naming its figure output and its data files.
All external data files are set in the DATA FILES block near the top (VAL_DIR = this folder).

Phase 0 (statistical foundation) needs no external data and is safe to
run first. Several later cells need external files (GSE101966's RNA-seq,
the ARID1A_RPE_KOs proteomics table, a BioGRID download) -- edit the path
variables in the DATA FILES block near the top (all default to this validation/ folder).
"""
import sys
import re
import warnings
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import t as t_dist, spearmanr, hypergeom, fisher_exact
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

SCRIPT_DIR = Path(__file__).resolve().parent
APP_ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(APP_ROOT))

# ---------------------------------------------------------------------------
# FIGURE CONFIGURATION -- every manuscript figure is controlled from here.
# Edit these, then re-run the final "ALL MANUSCRIPT FIGURES" cell: no other
# cell needs re-running (the figures are drawn from results already in memory).
# ---------------------------------------------------------------------------
# DATA FILES -- everything the validation reads or caches lives in validation/ (this script's own folder).
# (The four SimKO matrices -- abundance/expression/mutation/dependency.csv -- still come from the app's data/ folder.)
VAL_DIR = SCRIPT_DIR
ARID1A_RPE_PATH = str(VAL_DIR / 'ARID1A_RPE_KOs.csv')                      # RPE1 ARID1A-KO proteomics (Cells 3-5, 8)
HCT116_FPKM_PATH = str(VAL_DIR / 'GSE101966_HCT116_FPKM.txt.gz')           # HCT116 ARID1A-KO RNA-seq (Cells 6-8)
BIOGRID_PATH = str(VAL_DIR / 'BIOGRID-ORGANISM-Homo_sapiens-5.0.261.tab3.txt')  # BioGRID release (Cell 11)
BIOGRID_RESULTS_PATH = str(VAL_DIR / 'biogrid_at_scale_results.csv')       # BioGRID permutation cache -- KEEP this file (Cell 11)
BIOGRID_METHODS_PATH = str(VAL_DIR / 'biogrid_method_comparison.csv')       # BioGRID protein/RNA/correlation cache -- KEEP this file (Cell 11b)
PHASE0_TABLE_PATH = str(VAL_DIR / 'phase0_full_comparison.csv')            # Cell 2 output table

FIG_DIR = SCRIPT_DIR / 'figures_oct29'            # fig1..fig8 land here (fresh folder, so no stale files)
SUPP_DIR = FIG_DIR / 'supplementary'             # per-analysis scatter PNGs/CSVs land here
FIG_DIR.mkdir(exist_ok=True)
SUPP_DIR.mkdir(exist_ok=True)
FIG_DPI = 300
FIG_STYLE = {
    'font.size': 10, 'axes.titlesize': 11, 'axes.labelsize': 10,
    'xtick.labelsize': 9, 'ytick.labelsize': 9, 'legend.fontsize': 8.5,
    'axes.spines.top': False, 'axes.spines.right': False,
    'figure.max_open_warning': 0,
}
plt.rcParams.update(FIG_STYLE)

from simko_utils import (
    get_classes_by_joint_value, get_classes_by_mutation_status,
    get_differentials, get_mutation_differentials,
)


# ===========================================================================
# SECTION 0: Shared utilities (data loading, gene-list helpers, cohort sizing)
# ===========================================================================

def load_matrices(data_dir=None):
    """Load the four SimKO data matrices directly. Expects the same CSVs the
    live app uses (abundance.csv, expression.csv, mutation.csv,
    dependency.csv), each with a 'protein' column used as the index."""
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
    reference gene lists from essential_genes/."""
    if validation_dir is None:
        validation_dir = SCRIPT_DIR
    validation_dir = Path(validation_dir)
    ceg = pd.read_csv(validation_dir / 'essential_genes' / 'CEGv2.txt', sep='\t')
    neg = pd.read_csv(validation_dir / 'essential_genes' / 'NEGv1.txt', sep='\t')
    return list(dict.fromkeys(ceg['GENE'])), list(dict.fromkeys(neg['GENE']))


def find_covered_genes(candidate_genes, available_index, n=None):
    """Filter a candidate gene list down to genes actually present in a
    given matrix's index, preserving the candidate list's order."""
    covered = [g for g in candidate_genes if g in available_index]
    return covered[:n] if n else covered


def find_covered_genes_with_variance(candidate_genes, data_df, min_unique_frac=0.5, n=None):
    """Like find_covered_genes, but also excludes near-constant genes
    (heavy imputation artefact) and genes with duplicate index entries
    (multiple protein groups mapping to one symbol)."""
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
    """Cell lines belonging to a given tissue (name = <LINE>_<TISSUE>,
    tissue = everything after the first underscore)."""
    return [c for c in columns if '_' in c and c.split('_', 1)[1] == tissue]


def compute_cohort_n(n_cell_lines):
    """n=20 once 60+ cell lines are available, otherwise scaled down."""
    return ((n_cell_lines - 1) // 3) if n_cell_lines < 60 else 20


# ===========================================================================
# SECTION 1: Phase 0 -- permutation-based statistical validation
# ===========================================================================

def _vectorized_ttest_p(low_mean, low_std, median_mean, median_std, n_low, n_median):
    pooled_sd = np.sqrt((((n_low - 1) * (low_std ** 2)) + ((n_median - 1) * (median_std ** 2))) / (n_low + n_median - 2))
    with np.errstate(divide='ignore', invalid='ignore'):
        t_stat = (low_mean - median_mean) / (pooled_sd * np.sqrt(1 / n_low + 1 / n_median))
        p = 2 * (1 - t_dist.cdf(np.abs(t_stat), n_low + n_median - 2))
    p = np.where(pooled_sd == 0, 1.0, p)
    p = np.nan_to_num(p, nan=1.0)
    return p


def permutation_null_pvalues(data_df, low_class, median_class, n_permutations=200, seed=None):
    """Empirical null of per-feature p-values by repeatedly shuffling
    low/median cohort labels, preserving cohort sizes exactly."""
    rng = np.random.default_rng(seed)
    all_cls = list(low_class) + list(median_class)
    n_low = len(low_class)
    n_median = len(all_cls) - n_low

    working = data_df.filter(all_cls).to_numpy()
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
    """Like permutation_null_pvalues, but shuffles labels ONLY WITHIN each
    tissue stratum, holding tissue composition fixed across every
    permutation -- isolates tissue's contribution to apparent significance."""
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
    """Global vs tissue-stratified empirical FDR side by side, per target
    -- ratio near 1 means little tissue contribution; well above 1 means
    tissue composition, not target biology, is driving much of the signal."""
    tissue_of = {c: c.split('_', 1)[1] for c in layer_matrix.columns if '_' in c}
    rows = []
    for i, target in enumerate(target_genes):
        target_list = [target] if isinstance(target, str) else list(target)
        target_label = target if isinstance(target, str) else '+'.join(target)
        n_use = n if n is not None else compute_cohort_n(layer_matrix.shape[1])
        class_df = get_classes_by_joint_value(target_list, layer_matrix, n_use, mode=mode)
        low_class = list(class_df.loc[class_df['class'] == 'low'].index)
        median_class = list(class_df.loc[class_df['class'] == 'median'].index)

        observed = get_differentials(class_df, layer_matrix, len(low_class), len(median_class))
        observed_for_stats = observed.loc[~observed.index.isin(target_list)] if exclude_targets_from_stats else observed

        global_null = permutation_null_pvalues(layer_matrix, low_class, median_class,
                                                 n_permutations=n_permutations, seed=seed + i)
        strat_null = stratified_permutation_null_pvalues(layer_matrix, low_class, median_class, tissue_of,
                                                           n_permutations=n_permutations, seed=seed + i)

        global_curve = empirical_fdr_curve(observed_for_stats['p'].values, global_null, thresholds=thresholds)
        strat_curve = empirical_fdr_curve(observed_for_stats['p'].values, strat_null, thresholds=thresholds)
        strat_curve = strat_curve.drop(columns=['observed_hits'])

        merged = global_curve.merge(strat_curve, on='nominal_p_threshold', suffixes=('_global', '_stratified'))
        merged.insert(0, 'target', target_label)
        merged['tissue_bias_ratio'] = merged['empirical_fdr_stratified'] / merged['empirical_fdr_global']
        rows.append(merged)
    return pd.concat(rows, ignore_index=True)


def empirical_pvalue_from_permutations(observed_p, null_p_pool):
    """Fraction of pooled null p-values as small or smaller than each
    observed p-value -- borrows strength across the whole matrix."""
    pooled = np.sort(null_p_pool.ravel())
    ranks = np.searchsorted(pooled, observed_p, side='right')
    return ranks / len(pooled)


def empirical_fdr_curve(observed_p, null_p_pool, thresholds=None):
    """Real hits vs. average per-permutation hits at each threshold ->
    empirical FDR = expected_false_positives / observed_hits."""
    if thresholds is None:
        thresholds = [0.1, 0.05, 0.01, 0.005, 0.001, 0.0001]
    rows = []
    for thresh in thresholds:
        observed_hits = int(np.sum(observed_p < thresh))
        permuted_hits_per_run = (null_p_pool < thresh).sum(axis=1)
        expected_false_positives = float(permuted_hits_per_run.mean())
        empirical_fdr = expected_false_positives / observed_hits if observed_hits > 0 else np.nan
        rows.append({
            'nominal_p_threshold': thresh, 'observed_hits': observed_hits,
            'expected_false_positives_per_scan': round(expected_false_positives, 2),
            'empirical_fdr': round(empirical_fdr, 4) if not np.isnan(empirical_fdr) else np.nan,
        })
    return pd.DataFrame(rows)


def run_permutation_validation(target_genes, layer_matrix, n=None, mode='max',
                                n_permutations=200, seed=42, thresholds=None,
                                active_cell_lines=None):
    """End-to-end Phase 0 permutation validation for one KO target against
    one continuous data layer."""
    if active_cell_lines is not None:
        layer_matrix = layer_matrix.filter(active_cell_lines)
    if n is None:
        n = compute_cohort_n(layer_matrix.shape[1])
    else:
        total = layer_matrix.shape[1]
        median_start = max(round(total / 2 - n / 2), 0)
        if median_start < n:
            raise ValueError(
                f"n={n} is too large for {total} available cell lines -- the low and median "
                f"cohorts would overlap. Use compute_cohort_n({total})={compute_cohort_n(total)} "
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

    return {'class_df': class_df, 'observed': observed, 'null_p_pool': null_p, 'fdr_curve': fdr_curve}


def compare_correction_methods(observed_diff_df, target_label="", thresholds=None):
    """Nominal p vs Storey's q vs permutation-based empirical p, hit
    counts side by side at a few thresholds."""
    if thresholds is None:
        thresholds = [0.05, 0.01, 0.001]
    rows = []
    for thresh in thresholds:
        row = {
            'target': target_label, 'threshold': thresh,
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
    """run_permutation_validation() once per target, combined tables out."""
    fdr_rows, comparison_rows = [], []
    for i, target in enumerate(targets):
        target_list = [target] if isinstance(target, str) else list(target)
        target_label = target if isinstance(target, str) else '+'.join(target)
        result = run_permutation_validation(
            target_list, layer_matrix, n=n, mode=mode, n_permutations=n_permutations,
            seed=seed + i, thresholds=thresholds, active_cell_lines=active_cell_lines,
        )
        observed = result['observed']
        observed_for_stats = observed.loc[~observed.index.isin(target_list)] if exclude_targets_from_stats else observed
        fdr = empirical_fdr_curve(observed_for_stats['p'].values, result['null_p_pool'], thresholds=thresholds)
        fdr.insert(0, 'target', target_label)
        fdr_rows.append(fdr)
        comp = compare_correction_methods(observed_for_stats, target_label=target_label, thresholds=thresholds)
        comparison_rows.append(comp)
    return pd.concat(fdr_rows, ignore_index=True), pd.concat(comparison_rows, ignore_index=True)


def find_negative_controls_from_dependency(target_matrix, dependency_df, max_abs_mean_dependency=0.1,
                                            min_unique_frac=0.5, n=None):
    """Define candidate negative-control genes from THIS dataset's own
    dependency layer, sidestepping the NEGv1/proteomics coverage mismatch."""
    common_genes = target_matrix.index.intersection(dependency_df.index)
    mean_dep = dependency_df.loc[common_genes].mean(axis=1)
    candidates = list(mean_dep[mean_dep.abs() <= max_abs_mean_dependency].index)
    return find_covered_genes_with_variance(candidates, target_matrix, min_unique_frac=min_unique_frac, n=n)


def run_random_target_panel(layer_matrix, n_targets=10, min_unique_frac=0.5, n=None, mode='max',
                             n_permutations=200, seed=42, thresholds=None, active_cell_lines=None,
                             random_state=42, exclude_genes=None):
    """Permutation validation across genes sampled completely at random --
    no essentiality claim, isolates whether cohort construction itself
    (regardless of biology) tends to produce lineage-driven hits."""
    candidates = find_covered_genes_with_variance(list(layer_matrix.index), layer_matrix, min_unique_frac=min_unique_frac)
    if exclude_genes:
        candidates = [g for g in candidates if g not in exclude_genes]
    if len(candidates) < n_targets:
        raise ValueError(f"Only {len(candidates)} genes pass the variance filter; requested {n_targets}.")
    rng = np.random.default_rng(random_state)
    chosen = list(rng.choice(candidates, size=n_targets, replace=False))
    fdr_curves, comparisons = run_multi_target_validation(
        chosen, layer_matrix, n=n, mode=mode, n_permutations=n_permutations,
        seed=seed, thresholds=thresholds, active_cell_lines=active_cell_lines,
    )
    summary = (fdr_curves.groupby('nominal_p_threshold')['empirical_fdr']
               .agg(median='median', mean='mean', q25=lambda s: s.quantile(0.25),
                    q75=lambda s: s.quantile(0.75), min='min', max='max', n_genes='count')
               .reset_index())
    return {'targets': chosen, 'fdr_curves': fdr_curves, 'comparisons': comparisons, 'summary': summary}


def run_negative_control_panel(layer_matrix, n_targets=10, min_unique_frac=0.5, n=None, mode='max',
                                n_permutations=200, seed=42, thresholds=None, active_cell_lines=None,
                                validation_dir=None, random_state=42, backfill=True, dependency_df=None):
    """Permutation validation across a panel of NEGv1-derived non-essential
    genes, with dependency-derived backfill if NEGv1 coverage is too small."""
    _, neg_list = load_essential_gene_lists(validation_dir)
    strict_candidates = find_covered_genes_with_variance(neg_list, layer_matrix, min_unique_frac=min_unique_frac)
    rng = np.random.default_rng(random_state)

    chosen, passed_filter, source = [], {}, {}
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
            dep_candidates = find_negative_controls_from_dependency(layer_matrix, dependency_df, min_unique_frac=min_unique_frac)
            dep_pool = [g for g in dep_candidates if g not in chosen]
            n_from_dep = min(n_needed, len(dep_pool))
            if n_from_dep > 0:
                dep_chosen = list(rng.choice(dep_pool, size=n_from_dep, replace=False))
                print(f"Only {len(strict_candidates)} NEGv1 genes passed the variance filter; "
                      f"adding {n_from_dep} dependency-derived non-essential candidates: {dep_chosen}")
                chosen += dep_chosen
                passed_filter.update({g: True for g in dep_chosen})
                source.update({g: 'dependency_derived' for g in dep_chosen})
                n_needed = n_targets - len(chosen)

        if n_needed > 0:
            if not backfill:
                raise ValueError(f"Only {len(chosen)} candidates found; requested {n_targets} and backfill=False.")
            all_covered = find_covered_genes(neg_list, layer_matrix.index)
            remaining_pool = [g for g in all_covered if g not in chosen]
            if len(remaining_pool) < n_needed:
                raise ValueError(f"Only {len(chosen)} candidates found in total -- cannot reach n_targets={n_targets}.")
            backfilled = list(rng.choice(remaining_pool, size=n_needed, replace=False))
            print(f"Backfilling {n_needed} more from full NEGv1-covered pool WITHOUT variance screen: {backfilled}")
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
        return (df.groupby('nominal_p_threshold')['empirical_fdr']
                .agg(median='median', mean='mean', q25=lambda s: s.quantile(0.25),
                     q75=lambda s: s.quantile(0.75), min='min', max='max', n_genes='count')
                .reset_index())

    summary = _summarize(fdr_curves)
    clean_only = fdr_curves[fdr_curves['passed_variance_filter']]
    summary_clean_only = _summarize(clean_only) if len(clean_only) else None

    return {'targets': chosen, 'passed_variance_filter': passed_filter, 'source': source,
            'fdr_curves': fdr_curves, 'comparisons': comparisons,
            'summary': summary, 'summary_clean_only': summary_clean_only}


def run_phase0_tissue_comparison(layer_matrix, dependency_df, real_targets, n_permutations=200, seed=42,
                                  n_panel_targets=10, n_tissues=2, min_tissue_size=20):
    """Real targets vs essentiality-curated panel vs random panel, pan-cancer
    AND restricted to each of the n_tissues largest tissue lineages."""
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
    return pd.DataFrame(cols), results


_PANEL_COLORS = {
    'Real targets': '#1b9e77', 'Essentiality-curated panel': '#7570b3', 'Random panel': '#d95f02',
}


def plot_real_target_fdr_curves(fdr_curves, output_path, noise_baseline=(0.85, 0.96), title=None):
    """Empirical FDR vs. nominal p threshold, one line per target."""
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
    plt.savefig(output_path, dpi=FIG_DPI, bbox_inches='tight')
    plt.close(fig)
    return output_path


def plot_multipanel_multitissue_comparison(comparison_table, output_path, scope_labels=None, title=None):
    """Faceted line plot (one panel per tissue scope), one line per panel type."""
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
    plt.savefig(output_path, dpi=FIG_DPI, bbox_inches='tight')
    plt.close(fig)
    return output_path


def plot_tissue_bias_ratio(tissue_bias_df, target_categories, output_path, threshold=0.01, title=None):
    """Horizontal bar chart of per-gene tissue-bias ratio."""
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
    plt.savefig(output_path, dpi=FIG_DPI, bbox_inches='tight')
    plt.close(fig)
    return output_path


# ===========================================================================
# SECTION 2: Concordance with real perturbation data
# ===========================================================================

def convert_refseq_to_gene_symbol(fpkm_path, output_path, id_col='Transcript/RepeatID', sep='\t'):
    """One-time preprocessing: RefSeq transcript IDs -> gene symbols via
    mygene.info. Needs `pip install mygene` and outbound internet access --
    run locally if this environment is network-restricted, then point
    get_gene_symbol_fpkm_path at the cached output."""
    try:
        import mygene
    except ImportError:
        sys.exit("This step needs the mygene package: pip install mygene")

    df = pd.read_csv(fpkm_path, sep=sep)
    if id_col not in df.columns:
        sys.exit(f"Column {id_col!r} not found. Actual columns: {df.columns.tolist()}")

    is_refseq = df[id_col].astype(str).str.match(r'^[NX][MR]_')
    print(f"{is_refseq.sum()} of {len(df)} rows look like RefSeq transcript IDs; "
          f"{(~is_refseq).sum()} other rows will be dropped.")
    df = df[is_refseq].copy()

    mg = mygene.MyGeneInfo()
    unique_ids = df[id_col].unique().tolist()
    print(f"Querying mygene.info for {len(unique_ids)} unique RefSeq IDs...")
    results = mg.querymany(unique_ids, scopes='refseq', fields='symbol', species='human')

    id_to_symbol = {}
    n_no_hit, n_no_symbol = 0, 0
    for r in results:
        if r.get('notfound'):
            n_no_hit += 1
        elif 'symbol' in r:
            id_to_symbol[r['query']] = r['symbol']
        else:
            n_no_symbol += 1
    print(f"Mapped {len(id_to_symbol)} of {len(unique_ids)} IDs to a gene symbol "
          f"({n_no_hit} not found, {n_no_symbol} found but no symbol field).")

    df['gene_symbol'] = df[id_col].map(id_to_symbol)
    n_before = len(df)
    df = df.dropna(subset=['gene_symbol'])
    print(f"Dropped {n_before - len(df)} rows with no resolved gene symbol; {len(df)} remain.")

    fpkm_cols = [c for c in df.columns if c.strip().endswith('FPKM')]
    if not fpkm_cols:
        sys.exit(f"No columns ending in 'FPKM' found. Actual columns: {df.columns.tolist()}")

    gene_level = df.groupby('gene_symbol')[fpkm_cols].sum()
    print(f"Collapsed {len(df)} transcript rows into {gene_level.shape[0]} genes.")
    gene_level.to_csv(output_path, sep='\t')
    print(f"Saved {output_path} -- {gene_level.shape[0]} genes x {gene_level.shape[1]} sample columns.")
    return gene_level


def get_gene_symbol_fpkm_path(raw_fpkm_path, converted_cache_path=None, id_col='Transcript/RepeatID'):
    """Ensures a gene-symbol-indexed version of raw_fpkm_path exists on
    disk (converting + caching if needed) and returns its path."""
    if converted_cache_path is None:
        converted_cache_path = str(Path(raw_fpkm_path).with_suffix('')) + '_bysymbol.txt'
    if Path(converted_cache_path).exists():
        print(f"Using cached gene-symbol-converted file: {converted_cache_path}")
        return converted_cache_path
    print(f"No cached conversion found at {converted_cache_path} -- running RefSeq-to-symbol conversion now.")
    convert_refseq_to_gene_symbol(raw_fpkm_path, converted_cache_path, id_col=id_col)
    return converted_cache_path


def get_simko_differential(target_genes, layer_matrix, n=None, mode='max', active_cell_lines=None):
    """SimKO's own simulated-depletion comparison, genome-wide, no prefilter."""
    if active_cell_lines is not None:
        layer_matrix = layer_matrix.filter(active_cell_lines)
    if n is None:
        n = compute_cohort_n(layer_matrix.shape[1])
    class_df = get_classes_by_joint_value(target_genes, layer_matrix, n, mode=mode)
    n_low = (class_df['class'] == 'low').sum()
    n_median = (class_df['class'] == 'median').sum()
    diff_df = get_differentials(class_df, layer_matrix, n_low, n_median)
    target_list = target_genes if isinstance(target_genes, list) else [target_genes]
    return diff_df.loc[~diff_df.index.isin(target_list)]


def get_deficiency_composite(target_gene, abundance_df, expression_df, mutation_df, active_cell_lines=None):
    """
    A per-cell-line composite 'is this protein actually deficient' score for
    ONE target gene: composite = percentile_rank(abundance) *
    percentile_rank(expression) - mutation_count. Deliberately excludes
    Dependency -- dependency is the CONSEQUENCE of deficiency (a fitness
    cost), not evidence for it, so including it here would be circular
    when this score is later used to ask whether deficiency predicts a
    downstream dependency effect.

    Percentile rank runs 0 (lowest raw value) to 1 (highest), matching
    this codebase's existing get_classes_by_joint_value(mode='max')
    convention. Multiplying two percentiles rewards a cell line for being
    low in BOTH abundance and expression simultaneously (two low
    percentiles multiply to an even lower product) rather than averaging
    away a strong signal in one with a weak signal in the other, the way
    a summed z-score would.

    mutation_count (the number of deleterious mutations called for this
    gene in this cell line, unnormalised) is then SUBTRACTED directly.
    Because the rank-product term is bounded to [0, 1], even a single
    mutation can outweigh the entire abundance/expression signal and pull
    a cell line to the very bottom of the ranking regardless of its
    percentile there. This is a deliberate consequence of the formula as
    specified, not a bug -- a confirmed deleterious mutation is arguably
    stronger, more direct evidence of deficiency than sitting in a low
    percentile by chance -- but it means mutation-bearing cell lines can
    dominate cohort selection outright; the printed diagnostics below
    show directly how many cell lines that affects here.

    Lower composite = more deficient, matching
    get_classes_by_joint_value's fixed convention of taking the LOWEST n
    raw values as the 'low' cohort (confirmed directly from its source --
    it always sorts ascending and takes .head(n), regardless of mode).

    Returns a single-row DataFrame (index=[target_gene]) so it can be
    passed anywhere a normal layer_matrix is expected, e.g.
    get_classes_by_joint_value([target_gene], composite_df, n, mode='max').
    """
    shared_cls = [c for c in abundance_df.columns
                  if c in set(expression_df.columns) and c in set(mutation_df.columns)]
    if active_cell_lines is not None:
        shared_cls = [c for c in shared_cls if c in set(active_cell_lines)]
    n_dropped = abundance_df.shape[1] - len(shared_cls)
    if n_dropped:
        print(f"{n_dropped} of {abundance_df.shape[1]} Abundance cell lines lack full "
              f"Expression+Mutation coverage -- restricting to the {len(shared_cls)} cell lines common to all three.")

    abund_vals = abundance_df.loc[target_gene, shared_cls]
    expr_vals = expression_df.loc[target_gene, shared_cls]
    mut_vals = mutation_df.loc[target_gene, shared_cls]

    print(f"\n{target_gene} Mutation-layer values -- confirm this is a deleterious-mutation COUNT, "
          f"as assumed, before trusting the composite:")
    uniq = sorted(mut_vals.unique())
    print(f"  dtype={mut_vals.dtype}, {len(uniq)} unique value(s), "
          f"min={mut_vals.min()}, max={mut_vals.max()}, mean={mut_vals.mean():.3f}")
    print(f"  first few unique values: {uniq[:10]}")
    n_with_mutation = int((mut_vals > 0).sum())
    print(f"  {n_with_mutation} of {len(shared_cls)} cell lines have \u22651 deleterious mutation called "
          f"-- each is pulled toward the bottom of the ranking regardless of its abundance/expression "
          f"percentile, since the rank-product term is bounded to [0,1] and even a single mutation "
          f"subtracts more than that entire range.")

    abund_pct = abund_vals.rank(pct=True)  # 0 (lowest raw value) to 1 (highest)
    expr_pct = expr_vals.rank(pct=True)

    composite = (abund_pct * expr_pct) - mut_vals
    composite_df = pd.DataFrame([composite.values], index=[target_gene], columns=shared_cls)
    return composite_df, shared_cls


def load_real_differential_from_homer_fpkm(path, group_a_condition, group_b_condition,
                                            gene_col_index=0, sep='\t', pseudocount=1.0):
    """For GSE101966-style HOMER analyzeRepeats FPKM output. Matches
    columns to a condition by EXACT match (not substring -- substring
    matching silently pulls double-KO columns into single-KO comparisons)."""
    df = pd.read_csv(path, sep=sep)
    gene_col = df.columns[gene_col_index]

    def _extract_condition(col_name):
        return re.sub(r'[-_]Rep\d+.*$', '', col_name).strip()

    conditions = {c: _extract_condition(c) for c in df.columns}
    detected = sorted(set(conditions.values()) - {_extract_condition(gene_col)})
    print(f"Conditions detected in this file's column headers: {detected}")

    group_a_cols = [c for c, cond in conditions.items() if cond == group_a_condition]
    group_b_cols = [c for c, cond in conditions.items() if cond == group_b_condition]
    print(f"Group A ('{group_a_condition}'): {group_a_cols}")
    print(f"Group B ('{group_b_condition}'): {group_b_cols}")
    if not group_a_cols or not group_b_cols:
        raise ValueError(f"No columns matched one or both conditions exactly. Detected: {detected}")

    df = df.set_index(gene_col)
    df = df[~df.index.duplicated(keep='first')]
    mean_a = df[group_a_cols].mean(axis=1)
    mean_b = df[group_b_cols].mean(axis=1)
    logfc = np.log2((mean_b + pseudocount) / (mean_a + pseudocount))
    return pd.DataFrame({'logfc': logfc})


def load_real_differential_from_proteomics_table(path, gene_col='Gene Names (primary)',
                                                   logfc_col='mean log2 KO/WT', sep=None):
    """
    For an already-processed proteomics KO/WT differential table (e.g. a
    published per-protein log2FC table with a t-test already run) --
    unlike load_real_differential_from_homer_fpkm, this expects one row
    per gene with a ready-made fold-change column, not raw per-sample
    values needing collapsing. Column names are configurable since these
    vary a lot between papers/exports; defaults match a UniProt-style
    proteomics export (e.g. Blijlevens et al. 2021-style PBRM1 KO/WT data)
    with per-cell-line columns like 'median <CELL_LINE> log2 <GENE>
    KO/WT norm' alongside an aggregate 'mean log2 KO/WT' column -- pass
    logfc_col explicitly to use one specific cell line's column instead
    of the aggregate (e.g. for a tissue-matched comparison).

    'Gene Names (primary)' in UniProt-style exports occasionally contains
    more than one space-separated name for a single row; the first token
    is taken as the primary symbol, matching UniProt's own convention for
    which name in that field is "primary".
    """
    if sep is None:
        sep = ',' if str(path).lower().endswith('.csv') else '\t'
    df = pd.read_csv(path, sep=sep)
    if gene_col not in df.columns:
        raise ValueError(f"Column {gene_col!r} not found. Actual columns: {df.columns.tolist()}")
    if logfc_col not in df.columns:
        raise ValueError(f"Column {logfc_col!r} not found. Actual columns: {df.columns.tolist()}")

    genes = df[gene_col].astype(str).str.split().str[0]
    result = pd.DataFrame({'gene': genes, 'logfc': pd.to_numeric(df[logfc_col], errors='coerce')})
    n_before = len(result)
    result = result.dropna(subset=['gene', 'logfc'])
    result = result[result['gene'].str.lower() != 'nan']
    n_dupes = result['gene'].duplicated().sum()
    if n_dupes:
        print(f"{n_dupes} duplicate gene symbols found -- keeping the first occurrence of each.")
    result = result.drop_duplicates(subset='gene', keep='first')
    print(f"Loaded {len(result)} of {n_before} rows with a valid gene symbol and {logfc_col!r} value.")
    return result.set_index('gene')


def load_real_differential_from_raw_intensity_table(path, group_a_cols, group_b_cols, gene_col=2,
                                                      sep=None, pseudocount=1.0):
    """
    For a raw-intensity proteomics table, with or without a header row --
    columns identified purely by position (0-indexed). gene_col is the
    column holding the gene symbol (default 2, matching an
    accession / description / gene-symbol / replicate-intensities
    layout, e.g. ARID1A_RPE_KOs-style exports). group_a_cols/group_b_cols
    are 0-indexed positions of each condition's replicate columns.

    sep: field separator; if not given, auto-detected from the file
    extension (comma for .csv, tab otherwise) -- confirmed necessary
    directly (a real file arrived as .csv, not the tab-separated format
    originally assumed).

    Auto-detects a header row: if the values in group_a_cols/group_b_cols
    for the first row don't parse as numbers, that row is treated as a
    header and skipped -- also confirmed necessary directly (a real
    file's header row, read with no header handling, silently corrupted
    column dtype and broke the mean() call outright rather than
    producing a clear error).

    Computes log2((mean(group_b)+pseudocount)/(mean(group_a)+pseudocount))
    per gene -- group_b is the numerator, so pass the KO replicate columns
    as group_b for a conventional KO/WT log2FC (positive = higher in KO).

    GETTING group_a_cols/group_b_cols BACKWARDS SILENTLY FLIPS THE SIGN OF
    EVERY DOWNSTREAM RESULT. Confirm the actual WT-vs-KO column assignment
    from whatever pipeline/export generated the file before trusting
    anything built on this loader's output -- this function cannot verify
    that for you.
    """
    if sep is None:
        sep = ',' if str(path).lower().endswith('.csv') else '\t'
    probe = pd.read_csv(path, sep=sep, header=None, nrows=1)
    numeric_cols = list(group_a_cols) + list(group_b_cols)
    first_row_numeric = pd.to_numeric(probe.iloc[0, numeric_cols], errors='coerce').notna().all()
    skiprows = 0 if first_row_numeric else 1
    if skiprows:
        print(f"First row doesn't look numeric in the intensity columns -- treating it as a header "
              f"row and skipping it: {list(probe.iloc[0])}")

    df = pd.read_csv(path, sep=sep, header=None, skiprows=skiprows)
    genes = df.iloc[:, gene_col].astype(str).str.strip()
    group_a = df.iloc[:, list(group_a_cols)].apply(pd.to_numeric, errors='coerce').mean(axis=1)
    group_b = df.iloc[:, list(group_b_cols)].apply(pd.to_numeric, errors='coerce').mean(axis=1)
    logfc = np.log2((group_b + pseudocount) / (group_a + pseudocount))

    result = pd.DataFrame({'gene': genes, 'logfc': logfc})
    n_before = len(result)
    result = result.dropna(subset=['gene', 'logfc'])
    result = result[result['gene'].str.strip() != '']
    n_dupes = result['gene'].duplicated().sum()
    if n_dupes:
        print(f"{n_dupes} duplicate gene symbols found -- keeping the first occurrence of each.")
    result = result.drop_duplicates(subset='gene', keep='first')
    print(f"Loaded {len(result)} of {n_before} rows with a valid gene symbol and logfc value "
          f"(group_a cols={list(group_a_cols)}, group_b cols={list(group_b_cols)}).")
    return result.set_index('gene')


def run_concordance_analysis(simko_diff, real_diff, top_n=200, output_prefix='concordance',
                              simko_effect_col='diff', real_effect_col='logfc'):
    """Spearman correlation + top-N hit-list overlap (hypergeometric) --
    genome-wide, no hand-picked subset. Saves a scatter plot."""
    if Path(output_prefix).parent == Path('.'):
        output_prefix = str(SUPP_DIR / output_prefix)
    merged = simko_diff[[simko_effect_col]].join(real_diff[[real_effect_col]], how='inner').dropna()
    n_shared = len(merged)
    if n_shared == 0:
        raise ValueError("Zero genes shared -- check gene identifiers match on both sides.")
    if n_shared < 30:
        print(f"WARNING: only {n_shared} genes shared -- too few for a meaningful correlation.")

    rho, rho_p = spearmanr(merged[simko_effect_col], merged[real_effect_col])
    simko_up = set(merged[simko_effect_col].nlargest(top_n).index)
    simko_down = set(merged[simko_effect_col].nsmallest(top_n).index)
    real_up = set(merged[real_effect_col].nlargest(top_n).index)
    real_down = set(merged[real_effect_col].nsmallest(top_n).index)

    def _enrichment(set_a, set_b, background_n):
        overlap = len(set_a & set_b)
        pval = hypergeom.sf(overlap - 1, background_n, len(set_a), len(set_b))
        expected = len(set_a) * len(set_b) / background_n
        return overlap, expected, pval

    up_overlap, up_expected, up_p = _enrichment(simko_up, real_up, n_shared)
    down_overlap, down_expected, down_p = _enrichment(simko_down, real_down, n_shared)

    print(f"Shared genes: {n_shared}")
    print(f"Spearman rho: {rho:.4f} (p={rho_p:.2e})")
    print(f"Top-{top_n} UP overlap: {up_overlap} observed vs {up_expected:.1f} expected (p={up_p:.2e})")
    print(f"Top-{top_n} DOWN overlap: {down_overlap} observed vs {down_expected:.1f} expected (p={down_p:.2e})")

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(merged[real_effect_col], merged[simko_effect_col], s=6, alpha=0.35, color='#1b7a9c')
    ax.axhline(0, color='grey', linewidth=0.6)
    ax.axvline(0, color='grey', linewidth=0.6)
    ax.set_xlabel('Real experiment effect size (log2FC)')
    ax.set_ylabel("SimKO simulated-depletion effect size ('diff')")
    ax.set_title(f'Genome-wide concordance (n={n_shared})\nSpearman \u03c1={rho:.3f}, p={rho_p:.2e}')
    plt.tight_layout()
    scatter_path = f'{output_prefix}_scatter.png'
    plt.savefig(scatter_path, dpi=FIG_DPI, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved {scatter_path}")
    merged.to_csv(f'{output_prefix}_merged.csv')

    return {'n_shared': n_shared, 'spearman_rho': rho, 'spearman_p': rho_p, 'top_n': top_n,
            'up_overlap': up_overlap, 'up_expected': up_expected, 'up_p': up_p,
            'down_overlap': down_overlap, 'down_expected': down_expected, 'down_p': down_p, 'merged': merged}


def plot_cross_layer_grid(perm_abund_protein, perm_expr_protein,
                           output_path=FIG_DIR / 'fig5_cross_layer_grid.png'):
    """
    Figure 5: on REAL PROTEIN data (RPE1 ARID1A-knockout proteomics), which SimKO layer
    recovers the genes that change?

    Left column: SimKO's Abundance layer (cohorts split on protein abundance), i.e. protein vs protein.
    Right column: SimKO's Expression layer (cohorts split on RNA level), tested against the same real
    protein data. Rows: genes predicted UP / DOWN. For each panel, SimKO's 200 most-changed genes are
    compared with the real experiment's 200 most-changed genes and the genes in both lists are counted.
    Grey bars = that count when the cell-line groups are shuffled at random (200 shuffles), i.e. chance.
    Coloured line = SimKO's actual count. No new statistics: the nulls and p-values come from the
    permutation results already in memory. (The HCT116 RNA-seq comparison is in Figure 4, panel B.)
    """
    cols = [
        ('SimKO Abundance layer (protein)\nvs real RPE1 knockout proteomics\n[protein vs protein]', perm_abund_protein),
        ('SimKO Expression layer (RNA)\nvs real RPE1 knockout proteomics\n[RNA vs protein]', perm_expr_protein),
    ]
    rows = [('Genes predicted UP\n(higher after knockout)', 'observed_up', 'null_ups', 'empirical_p_up'),
            ('Genes predicted DOWN\n(lower after knockout)', 'observed_down', 'null_downs', 'empirical_p_down')]
    GREEN, ORANGE = '#1b7a3c', '#d9822b'
    fig, axes = plt.subplots(2, 2, figsize=(9.5, 7.8))
    for r, (label, obs_key, null_key, p_key) in enumerate(rows):
        for c, (title, perm) in enumerate(cols):
            ax = axes[r, c]
            null = np.asarray(perm[null_key]); obs = perm[obs_key]; p = perm[p_key]
            sig = p < 0.05
            hi = max(null.max(), obs) + 2
            bins = np.arange(min(null.min(), obs) - 1, hi + 1, 1)
            ax.hist(null, bins=bins, color='#c8c8c8', edgecolor='white')
            ax.axvline(obs, color=GREEN if sig else ORANGE, linewidth=3)
            ax.set_xlim(bins[0], bins[-1])
            ax.set_ylim(0, ax.get_ylim()[1] * 1.75)   # headroom so the text never sits on the bars
            on_right = obs > (bins[0] + bins[-1]) / 2
            verdict = 'REAL SIGNAL' if sig else 'not above chance'
            ax.text(0.04 if on_right else 0.96, 0.96,
                    f"SimKO: {obs} shared genes\nChance: {null.mean():.0f} on average\n"
                    f"p = {p:.3f}\n{verdict}",
                    transform=ax.transAxes, ha='left' if on_right else 'right', va='top', fontsize=10,
                    color=GREEN if sig else '#555555', fontweight='bold' if sig else 'normal',
                    linespacing=1.35)
            if r == 0:
                ax.set_title(title, fontsize=10.5, fontweight='bold')
            if c == 0:
                ax.set_ylabel(f'{label}\n\nNumber of random shuffles', fontsize=10)
            ax.set_xlabel('Genes in common with the real experiment\n(out of the top 200 in each list)', fontsize=9)
            ax.spines[['top', 'right']].set_visible(False)
    fig.suptitle('On real protein data, which SimKO layer recovers the genes that change after ARID1A loss?',
                 fontsize=12.5, fontweight='bold', y=0.995)
    handles = [Patch(facecolor='#c8c8c8', label='Chance: overlap from 200 random shuffles of the cell-line groups'),
               Line2D([0], [0], color=GREEN, lw=3, label='SimKO\'s actual overlap: above chance (p < 0.05)'),
               Line2D([0], [0], color=ORANGE, lw=3, label='SimKO\'s actual overlap: within chance')]
    fig.legend(handles=handles, loc='lower center', ncol=1, fontsize=10, frameon=False, bbox_to_anchor=(0.5, 0.0))
    fig.text(0.5, 0.115, 'How to read: the further the coloured line sits to the right of the grey bars, the more '
             'SimKO\'s prediction agrees with the real experiment.', ha='center', fontsize=10, style='italic')
    plt.tight_layout(rect=(0, 0.17, 1, 0.97))
    plt.savefig(output_path, dpi=FIG_DPI, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved {output_path}")
    return output_path


def get_cohort_dependency_values(target_gene, downstream_gene, cohort_matrix, dependency_df,
                                  n=None, mode='max', active_cell_lines=None):
    """
    Companion to test_single_gene_dependency_shift: returns the raw
    per-cell-line downstream_gene Dependency values for the SAME low/median
    cohort split, for plotting (e.g. a box/strip plot) rather than testing.
    Recomputes only the cohort split (cheap, deterministic) -- does not
    repeat any permutation work.
    """
    if active_cell_lines is not None:
        cohort_matrix = cohort_matrix.filter(active_cell_lines)
    if n is None:
        n = compute_cohort_n(cohort_matrix.shape[1])
    target_list = [target_gene] if isinstance(target_gene, str) else list(target_gene)
    class_df = get_classes_by_joint_value(target_list, cohort_matrix, n, mode=mode)
    low_cls = list(class_df.loc[class_df['class'] == 'low'].index)
    median_cls = list(class_df.loc[class_df['class'] == 'median'].index)
    dep_vals = dependency_df.loc[downstream_gene]
    return pd.DataFrame({
        'cell_line': low_cls + median_cls,
        'cohort': ['low'] * len(low_cls) + ['median'] * len(median_cls),
        'dependency': [dep_vals.get(c, np.nan) for c in low_cls + median_cls],
    })


def plot_benchmark_figure(panels, output_path=FIG_DIR / 'fig4_benchmark_vs_correlation.png', title=None):
    """
    Figure 4: SimKO's cohort method vs a plain correlation, benchmarked against
    real knockout data. Bars show -log10(empirical permutation p-value) for the
    overall correlation, the top-N UP overlap and the top-N DOWN overlap, plus
    the combined best-of-both result; the dotted line marks p=0.05.

    panels: list of dicts, one per panel, each with keys
        'title'        -- panel title
        'simko'        -- dict returned by permutation_null_for_concordance
        'correlation'  -- dict returned by permutation_null_for_correlation_concordance
        'combined'     -- dict returned by test_combined_strategy

    A permutation p-value cannot be smaller than 1/n_permutations, so any bar
    that has hit that floor is drawn at the floor, hatched, and labelled
    '<floor' rather than implying a precision the permutations cannot give.
    """
    colors = {'simko': '#1b9e77', 'correlation': '#d95f02', 'combined': '#7570b3'}
    fig, axes = plt.subplots(1, len(panels), figsize=(6.2 * len(panels), 4.8), squeeze=False)
    for ax, panel in zip(axes[0], panels):
        simko, corr, comb = panel['simko'], panel['correlation'], panel['combined']
        n_perm = len(simko['null_rhos'])
        floor = 1.0 / n_perm
        comb_floor = 1.0 / len(comb['null_combined'])

        def _bar(x, p, width, color, label=None, p_floor=floor):
            hit_floor = p <= p_floor
            value = -np.log10(max(p, p_floor))
            ax.bar(x, value, width, color=color, label=label,
                   hatch='//' if hit_floor else None, edgecolor='white' if hit_floor else color)
            text = '<%g' % p_floor if hit_floor else '%.2f' % value
            ax.text(x, value + 0.04, text, ha='center', va='bottom', fontsize=8)

        w = 0.36
        groups = [('Overall\ncorrelation', 'empirical_p_rho'),
                  ('Top-200\nUP overlap', 'empirical_p_up'),
                  ('Top-200\nDOWN overlap', 'empirical_p_down')]
        for i, (name, key) in enumerate(groups):
            _bar(i - w / 2, simko[key], w, colors['simko'])
            _bar(i + w / 2, corr[key], w, colors['correlation'])
        _bar(3, comb['empirical_p_combined'], 0.5, colors['combined'], p_floor=comb_floor)

        ax.axhline(-np.log10(0.05), color='grey', linestyle=':', linewidth=1)
        ax.set_xticks(range(4))
        ax.set_xticklabels([g[0] for g in groups] + ['Combined\n(best of both)'])
        ax.set_ylabel('-log10(empirical permutation p-value)')
        ax.set_title(panel['title'])
        ax.set_ylim(0, max(-np.log10(floor), -np.log10(comb_floor)) * 1.18)
        handles = [Patch(facecolor=colors['simko'], label='SimKO cohort method'),
                   Patch(facecolor=colors['correlation'], label='Plain correlation'),
                   Patch(facecolor=colors['combined'], label='Combined (best of both)'),
                   Line2D([0], [0], color='grey', linestyle=':', label='p = 0.05'),
                   Patch(facecolor='white', edgecolor='grey', hatch='//',
                         label='at permutation floor (p < 1/%d)' % n_perm)]
        ax.legend(handles=handles, loc='upper left', frameon=False)
    if title:
        fig.suptitle(title)
    plt.tight_layout()
    plt.savefig(output_path, dpi=FIG_DPI, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved {output_path}")
    return output_path


def plot_brca_parp1_boxplot(brca1_df, brca2_df, output_path=FIG_DIR / 'fig6_brca_parp1_boxplot.png'):
    """
    Box + strip plot of PARP1 dependency by BRCA1/BRCA2 cohort (low vs
    median), pan-cancer -- makes the BRCA1/BRCA2->PARP1 case study's
    numeric result (a diff and a p-value) visually concrete: the actual
    per-cell-line spread, not just a summary statistic.
    """
    fig, axes = plt.subplots(1, 2, figsize=(8, 5), sharey=True)
    rng = np.random.default_rng(0)
    for ax, df, gene in zip(axes, [brca1_df, brca2_df], ['BRCA1', 'BRCA2']):
        order = ['low', 'median']
        data = [df.loc[df['cohort'] == c, 'dependency'].dropna().values for c in order]
        ax.boxplot(data, positions=[0, 1], widths=0.5, showfliers=False,
                   patch_artist=True, boxprops=dict(facecolor='#e8e8e8'))
        for i, vals in enumerate(data):
            jitter = rng.normal(0, 0.04, len(vals))
            ax.scatter(np.full(len(vals), i) + jitter, vals, s=20, alpha=0.7, color='#1b7a9c', zorder=3)
        ax.set_xticks([0, 1])
        ax.set_xticklabels([f'Low {gene}\n(n={len(data[0])})', f'Median {gene}\n(n={len(data[1])})'])
        ax.set_title(gene)
        ax.grid(axis='y', alpha=0.25, linewidth=0.5)
    axes[0].set_ylabel('PARP1 dependency score (0\u20131, higher = more essential)')
    plt.tight_layout()
    plt.savefig(output_path, dpi=FIG_DPI, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved {output_path}")
    return output_path


def permutation_null_for_concordance(target_genes, layer_matrix, real_diff, n, active_cell_lines=None,
                                      n_permutations=200, seed=42, mode='max', top_n=200,
                                      simko_effect_col='diff', real_effect_col='logfc'):
    """Whether SimKO's concordance with real data is distinguishable from
    chance -- shuffles cohort labels (NOT based on the target's real
    ranking), recomputes the genome-wide differential each time."""
    if active_cell_lines is not None:
        layer_matrix = layer_matrix.filter(active_cell_lines)

    class_df = get_classes_by_joint_value(target_genes, layer_matrix, n, mode=mode)
    n_low = (class_df['class'] == 'low').sum()
    n_median = (class_df['class'] == 'median').sum()
    observed_diff = get_differentials(class_df, layer_matrix, n_low, n_median)
    target_list = target_genes if isinstance(target_genes, list) else [target_genes]
    observed_diff = observed_diff.loc[~observed_diff.index.isin(target_list)]
    merged_obs = observed_diff[[simko_effect_col]].join(real_diff[[real_effect_col]], how='inner').dropna()
    observed_rho, _ = spearmanr(merged_obs[simko_effect_col], merged_obs[real_effect_col])

    real_up_genes = set(merged_obs[real_effect_col].nlargest(top_n).index)
    real_down_genes = set(merged_obs[real_effect_col].nsmallest(top_n).index)
    observed_up = len(set(merged_obs[simko_effect_col].nlargest(top_n).index) & real_up_genes)
    observed_down = len(set(merged_obs[simko_effect_col].nsmallest(top_n).index) & real_down_genes)

    rng = np.random.default_rng(seed)
    all_cls = list(layer_matrix.columns)
    gene_index = layer_matrix.index
    working = layer_matrix.to_numpy()

    aligned_real = real_diff[real_effect_col].reindex(gene_index)
    valid_mask = aligned_real.notna().values & (~gene_index.isin(target_list))
    aligned_real_vals = aligned_real.values
    valid_gene_index = gene_index[valid_mask]

    null_rhos, null_ups, null_downs = [], [], []
    for i in range(n_permutations):
        shuffled = rng.permutation(len(all_cls))
        low_idx = shuffled[:n]
        median_idx = shuffled[n:2 * n]
        diff = working[:, low_idx].mean(axis=1) - working[:, median_idx].mean(axis=1)
        diff_valid = diff[valid_mask]
        rho, _ = spearmanr(diff_valid, aligned_real_vals[valid_mask])
        null_rhos.append(rho)
        diff_series = pd.Series(diff_valid, index=valid_gene_index)
        perm_up_genes = set(diff_series.nlargest(top_n).index)
        perm_down_genes = set(diff_series.nsmallest(top_n).index)
        null_ups.append(len(perm_up_genes & real_up_genes))
        null_downs.append(len(perm_down_genes & real_down_genes))

    null_rhos, null_ups, null_downs = np.array(null_rhos), np.array(null_ups), np.array(null_downs)
    emp_p_rho = float(np.mean(np.abs(null_rhos) >= abs(observed_rho)))
    emp_p_up = float(np.mean(null_ups >= observed_up))
    emp_p_down = float(np.mean(null_downs >= observed_down))

    print(f"Observed rho (n={n}): {observed_rho:.4f}")
    print(f"  Null rho: mean={null_rhos.mean():.4f}, sd={null_rhos.std():.4f}  Empirical two-sided p={emp_p_rho:.4f}")
    print(f"Observed top-{top_n} UP overlap: {observed_up}  Null mean={null_ups.mean():.2f}  Empirical p={emp_p_up:.4f}")
    print(f"Observed top-{top_n} DOWN overlap: {observed_down}  Null mean={null_downs.mean():.2f}  Empirical p={emp_p_down:.4f}")

    return {'n': n, 'observed_rho': observed_rho, 'null_rhos': null_rhos, 'empirical_p_rho': emp_p_rho,
            'observed_up': observed_up, 'null_ups': null_ups, 'empirical_p_up': emp_p_up,
            'observed_down': observed_down, 'null_downs': null_downs, 'empirical_p_down': emp_p_down}


# ===========================================================================
# SECTION 3: Benchmark against the simplest alternative (plain correlation)
# ===========================================================================

def get_correlation_differential(target_gene, layer_matrix, active_cell_lines=None, method='spearman',
                                  target_matrix=None):
    """Plain correlation, target vs every other gene, no cohort split.

    target_matrix: optional. If given, the target's values come from THIS matrix (e.g. the
    Abundance matrix) while the genes it is correlated against come from layer_matrix (e.g.
    the Expression matrix), using the cell lines the two share. Default None = both from
    layer_matrix, as before."""
    if active_cell_lines is not None:
        layer_matrix = layer_matrix.filter(active_cell_lines)
    if target_matrix is not None:
        common = [c for c in layer_matrix.columns if c in target_matrix.columns]
        layer_matrix = layer_matrix[common]
        target_vals = target_matrix.loc[target_gene, common].values
        other = layer_matrix.drop(index=target_gene) if target_gene in layer_matrix.index else layer_matrix
    else:
        target_vals = layer_matrix.loc[target_gene].values
        other = layer_matrix.drop(index=target_gene)

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
    corr = cov / (np.sqrt((y ** 2).sum(axis=1)) * np.sqrt((x ** 2).sum()))
    with np.errstate(divide='ignore', invalid='ignore'):
        t_stat = corr * np.sqrt((n - 2) / (1 - corr ** 2))
    pvals = 2 * t_dist.sf(np.abs(t_stat), df=n - 2)
    # SIGN CONVENTION: SimKO's 'diff' is (low-target cohort) minus (median cohort), so positive = the gene is HIGHER
    # when the target is low = predicted to RISE after a knockout. A gene positively correlated with the target does
    # the opposite (it falls when the target is low), so the correlation must be negated to mean the same thing.
    # Without this, the correlation method's UP and DOWN overlaps and its signed rho are scored upside down.
    return pd.DataFrame({'diff': -corr, 'p': pvals}, index=other.index)


def permutation_null_for_correlation_concordance(target_gene, layer_matrix, real_diff, active_cell_lines=None,
                                                   n_permutations=200, seed=42, top_n=200, method='spearman',
                                                   real_effect_col='logfc', target_matrix=None):
    """Correlation-method analog of permutation_null_for_concordance --
    shuffles which cell line's target value pairs with which cell line's
    downstream expression, not cohort labels (there's no cohort here)."""
    if active_cell_lines is not None:
        layer_matrix = layer_matrix.filter(active_cell_lines)
    if target_matrix is not None:      # target value from a different layer (e.g. protein abundance vs gene RNA)
        common = [c for c in layer_matrix.columns if c in target_matrix.columns]
        layer_matrix = layer_matrix[common]
        target_matrix = target_matrix[common]

    observed_diff = get_correlation_differential(target_gene, layer_matrix, method=method,
                                                 target_matrix=target_matrix)
    merged_obs = observed_diff[['diff']].join(real_diff[[real_effect_col]], how='inner').dropna()
    observed_rho, _ = spearmanr(merged_obs['diff'], merged_obs[real_effect_col])

    real_up_genes = set(merged_obs[real_effect_col].nlargest(top_n).index)
    real_down_genes = set(merged_obs[real_effect_col].nsmallest(top_n).index)
    observed_up = len(set(merged_obs['diff'].nlargest(top_n).index) & real_up_genes)
    observed_down = len(set(merged_obs['diff'].nsmallest(top_n).index) & real_down_genes)

    rng = np.random.default_rng(seed)
    all_cls = list(layer_matrix.columns)
    other = layer_matrix.drop(index=target_gene) if target_gene in layer_matrix.index else layer_matrix
    other_ranks = other.rank(axis=1).to_numpy() if method == 'spearman' else other.to_numpy()
    gene_index = other.index

    aligned_real = real_diff[real_effect_col].reindex(gene_index)
    valid_mask = aligned_real.notna().values
    aligned_real_vals = aligned_real.values
    valid_gene_index = gene_index[valid_mask]

    target_vals = (target_matrix if target_matrix is not None else layer_matrix).loc[target_gene].values
    target_ranks_base = pd.Series(target_vals).rank().values if method == 'spearman' else target_vals

    null_rhos, null_ups, null_downs = [], [], []
    for i in range(n_permutations):
        shuffled_idx = rng.permutation(len(all_cls))
        target_ranks = target_ranks_base[shuffled_idx]
        x = target_ranks - target_ranks.mean()
        y = other_ranks - other_ranks.mean(axis=1, keepdims=True)
        cov = (y * x).sum(axis=1)
        corr = cov / (np.sqrt((y ** 2).sum(axis=1)) * np.sqrt((x ** 2).sum()))
        corr_valid = corr[valid_mask]
        corr_valid = -corr_valid   # same sign convention as get_correlation_differential (positive = predicted to rise)
        rho, _ = spearmanr(corr_valid, aligned_real_vals[valid_mask])
        null_rhos.append(rho)
        corr_series = pd.Series(corr_valid, index=valid_gene_index)
        null_ups.append(len(set(corr_series.nlargest(top_n).index) & real_up_genes))
        null_downs.append(len(set(corr_series.nsmallest(top_n).index) & real_down_genes))

    null_rhos, null_ups, null_downs = np.array(null_rhos), np.array(null_ups), np.array(null_downs)
    emp_p_rho = float(np.mean(np.abs(null_rhos) >= abs(observed_rho)))
    emp_p_up = float(np.mean(null_ups >= observed_up))
    emp_p_down = float(np.mean(null_downs >= observed_down))

    print(f"[Correlation] Observed rho: {observed_rho:.4f}  Empirical two-sided p={emp_p_rho:.4f}")
    print(f"[Correlation] Observed UP overlap: {observed_up}  Empirical p={emp_p_up:.4f}")
    print(f"[Correlation] Observed DOWN overlap: {observed_down}  Empirical p={emp_p_down:.4f}")

    return {'observed_rho': observed_rho, 'null_rhos': null_rhos, 'empirical_p_rho': emp_p_rho,
            'observed_up': observed_up, 'null_ups': null_ups, 'empirical_p_up': emp_p_up,
            'observed_down': observed_down, 'null_downs': null_downs, 'empirical_p_down': emp_p_down}


def test_combined_strategy(perm_simko, perm_corr):
    """Best-of-both-directions combined significance test, reusing the
    already-computed null arrays rather than running new permutations."""
    up_choice = 'simko' if perm_simko['empirical_p_up'] <= perm_corr['empirical_p_up'] else 'correlation'
    down_choice = 'simko' if perm_simko['empirical_p_down'] <= perm_corr['empirical_p_down'] else 'correlation'
    up_source = perm_simko if up_choice == 'simko' else perm_corr
    down_source = perm_simko if down_choice == 'simko' else perm_corr

    observed_combined = up_source['observed_up'] + down_source['observed_down']
    null_combined = up_source['null_ups'] + down_source['null_downs']
    emp_p_combined = float(np.mean(null_combined >= observed_combined))

    print(f"Best for UP: {up_choice} (p={up_source['empirical_p_up']:.4f}, {up_source['observed_up']} hits)")
    print(f"Best for DOWN: {down_choice} (p={down_source['empirical_p_down']:.4f}, {down_source['observed_down']} hits)")
    print(f"Combined observed: {observed_combined}  Combined null mean: {null_combined.mean():.2f}")
    print(f"Combined empirical p: {emp_p_combined:.4f}")

    return {'up_choice': up_choice, 'down_choice': down_choice, 'observed_combined': observed_combined,
            'null_combined': null_combined, 'empirical_p_combined': emp_p_combined}


# ===========================================================================
# SECTION 4: Power check (does a known-real signal survive a small cohort?)
# ===========================================================================

def run_power_check(expression, target='PCNA', n_small=8, n_large=20, n_permutations=200, seed=42):
    """
    Compares empirical FDR for a known-strong real signal (PCNA, Expression
    layer) at a small cohort size (n_small, e.g. matching a tissue-
    restricted run) against the standard cohort size (n_large) -- a
    reference point for whether a null/weak result at small n is more
    likely a genuine null or just a power problem. Not specific to PCNA or
    ARID1B -- pass any target/comparison n's relevant to whatever small-
    cohort result you're trying to contextualise.
    """
    print(f"{'=' * 70}\n{target}, Expression layer, standard n={n_large}\n{'=' * 70}")
    result_large = run_permutation_validation([target], expression, n=n_large, n_permutations=n_permutations, seed=seed)
    print(result_large['fdr_curve'].to_string(index=False))

    print(f"\n{'=' * 70}\n{target}, Expression layer, forced n={n_small}\n{'=' * 70}")
    result_small = run_permutation_validation([target], expression, n=n_small, n_permutations=n_permutations, seed=seed)
    print(result_small['fdr_curve'].to_string(index=False))

    print(f"\n{'=' * 70}\nComparison at p<0.01\n{'=' * 70}")
    row_large = result_large['fdr_curve'].query('nominal_p_threshold == 0.01').iloc[0]
    row_small = result_small['fdr_curve'].query('nominal_p_threshold == 0.01').iloc[0]
    print(f"n={n_large}: empirical FDR = {row_large['empirical_fdr']:.4f} ({row_large['observed_hits']} hits)")
    print(f"n={n_small}:  empirical FDR = {row_small['empirical_fdr']:.4f} ({row_small['observed_hits']} hits)")
    print(
        f"\nIf n={n_small}'s empirical FDR is still clearly well below the ~0.85-0.96 pure-noise "
        f"ceiling (Phase 0's synthetic testing), the known-real signal survives the reduced cohort "
        f"reasonably well -- low power alone probably isn't the explanation for a null result at "
        f"this n. If it's risen sharply toward that ceiling, low power is a genuine, mundane "
        f"explanation worth taking seriously before concluding anything else."
    )
    return result_large, result_small


# ===========================================================================
# SECTION 5: Synthetic lethality case study (BRCA1/BRCA2 -> PARP1)
# ===========================================================================

def test_single_gene_dependency_shift(target_gene, downstream_gene, cohort_matrix, dependency_df,
                                       n=None, mode='max', n_permutations=200, seed=42,
                                       active_cell_lines=None):
    """
    Whether downstream_gene's Dependency-layer score differs between a
    natural low-target_gene cohort (cohort_matrix, typically Abundance)
    and a median-target_gene cohort, permutation-checked. One-sided
    empirical p tests specifically for HIGHER dependency score in the low
    cohort (this codebase's confirmed 0-1 probability-of-dependency
    convention: higher = more essential).

    target_gene accepts either a single gene (str) or a list of genes,
    for a joint low-cohort across several targets at once (e.g. "low in
    both HSPA4 and HSPA4L") -- str and list are handled identically
    throughout except for the target_label used in printed output.
    """
    target_list = [target_gene] if isinstance(target_gene, str) else list(target_gene)
    target_label = target_gene if isinstance(target_gene, str) else '+'.join(target_list)

    if active_cell_lines is not None:
        cohort_matrix = cohort_matrix.filter(active_cell_lines)
    if n is None:
        n = compute_cohort_n(cohort_matrix.shape[1])

    # Guard against too few cell lines to form two cohorts of size n at
    # all (e.g. active_cell_lines resolves to 0 or very few columns for a
    # given tissue) -- caught directly from a real failure mode: with 0
    # cell lines, observed_diff and every permutation's null value become
    # NaN, and numpy's NaN >= NaN silently evaluates to False, making
    # np.mean(null_diffs >= observed_diff) compute to 0.0 -- which then
    # reads as "empirical p=0.0000, maximally significant" instead of the
    # true answer, "undefined". A loud, explicit error here is much safer
    # than that silent false-positive.
    if n < 1 or cohort_matrix.shape[1] < 2 * n:
        raise ValueError(
            f"Not enough cell lines to form two cohorts of size n={n} ({cohort_matrix.shape[1]} available) "
            f"-- likely an empty or near-empty active_cell_lines filter (e.g. a tissue with 0 matches). "
            f"Fix the filter rather than let this silently produce a NaN-based, spuriously 'significant' result."
        )

    # Guard against duplicate index entries -- a duplicated downstream_gene
    # row makes dependency_df.loc[downstream_gene] return a DataFrame
    # instead of a Series, silently breaking the reindex/to_numpy step
    # below rather than raising (this produced a near-zero null SD in an
    # earlier PCNA self-check before being traced to this).
    for g in target_list:
        n_target_rows = (cohort_matrix.index == g).sum()
        if n_target_rows != 1:
            raise ValueError(f"{g!r} appears {n_target_rows} times in cohort_matrix's index (expected 1).")
    n_downstream_rows = (dependency_df.index == downstream_gene).sum()
    if n_downstream_rows != 1:
        raise ValueError(f"{downstream_gene!r} appears {n_downstream_rows} times in dependency_df's index (expected 1).")
    if downstream_gene in target_list:
        raise ValueError(f"{downstream_gene!r} is both a target and the downstream readout gene.")

    class_df = get_classes_by_joint_value(target_list, cohort_matrix, n, mode=mode)
    n_low = (class_df['class'] == 'low').sum()
    n_median = (class_df['class'] == 'median').sum()
    diff_df = get_differentials(class_df, dependency_df, n_low, n_median)
    observed_diff = diff_df.loc[downstream_gene, 'diff']
    observed_p = diff_df.loc[downstream_gene, 'p']

    rng = np.random.default_rng(seed)
    all_cls = list(cohort_matrix.columns)
    dep_vals = dependency_df.loc[downstream_gene].reindex(all_cls).to_numpy()
    valid = ~np.isnan(dep_vals)

    null_diffs = []
    for i in range(n_permutations):
        shuffled = rng.permutation(np.where(valid)[0])
        if len(shuffled) < 2 * n:
            continue
        low_idx, median_idx = shuffled[:n], shuffled[n:2 * n]
        null_diffs.append(dep_vals[low_idx].mean() - dep_vals[median_idx].mean())
    null_diffs = np.array(null_diffs)

    emp_p_one_sided = float(np.mean(null_diffs >= observed_diff))
    emp_p_two_sided = float(np.mean(np.abs(null_diffs) >= abs(observed_diff)))

    print(f"{target_label} low (n={n_low}) vs median (n={n_median}) {downstream_gene} dependency:")
    print(f"  Observed diff: {observed_diff:.4f} (positive = more essential in low-{target_label} cells)")
    print(f"  Nominal p: {observed_p:.2e}")
    print(f"  Null diff ({len(null_diffs)} perms): mean={null_diffs.mean():.4f}, sd={null_diffs.std():.4f}")
    print(f"  Empirical one-sided p (predicted direction): {emp_p_one_sided:.4f}")
    print(f"  Empirical two-sided p: {emp_p_two_sided:.4f}")

    return {'target': target_label, 'downstream': downstream_gene, 'observed_diff': observed_diff,
            'observed_p': observed_p, 'null_diffs': null_diffs, 'empirical_p_one_sided': emp_p_one_sided,
            'empirical_p_two_sided': emp_p_two_sided, 'n_low': n_low, 'n_median': n_median}


def joint_permutation_dependency_battery(tests, cohort_matrix, dependency_df, n=None, n_permutations=200, seed=42):
    """
    Family-wise version of test_single_gene_dependency_shift for a BATTERY
    of related (cohort_targets, readout_gene) tests -- e.g. every
    directional pair within a small gene family. Running each test with
    its own independent permutation draw and then correcting the resulting
    p-values (e.g. Benjamini-Hochberg) treats the tests as independent,
    which they are not: tests sharing a cohort definition share that
    cohort's exact low/median split, and tests sharing a readout gene
    share that gene's real dependency values -- both induce correlation a
    post-hoc correction on separately-drawn nulls cannot see or account
    for.

    This instead runs ONE combined permutation loop: each iteration draws
    ONE shuffled low/median split per DISTINCT cohort definition in the
    battery (shared by every test that uses it), and reuses each gene's
    real dependency values throughout -- so within an iteration, the same
    correlation structure the real, unpermuted data has is preserved.

    tests: list of (cohort_targets, readout_gene) tuples, where
    cohort_targets is itself a tuple of one or more gene names (e.g.
    [(('A',), 'B'), (('A','C'), 'B')] for "low-A -> dep(B)" and
    "low-A-and-C -> dep(B)").

    Returns per-test observed diffs and marginal empirical p-values (each
    computed from this same joint null, so already internally consistent
    with each other), plus the direct family-wise answer: how many
    permutations produced as many marginally-significant tests (at
    one-sided p<0.05, the same criterion used for each test individually)
    as were actually observed, purely from this joint null.
    """
    # cohort_matrix (typically Abundance) and dependency_df do not
    # necessarily cover the exact same cell lines. Confirmed necessary
    # directly: without this, a missing cell line's dependency value
    # reindexes to NaN, numpy's mean() of an array containing NaN is
    # itself NaN, and NaN >= observed_diff silently evaluates to False
    # for every permutation -- making every test read as "0% of nulls
    # exceeded observed", i.e. maximally, spuriously significant, on
    # what was actually pure noise.
    shared_cls = [c for c in cohort_matrix.columns if c in set(dependency_df.columns)]
    n_dropped = cohort_matrix.shape[1] - len(shared_cls)
    if n_dropped:
        print(f"{n_dropped} of {cohort_matrix.shape[1]} cohort_matrix cell lines are absent from "
              f"dependency_df -- restricting this battery to the {len(shared_cls)} cell lines common to both.")
    cohort_matrix = cohort_matrix[shared_cls]
    dependency_df = dependency_df[shared_cls]

    if n is None:
        n = compute_cohort_n(cohort_matrix.shape[1])
    all_cls = list(cohort_matrix.columns)
    if cohort_matrix.shape[1] < 2 * n:
        raise ValueError(f"n={n} needs at least {2 * n} cell lines; only {cohort_matrix.shape[1]} available.")

    distinct_cohorts = sorted(set(cohort_targets for cohort_targets, _ in tests))
    for cohort_targets, readout in tests:
        for g in cohort_targets:
            if (cohort_matrix.index == g).sum() != 1:
                raise ValueError(f"{g!r} does not appear exactly once in cohort_matrix's index.")
        if (dependency_df.index == readout).sum() != 1:
            raise ValueError(f"{readout!r} does not appear exactly once in dependency_df's index.")
        if readout in cohort_targets:
            raise ValueError(f"{readout!r} is both a cohort target and the readout gene in the same test.")

    # Real cohort split and observed diff for each distinct cohort / test
    real_splits, observed_diffs = {}, {}
    for cohort_targets in distinct_cohorts:
        class_df = get_classes_by_joint_value(list(cohort_targets), cohort_matrix, n, mode='max')
        real_splits[cohort_targets] = (
            list(class_df.loc[class_df['class'] == 'low'].index),
            list(class_df.loc[class_df['class'] == 'median'].index),
        )
    dep_cache = {readout: dependency_df.loc[readout].reindex(all_cls) for _, readout in tests}
    for cohort_targets, readout in tests:
        low_cls, median_cls = real_splits[cohort_targets]
        observed_diffs[(cohort_targets, readout)] = (
            dep_cache[readout].reindex(low_cls).mean() - dep_cache[readout].reindex(median_cls).mean()
        )

    rng = np.random.default_rng(seed)
    n_cls = len(all_cls)
    null_diffs = {test: np.empty(n_permutations) for test in tests}

    for i in range(n_permutations):
        # One shuffle per DISTINCT cohort, shared by every test using it
        shuffles = {}
        for cohort_targets in distinct_cohorts:
            shuffled = rng.permutation(n_cls)
            shuffles[cohort_targets] = (shuffled[:n], shuffled[n:2 * n])
        for cohort_targets, readout in tests:
            low_idx, median_idx = shuffles[cohort_targets]
            dep_vals = dep_cache[readout].to_numpy()
            null_diffs[(cohort_targets, readout)][i] = dep_vals[low_idx].mean() - dep_vals[median_idx].mean()

    marginal_p = {test: float(np.mean(null_diffs[test] >= observed_diffs[test])) for test in tests}
    observed_hit_count = sum(1 for test in tests if marginal_p[test] < 0.05)

    # For each permutation, how many tests would have been called a "hit"
    # (null diff at least as extreme as the REAL observed diff) purely by
    # chance, in that one internally-consistent draw
    null_hit_counts = np.zeros(n_permutations, dtype=int)
    for test in tests:
        null_hit_counts += (null_diffs[test] >= observed_diffs[test]).astype(int)
    emp_p_family = float(np.mean(null_hit_counts >= observed_hit_count))

    print(f"Joint permutation battery: {len(tests)} tests across {len(distinct_cohorts)} distinct cohorts, "
          f"{n_permutations} joint permutations.")
    print(f"\n{'Cohort':<20}{'Readout':<10}{'diff':>10}{'marginal p':>13}")
    for cohort_targets, readout in tests:
        label = '+'.join(cohort_targets)
        flag = '  *' if marginal_p[(cohort_targets, readout)] < 0.05 else ''
        print(f"{label:<20}{readout:<10}{observed_diffs[(cohort_targets, readout)]:>+10.4f}"
              f"{marginal_p[(cohort_targets, readout)]:>13.4f}{flag}")
    print(f"\nObserved: {observed_hit_count} of {len(tests)} tests at marginal one-sided p<0.05")
    print(f"Null hit count across {n_permutations} joint permutations: mean={null_hit_counts.mean():.2f}, "
          f"max={null_hit_counts.max()}")
    print(f"Family-wise empirical p (>= {observed_hit_count} hits, this correlation structure): {emp_p_family:.4f}")

    return {'tests': tests, 'observed_diffs': observed_diffs, 'null_diffs': null_diffs,
            'marginal_p': marginal_p, 'observed_hit_count': observed_hit_count,
            'null_hit_counts': null_hit_counts, 'empirical_p_family': emp_p_family}


def verify_dependency_convention(dependency_df, essential_genes=('PCNA', 'RPL3', 'POLR2A'),
                                  nonessential_genes=('OR4A16', 'TAS2R38', 'KRTAP4-1')):
    """
    Directly confirms the dependency matrix's scale/sign convention (no
    cohort split needed): essential genes' mean scores should be clearly
    higher (this matrix: 0-1, higher = more essential) than non-essential
    genes'. Run this before trusting the direction of any dependency-layer
    result -- an incorrect but internally consistent sign assumption does
    not announce itself.
    """
    def _panel_means(genes):
        found = [g for g in genes if g in dependency_df.index]
        missing = [g for g in genes if g not in dependency_df.index]
        if missing:
            print(f"  (not found, skipped: {missing})")
        return dependency_df.loc[found].mean(axis=1) if found else pd.Series(dtype=float)

    print("Essential-gene panel mean dependency (expect strongly positive, near 1):")
    ess = _panel_means(essential_genes)
    for g, v in ess.items():
        print(f"  {g}: {v:.3f}")
    print("Non-essential-gene panel mean dependency (expect near zero):")
    non = _panel_means(nonessential_genes)
    for g, v in non.items():
        print(f"  {g}: {v:.3f}")

    if len(ess) == 0 or len(non) == 0:
        print("\nCould not complete the check -- too few of either panel's genes found.")
        return None

    print(f"\nEssential mean: {ess.mean():.3f}   Non-essential mean: {non.mean():.3f}")
    convention_ok = ess.mean() > non.mean() + 0.2
    if convention_ok:
        print("Convention confirmed: essential genes score clearly higher than non-essential genes.")
    else:
        print("WARNING: not clearly separated in the expected direction -- treat other results with caution.")
    return convention_ok


# ===========================================================================
# SECTION 6: At-scale validation against real BioGRID genetic-interaction data
# ===========================================================================

SL_CATEGORIES = {'Synthetic Lethality', 'Synthetic Growth Defect', 'Negative Genetic'}


def load_biogrid_sl_pairs(path, categories=SL_CATEGORIES, sep='\t'):
    """Loads a raw BioGRID *.tab3.txt file, filters to human-human genetic
    interactions in the given categories. Positional column access
    (confirmed against real rows), auto-skips a stray non-tab-delimited
    line before the header if present."""
    with open(path, 'r', encoding='utf-8', errors='replace') as f:
        first_line = f.readline()
    skiprows = 0
    if first_line.count('\t') < 10:
        print(f"First line doesn't look like a real BioGRID header ({first_line.count(chr(9))} tabs, "
              f"expected 36) -- treating as a stray line and skipping: {first_line.strip()!r}")
        skiprows = 1

    df = pd.read_csv(path, sep=sep, low_memory=False, dtype=str, skiprows=skiprows)
    print(f"Loaded {len(df)} rows, {len(df.columns)} columns.")

    expected = {7: 'Official Symbol Interactor A', 8: 'Official Symbol Interactor B',
                11: 'Experimental System', 12: 'Experimental System Type',
                15: 'Organism ID Interactor A', 16: 'Organism ID Interactor B', 17: 'Throughput'}
    print("Confirming column positions against this file's actual header row:")
    for idx, exp_name in expected.items():
        actual = df.columns[idx] if idx < len(df.columns) else '<column missing>'
        flag = '' if exp_name.split()[0].lower() in str(actual).lower() else '  <-- CHECK THIS'
        print(f"  Column {idx}: expected ~{exp_name!r}, actual header {actual!r}{flag}")

    sym_a, sym_b = df.iloc[:, 7], df.iloc[:, 8]
    exp_system = df.iloc[:, 11]
    exp_type = df.iloc[:, 12].str.strip().str.lower()
    org_a = df.iloc[:, 15].astype(str).str.strip()
    org_b = df.iloc[:, 16].astype(str).str.strip()
    throughput = df.iloc[:, 17]

    mask = (exp_type == 'genetic') & (org_a == '9606') & (org_b == '9606') & (exp_system.isin(categories))
    print(f"\n{mask.sum()} of {len(df)} rows match: genetic, human-human, category in {sorted(categories)}")

    filtered = pd.DataFrame({'gene_a': sym_a[mask].values, 'gene_b': sym_b[mask].values,
                              'category': exp_system[mask].values, 'throughput': throughput[mask].values})
    filtered = filtered[filtered['gene_a'] != filtered['gene_b']]
    n_before = len(filtered)
    filtered = filtered.drop_duplicates(subset=['gene_a', 'gene_b', 'category'])
    print(f"{len(filtered)} unique rows after dedup ({n_before - len(filtered)} duplicates removed)")
    print(f"{filtered[['gene_a', 'gene_b']].drop_duplicates().shape[0]} unique pairs, "
          f"{pd.concat([filtered['gene_a'], filtered['gene_b']]).nunique()} unique genes")
    return filtered


def compute_target_null(target_gene, abundance_df, dependency_df, n=None, mode='max',
                         n_permutations=200, seed=42):
    """One Abundance-based permutation null (full genome-wide diff vector),
    computed once per target and reused for every partner/random control."""
    if target_gene not in abundance_df.index:
        return None
    if n is None:
        n = compute_cohort_n(abundance_df.shape[1])

    class_df = get_classes_by_joint_value([target_gene], abundance_df, n, mode=mode)
    n_low = (class_df['class'] == 'low').sum()
    n_median = (class_df['class'] == 'median').sum()
    observed = get_differentials(class_df, dependency_df, n_low, n_median)['diff']

    rng = np.random.default_rng(seed)
    all_cls = list(abundance_df.columns)
    dep_vals = dependency_df.reindex(columns=all_cls).to_numpy()
    gene_index = dependency_df.index

    null_diffs = np.empty((n_permutations, len(gene_index)))
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', message='Mean of empty slice')
        for i in range(n_permutations):
            shuffled = rng.permutation(len(all_cls))
            low_idx, median_idx = shuffled[:n], shuffled[n:2 * n]
            null_diffs[i] = (np.nanmean(dep_vals[:, low_idx], axis=1)
                              - np.nanmean(dep_vals[:, median_idx], axis=1))

    return {'observed': observed, 'null_diffs': null_diffs, 'gene_index': gene_index,
            'n_low': n_low, 'n_median': n_median}


def empirical_p_for_gene(target_null, readout_gene):
    """One-sided empirical p (higher dependency = predicted direction).
    NaN permutation entries excluded from both numerator and denominator."""
    gene_index = target_null['gene_index']
    if readout_gene not in gene_index:
        return None
    pos = gene_index.get_loc(readout_gene)
    if not np.isscalar(pos):
        return None
    observed = target_null['observed'].loc[readout_gene]
    if pd.isna(observed):
        return None
    null_col = target_null['null_diffs'][:, pos]
    valid = ~np.isnan(null_col)
    if valid.sum() < 20:
        return None
    return float(np.mean(null_col[valid] >= observed)), float(observed)


def run_at_scale_validation(pairs_df, abundance_df, dependency_df, n_random_controls=1,
                             n_permutations=200, seed=42, max_targets=None):
    """Every documented pair, both directions, against a matched
    random-control set per target, from the same target-specific null."""
    genes_as_target = pd.concat([pairs_df['gene_a'], pairs_df['gene_b']]).unique()
    rng = np.random.default_rng(seed)
    if max_targets is not None and len(genes_as_target) > max_targets:
        genes_as_target = rng.choice(genes_as_target, size=max_targets, replace=False)
        print(f"Subsampled to {max_targets} of the full target gene set.")

    partners_of = {}
    tests = []
    for _, row in pairs_df.iterrows():
        a, b = row['gene_a'], row['gene_b']
        partners_of.setdefault(a, set()).add(b)
        partners_of.setdefault(b, set()).add(a)
    for _, row in pairs_df.iterrows():
        tests.append((row['gene_a'], row['gene_b']))
        tests.append((row['gene_b'], row['gene_a']))
    tests = list(set(tests))

    target_cache, results, n_targets_done = {}, [], 0
    unique_targets_needed = set(t for t, r in tests if t in genes_as_target)
    print(f"\nComputing permutation nulls for {len(unique_targets_needed)} unique target genes "
          f"({n_permutations} permutations each)...")

    for target, readout in tests:
        if target not in genes_as_target:
            continue
        if target not in target_cache:
            target_cache[target] = compute_target_null(target, abundance_df, dependency_df,
                                                         n_permutations=n_permutations, seed=seed)
            n_targets_done += 1
            if n_targets_done % 25 == 0:
                print(f"  ...{n_targets_done}/{len(unique_targets_needed)} target nulls computed")
        tnull = target_cache[target]
        if tnull is None:
            continue
        real_result = empirical_p_for_gene(tnull, readout)
        if real_result is not None:
            p, diff = real_result
            results.append({'target': target, 'readout': readout, 'is_real': True, 'empirical_p': p, 'diff': diff})

        candidates = [g for g in tnull['gene_index'] if g != target and g not in partners_of.get(target, set())]
        if len(candidates) >= n_random_controls:
            controls = rng.choice(candidates, size=n_random_controls, replace=False)
            for c in controls:
                ctrl_result = empirical_p_for_gene(tnull, c)
                if ctrl_result is not None:
                    p, diff = ctrl_result
                    results.append({'target': target, 'readout': c, 'is_real': False, 'empirical_p': p, 'diff': diff})

    return pd.DataFrame(results)


def attach_category_and_throughput(results_df, pairs_df):
    """Re-attaches each real pair's category/throughput (not carried in
    the results CSV itself), matching either test direction."""
    fwd = pairs_df.rename(columns={'gene_a': 'target', 'gene_b': 'readout'})
    rev = pairs_df.rename(columns={'gene_b': 'target', 'gene_a': 'readout'})
    lookup = pd.concat([fwd, rev], ignore_index=True).drop_duplicates(subset=['target', 'readout'])

    real = results_df[results_df['is_real']].merge(lookup, on=['target', 'readout'], how='left')
    random_ = results_df[~results_df['is_real']].copy()
    random_['category'] = 'random control'
    random_['throughput'] = 'n/a'

    n_unmatched = real['category'].isna().sum()
    if n_unmatched:
        print(f"Note: {n_unmatched} of {len(real)} real-pair rows couldn't be matched to a category -- excluded.")
        real = real.dropna(subset=['category'])
    return pd.concat([real, random_], ignore_index=True)


def breakdown_by(combined_df, column, thresholds=(0.05, 0.01)):
    print(f"\n{'=' * 70}\nHit rate by {column}\n{'=' * 70}")
    for value, group in combined_df.groupby(column):
        row = f"{value:<25} n={len(group):<8}"
        for t in thresholds:
            row += f"  p<{t}: {(group['empirical_p'] < t).mean():.3f}"
        print(row)


def test_slice_vs_random(combined_df, slice_mask, label, thresholds=(0.05, 0.01)):
    """One-sided Fisher's exact test (slice rate > random rate)."""
    random_df = combined_df[combined_df['category'] == 'random control']
    slice_df = combined_df[slice_mask]
    print(f"\n{label} (n={len(slice_df)}) vs random control (n={len(random_df)}):")
    if len(slice_df) == 0:
        print("  (empty slice -- nothing to test)")
        return
    for t in thresholds:
        slice_hits = int((slice_df['empirical_p'] < t).sum())
        random_hits = int((random_df['empirical_p'] < t).sum())
        table = [[slice_hits, len(slice_df) - slice_hits], [random_hits, len(random_df) - random_hits]]
        odds_ratio, p_value = fisher_exact(table, alternative='greater')
        slice_rate = slice_hits / len(slice_df)
        random_rate = random_hits / len(random_df) if len(random_df) else float('nan')
        print(f"  p<{t}: slice rate={slice_rate:.4f} ({slice_hits}/{len(slice_df)}), "
              f"random rate={random_rate:.4f} ({random_hits}/{len(random_df)}), "
              f"odds ratio={odds_ratio:.2f}, Fisher's exact p={p_value:.4f}")


def plot_category_throughput_breakdown(combined, output_path=FIG_DIR / 'fig7_biogrid_breakdown.png', thresholds=(0.05, 0.01)):
    random_rate = {t: (combined.loc[combined['category'] == 'random control', 'empirical_p'] < t).mean()
                   for t in thresholds}
    rows = []
    for (cat, tp), group in combined.groupby(['category', 'throughput']):
        if cat == 'random control' or len(group) < 30:
            continue
        rows.append({'label': f"{cat}\n({tp})", 'n': len(group),
                     **{f'p<{t}': (group['empirical_p'] < t).mean() for t in thresholds}})
    if not rows:
        print(f"No category/throughput slice has >=30 pairs -- skipping {output_path}.")
        return
    df = pd.DataFrame(rows).sort_values(f'p<{thresholds[-1]}', ascending=False)

    fig, ax = plt.subplots(figsize=(8, 5))
    x = np.arange(len(df))
    width = 0.35
    for i, t in enumerate(thresholds):
        offset = (i - (len(thresholds) - 1) / 2) * width
        ax.bar(x + offset, df[f'p<{t}'], width, label=f'empirical p<{t}')
    for t in thresholds:
        ax.axhline(random_rate[t], linestyle=':', linewidth=1, color='grey')
    ax.set_xticks(x)
    ax.set_xticklabels([f"{row['label']}\nn={row['n']}" for _, row in df.iterrows()], fontsize=8)
    ax.set_ylabel('Hit rate (fraction of pairs)')
    ax.set_title('BioGRID pair hit rate by evidence category and throughput\n(dotted lines = random-control baseline)')
    ax.legend(fontsize=8.5)
    plt.tight_layout()
    plt.savefig(output_path, dpi=FIG_DPI, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved {output_path}")


def plot_enrichment_curve(combined, output_path=FIG_DIR / 'fig8_biogrid_enrichment.png'):
    """
    PRIMARY summary figure: fold-enrichment (real hit rate / random hit
    rate) vs significance threshold. Built because whole-curve ROC-AUC
    understates this effect -- it integrates across the entire score
    range, but the real signal here is concentrated specifically in the
    significant tail (confirmed: Synthetic Lethality/Low Throughput gave
    AUC=0.571 despite Fisher's exact p=0.0020 at p<0.01 on the same data).
    This curve shows that tail concentration directly instead.
    """
    thresholds = np.logspace(np.log10(0.005), 0, 60)
    random_df = combined[combined['category'] == 'random control']
    slices = [
        ('Synthetic Lethality, Low Throughput',
         (combined['category'] == 'Synthetic Lethality') & (combined['throughput'] == 'Low Throughput'), '#d95f02'),
        ('All Low Throughput', combined['throughput'] == 'Low Throughput', '#1b9e77'),
        ('Negative Genetic, High Throughput',
         (combined['category'] == 'Negative Genetic') & (combined['throughput'] == 'High Throughput'), '#7570b3'),
    ]
    fig, ax = plt.subplots(figsize=(7, 5))
    for label, mask, color in slices:
        slice_df = combined[mask]
        n = len(slice_df)
        fold = []
        for t in thresholds:
            slice_rate = (slice_df['empirical_p'] < t).mean()
            random_rate = (random_df['empirical_p'] < t).mean()
            fold.append(slice_rate / random_rate if random_rate > 0 else np.nan)
        ax.plot(thresholds, fold, label=f'{label} (n={n})', color=color, linewidth=1.8)

    ax.axhline(1.0, color='grey', linestyle=':', linewidth=1, label='No enrichment (fold=1)')
    ax.axvline(0.05, color='lightgrey', linestyle='--', linewidth=0.8)
    ax.axvline(0.01, color='lightgrey', linestyle='--', linewidth=0.8)
    ax.set_xscale('log')
    ax.invert_xaxis()
    ax.set_xlabel('Empirical p-value threshold (stricter \u2192)')
    ax.set_ylabel('Fold enrichment over random-control hit rate')
    ax.set_title('Real BioGRID pairs are enriched specifically in the significant tail,\nnot uniformly across the score range')
    ax.legend(fontsize=8, loc='upper left')
    plt.tight_layout()
    plt.savefig(output_path, dpi=FIG_DPI, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved {output_path}")


def plot_biogrid_method_comparison(summary, output_path=None):
    """Figure 9: fold enrichment over the random control for each method, by BioGRID slice, at p<0.05 and p<0.01.
    A star marks a one-sided Fisher p<0.05 against the random control. The dotted line (fold=1) is no enrichment."""
    output_path = output_path or (FIG_DIR / 'fig9_biogrid_method_comparison.png')
    methods = list(METHOD_LABELS.values())
    colors = ['#d95f02', '#1b9e77', '#7570b3', '#e7298a']
    slices = [sname for sname in summary['slice'].unique() if sname != 'All real pairs']
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8), sharey=True)
    width = 0.2
    for ax, t in zip(axes, (0.05, 0.01)):
        for i, (m, c) in enumerate(zip(methods, colors)):
            vals, stars = [], []
            for sname in slices:
                row = summary[(summary['slice'] == sname) & (summary['method'] == m) & (summary['threshold'] == t)].iloc[0]
                vals.append(row['fold']); stars.append(row['fisher_p'] < 0.05)
            xs = np.arange(len(slices)) + (i - 1.5) * width
            ax.bar(xs, vals, width, color=c, label=m)
            for x_, v_, st in zip(xs, vals, stars):
                if st:
                    ax.text(x_, v_, '*', ha='center', va='bottom', fontsize=13)
        ax.axhline(1.0, color='grey', linestyle=':', linewidth=1)
        ax.set_xticks(np.arange(len(slices)))
        ax.set_xticklabels([sn.replace(', ', ',\n').replace(' [primary]', '\n[primary]') for sn in slices], fontsize=8)
        ax.set_title('Hit threshold: empirical p<%g' % t, fontsize=10)
    axes[0].set_ylabel('Fold enrichment over the random control')
    axes[0].legend(fontsize=8, loc='upper right')
    fig.suptitle('Which method best finds known synthetic lethal pairs?\n(* = Fisher p<0.05 against the random control)', fontsize=11)
    plt.tight_layout()
    plt.savefig(output_path, dpi=FIG_DPI, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved {output_path}")


def compute_roc(y_true, y_score):
    """Dependency-free ROC/AUC (verified to match sklearn exactly during
    development). Secondary figure only -- plot_enrichment_curve above is
    the primary summary for this analysis."""
    y_true, y_score = np.asarray(y_true), np.asarray(y_score)
    order = np.argsort(-y_score)
    y_true_sorted = y_true[order]
    n_pos, n_neg = y_true.sum(), len(y_true) - y_true.sum()
    if n_pos == 0 or n_neg == 0:
        raise ValueError("Need at least one positive and one negative example.")
    tps = np.cumsum(y_true_sorted)
    fps = np.cumsum(1 - y_true_sorted)
    tpr = np.concatenate([[0], tps / n_pos, [1]])
    fpr = np.concatenate([[0], fps / n_neg, [1]])
    _trapezoid = getattr(np, 'trapezoid', None) or np.trapz
    return fpr, tpr, _trapezoid(tpr, fpr)


def plot_roc_curves_secondary(combined, output_path=FIG_DIR / 'fig8b_biogrid_roc_secondary.png'):
    fig, ax = plt.subplots(figsize=(5.5, 5.5))
    random_mask = combined['category'] == 'random control'
    slices = [
        ('All real pairs', combined['category'].isin(['Synthetic Lethality', 'Negative Genetic', 'Synthetic Growth Defect']), '#1b9e77'),
        ('Synthetic Lethality, Low Throughput', (combined['category'] == 'Synthetic Lethality') & (combined['throughput'] == 'Low Throughput'), '#d95f02'),
        ('Negative Genetic, High Throughput', (combined['category'] == 'Negative Genetic') & (combined['throughput'] == 'High Throughput'), '#7570b3'),
    ]
    aucs = {}
    for label, mask, color in slices:
        subset = combined[mask | random_mask]
        y_true = (~subset['category'].eq('random control')).astype(int).values
        y_score = 1 - subset['empirical_p'].values
        fpr, tpr, roc_auc = compute_roc(y_true, y_score)
        ax.plot(fpr, tpr, label=f'{label} (AUC={roc_auc:.3f})', color=color)
        aucs[label] = roc_auc
    ax.plot([0, 1], [0, 1], 'k--', linewidth=0.8, label='Chance (AUC=0.500)')
    ax.set_xlabel('False positive rate')
    ax.set_ylabel('True positive rate')
    ax.set_title('Secondary figure: whole-curve ROC-AUC\n(understates a tail-concentrated effect -- see enrichment curve)')
    ax.legend(loc='lower right', fontsize=7.5)
    plt.tight_layout()
    plt.savefig(output_path, dpi=FIG_DPI, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved {output_path}")
    return aucs



# ===========================================================================
# RUN PIPELINE -- edit the config in each cell below, then run cells with
# Shift+Enter (VS Code Interactive Window), or run this whole file as a
# script. Cells are independent; run only the ones you need.
# ===========================================================================

# %% [CELL 1 of 17] Load the four SimKO matrices
# #####################################################################
# CELL 1 of 17: Load the four SimKO matrices
#   FIGURE OUTPUT: (no figure) -- needed by every later cell
#   DATA / DEPENDS ON: app data/ folder
# #####################################################################
abundance, expression, mutation, dependency = load_matrices()

# %% [CELL 2 of 17] Phase 0: permutation FDR, tissue comparison, tissue bias
# #####################################################################
# CELL 2 of 17: Phase 0: permutation FDR, tissue comparison, tissue bias
#   FIGURE OUTPUT: FIGURES 1, 2, 3 (drawn in Cell 17)
#   DATA / DEPENDS ON: none external
# #####################################################################
# --- Phase 0: statistical foundation (no external data needed) ---
TARGETS = ['PCNA', 'ARID1B', 'WRN', 'SMARCA4']
targets_covered = [t for t in TARGETS if t in abundance.index]

comparison_table, phase0_results = run_phase0_tissue_comparison(
    abundance, dependency, targets_covered, n_permutations=200, seed=42, n_panel_targets=10, n_tissues=2,
)
print(comparison_table.to_string())
comparison_table.to_csv(PHASE0_TABLE_PATH)

# (Figures 1-3 are drawn from these results in the final 'ALL MANUSCRIPT FIGURES' cell.)
tissue_bias_targets = targets_covered + phase0_results['pancancer']['raw']['neg_control']['targets']
tissue_bias_df = quantify_tissue_bias(tissue_bias_targets, abundance, n_permutations=200, seed=42)
target_categories = {t: 'Real targets' for t in targets_covered}
target_categories.update({t: 'Essentiality-curated panel' for t in phase0_results['pancancer']['raw']['neg_control']['targets']})

# %% [CELL 3 of 17] ARID1A vs real RPE1 knockout proteomics (Abundance layer)
# #####################################################################
# CELL 3 of 17: ARID1A vs real RPE1 knockout proteomics (Abundance layer)
#   FIGURE OUTPUT: FIGURE 4 panel A and FIGURE 5 left column (Abundance vs real protein; drawn in Cell 17)
#   DATA / DEPENDS ON: ARID1A_RPE_PATH
# #####################################################################
# --- Concordance: ARID1A, real RPE1 KO proteomics (Abundance layer) ---
# Real header: Accession, Description, Gene name, R1_WT, R1_KO27, R1_KO32,
# R2_WT, R2_KO27, R2_KO32 -- two replicates, each with one WT and two
# independent KO clones. Pooling both clones as "KO" for the primary test.
# (ARID1A_RPE_PATH is set in the DATA FILES block at the top of the script.)
WT_COLS = [3, 6]              # R1_WT, R2_WT
KO_COLS = [4, 5, 7, 8]        # both clones, both reps, pooled

real_diff_arid1a = load_real_differential_from_raw_intensity_table(
    ARID1A_RPE_PATH, group_a_cols=WT_COLS, group_b_cols=KO_COLS, gene_col=2
)
simko_diff_arid1a = get_simko_differential(['ARID1A'], abundance)
result_arid1a = run_concordance_analysis(simko_diff_arid1a, real_diff_arid1a, top_n=200,
                                          output_prefix='concordance_ARID1A_abundance')
n_pc = compute_cohort_n(abundance.shape[1])
perm_arid1a = permutation_null_for_concordance(['ARID1A'], abundance, real_diff_arid1a, n=n_pc,
                                                n_permutations=200, seed=42)

# %% [CELL 4 of 17] Robustness: do the two RPE1 KO clones agree?
# #####################################################################
# CELL 4 of 17: Robustness: do the two RPE1 KO clones agree?
#   FIGURE OUTPUT: (no figure; supplementary scatter files only)
#   DATA / DEPENDS ON: ARID1A_RPE_PATH
# #####################################################################
# --- Optional robustness check: do the two independent KO clones agree? ---
# Two clones both showing concordance with SimKO is stronger evidence than
# the pooled result alone; if they disagree, that's worth knowing too.
real_diff_arid1a_clone27 = load_real_differential_from_raw_intensity_table(
    ARID1A_RPE_PATH, group_a_cols=WT_COLS, group_b_cols=[4, 7], gene_col=2  # R1_KO27, R2_KO27
)
real_diff_arid1a_clone32 = load_real_differential_from_raw_intensity_table(
    ARID1A_RPE_PATH, group_a_cols=WT_COLS, group_b_cols=[5, 8], gene_col=2  # R1_KO32, R2_KO32
)
result_arid1a_clone27 = run_concordance_analysis(simko_diff_arid1a, real_diff_arid1a_clone27, top_n=200,
                                                  output_prefix='concordance_ARID1A_clone27')
result_arid1a_clone32 = run_concordance_analysis(simko_diff_arid1a, real_diff_arid1a_clone32, top_n=200,
                                                  output_prefix='concordance_ARID1A_clone32')

# %% [CELL 5 of 17] Benchmark: cohort method vs plain correlation, RPE1 protein data
# #####################################################################
# CELL 5 of 17: Benchmark: cohort method vs plain correlation, RPE1 protein data
#   FIGURE OUTPUT: FIGURE 4, panel A (drawn in Cell 17)
#   DATA / DEPENDS ON: needs Cell 3
# #####################################################################
# --- Benchmark: SimKO's cohort method vs. a plain correlation (ARID1A RPE data) ---
corr_diff_arid1a = get_correlation_differential('ARID1A', abundance)
result_corr_arid1a = run_concordance_analysis(corr_diff_arid1a, real_diff_arid1a, top_n=200,
                                               output_prefix='benchmark_correlation_ARID1A')
perm_corr_arid1a = permutation_null_for_correlation_concordance('ARID1A', abundance, real_diff_arid1a,
                                                                  n_permutations=200, seed=42)
combined_arid1a = test_combined_strategy(perm_arid1a, perm_corr_arid1a)

# %% [CELL 6 of 17] ARID1A vs real HCT116 knockout RNA-seq (Expression layer)
# #####################################################################
# CELL 6 of 17: ARID1A vs real HCT116 knockout RNA-seq (Expression layer)
#   FIGURE OUTPUT: (no figure directly; Cell 7 reuses this to make Figure 4 panel B)
#   DATA / DEPENDS ON: HCT116_FPKM_PATH
# #####################################################################
# --- Concordance: ARID1A, real HCT116 RNA-seq (Expression layer) ---
# Independent ground truth, own-layer test: isogenic HCT116 cells, CRISPR
# knockout of ARID1A (GSE101966, HOMER analyzeRepeats FPKM output). Needs
# the gene-symbol-converted cache; if it doesn't exist yet, this runs a
# one-time RefSeq->symbol conversion (needs `pip install mygene` and
# outbound internet access -- run that conversion step locally first if
# this environment is network-restricted, then point HCT116_FPKM_PATH at
# the raw file as usual and the cached _bysymbol.txt will be picked up).
# (HCT116_FPKM_PATH is set in the DATA FILES block at the top of the script.)

hct116_symbol_path = get_gene_symbol_fpkm_path(HCT116_FPKM_PATH)
real_diff_arid1a_transcript = load_real_differential_from_homer_fpkm(
    hct116_symbol_path, gene_col_index=0, group_a_condition='WT', group_b_condition='ARID1AKO',
)
simko_diff_arid1a_expression = get_simko_differential(['ARID1A'], expression)
result_arid1a_expression = run_concordance_analysis(simko_diff_arid1a_expression, real_diff_arid1a_transcript,
                                                      top_n=200, output_prefix='concordance_ARID1A_expression')
n_pc_expr = compute_cohort_n(expression.shape[1])
perm_arid1a_expression = permutation_null_for_concordance(['ARID1A'], expression, real_diff_arid1a_transcript,
                                                            n=n_pc_expr, n_permutations=200, seed=42)

# %% [CELL 7 of 17] Benchmark: cohort method vs plain correlation, HCT116 RNA-seq (Expression AND Abundance layers)
# #####################################################################
# CELL 7 of 17: Benchmark: cohort method vs plain correlation, HCT116 RNA-seq (Expression AND Abundance layers)
#   FIGURE OUTPUT: FIGURE 4, panel B (Abundance layer vs HCT116 RNA); the Expression-layer comparison is optional (drawn in Cell 17)
#   DATA / DEPENDS ON: needs Cells 3 and 6
# #####################################################################
# --- Benchmark: SimKO cohort method vs plain correlation, HCT116 RNA-seq (Expression layer) ---
# This is the comparison Figure 4 shows. Both methods get the SAME number of
# permutations (the combined best-of-both test sums their null arrays
# element-wise, so the lengths must match). A permutation p-value cannot be
# smaller than 1/N_PERM_BENCH, so raise this if you want to resolve smaller
# p-values than the floor printed on the figure. Needs the HCT116 cell above.
N_PERM_BENCH = 1000

perm_simko_expr = permutation_null_for_concordance(['ARID1A'], expression, real_diff_arid1a_transcript,
                                                     n=n_pc_expr, n_permutations=N_PERM_BENCH, seed=42)
corr_diff_arid1a_expression = get_correlation_differential('ARID1A', expression)
result_corr_arid1a_expression = run_concordance_analysis(corr_diff_arid1a_expression, real_diff_arid1a_transcript,
                                                           top_n=200, output_prefix='benchmark_correlation_ARID1A_expression')
perm_corr_expr = permutation_null_for_correlation_concordance('ARID1A', expression, real_diff_arid1a_transcript,
                                                               n_permutations=N_PERM_BENCH, seed=42)
combined_expr = test_combined_strategy(perm_simko_expr, perm_corr_expr)

# --- Same benchmark with SimKO's DEFAULT layer (Abundance) against the same HCT116 RNA-seq ---
# SimKO side: cohorts formed on the target's PROTEIN abundance (exactly the Cell 8 cross-check, but with
# N_PERM_BENCH permutations so the two methods can be combined). Correlation side: ARID1A protein abundance
# correlated, across cell lines, with every gene's RNA expression -- restricted to genes that also exist in the
# Abundance matrix so both methods rank the same gene universe. Needs Cells 3 and 6.
perm_simko_abund_rna = permutation_null_for_concordance(['ARID1A'], abundance, real_diff_arid1a_transcript,
                                                          n=n_pc, n_permutations=N_PERM_BENCH, seed=42)
_expr_same_genes = expression.loc[expression.index.isin(abundance.index) | (expression.index == 'ARID1A')]
corr_diff_abund_rna = get_correlation_differential('ARID1A', _expr_same_genes, target_matrix=abundance)
result_corr_abund_rna = run_concordance_analysis(corr_diff_abund_rna, real_diff_arid1a_transcript, top_n=200,
                                                   output_prefix='benchmark_correlation_ARID1A_abundance_vs_rna')
perm_corr_abund_rna = permutation_null_for_correlation_concordance('ARID1A', _expr_same_genes, real_diff_arid1a_transcript,
                                                                    n_permutations=N_PERM_BENCH, seed=42,
                                                                    target_matrix=abundance)
combined_abund_rna = test_combined_strategy(perm_simko_abund_rna, perm_corr_abund_rna)

# %% [CELL 7b of 17 -- optional] Quick look: SimKO vs plain correlation on the strongest genes
# #####################################################################
# CELL 7b (optional, exploratory): SimKO vs plain correlation on the top third / top 10% of genes
#   FIGURE OUTPUT: none (prints a table only; descriptive, no permutation test)
#   DATA / DEPENDS ON: needs Cells 3, 5 and 7
# #####################################################################
# The overall rho is dragged down by the ~6,000 genes that barely change. This asks the same question on subsets.
#   REAL = the genes with the largest |change| in the REAL experiment (does the method rank the genes that really moved?)
#   PRED = the genes with the largest |score| from the METHOD (are its confident calls right?)
# Per subset: Spearman rho against the real change, and % of genes whose direction (up/down) is right (50% = chance).
from scipy.stats import spearmanr

def _top_frac_table(simko_diff, corr_diff, real_diff, label, fracs=(1.0, 1/3, 0.10)):
    rows = []
    # same gene universe for both methods, so the n_genes and the subsets match exactly
    common = simko_diff.index.intersection(corr_diff.index).intersection(real_diff.index)
    for name, d in (('SimKO cohort', simko_diff), ('Plain correlation', corr_diff)):
        m = d.loc[d.index.isin(common), ['diff']].join(real_diff[['logfc']], how='inner').dropna()
        for pick in ('REAL', 'PRED'):
            key = m['logfc'].abs() if pick == 'REAL' else m['diff'].abs()
            for f in fracs:
                if f == 1.0 and pick == 'PRED':
                    continue                       # all genes is the same for both picks
                k = int(round(len(m) * f))
                sub = m.loc[key.nlargest(k).index]
                rho = spearmanr(sub['diff'], sub['logfc'])[0]
                agree = float((np.sign(sub['diff']) == np.sign(sub['logfc'])).mean())
                rows.append({'comparison': label, 'method': name,
                             'subset': 'all genes' if f == 1.0 else f'top {f:.0%} by {pick}',
                             'n_genes': k, 'rho': round(rho, 3), 'direction_agree': round(agree, 3)})
    return pd.DataFrame(rows)

_a = _top_frac_table(simko_diff_arid1a, corr_diff_arid1a, real_diff_arid1a,
                     'A: Abundance vs RPE1 protein')
_b = _top_frac_table(simko_diff_arid1a, corr_diff_abund_rna, real_diff_arid1a_transcript,
                     'B: Abundance vs HCT116 RNA')
_quick = pd.concat([_a, _b], ignore_index=True)
print(_quick.pivot_table(index=['comparison', 'subset', 'n_genes'], columns='method',
                         values=['rho', 'direction_agree']).to_string())

# %% [CELL 7c of 17 -- optional] Quick look: does a parameter-free SimKO + correlation combination help?
# #####################################################################
# CELL 7c (optional, exploratory): SimKO + plain correlation combined into one score
#   FIGURE OUTPUT: none (prints a table only; descriptive, no permutation test)
#   DATA / DEPENDS ON: needs Cells 3, 5 and 7
# #####################################################################
# ONE combination, fixed before looking at the results, with NO tuned parameters:
#   combined score = the average of the two methods' rank-scaled scores (each rescaled to -1..+1 by rank).
# Both scores already share one sign convention (positive = predicted to RISE after knockout).
# Do not add weights, thresholds or other variants and re-run until it looks good: that tunes the method on the very
# data it is then scored on. If this helps, the honest next step is a joint permutation null and a SECOND target.
def _rank_scale(s):
    return (s.rank() - 0.5) / len(s) * 2 - 1

def _combine_scores(simko_diff, corr_diff):
    common = simko_diff.index.intersection(corr_diff.index)
    comb = (_rank_scale(simko_diff.loc[common, 'diff']) + _rank_scale(corr_diff.loc[common, 'diff'])) / 2
    return pd.DataFrame({'diff': comb})

def _three_way_table(simko_diff, corr_diff, real_diff, label, fracs=(1.0, 1/3, 0.10), top_n=200):
    common = simko_diff.index.intersection(corr_diff.index).intersection(real_diff.index)
    s, c = simko_diff.loc[common, ['diff']], corr_diff.loc[common, ['diff']]
    methods = {'SimKO cohort': s, 'Plain correlation': c, 'Combined': _combine_scores(s, c)}
    real = real_diff.loc[common, 'logfc'].dropna()
    rows = []
    for name, d in methods.items():
        m = d[['diff']].join(real.rename('logfc'), how='inner').dropna()
        for pick in ('REAL', 'PRED'):
            key = m['logfc'].abs() if pick == 'REAL' else m['diff'].abs()
            for f in fracs:
                if f == 1.0 and pick == 'PRED':
                    continue
                k = int(round(len(m) * f)); sub = m.loc[key.nlargest(k).index]
                rows.append({'comparison': label, 'method': name,
                             'subset': 'all genes' if f == 1.0 else f'top {f:.0%} by {pick}', 'n_genes': k,
                             'rho': round(spearmanr(sub['diff'], sub['logfc'])[0], 3),
                             'direction_agree': round(float((np.sign(sub['diff']) == np.sign(sub['logfc'])).mean()), 3)})
        up_real = set(m['logfc'].nlargest(top_n).index); dn_real = set(m['logfc'].nsmallest(top_n).index)
        rows.append({'comparison': label, 'method': name, 'subset': f'top-{top_n} UP overlap', 'n_genes': top_n,
                     'rho': len(set(m['diff'].nlargest(top_n).index) & up_real), 'direction_agree': np.nan})
        rows.append({'comparison': label, 'method': name, 'subset': f'top-{top_n} DOWN overlap', 'n_genes': top_n,
                     'rho': len(set(m['diff'].nsmallest(top_n).index) & dn_real), 'direction_agree': np.nan})
    return pd.DataFrame(rows)

_t = pd.concat([
    _three_way_table(simko_diff_arid1a, corr_diff_arid1a, real_diff_arid1a, 'A: Abundance vs RPE1 protein'),
    _three_way_table(simko_diff_arid1a, corr_diff_abund_rna, real_diff_arid1a_transcript, 'B: Abundance vs HCT116 RNA')],
    ignore_index=True)
print("(rho column = Spearman rho for the gene subsets; for the 'overlap' rows it is the COUNT of shared genes, ~6 expected by chance)")
print(_t.pivot_table(index=['comparison', 'subset', 'n_genes'], columns='method',
                     values=['rho', 'direction_agree'], sort=False).to_string())

# %% [CELL 7d of 17 -- optional] Direct test: is SimKO's difference from the plain correlation more than chance?
# #####################################################################
# CELL 7d (optional): direct SimKO-vs-correlation difference test (supports, or removes, any "better than" wording)
#   FIGURE OUTPUT: none (prints a table only)
#   DATA / DEPENDS ON: needs Cells 3, 5 and 7 (uses their stored permutation nulls; runs no new permutations)
# #####################################################################
# Until now each method was tested against chance separately. Here the DIFFERENCE is tested. For each metric
# (overall rho, top-200 UP overlap, top-200 DOWN overlap) the observed difference (SimKO minus correlation) is compared
# with the differences you get when neither method has any real signal (every SimKO null value minus every correlation
# null value). Two-sided, so it reports a correlation advantage as readily as a SimKO advantage.
# The pairing between the two methods is ignored, which makes the test conservative. The smallest possible p is
# 1/(number of SimKO nulls x number of correlation nulls).
def _direct_difference_table(label, perm_s, perm_c):
    rows = []
    for metric, k_obs, k_null in (('overall rho', 'observed_rho', 'null_rhos'),
                                  ('top-200 UP overlap', 'observed_up', 'null_ups'),
                                  ('top-200 DOWN overlap', 'observed_down', 'null_downs')):
        obs = perm_s[k_obs] - perm_c[k_obs]
        null = np.subtract.outer(np.asarray(perm_s[k_null], float), np.asarray(perm_c[k_null], float))
        p_ge, p_le = float(np.mean(null >= obs)), float(np.mean(null <= obs))
        rows.append({'comparison': label, 'metric': metric,
                     'SimKO': round(float(perm_s[k_obs]), 3), 'Plain correlation': round(float(perm_c[k_obs]), 3),
                     'difference (SimKO - corr)': round(float(obs), 3), 'null difference mean': round(float(null.mean()), 3),
                     'two-sided p': round(min(1.0, 2 * min(p_ge, p_le)), 4)})
    return pd.DataFrame(rows)

_direct = pd.concat([
    _direct_difference_table('A: Abundance vs RPE1 protein', perm_arid1a, perm_corr_arid1a),
    _direct_difference_table('B: Abundance vs HCT116 RNA', perm_simko_abund_rna, perm_corr_abund_rna)], ignore_index=True)
print(_direct.to_string(index=False))

# %% [CELL 8 of 17] Headline 2x2: Abundance vs Expression, each against both real datasets
# #####################################################################
# CELL 8 of 17: Headline 2x2: Abundance vs Expression, each against both real datasets
#   FIGURE OUTPUT: FIGURE 5 right column (Expression layer vs real protein; drawn in Cell 17)
#   DATA / DEPENDS ON: needs Cells 3 and 6
# #####################################################################
# --- Headline: does Abundance or Expression concord better, tested fairly? ---
# The two results above each use a DIFFERENT ground-truth dataset (RPE1
# protein vs HCT116 transcript) -- not a fair comparison, since any
# difference could just reflect which real dataset is less noisy, not
# which SimKO layer is more predictive. This crosses each SimKO layer's
# output against the OTHER ground truth as well, giving a real 2x2: for
# each ground truth, which layer concords with it better. Reuses exactly
# the same concordance/permutation machinery as every result above --
# only the layer/ground-truth pairing is new.
print(f"\n{'=' * 70}\nCross-check 1: Abundance-layer SimKO output vs HCT116 RNA-seq ground truth\n{'=' * 70}")
result_abund_vs_transcript = run_concordance_analysis(simko_diff_arid1a, real_diff_arid1a_transcript, top_n=200,
                                                        output_prefix='cross_abundance_vs_transcript')
perm_abund_vs_transcript = permutation_null_for_concordance(['ARID1A'], abundance, real_diff_arid1a_transcript,
                                                              n=n_pc, n_permutations=200, seed=42)

print(f"\n{'=' * 70}\nCross-check 2: Expression-layer SimKO output vs RPE1 proteomics ground truth\n{'=' * 70}")
result_expr_vs_protein = run_concordance_analysis(simko_diff_arid1a_expression, real_diff_arid1a, top_n=200,
                                                    output_prefix='cross_expression_vs_protein')
perm_expr_vs_protein = permutation_null_for_concordance(['ARID1A'], expression, real_diff_arid1a,
                                                          n=n_pc_expr, n_permutations=200, seed=42)

print(f"\n{'=' * 70}\n2x2 summary: empirical p for top-200 UP / DOWN overlap\n{'=' * 70}")
print(f"{'SimKO layer':<14}{'vs RPE1 protein (UP/DOWN)':<28}{'vs HCT116 transcript (UP/DOWN)':<32}")
abund_vs_protein = f"{perm_arid1a['empirical_p_up']:.4f} / {perm_arid1a['empirical_p_down']:.4f}"
abund_vs_transcript = f"{perm_abund_vs_transcript['empirical_p_up']:.4f} / {perm_abund_vs_transcript['empirical_p_down']:.4f}"
expr_vs_protein = f"{perm_expr_vs_protein['empirical_p_up']:.4f} / {perm_expr_vs_protein['empirical_p_down']:.4f}"
expr_vs_transcript = f"{perm_arid1a_expression['empirical_p_up']:.4f} / {perm_arid1a_expression['empirical_p_down']:.4f}"
print(f"{'Abundance':<14}{abund_vs_protein:<28}{abund_vs_transcript:<32}")
print(f"{'Expression':<14}{expr_vs_protein:<28}{expr_vs_transcript:<32}")
print("\n(Diagonal cells -- Abundance/RPE1-protein and Expression/HCT116-transcript -- are each layer's")
print(" 'home' ground truth, already reported above. Off-diagonal cells are the new cross-checks: do")
print(" they still concord, or was the original result partly a property of that specific pairing?)")

# %% [CELL 9 of 17] Power check: does a real signal survive a small cohort?
# #####################################################################
# CELL 9 of 17: Power check: does a real signal survive a small cohort?
#   FIGURE OUTPUT: (no figure; quoted in text)
#   DATA / DEPENDS ON: none external
# #####################################################################
# --- Power check: does a known-real signal survive a small cohort? ---
power_large, power_small = run_power_check(expression, target='PCNA', n_small=8, n_large=20, n_permutations=200)

# %% [CELL 10 of 17] Synthetic lethality case study: BRCA1/BRCA2 -> PARP1
# #####################################################################
# CELL 10 of 17: Synthetic lethality case study: BRCA1/BRCA2 -> PARP1
#   FIGURE OUTPUT: FIGURE 6 (drawn in Cell 17)
#   DATA / DEPENDS ON: none external
# #####################################################################
# --- Synthetic lethality case study: BRCA1/BRCA2 -> PARP1 ---
verify_dependency_convention(dependency)

sl_results = {}
for target in ['BRCA1', 'BRCA2']:
    if target in abundance.index and 'PARP1' in dependency.index:
        sl_results[target] = test_single_gene_dependency_shift(target, 'PARP1', abundance, dependency, n_permutations=200)

sl_tissue_results = {}
for tissue in ['BREAST', 'OVARY']:
    tissue_cls = get_tissue_cell_lines(abundance.columns, tissue)
    if len(tissue_cls) < 2:
        print(f"Skipping {tissue} -- too few cell lines.")
        continue
    for target in ['BRCA1', 'BRCA2']:
        if target in abundance.index and 'PARP1' in dependency.index:
            sl_tissue_results[(target, tissue)] = test_single_gene_dependency_shift(
                target, 'PARP1', abundance, dependency, active_cell_lines=tissue_cls, n_permutations=200
            )

# %% [CELL 11 of 17] At-scale test against BioGRID synthetic lethal / genetic interaction pairs
# #####################################################################
# CELL 11 of 17: At-scale test against BioGRID synthetic lethal / genetic interaction pairs
#   FIGURE OUTPUT: FIGURES 7, 8 and supplementary ROC (drawn in Cell 17)
#   DATA / DEPENDS ON: BIOGRID_PATH; caches to BIOGRID_RESULTS_PATH
# #####################################################################
# --- At-scale validation: real BioGRID synthetic lethal / genetic interaction pairs ---
# (BIOGRID_PATH and BIOGRID_RESULTS_PATH are set in the DATA FILES block at the top of the script.)

biogrid_pairs = load_biogrid_sl_pairs(BIOGRID_PATH)

biogrid_results_path = BIOGRID_RESULTS_PATH
if Path(biogrid_results_path).exists():
    print(f"Reusing cached {biogrid_results_path} -- delete it to force a fresh permutation run.")
    biogrid_results = pd.read_csv(biogrid_results_path)
else:
    genes_present = set(abundance.index) & set(dependency.index)
    biogrid_pairs_testable = biogrid_pairs[biogrid_pairs['gene_a'].isin(genes_present) | biogrid_pairs['gene_b'].isin(genes_present)]
    biogrid_results = run_at_scale_validation(biogrid_pairs_testable, abundance, dependency, max_targets=None)
    biogrid_results.to_csv(biogrid_results_path, index=False)

biogrid_combined = attach_category_and_throughput(biogrid_results, biogrid_pairs)
breakdown_by(biogrid_combined, 'category')
breakdown_by(biogrid_combined, 'throughput')

# Pooled real-vs-random rates (quoted in the Results text) and the per-slice rates behind Figure 7's caption
is_real = biogrid_combined['category'] != 'random control'
print("\nPooled hit rate, all real pairs vs random control:")
for t in (0.05, 0.01):
    real_rate = (biogrid_combined.loc[is_real, 'empirical_p'] < t).mean()
    rand_rate = (biogrid_combined.loc[~is_real, 'empirical_p'] < t).mean()
    print("  p<%g: real %.4f (n=%d) vs random %.4f (n=%d)" % (t, real_rate, is_real.sum(), rand_rate, (~is_real).sum()))
print("\nHit rate by category AND throughput (slices with >=30 pairs; compare with random control above):")
for (cat, tp), grp in biogrid_combined[is_real].groupby(['category', 'throughput']):
    if len(grp) >= 30:
        print("  %-24s %-18s n=%-6d p<0.05: %.3f   p<0.01: %.3f" % (
            cat, tp, len(grp), (grp['empirical_p'] < 0.05).mean(), (grp['empirical_p'] < 0.01).mean()))

test_slice_vs_random(biogrid_combined, (biogrid_combined['category'] == 'Synthetic Lethality') & (biogrid_combined['throughput'] == 'Low Throughput'),
                      "Synthetic Lethality, Low Throughput")
test_slice_vs_random(biogrid_combined, biogrid_combined['throughput'] == 'Low Throughput', "All Low Throughput")
test_slice_vs_random(biogrid_combined, biogrid_combined['category'] == 'Synthetic Lethality', "All Synthetic Lethality")
test_slice_vs_random(biogrid_combined, (biogrid_combined['category'] == 'Negative Genetic') & (biogrid_combined['throughput'] == 'High Throughput'),
                      "Negative Genetic, High Throughput (sanity check)")

# (Figures 7, 8 and 8b are drawn from biogrid_combined in the final 'ALL MANUSCRIPT FIGURES' cell.)

# %% [CELL 11b of 17 -- optional] BioGRID: protein vs RNA vs plain correlation for predicting synthetic lethal pairs
# #####################################################################
# CELL 11b (optional): which layer / method best predicts BioGRID synthetic lethal pairs?
#   FIGURE OUTPUT: FIGURE 9 (drawn in Cell 17)
#   DATA / DEPENDS ON: needs Cell 11 (biogrid_results, biogrid_pairs) and Cell 1 (matrices); caches to BIOGRID_METHODS_PATH
# #####################################################################
# DESIGN (fixed before running -- do not add methods, slices or thresholds after seeing the table):
#   Four methods, each asked the same question for every BioGRID test pair (target A, readout B): "when A is low, is B's
#   CRISPR dependency higher?"
#     1. Protein cohort     : SimKO default. 20 lowest vs 20 median cell lines by A's protein abundance.
#     2. RNA cohort         : the same, ranked by A's RNA expression.
#     3. Protein correlation: rank correlation of A's protein abundance with B's dependency across cell lines (sign
#                             flipped, so positive = B more essential when A is low).
#     4. RNA correlation    : the same, with A's RNA expression.
#   Every method uses the SAME cell lines (those in all of abundance, expression and dependency), the SAME real pairs and
#   the SAME random-control pairs as Cell 11, and its OWN shuffle test (N_PERM_METHODS shuffles of the cell lines).
#   Slices, fixed in advance: Synthetic Lethality / Low Throughput (PRIMARY); Synthetic Lethality / High Throughput and
#   Negative Genetic / High Throughput (negative controls); all real pairs.
#   Metrics: hit rate vs the random control at p<0.05 and p<0.01 (one-sided Fisher), fold enrichment, and ROC-AUC with a
#   bootstrap 95% CI for each method's difference from the protein cohort.
# Only the readout genes that are actually tested are computed for each target (not the whole genome), so this is fast.
N_PERM_METHODS = 200
N_BOOT_METHODS = 1000
if 'BIOGRID_METHODS_PATH' not in globals():   # defined in the DATA FILES cell; repeated here so this cell runs without re-running that one
    BIOGRID_METHODS_PATH = str(VAL_DIR / 'biogrid_method_comparison.csv')

def biogrid_method_pvalues(tests_df, abundance_df, expression_df, dependency_df, n_permutations=200, seed=42):
    """One empirical p per (target, readout) test for each of the four methods. One-sided: the predicted direction
    (readout more essential when the target is low). NaN where a method cannot be computed for that test."""
    common = [c for c in abundance_df.columns if c in expression_df.columns and c in dependency_df.columns]
    A, E = abundance_df[common], expression_df[common]
    D = dependency_df[common]
    D = D[~D.index.duplicated()]
    C = len(common)
    n = compute_cohort_n(C)
    print(f"Cell lines shared by Abundance, Expression and Dependency: {C} (cohort size n={n}).")
    D_vals = D.to_numpy(float)
    D_rank = D.rank(axis=1).to_numpy()                 # NaN stays NaN
    M = (~np.isnan(D_rank)).astype(float)
    Y0 = np.nan_to_num(D_rank, nan=0.0)
    dep_pos = {g: i for i, g in enumerate(D.index)}
    layers = {'protein': A, 'rna': E}

    def corr_stat(gi, xvecs):
        """-(rank correlation) of each readout (rows gi) with each column of xvecs (C x Q), over the cell lines where the
        readout is observed. Positive = readout more essential when the target is low."""
        Ms, Ys = M[gi], Y0[gi]
        n_g = Ms.sum(axis=1)[:, None]
        sx, sxx, sxy = Ms @ xvecs, Ms @ (xvecs ** 2), Ys @ xvecs
        sy, syy = Ys.sum(axis=1)[:, None], (Ys ** 2).sum(axis=1)[:, None]
        with np.errstate(divide='ignore', invalid='ignore'):
            r = (sxy - sx * sy / n_g) / np.sqrt((sxx - sx ** 2 / n_g) * (syy - sy ** 2 / n_g))
        r[(n_g < 20).ravel(), :] = np.nan
        return -r

    def emp_p(obs, null):
        out = np.full(len(obs), np.nan)
        for j in range(len(obs)):
            v = null[j][~np.isnan(null[j])]
            if len(v) >= 20 and not np.isnan(obs[j]):
                out[j] = float(np.mean(v >= obs[j]))
        return out

    tests = tests_df.reset_index(drop=True)
    cols = ('p_protein_cohort', 'p_rna_cohort', 'p_protein_corr', 'p_rna_corr')
    res = {c: np.full(len(tests), np.nan) for c in cols}
    rng = np.random.default_rng(seed)
    groups = tests.groupby('target').indices
    for ti, (target, rows) in enumerate(groups.items()):
        readouts = tests['readout'].to_numpy()[rows]
        ok = np.array([r in dep_pos for r in readouts])
        if not ok.any():
            continue
        gi = np.array([dep_pos[r] for r in readouts[ok]])
        rows_ok = rows[ok]
        perms = np.stack([rng.permutation(C) for _ in range(n_permutations)])      # P x C, shared by all four methods
        Dsub = D_vals[gi]
        for lname, L in layers.items():
            if target not in L.index:
                continue
            x = L.loc[target].to_numpy(float)
            if np.isnan(x).any():
                continue
            class_df = get_classes_by_joint_value([target], L, n, mode='max')
            low_pos = L.columns.get_indexer(class_df.index[class_df['class'] == 'low'])
            med_pos = L.columns.get_indexer(class_df.index[class_df['class'] == 'median'])
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore', message='Mean of empty slice')
                obs = np.nanmean(Dsub[:, low_pos], axis=1) - np.nanmean(Dsub[:, med_pos], axis=1)
                null = np.nanmean(Dsub[:, perms[:, :n]], axis=2) - np.nanmean(Dsub[:, perms[:, n:2 * n]], axis=2)
            res[f'p_{lname}_cohort'][rows_ok] = emp_p(obs, null)
            xr = pd.Series(x).rank().to_numpy()
            obs_c = corr_stat(gi, xr[:, None])[:, 0]
            null_c = corr_stat(gi, xr[perms.T])
            res[f'p_{lname}_corr'][rows_ok] = emp_p(obs_c, null_c)
        if (ti + 1) % 500 == 0:
            print(f"  ...{ti + 1}/{len(groups)} targets done")
    out = tests.copy()
    for c in cols:
        out[c] = res[c]
    return out

def _auc_low_p_is_real(p_real, p_random):
    from scipy.stats import rankdata
    r = rankdata(-np.concatenate([p_real, p_random]))
    npos, nneg = len(p_real), len(p_random)
    return (r[:npos].sum() - npos * (npos + 1) / 2) / (npos * nneg)

METHOD_LABELS = {'p_protein_cohort': 'Protein cohort (SimKO default)', 'p_rna_cohort': 'RNA cohort',
                 'p_protein_corr': 'Protein correlation', 'p_rna_corr': 'RNA correlation'}

if Path(BIOGRID_METHODS_PATH).exists():
    print(f"Reusing cached {BIOGRID_METHODS_PATH} -- delete it to force a fresh run.")
    _methods_raw = pd.read_csv(BIOGRID_METHODS_PATH)
else:
    _tests = biogrid_results[['target', 'readout', 'is_real']].drop_duplicates().reset_index(drop=True)
    print(f"{len(_tests)} tests ({int(_tests['is_real'].sum())} real, {int((~_tests['is_real']).sum())} random controls).")
    _methods_raw = biogrid_method_pvalues(_tests, abundance, expression, dependency, n_permutations=N_PERM_METHODS)
    _methods_raw.to_csv(BIOGRID_METHODS_PATH, index=False)

biogrid_methods_combined = attach_category_and_throughput(_methods_raw, biogrid_pairs)
_before = len(biogrid_methods_combined)
biogrid_methods_combined = biogrid_methods_combined.dropna(subset=list(METHOD_LABELS)).reset_index(drop=True)
print(f"Kept {len(biogrid_methods_combined)} of {_before} tests where all four methods could be computed.")

_cat, _tp = biogrid_methods_combined['category'], biogrid_methods_combined['throughput']
_is_real = _cat != 'random control'
# UNIT OF ANALYSIS: the two directions of one pair (A low -> B, B low -> A) are NOT independent tests. Bootstrap CIs below
# resample whole PAIRS (both directions together) so that they are not too narrow.
biogrid_methods_combined['pair_id'] = ['|'.join(sorted((a, b))) for a, b in zip(biogrid_methods_combined['target'], biogrid_methods_combined['readout'])]
def make_cluster_sampler(pair_ids):
    """Returns f(rng) -> row indices of a bootstrap resample that draws whole pairs (all their rows) with replacement."""
    codes, uniq = pd.factorize(pd.Series(pair_ids))
    order = np.argsort(codes, kind='stable')
    starts = np.r_[0, np.cumsum(np.bincount(codes))]
    groups = [order[starts[k]:starts[k + 1]] for k in range(len(uniq))]
    def sample(rng_):
        return np.concatenate([groups[k] for k in rng_.integers(0, len(groups), len(groups))])
    return sample
_SLICES = [('Synthetic Lethality, Low Throughput [primary]', (_cat == 'Synthetic Lethality') & (_tp == 'Low Throughput')),
           ('Synthetic Lethality, High Throughput', (_cat == 'Synthetic Lethality') & (_tp == 'High Throughput')),
           ('Negative Genetic, High Throughput', (_cat == 'Negative Genetic') & (_tp == 'High Throughput')),
           ('All real pairs', _is_real)]
_rand = biogrid_methods_combined[~_is_real]
_rows = []
for _slice_name, _mask in _SLICES:
    _sl = biogrid_methods_combined[_mask]
    for _col, _mlabel in METHOD_LABELS.items():
        for _t in (0.05, 0.01):
            _sh, _rh = int((_sl[_col] < _t).sum()), int((_rand[_col] < _t).sum())
            _fp = fisher_exact([[_sh, len(_sl) - _sh], [_rh, len(_rand) - _rh]], alternative='greater')[1]
            _sr, _rr = _sh / len(_sl), _rh / len(_rand)
            _rows.append({'slice': _slice_name, 'method': _mlabel, 'threshold': _t, 'n_real': len(_sl),
                          'real_rate': _sr, 'random_rate': _rr, 'fold': _sr / _rr if _rr > 0 else np.nan,
                          'fisher_p': _fp, 'auc': _auc_low_p_is_real(_sl[_col].to_numpy(), _rand[_col].to_numpy())})
biogrid_method_summary = pd.DataFrame(_rows)

print("\n" + "=" * 100 + "\nBioGRID: hit rate and fold enrichment over the random control, by method\n" + "=" * 100)
for _slice_name, _ in _SLICES:
    print(f"\n{_slice_name}  (n={int(biogrid_method_summary.loc[biogrid_method_summary['slice'] == _slice_name, 'n_real'].iloc[0])})")
    _s = biogrid_method_summary[biogrid_method_summary['slice'] == _slice_name]
    print(f"  {'method':<32}{'p<0.05 rate':>12}{'(random)':>10}{'fold':>7}{'Fisher p':>10}   {'p<0.01 rate':>12}{'(random)':>10}{'fold':>7}{'Fisher p':>10}   {'ROC-AUC':>8}")
    for _mlabel in METHOD_LABELS.values():
        _a = _s[(_s['method'] == _mlabel) & (_s['threshold'] == 0.05)].iloc[0]
        _b = _s[(_s['method'] == _mlabel) & (_s['threshold'] == 0.01)].iloc[0]
        print(f"  {_mlabel:<32}{_a['real_rate']:>12.4f}{_a['random_rate']:>10.4f}{_a['fold']:>7.2f}{_a['fisher_p']:>10.4f}   "
              f"{_b['real_rate']:>12.4f}{_b['random_rate']:>10.4f}{_b['fold']:>7.2f}{_b['fisher_p']:>10.4f}   {_a['auc']:>8.3f}")

# Direct comparison on the PRIMARY slice: ROC-AUC difference from the protein cohort, bootstrap over pairs (the same resampled
# pairs for every method, so the methods stay paired).
_prim = biogrid_methods_combined[_SLICES[0][1]]
_brng = np.random.default_rng(7)
_boot = {c: [] for c in METHOD_LABELS}
_prim_sampler = make_cluster_sampler(_prim['pair_id'].to_numpy())      # resample whole pairs (both directions together)
_boot_hit = {c: [] for c in METHOD_LABELS}
for _ in range(N_BOOT_METHODS):
    _ri, _qi = _prim_sampler(_brng), _brng.integers(0, len(_rand), len(_rand))
    for _col in METHOD_LABELS:
        _boot[_col].append(_auc_low_p_is_real(_prim[_col].to_numpy()[_ri], _rand[_col].to_numpy()[_qi]))
        _boot_hit[_col].append((_prim[_col].to_numpy()[_ri] < 0.05).mean() - (_rand[_col].to_numpy()[_qi] < 0.05).mean())
print(f"\nPRIMARY slice: {len(_prim)} directional tests from {_prim['pair_id'].nunique()} unique pairs. The Fisher p-values above treat every")
print("directional test as independent, which overstates the evidence. Pair-clustered version (resamples whole pairs):")
print(f"  {'method':<32}{'hit rate minus random rate at p<0.05':>40}{'95% CI':>22}{'one-sided p':>14}")
for _col, _mlabel in METHOD_LABELS.items():
    _bh = np.array(_boot_hit[_col])
    _lo, _hi = np.percentile(_bh, [2.5, 97.5])
    print(f"  {_mlabel:<32}{(_prim[_col] < 0.05).mean() - (_rand[_col] < 0.05).mean():>40.4f}{'[' + format(_lo, '+.3f') + ', ' + format(_hi, '+.3f') + ']':>22}{(np.sum(_bh <= 0) + 1) / (len(_bh) + 1):>14.4f}")
print("\nPRIMARY slice (Synthetic Lethality, Low Throughput): ROC-AUC, and difference from the protein cohort (bootstrap 95% CI)")
for _col, _mlabel in METHOD_LABELS.items():
    _auc = _auc_low_p_is_real(_prim[_col].to_numpy(), _rand[_col].to_numpy())
    _d = np.array(_boot[_col]) - np.array(_boot['p_protein_cohort'])
    if _col == 'p_protein_cohort':
        print(f"  {_mlabel:<32}AUC={_auc:.3f}")
    else:
        _lo, _hi = np.percentile(_d, [2.5, 97.5])
        print(f"  {_mlabel:<32}AUC={_auc:.3f}   difference from protein cohort={_auc - _auc_low_p_is_real(_prim['p_protein_cohort'].to_numpy(), _rand['p_protein_cohort'].to_numpy()):+.3f}  95% CI [{_lo:+.3f}, {_hi:+.3f}]")
print("(A CI that excludes 0 is a real difference between the methods. If it includes 0, say the methods performed similarly.)")

# Sensitivity and specificity. Positives = real BioGRID pairs in the slice; negatives = the random-control pairs.
#   Sensitivity = share of real pairs called a hit.  Specificity = share of random pairs NOT called a hit.
# Part 1: at the fixed p thresholds (specificity then differs by method, so it is not a like-for-like comparison).
# Part 2: at MATCHED specificity (95% and 90%): each method's own cut-off is set so that 5% / 10% of random pairs are called,
#         then sensitivity is compared. This is the fair comparison, because the correlation methods call more random pairs.
def _sens_at_spec(p_real, p_random, spec):
    """Sensitivity at the most permissive cut-off whose specificity is still >= spec. Returns (sensitivity, achieved specificity)."""
    best = None
    for t in np.unique(p_random):
        if (p_random <= t).mean() <= (1 - spec) + 1e-12:
            best = t
    if best is None:
        return 0.0, 1.0
    return float((p_real <= best).mean()), float((p_random > best).mean())

_ss_rows = []
for _slice_name, _mask in _SLICES:
    _sl = biogrid_methods_combined[_mask]
    for _col, _mlabel in METHOD_LABELS.items():
        _pr, _pn = _sl[_col].to_numpy(), _rand[_col].to_numpy()
        _r = {'slice': _slice_name, 'method': _mlabel, 'n_real': len(_sl)}
        for _t in (0.05, 0.01):
            _r[f'sens_p{_t}'], _r[f'spec_p{_t}'] = float((_pr < _t).mean()), float((_pn >= _t).mean())
        for _sp in (0.95, 0.90):
            _r[f'sens_at_spec{int(_sp * 100)}'], _r[f'achieved_spec{int(_sp * 100)}'] = _sens_at_spec(_pr, _pn, _sp)
        _ss_rows.append(_r)
biogrid_method_sensspec = pd.DataFrame(_ss_rows)

print("\n" + "=" * 100 + "\nBioGRID: sensitivity and specificity by method (negatives = random-control pairs)\n" + "=" * 100)
for _slice_name, _ in _SLICES:
    _s = biogrid_method_sensspec[biogrid_method_sensspec['slice'] == _slice_name]
    print(f"\n{_slice_name}  (n={int(_s['n_real'].iloc[0])} real vs {len(_rand)} random)")
    print(f"  {'method':<32}{'sens p<.05':>11}{'spec p<.05':>11}{'sens p<.01':>11}{'spec p<.01':>11}   {'sens@95%spec':>13}{'(achieved)':>11}  {'sens@90%spec':>13}{'(achieved)':>11}")
    for _, _x in _s.iterrows():
        print(f"  {_x['method']:<32}{_x['sens_p0.05']:>11.3f}{_x['spec_p0.05']:>11.3f}{_x['sens_p0.01']:>11.3f}{_x['spec_p0.01']:>11.3f}   "
              f"{_x['sens_at_spec95']:>13.3f}{_x['achieved_spec95']:>11.3f}  {_x['sens_at_spec90']:>13.3f}{_x['achieved_spec90']:>11.3f}")

# Primary slice: sensitivity at 95% specificity, difference from the protein cohort (same paired bootstrap design as the AUC).
_bs = {c: [] for c in METHOD_LABELS}
_brng2 = np.random.default_rng(11)
for _ in range(N_BOOT_METHODS):
    _ri, _qi = _prim_sampler(_brng2), _brng2.integers(0, len(_rand), len(_rand))
    for _col in METHOD_LABELS:
        _bs[_col].append(_sens_at_spec(_prim[_col].to_numpy()[_ri], _rand[_col].to_numpy()[_qi], 0.95)[0])
print("\nPRIMARY slice: sensitivity at 95% specificity, and difference from the protein cohort (bootstrap 95% CI)")
_base = _sens_at_spec(_prim['p_protein_cohort'].to_numpy(), _rand['p_protein_cohort'].to_numpy(), 0.95)[0]
for _col, _mlabel in METHOD_LABELS.items():
    _sv = _sens_at_spec(_prim[_col].to_numpy(), _rand[_col].to_numpy(), 0.95)[0]
    if _col == 'p_protein_cohort':
        print(f"  {_mlabel:<32}sensitivity={_sv:.3f}")
    else:
        _lo, _hi = np.percentile(np.array(_bs[_col]) - np.array(_bs['p_protein_cohort']), [2.5, 97.5])
        print(f"  {_mlabel:<32}sensitivity={_sv:.3f}   difference from protein cohort={_sv - _base:+.3f}  95% CI [{_lo:+.3f}, {_hi:+.3f}]")
print("(Permutation p-values come in steps of 1/N_PERM_METHODS, so 'achieved' specificity can sit slightly above the target.)")

# Precision at the top of the ranking (defined BEFORE running; report whichever way it comes out).
#   Pool = the slice's real pairs + an EQUAL number of random-control pairs, so chance = 50%.
#   Each method ranks the pool by its own p-value (ties broken at random). Precision@top-F = share of real pairs in the top F.
#   F = 10% and 25%. Repeated over many random-control subsamples. Primary slice: real pairs are also resampled (bootstrap)
#   so the 95% interval and the difference from the protein cohort include uncertainty about the real pairs.
def _precision_top(p_real, p_rand, fracs, rng):
    allp = np.concatenate([p_real, p_rand])
    lab = np.r_[np.ones(len(p_real)), np.zeros(len(p_rand))]
    order = np.lexsort((rng.random(len(allp)), allp))
    return [float(lab[order[:max(1, int(round(f * len(allp))))]].mean()) for f in fracs]

_FRACS = (0.10, 0.25)
_prng = np.random.default_rng(23)
_prim_sampler_p = make_cluster_sampler(biogrid_methods_combined[_SLICES[0][1]]['pair_id'].to_numpy())
biogrid_method_precision = {}
print("\n" + "=" * 100 + "\nBioGRID: precision at the top of the ranking (real pairs vs an equal number of random pairs; chance = 0.50)\n" + "=" * 100)
for _si, (_slice_name, _mask) in enumerate(_SLICES):
    _sl = biogrid_methods_combined[_mask]
    _is_primary = (_si == 0)
    _nit = 1000 if _is_primary else 200
    _res = {c: np.zeros((_nit, len(_FRACS))) for c in METHOD_LABELS}
    for _it in range(_nit):
        _ri = _prim_sampler_p(_prng) if _is_primary else np.arange(len(_sl))
        _qi = _prng.choice(len(_rand), size=len(_ri), replace=False)
        for _col in METHOD_LABELS:
            _res[_col][_it] = _precision_top(_sl[_col].to_numpy()[_ri], _rand[_col].to_numpy()[_qi], _FRACS, _prng)
    biogrid_method_precision[_slice_name] = _res
    print(f"\n{_slice_name}  (n={len(_sl)} real + {len(_sl)} random per draw)")
    print(f"  {'method':<32}" + ''.join(f"{'top ' + str(int(f * 100)) + '%':>10}" + ("   95% interval        " if _is_primary else '') for f in _FRACS))
    for _col, _mlabel in METHOD_LABELS.items():
        _line = f"  {_mlabel:<32}"
        for _j, _f in enumerate(_FRACS):
            _m = _res[_col][:, _j].mean()
            _line += f"{_m:>10.3f}"
            if _is_primary:
                _lo, _hi = np.percentile(_res[_col][:, _j], [2.5, 97.5])
                _line += f"   [{_lo:.3f}, {_hi:.3f}]   "
        print(_line)
    if _is_primary:
        print("  Difference from the protein cohort (95% CI):")
        for _col, _mlabel in METHOD_LABELS.items():
            if _col == 'p_protein_cohort':
                continue
            _line = f"    {_mlabel:<30}"
            for _j, _f in enumerate(_FRACS):
                _d = _res[_col][:, _j] - _res['p_protein_cohort'][:, _j]
                _lo, _hi = np.percentile(_d, [2.5, 97.5])
                _line += f"top {int(_f * 100)}%: {_d.mean():+.3f} [{_lo:+.3f}, {_hi:+.3f}]   "
            print(_line)
print("(Chance is 0.50. A CI for a difference that includes 0 means the methods performed similarly.)")

# Coverage: how much of BioGRID can SimKO's Abundance layer actually test?
#   A directional test (target low -> readout dependency) needs the TARGET in the Abundance matrix and the READOUT in the Dependency matrix.
#   Each pair can be tested in two directions. Pairs that cannot be tested in either direction are silently absent from every table above.
_ab, _dp = set(abundance.index), set(dependency.index)
_bp = biogrid_pairs.copy()
_bp['fwd'] = _bp['gene_a'].isin(_ab) & _bp['gene_b'].isin(_dp)
_bp['rev'] = _bp['gene_b'].isin(_ab) & _bp['gene_a'].isin(_dp)
_bp['any'] = _bp['fwd'] | _bp['rev']
_bp['both'] = _bp['fwd'] & _bp['rev']
_bp['a_in_ab'] = _bp['gene_a'].isin(_ab)
_bp['b_in_ab'] = _bp['gene_b'].isin(_ab)
_bp['neither_in_ab'] = ~_bp['a_in_ab'] & ~_bp['b_in_ab']
_used = biogrid_methods_combined[_is_real].copy()
_used['pair'] = [tuple(sorted(x)) for x in zip(_used['target'], _used['readout'])]
_bp['pair'] = [tuple(sorted(x)) for x in zip(_bp['gene_a'], _bp['gene_b'])]
print("\n" + "=" * 100 + "\nBioGRID coverage: how many pairs can the Abundance layer test?\n" + "=" * 100)
print(f"{'slice':<46}{'pairs':>7}{'testable':>10}{'%':>6}{'both dirs':>11}{'neither gene':>14}{'pairs used in 11b':>19}")
_cov_slices = [('Synthetic Lethality, Low Throughput', (_bp['category'] == 'Synthetic Lethality') & (_bp['throughput'] == 'Low Throughput')),
               ('Synthetic Lethality, High Throughput', (_bp['category'] == 'Synthetic Lethality') & (_bp['throughput'] == 'High Throughput')),
               ('Negative Genetic, High Throughput', (_bp['category'] == 'Negative Genetic') & (_bp['throughput'] == 'High Throughput')),
               ('All pairs', pd.Series(True, index=_bp.index))]
biogrid_coverage = []
for _nm, _m in _cov_slices:
    _g = _bp[_m].drop_duplicates(subset=['pair'])
    _u = _used[_used['pair'].isin(set(_g['pair']))]['pair'].nunique() if _nm != 'All pairs' else _used['pair'].nunique()
    _row = (_nm, len(_g), int(_g['any'].sum()), 100 * _g['any'].mean(), int(_g['both'].sum()), int(_g['neither_in_ab'].sum()), _u)
    biogrid_coverage.append(_row)
    print(f"{_row[0]:<46}{_row[1]:>7}{_row[2]:>10}{_row[3]:>6.1f}{_row[4]:>11}{_row[5]:>14}{_row[6]:>19}")
print("(testable = at least one direction can be tested; neither gene = neither partner is in the Abundance matrix.)")
print("Note: the tables above count directional tests. A pair tested in both directions appears twice.")

# %% [CELL 11c of 17 -- optional] SynLethDB: the same four-method comparison on a second synthetic lethality database
# #####################################################################
# CELL 11c (optional): SynLethDB replication of the four-method comparison
#   FIGURE OUTPUT: (no figure; supplementary tables / Results text)
#   DATA / DEPENDS ON: SynLethDB SL edge file (tab-separated, set SYNLETHDB_PATH below); needs Cell 11b run first
#                      (it re-uses biogrid_method_pvalues, METHOD_LABELS, _auc_low_p_is_real, _sens_at_spec, _precision_top)
# #####################################################################
# DESIGN (fixed before running -- do not add sources, slices or metrics after seeing the output):
#   Input: SynLethDB rows with relation == 'SL'. Evidence types kept: CRISPR/CRISPRi, GenomeRNAi, Decipher, Low Throughput,
#   High Throughput. 'Text Mining' and 'Computational Prediction' are EXCLUDED always: they are not experimental evidence, and
#   computational predictions can be derived from the same kind of data we are testing (circular).
#   Tests: every pair in both directions (A low -> B more essential), gene symbols matched to the matrices; untestable pairs dropped.
#   Random controls: for each real test, 3 readouts drawn from the dependency genes, excluding the target and every known partner
#   of that target in SynLethDB or BioGRID (same idea as Cell 11, matched per target).
#   Slices: (1) CRISPR/CRISPRi [PRIMARY]; (2) CRISPR/CRISPRi pairs NOT also in BioGRID (checks independence from Cell 11b);
#           (3) RNAi (GenomeRNAi + Decipher); (4) all experimental evidence.
#   Metrics: as Cell 11b -- fold enrichment and Fisher p vs random at p<0.05, ROC-AUC, sensitivity at 95% specificity,
#   precision in the top 10% / 25% of a 1:1 real:random pool -- with bootstrap CIs for the differences on the PRIMARY slice.
#   Note: pairs measured in ONE cell line (e.g. a K562 double-knockdown screen) are induced synthetic lethality; the cohort
#   method asks about naturally low abundance across cell lines. This is the same assumption as everywhere else in the paper.
if 'SYNLETHDB_PATH' not in globals():
    SYNLETHDB_PATH = str(VAL_DIR / 'Human.SL.detailed.tsv')            # SynLethDB human SL edge file (tab-separated)
SYNLETHDB_METHODS_PATH = str(VAL_DIR / 'synlethdb_method_comparison.csv')   # cache -- KEEP this file
N_RAND_PER_TEST_SLDB = 3

def _find_synlethdb(configured):
    """Use the configured path if it exists; otherwise look for Human.SL.detailed.tsv (or any Human.SL*.tsv) near the script/notebook."""
    if Path(configured).exists():
        return str(configured)
    roots = [Path(VAL_DIR), Path(VAL_DIR).parent, Path.cwd(), Path.cwd().parent]
    for root in roots:
        for pat in ('Human.SL.detailed.tsv', 'Human.SL*.tsv', '**/Human.SL.detailed.tsv'):
            hits = sorted(root.glob(pat)) if root.exists() else []
            if hits:
                print(f"Configured SYNLETHDB_PATH not found ({configured}); using {hits[0]} instead.")
                return str(hits[0])
    raise FileNotFoundError(f"SynLethDB file not found. Tried: {configured}\nand searched {[str(r) for r in roots]} for Human.SL*.tsv.\n"
                            f"Set SYNLETHDB_PATH = '/full/path/to/Human.SL.detailed.tsv' above this line and re-run.")
SYNLETHDB_PATH = _find_synlethdb(SYNLETHDB_PATH)

import csv as _csv
_sraw = pd.read_csv(SYNLETHDB_PATH, sep='\t', dtype=str, low_memory=False, quoting=_csv.QUOTE_NONE)   # QUOTE_NONE: a stray quote in a text field cannot swallow rows
print(f"Read {len(_sraw)} rows from {SYNLETHDB_PATH}. Columns: {list(_sraw.columns)}")
for _c in _sraw.columns:
    _sraw[_c] = _sraw[_c].str.strip().str.strip('"')
_need = ['x_name', 'y_name', 'relation', 'rel_source']
_miss = [c for c in _need if c not in _sraw.columns]
if _miss:
    raise ValueError(f"SynLethDB file is missing columns {_miss}. Columns found: {list(_sraw.columns)}")
_sraw = _sraw[_sraw['relation'].str.strip() == 'SL'].copy()
print(f"SynLethDB: {len(_sraw)} SL rows. By evidence source:")
print(_sraw['rel_source'].value_counts().to_string())

_CRISPR, _RNAI, _OTHER = {'CRISPR/CRISPRi'}, {'GenomeRNAi', 'Decipher'}, {'Low Throughput', 'High Throughput'}
_sall = _sraw.rename(columns={'x_name': 'gene_a', 'y_name': 'gene_b'}).dropna(subset=['gene_a', 'gene_b'])
_sall = _sall[_sall['gene_a'] != _sall['gene_b']]
_sexp = _sall[_sall['rel_source'].isin(_CRISPR | _RNAI | _OTHER)].copy()
print(f"\nKept {len(_sexp)} experimental rows ({len(_sall) - len(_sexp)} Text Mining / Computational Prediction rows excluded).")
_sexp['pair'] = [tuple(sorted(x)) for x in zip(_sexp['gene_a'], _sexp['gene_b'])]
_sexp['crispr'] = _sexp['rel_source'].isin(_CRISPR)
_sexp['rnai'] = _sexp['rel_source'].isin(_RNAI)
print(f"{_sexp['pair'].nunique()} unique experimental pairs; CRISPR/CRISPRi: {_sexp.loc[_sexp['crispr'], 'pair'].nunique()}, RNAi: {_sexp.loc[_sexp['rnai'], 'pair'].nunique()}")

_bg_pairs = set(tuple(sorted(x)) for x in zip(biogrid_pairs['gene_a'], biogrid_pairs['gene_b']))
_pf = _sexp.groupby('pair').agg(crispr=('crispr', 'max'), rnai=('rnai', 'max')).reset_index()
_pf['in_biogrid'] = _pf['pair'].isin(_bg_pairs)
print(f"{int(_pf['in_biogrid'].sum())} of {len(_pf)} experimental SynLethDB pairs are also in BioGRID "
      f"({int((_pf['in_biogrid'] & _pf['crispr']).sum())} of {int(_pf['crispr'].sum())} CRISPR/CRISPRi pairs).")

_ab_set, _dp_set = set(abundance.index), set(dependency.index)
_real = []
for _, _r in _pf.iterrows():
    a, b = _r['pair']
    for t, u in ((a, b), (b, a)):
        if t in _ab_set and u in _dp_set:
            _real.append({'target': t, 'readout': u, 'is_real': True, 'crispr': bool(_r['crispr']), 'rnai': bool(_r['rnai']), 'in_biogrid': bool(_r['in_biogrid'])})
_real = pd.DataFrame(_real).groupby(['target', 'readout', 'is_real'], as_index=False).agg(crispr=('crispr', 'max'), rnai=('rnai', 'max'), in_biogrid=('in_biogrid', 'max'))
_testable_pairs = set(tuple(sorted(x)) for x in zip(_real['target'], _real['readout']))
print(f"Testable: {len(_testable_pairs)} of {len(_pf)} pairs ({100 * len(_testable_pairs) / len(_pf):.1f}%), {len(_real)} directional tests "
      f"(CRISPR/CRISPRi: {int(_real['crispr'].sum())}).")

_partners = {}
for _a, _b in list(_sall[['gene_a', 'gene_b']].itertuples(index=False, name=None)) + list(_bg_pairs):
    _partners.setdefault(_a, set()).add(_b)
    _partners.setdefault(_b, set()).add(_a)
_dep_genes = np.array(sorted(_dp_set))
_srng = np.random.default_rng(42)
_rand_rows = []
for _tgt, _grp in _real.groupby('target'):
    _cand = _dep_genes[~np.isin(_dep_genes, list(_partners.get(_tgt, set()) | {_tgt}))]
    _k = min(len(_grp) * N_RAND_PER_TEST_SLDB, len(_cand))
    for _g in _srng.choice(_cand, size=_k, replace=False):
        _rand_rows.append({'target': _tgt, 'readout': _g, 'is_real': False, 'crispr': False, 'rnai': False, 'in_biogrid': False})
_stests = pd.concat([_real, pd.DataFrame(_rand_rows)], ignore_index=True)

if Path(SYNLETHDB_METHODS_PATH).exists():
    print(f"Reusing cached {SYNLETHDB_METHODS_PATH} -- delete it to force a fresh run.")
    _sres = pd.read_csv(SYNLETHDB_METHODS_PATH)
else:
    print(f"{len(_stests)} tests ({int(_stests['is_real'].sum())} real, {int((~_stests['is_real']).sum())} random controls).")
    _sres = biogrid_method_pvalues(_stests, abundance, expression, dependency, n_permutations=N_PERM_METHODS)
    _sres.to_csv(SYNLETHDB_METHODS_PATH, index=False)
_before = len(_sres)
synlethdb_methods = _sres.dropna(subset=list(METHOD_LABELS)).reset_index(drop=True)
print(f"Kept {len(synlethdb_methods)} of {_before} tests where all four methods could be computed.")

synlethdb_methods['pair_id'] = ['|'.join(sorted((x, y))) for x, y in zip(synlethdb_methods['target'], synlethdb_methods['readout'])]
_s_real = synlethdb_methods['is_real'].astype(bool)
_SL_SLICES = {'CRISPR/CRISPRi [primary]': _s_real & synlethdb_methods['crispr'].astype(bool),
              'CRISPR/CRISPRi, not also in BioGRID': _s_real & synlethdb_methods['crispr'].astype(bool) & ~synlethdb_methods['in_biogrid'].astype(bool),
              'RNAi (GenomeRNAi + Decipher)': _s_real & synlethdb_methods['rnai'].astype(bool),
              'All experimental evidence': _s_real}
_s_rand = synlethdb_methods[~_s_real]
_primary_name = 'CRISPR/CRISPRi [primary]'
_srng2 = np.random.default_rng(31)
print("\n" + "=" * 112 + "\nSynLethDB: four-method comparison (negatives = random-control pairs)\n" + "=" * 112)
for _name, _mask in _SL_SLICES.items():
    _sl_df = synlethdb_methods[_mask]
    print(f"\n{_name}  (n={len(_sl_df)} real tests vs {len(_s_rand)} random)")
    if len(_sl_df) < 30:
        print("  (fewer than 30 tests -- not evaluated)")
        continue
    print(f"  {'method':<32}{'hit .05':>8}{'(rand)':>8}{'fold':>6}{'Fisher p':>10}{'AUC':>7}{'sens@95%':>10}{'prec top10%':>13}{'prec top25%':>13}")
    for _col, _lab in METHOD_LABELS.items():
        _pr, _pn = _sl_df[_col].to_numpy(), _s_rand[_col].to_numpy()
        _h, _rh = int((_pr < 0.05).sum()), int((_pn < 0.05).sum())
        _fp = fisher_exact([[_h, len(_pr) - _h], [_rh, len(_pn) - _rh]], alternative='greater')[1]
        _pre = np.mean([_precision_top(_pr, _pn[_srng2.choice(len(_pn), size=len(_pr), replace=False)], (0.10, 0.25), _srng2) for _ in range(200)], axis=0)
        print(f"  {_lab:<32}{_h / len(_pr):>8.3f}{_rh / len(_pn):>8.3f}{(_h / len(_pr)) / (_rh / len(_pn)) if _rh else float('nan'):>6.2f}{_fp:>10.4f}"
              f"{_auc_low_p_is_real(_pr, _pn):>7.3f}{_sens_at_spec(_pr, _pn, 0.95)[0]:>10.3f}{_pre[0]:>13.3f}{_pre[1]:>13.3f}")

_prim_df = synlethdb_methods[_SL_SLICES[_primary_name]]
if len(_prim_df) >= 30:
    _bt = {c: {'auc': [], 'sens': [], 'prec': []} for c in METHOD_LABELS}
    _sl_sampler = make_cluster_sampler(_prim_df['pair_id'].to_numpy())      # resample whole pairs (both directions together)
    for _ in range(N_BOOT_METHODS):
        _ri, _qi = _sl_sampler(_srng2), _srng2.integers(0, len(_s_rand), len(_s_rand))
        _q2 = _srng2.choice(len(_s_rand), size=len(_ri), replace=False)
        for _col in METHOD_LABELS:
            _pr, _pn = _prim_df[_col].to_numpy()[_ri], _s_rand[_col].to_numpy()
            _bt[_col]['auc'].append(_auc_low_p_is_real(_pr, _pn[_qi]))
            _bt[_col]['sens'].append(_sens_at_spec(_pr, _pn[_qi], 0.95)[0])
            _bt[_col]['prec'].append(_precision_top(_pr, _pn[_q2], (0.10,), _srng2)[0])
    print(f"\nPRIMARY slice ({_primary_name}): difference from the protein cohort (bootstrap 95% CI)")
    print(f"  {'method':<32}{'AUC diff':>30}{'sens@95% diff':>30}{'precision top10% diff':>30}")
    for _col, _lab in METHOD_LABELS.items():
        if _col == 'p_protein_cohort':
            continue
        _cells = []
        for _k in ('auc', 'sens', 'prec'):
            _d = np.array(_bt[_col][_k]) - np.array(_bt['p_protein_cohort'][_k])
            _lo, _hi = np.percentile(_d, [2.5, 97.5])
            _cells.append(f"{_d.mean():+.3f} [{_lo:+.3f}, {_hi:+.3f}]")
        print(f"  {_lab:<32}{_cells[0]:>30}{_cells[1]:>30}{_cells[2]:>30}")
    print("(A CI that includes 0 means the methods performed similarly.)")

# %% [CELL 11d of 17 -- optional] BioGRID: reference-group width and tissue-stratified null (protein layer)
# #####################################################################
# CELL 11d (optional): does a wider reference group and/or a tissue-stratified null change the BioGRID result?
#   FIGURE OUTPUT: (no figure; supplementary table / Results text)
#   DATA / DEPENDS ON: needs Cell 11b run first (biogrid_methods_combined, METHOD_LABELS and its helper functions);
#                      cell-line names must end in the tissue (e.g. HGC27_STOMACH); caches to VARIANT_PATH
# #####################################################################
# DESIGN (fixed before running -- do not add variants, cut-offs or slices after seeing the table):
#   Protein layer only. The same real pairs, random controls, cell lines and slices as Cell 11b.
#   Low group = the 20 lowest lines by the target's protein abundance (as before). What varies is (a) the REFERENCE group and
#   (b) the NULL:
#     1. Median 20,  pooled null            -- the current method (Cell 11b, 'p_protein_cohort')
#     2. Median 20,  tissue-stratified null -- isolates the effect of stratification
#     3. Central 40% of lines (ranks 30-70%), tissue-stratified null
#     4. Top 50% of lines (ranks 50-100%), tissue-stratified null
#     5. Protein correlation, pooled null          (Cell 11b, 'p_protein_corr')
#     6. Protein correlation, tissue-stratified null
#   Tissue-stratified null: the target's abundance values are shuffled ONLY among cell lines of the same tissue, then the low
#   and reference groups are re-chosen from the shuffled values. Tissue structure is therefore preserved under the null, so
#   a lineage that drives both abundance and dependency cannot, by itself, produce a hit.
#   Metrics and slices as Cell 11b; differences from variant 1 are reported with bootstrap 95% CIs on the PRIMARY slice.
N_PERM_VARIANTS = 500
VARIANT_PATH = str(VAL_DIR / 'biogrid_cohort_variants.csv')      # cache -- KEEP this file

VARIANT_LABELS = {'p_protein_cohort': '1 Median 20, pooled null (current)',
                  'q_med20_strat': '2 Median 20, tissue-stratified',
                  'q_central40_strat': '3 Central 40%, tissue-stratified',
                  'q_top50_strat': '4 Top 50%, tissue-stratified',
                  'p_protein_corr': '5 Correlation, pooled null',
                  'q_corr_strat': '6 Correlation, tissue-stratified'}

def cohort_variant_pvalues(tests_df, abundance_df, expression_df, dependency_df, n_permutations=500, seed=42):
    """Empirical one-sided p (readout more essential when the target is low) for the tissue-stratified variants.
    Uses the same cell lines as Cell 11b (those in abundance, expression and dependency)."""
    from scipy.stats import rankdata
    common = [c for c in abundance_df.columns if c in expression_df.columns and c in dependency_df.columns]
    A = abundance_df[common]
    D = dependency_df[common]
    D = D[~D.index.duplicated()]
    C = len(common)
    n = compute_cohort_n(C)
    tissue = np.array([c.split('_', 1)[1] if '_' in c else 'UNKNOWN' for c in common])
    tgroups = [np.where(tissue == t)[0] for t in np.unique(tissue)]
    print(f"{C} cell lines, {len(tgroups)} tissues (largest {max(len(g) for g in tgroups)}, smallest {min(len(g) for g in tgroups)}); cohort size n={n}.")
    if (tissue == 'UNKNOWN').any():
        print(f"WARNING: {int((tissue == 'UNKNOWN').sum())} cell-line names have no tissue suffix and are pooled together.")
    rng = np.random.default_rng(seed)
    P = n_permutations
    perm_idx = np.tile(np.arange(C), (P, 1))                 # P x C: shuffles WITHIN tissue
    for g in tgroups:
        if len(g) > 1:
            for i in range(P):
                perm_idx[i, g] = rng.permutation(g)
    D_vals = D.to_numpy(float)
    D_rank = D.rank(axis=1).to_numpy()
    M = (~np.isnan(D_rank)).astype(float)
    Y0 = np.nan_to_num(D_rank, nan=0.0)
    dep_pos = {g: i for i, g in enumerate(D.index)}

    def corr_stat(gi, xvecs):
        Ms, Ys = M[gi], Y0[gi]
        n_g = Ms.sum(axis=1)[:, None]
        sx, sxx, sxy = Ms @ xvecs, Ms @ (xvecs ** 2), Ys @ xvecs
        sy, syy = Ys.sum(axis=1)[:, None], (Ys ** 2).sum(axis=1)[:, None]
        with np.errstate(divide='ignore', invalid='ignore'):
            r = (sxy - sx * sy / n_g) / np.sqrt((sxx - sx ** 2 / n_g) * (syy - sy ** 2 / n_g))
        r[(n_g < 20).ravel(), :] = np.nan
        return -r

    def emp_p(obs, null):
        out = np.full(len(obs), np.nan)
        for j in range(len(obs)):
            v = null[j][~np.isnan(null[j])]
            if len(v) >= 20 and not np.isnan(obs[j]):
                out[j] = float(np.mean(v >= obs[j]))
        return out

    ref_slices = {'q_med20_strat': slice(C // 2 - n // 2, C // 2 - n // 2 + n),
                  'q_central40_strat': slice(int(C * 0.3), int(C * 0.7)),
                  'q_top50_strat': slice(C // 2, C)}
    tests = tests_df.reset_index(drop=True)
    cols = list(ref_slices) + ['q_corr_strat']
    res = {c: np.full(len(tests), np.nan) for c in cols}
    by_target = tests.groupby('target').indices
    for ti, (target, rows) in enumerate(by_target.items()):
        if target not in A.index:
            continue
        x = A.loc[target].to_numpy(float)
        if np.isnan(x).any():
            continue
        readouts = tests['readout'].to_numpy()[rows]
        ok = np.array([r in dep_pos for r in readouts])
        if not ok.any():
            continue
        gi_all = np.array([dep_pos[r] for r in readouts[ok]])
        rows_ok = rows[ok]
        o = np.argsort(x, kind='stable')
        O = np.argsort(x[perm_idx], axis=1, kind='stable')   # P x C: group order under each within-tissue shuffle
        xr = rankdata(x)
        for s in range(0, len(gi_all), 100):                 # chunk readouts to bound memory
            gi, rws = gi_all[s:s + 100], rows_ok[s:s + 100]
            Dsub = D_vals[gi]
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore', message='Mean of empty slice')
                low_obs = np.nanmean(Dsub[:, o[:n]], axis=1)
                low_null = np.nanmean(Dsub[:, O[:, :n]], axis=2)
                for name, sl in ref_slices.items():
                    obs = low_obs - np.nanmean(Dsub[:, o[sl]], axis=1)
                    null = low_null - np.nanmean(Dsub[:, O[:, sl]], axis=2)
                    res[name][rws] = emp_p(obs, null)
            obs_c = corr_stat(gi, xr[:, None])[:, 0]
            null_c = corr_stat(gi, xr[perm_idx].T)
            res['q_corr_strat'][rws] = emp_p(obs_c, null_c)
        if (ti + 1) % 500 == 0:
            print(f"  ...{ti + 1}/{len(by_target)} targets done")
    out = tests.copy()
    for c in cols:
        out[c] = res[c]
    return out

def evaluate_columns(df, slices, labels, base_col, primary_name, n_boot=1000, seed=19):
    """Hit rate / fold / Fisher / AUC / sensitivity at 95% specificity / precision at the top, for each p-value column in
    `labels`, on each slice (negatives = the random controls in df). Bootstrap CIs for differences from base_col on the primary slice."""
    rand = df[~df['is_real'].astype(bool)]
    rng_ = np.random.default_rng(seed)
    for name, mask in slices.items():
        sl = df[mask]
        print(f"\n{name}  (n={len(sl)} real tests vs {len(rand)} random)")
        if len(sl) < 30:
            print("  (fewer than 30 tests -- not evaluated)")
            continue
        print(f"  {'variant':<40}{'hit .05':>8}{'(rand)':>8}{'fold':>6}{'Fisher p':>10}{'AUC':>7}{'sens@95%':>10}{'prec top10%':>13}{'prec top25%':>13}")
        for col, lab in labels.items():
            pr, pn = sl[col].to_numpy(), rand[col].to_numpy()
            h, rh = int((pr < 0.05).sum()), int((pn < 0.05).sum())
            fp = fisher_exact([[h, len(pr) - h], [rh, len(pn) - rh]], alternative='greater')[1]
            pre = np.mean([_precision_top(pr, pn[rng_.choice(len(pn), size=len(pr), replace=False)], (0.10, 0.25), rng_) for _ in range(200)], axis=0)
            print(f"  {lab:<40}{h / len(pr):>8.3f}{rh / len(pn):>8.3f}{(h / len(pr)) / (rh / len(pn)) if rh else float('nan'):>6.2f}{fp:>10.4f}"
                  f"{_auc_low_p_is_real(pr, pn):>7.3f}{_sens_at_spec(pr, pn, 0.95)[0]:>10.3f}{pre[0]:>13.3f}{pre[1]:>13.3f}")
    prim = df[slices[primary_name]]
    if len(prim) >= 30:
        bt = {c: {'auc': [], 'sens': [], 'prec': []} for c in labels}
        sampler = make_cluster_sampler(prim['pair_id'].to_numpy()) if 'pair_id' in prim.columns else None   # whole pairs, both directions
        print(f"\n(Bootstrap resamples whole pairs: {len(prim)} directional tests from {prim['pair_id'].nunique() if sampler else '?'} unique pairs.)")
        for _ in range(n_boot):
            ri = sampler(rng_) if sampler else rng_.integers(0, len(prim), len(prim))
            qi = rng_.integers(0, len(rand), len(rand))
            q2 = rng_.choice(len(rand), size=len(ri), replace=False)
            for c in labels:
                pr, pn = prim[c].to_numpy()[ri], rand[c].to_numpy()
                bt[c]['auc'].append(_auc_low_p_is_real(pr, pn[qi]))
                bt[c]['sens'].append(_sens_at_spec(pr, pn[qi], 0.95)[0])
                bt[c]['prec'].append(_precision_top(pr, pn[q2], (0.10,), rng_)[0])
        print(f"\nPRIMARY slice ({primary_name}): difference from '{labels[base_col]}' (bootstrap 95% CI)")
        print(f"  {'variant':<40}{'AUC diff':>28}{'sens@95% diff':>28}{'precision top10% diff':>28}")
        for c, lab in labels.items():
            if c == base_col:
                continue
            cells = []
            for k in ('auc', 'sens', 'prec'):
                d = np.array(bt[c][k]) - np.array(bt[base_col][k])
                lo, hi = np.percentile(d, [2.5, 97.5])
                cells.append(f"{d.mean():+.3f} [{lo:+.3f}, {hi:+.3f}]")
            print(f"  {lab:<40}{cells[0]:>28}{cells[1]:>28}{cells[2]:>28}")
        print("(A CI that includes 0 means the variants performed similarly.)")

if Path(VARIANT_PATH).exists():
    print(f"Reusing cached {VARIANT_PATH} -- delete it to force a fresh run.")
    _vraw = pd.read_csv(VARIANT_PATH)
else:
    _vt = biogrid_methods_combined[['target', 'readout', 'is_real']].drop_duplicates().reset_index(drop=True)
    print(f"{len(_vt)} tests ({int(_vt['is_real'].sum())} real).")
    _vraw = cohort_variant_pvalues(_vt, abundance, expression, dependency, n_permutations=N_PERM_VARIANTS)
    _vraw.to_csv(VARIANT_PATH, index=False)
biogrid_variants = biogrid_methods_combined.merge(_vraw[['target', 'readout', 'is_real', 'q_med20_strat', 'q_central40_strat', 'q_top50_strat', 'q_corr_strat']],
                                                  on=['target', 'readout', 'is_real'], how='left')
_vb = len(biogrid_variants)
biogrid_variants = biogrid_variants.dropna(subset=list(VARIANT_LABELS)).reset_index(drop=True)
print(f"Kept {len(biogrid_variants)} of {_vb} tests where every variant could be computed.")
_vc, _vt_ = biogrid_variants['category'], biogrid_variants['throughput']
_V_SLICES = {'Synthetic Lethality, Low Throughput [primary]': (_vc == 'Synthetic Lethality') & (_vt_ == 'Low Throughput'),
             'Synthetic Lethality, High Throughput': (_vc == 'Synthetic Lethality') & (_vt_ == 'High Throughput'),
             'Negative Genetic, High Throughput': (_vc == 'Negative Genetic') & (_vt_ == 'High Throughput'),
             'All real pairs': biogrid_variants['is_real'].astype(bool)}
print("\n" + "=" * 130 + "\nBioGRID: reference-group width and tissue-stratified null (protein layer)\n" + "=" * 130)
evaluate_columns(biogrid_variants, _V_SLICES, VARIANT_LABELS, 'p_protein_cohort', 'Synthetic Lethality, Low Throughput [primary]', n_boot=N_BOOT_METHODS)

# %% [CELL 11e of 17 -- optional] Readout-matched random controls; SynLethDB high-confidence slice
# #####################################################################
# CELL 11e (optional): do the BioGRID / SynLethDB results survive READOUT-MATCHED random controls?
#   FIGURE OUTPUT: (no figure; supplementary table / Results text)
#   DATA / DEPENDS ON: needs Cells 11b and 11c run first (it re-uses biogrid_methods_combined, synlethdb_methods, _partners and
#                      the helper functions); caches to MATCHED_PATH
# #####################################################################
# DESIGN (fixed before running -- no further slices, controls or metrics after seeing the table):
#   Why: the original random controls are matched to the real test on the TARGET only. The AUCs below 0.5 on SynLethDB pairs suggest
#   the real READOUTS differ from random readouts (e.g. broadly essential genes with little dependency variance). So:
#   (1) Print the dependency mean / SD of the real readouts vs the original random readouts.
#   (2) Matched controls: for each real test (A -> B), draw MATCH_PER_TEST random readouts for the SAME target A from genes whose
#       dependency mean (8 quantile bins) and SD (6 quantile bins within each mean bin) match B, excluding A, B and every known partner of A.
#   (3) Slices: BioGRID Synthetic Lethality / Low Throughput [primary]; SynLethDB CRISPR/CRISPRi [pre-specified primary];
#       SynLethDB HIGH-CONFIDENCE [POST HOC -- added after the CRISPR result was seen]. High-confidence = a pair supported by
#       >= 2 distinct PubMed IDs among experimental rows, or by the 'Low Throughput' source. Multi-source rows ('A;B') are
#       split on ';' and kept if ANY source is experimental (this recovers rows the exact-match filter in Cell 11c dropped).
#   (4) Report BOTH control types for every slice: fold, ROC-AUC, and the pair-clustered hit-rate excess over the matched controls.
MATCH_PER_TEST = 5
MATCHED_PATH = str(VAL_DIR / 'matched_controls_methods.csv')          # cache -- KEEP this file
for _need_name in ('biogrid_methods_combined', 'synlethdb_methods', '_partners', '_sall', '_CRISPR', '_RNAI', '_OTHER', 'make_cluster_sampler'):
    if _need_name not in globals():
        raise RuntimeError(f"Cell 11e needs '{_need_name}' -- run Cells 11b and 11c first.")

_common = [c for c in abundance.columns if c in expression.columns and c in dependency.columns]
_Dm = dependency[_common]
_Dm = _Dm[~_Dm.index.duplicated()]
_gdf = pd.DataFrame({'m': _Dm.mean(axis=1), 's': _Dm.std(axis=1)}).dropna()
_gdf['mb'] = pd.qcut(_gdf['m'], 8, labels=False, duplicates='drop')
_gdf['sb'] = _gdf.groupby('mb')['s'].transform(lambda x: pd.qcut(x, 6, labels=False, duplicates='drop'))
_cell_members = {k: np.array(list(v)) for k, v in _gdf.groupby(['mb', 'sb']).groups.items()}

# ---- real tests for the three slices -------------------------------------------------------------------------------------------
_bg_real = biogrid_methods_combined[(biogrid_methods_combined['category'] == 'Synthetic Lethality') & (biogrid_methods_combined['throughput'] == 'Low Throughput')].copy()
_bg_ctrl_orig = biogrid_methods_combined[biogrid_methods_combined['category'] == 'random control'].copy()
_sl_real = synlethdb_methods[synlethdb_methods['is_real'].astype(bool) & synlethdb_methods['crispr'].astype(bool)].copy()
_sl_ctrl_orig = synlethdb_methods[~synlethdb_methods['is_real'].astype(bool)].copy()

# ---- SynLethDB high-confidence pairs -------------------------------------------------------------------------------------------
_x = _sall.copy()
_x['srcs'] = _x['rel_source'].fillna('').str.split(r'[;|]')
_exp_set = _CRISPR | _RNAI | _OTHER
_x['exp'] = _x['srcs'].apply(lambda L: any(t.strip() in _exp_set for t in L))
_x['lt'] = _x['srcs'].apply(lambda L: any(t.strip() == 'Low Throughput' for t in L))
_x = _x[_x['exp']].copy()
_x['pair_id'] = ['|'.join(sorted((a, b))) for a, b in zip(_x['gene_a'], _x['gene_b'])]
_pairs_hc = _x.groupby('pair_id').agg(n_pmid=('pubmed_id', lambda v: v.dropna().nunique()), lt=('lt', 'max'), a=('gene_a', 'first'), b=('gene_b', 'first')).reset_index()
_pairs_hc['hc'] = (_pairs_hc['n_pmid'] >= 2) | _pairs_hc['lt']
_ab_set, _dp_set = set(abundance.index), set(dependency.index)
_hc_rows = []
for _, r in _pairs_hc[_pairs_hc['hc']].iterrows():
    for t, u in ((r['a'], r['b']), (r['b'], r['a'])):
        if t in _ab_set and u in _dp_set:
            _hc_rows.append({'target': t, 'readout': u, 'is_real': True, 'pair_id': r['pair_id']})
_hc_real = pd.DataFrame(_hc_rows)
_bg_all_ids = {'|'.join(sorted(t)) for t in zip(biogrid_pairs['gene_a'], biogrid_pairs['gene_b'])}
_bg_lt_ids = {'|'.join(sorted(t)) for t in zip(biogrid_pairs.loc[(biogrid_pairs['category'] == 'Synthetic Lethality') & (biogrid_pairs['throughput'] == 'Low Throughput'), 'gene_a'],
                                              biogrid_pairs.loc[(biogrid_pairs['category'] == 'Synthetic Lethality') & (biogrid_pairs['throughput'] == 'Low Throughput'), 'gene_b'])}
_in_bg_map = {pid: (pid in _bg_all_ids) for pid in _pairs_hc['pair_id']}
if len(_hc_real):
    _hc_real['in_bg'] = _hc_real['pair_id'].map(_in_bg_map)
    _hc_ids = set(_hc_real['pair_id'])
    print(f"Overlap of the testable high-confidence pairs with BioGRID: {sum(_in_bg_map[i] for i in _hc_ids)} of {len(_hc_ids)} are in BioGRID at all; "
          f"{len(_hc_ids & _bg_lt_ids)} are in the BioGRID Synthetic Lethality / Low Throughput slice. "
          f"Only the pairs NOT in BioGRID are an independent check of the BioGRID result.")
print(f"SynLethDB experimental pairs (multi-source rows split): {len(_pairs_hc)}; high-confidence: {int(_pairs_hc['hc'].sum())} "
      f"({int(_pairs_hc['lt'].sum())} via 'Low Throughput', {int((_pairs_hc['n_pmid'] >= 2).sum())} via >=2 PubMed IDs); "
      f"testable: {_hc_real['pair_id'].nunique() if len(_hc_real) else 0} pairs, {len(_hc_real)} directional tests.")

# ---- (1) do real readouts differ from the original random readouts? -------------------------------------------------------------
def _readout_profile(df, label):
    g = _gdf.reindex(df['readout'].unique())
    d = df[['readout']].drop_duplicates()
    g = _gdf.reindex(d['readout']).dropna()
    print(f"  {label:<46}n={len(g):>6}  median dependency mean={g['m'].median():+.3f}  median dependency SD={g['s'].median():.3f}")
print("\nDependency profile of the READOUT genes (mean / SD across cell lines):")
_readout_profile(_bg_real, 'BioGRID SL-LT real readouts')
_readout_profile(_bg_ctrl_orig, 'BioGRID original random readouts')
_readout_profile(_sl_real, 'SynLethDB CRISPR real readouts')
_readout_profile(_sl_ctrl_orig, 'SynLethDB original random readouts')
if len(_hc_real):
    _readout_profile(_hc_real, 'SynLethDB high-confidence real readouts')

# ---- (2) matched controls ------------------------------------------------------------------------------------------------------
_mrng = np.random.default_rng(101)
def _matched_controls(real_df, tag):
    rows = []
    for t, r in real_df[['target', 'readout']].drop_duplicates().itertuples(index=False, name=None):
        if r not in _gdf.index:
            continue
        key = (int(_gdf.at[r, 'mb']), int(_gdf.at[r, 'sb']))
        cand = _cell_members.get(key, np.array([]))
        excl = _partners.get(t, set()) | {t, r}
        cand = cand[~np.isin(cand, list(excl))]
        if len(cand) == 0:
            continue
        for g in _mrng.choice(cand, size=min(MATCH_PER_TEST, len(cand)), replace=False):
            rows.append({'tag': tag, 'target': t, 'readout': g, 'is_real': False, 'for_readout': r})
    return pd.DataFrame(rows)
_ctrl_specs = pd.concat([_matched_controls(_bg_real, 'bg'), _matched_controls(_sl_real, 'sl'),
                         _matched_controls(_hc_real, 'hc') if len(_hc_real) else pd.DataFrame()], ignore_index=True)
_calc = pd.concat([_ctrl_specs[['target', 'readout', 'is_real']], _hc_real[['target', 'readout', 'is_real']] if len(_hc_real) else pd.DataFrame()],
                  ignore_index=True).drop_duplicates().reset_index(drop=True)
if Path(MATCHED_PATH).exists():
    print(f"Reusing cached {MATCHED_PATH} -- delete it to force a fresh run.")
    _mp = pd.read_csv(MATCHED_PATH)
else:
    print(f"{len(_calc)} tests to compute ({len(_ctrl_specs)} matched controls across slices).")
    _mp = biogrid_method_pvalues(_calc, abundance, expression, dependency, n_permutations=N_PERM_METHODS)
    _mp.to_csv(MATCHED_PATH, index=False)
_mp = _mp.dropna(subset=list(METHOD_LABELS)).drop_duplicates(subset=['target', 'readout', 'is_real'])
_ctrl_p = _ctrl_specs.merge(_mp, on=['target', 'readout', 'is_real'], how='inner')
_hc_p = _hc_real.merge(_mp, on=['target', 'readout', 'is_real'], how='inner') if len(_hc_real) else pd.DataFrame()
print(f"Matched controls with p-values: {len(_ctrl_p)} of {len(_ctrl_specs)}.")
print("Matched-control readout profile:")
for _tag, _lab in (('bg', 'BioGRID SL-LT matched'), ('sl', 'SynLethDB CRISPR matched'), ('hc', 'SynLethDB high-confidence matched')):
    _sub = _ctrl_p[_ctrl_p['tag'] == _tag]
    if len(_sub):
        _readout_profile(_sub, _lab)

# ---- (3) evaluation with both control types ------------------------------------------------------------------------------------
def evaluate_two_controls(name, real, ctrl_orig, ctrl_matched, n_boot=500, seed=77):
    rng_ = np.random.default_rng(seed)
    print(f"\n{name}\n  real: {len(real)} directional tests from {real['pair_id'].nunique()} pairs; original controls {len(ctrl_orig)}; matched controls {len(ctrl_matched)}")
    if len(real) < 30 or len(ctrl_matched) < 30:
        print("  (too few tests -- not evaluated)")
        return []
    samp = make_cluster_sampler(real['pair_id'].to_numpy())
    print(f"  {'method':<32}{'hit':>7}{'orig ctrl':>10}{'fold':>6}{'AUC':>7}   {'matched ctrl':>12}{'fold':>6}{'AUC':>7}  {'excess vs matched (pair-clustered 95% CI)':>44}{'p':>8}")
    out = []
    boots = {c: [] for c in METHOD_LABELS}
    for _ in range(n_boot):
        ri, qi = samp(rng_), rng_.integers(0, len(ctrl_matched), len(ctrl_matched))
        for c in METHOD_LABELS:
            boots[c].append((real[c].to_numpy()[ri] < 0.05).mean() - (ctrl_matched[c].to_numpy()[qi] < 0.05).mean())
    for c, lab in METHOD_LABELS.items():
        pr, po, pm = real[c].to_numpy(), ctrl_orig[c].to_numpy(), ctrl_matched[c].to_numpy()
        hit, ro, rm = (pr < 0.05).mean(), (po < 0.05).mean(), (pm < 0.05).mean()
        b = np.array(boots[c]); lo, hi = np.percentile(b, [2.5, 97.5])
        pcl = (np.sum(b <= 0) + 1) / (len(b) + 1)
        print(f"  {lab:<32}{hit:>7.3f}{ro:>10.3f}{hit / ro if ro else float('nan'):>6.2f}{_auc_low_p_is_real(pr, po):>7.3f}   {rm:>12.3f}{hit / rm if rm else float('nan'):>6.2f}"
              f"{_auc_low_p_is_real(pr, pm):>7.3f}  {hit - rm:>+14.3f} [{lo:+.3f}, {hi:+.3f}]{'':>9}{pcl:>8.4f}")
        out.append({'slice': name, 'method': lab, 'hit_rate': hit, 'orig_rand': ro, 'matched_rand': rm, 'auc_orig': _auc_low_p_is_real(pr, po),
                    'auc_matched': _auc_low_p_is_real(pr, pm), 'excess_matched': hit - rm, 'ci_lo': lo, 'ci_hi': hi, 'p_clustered': pcl})
    return out

print("\n" + "=" * 150 + "\nReadout-matched controls vs original controls (protein / RNA, cohort / correlation)\n" + "=" * 150)
_ms = []
_ms += evaluate_two_controls('BioGRID Synthetic Lethality, Low Throughput [primary]', _bg_real, _bg_ctrl_orig, _ctrl_p[_ctrl_p['tag'] == 'bg'])
_ms += evaluate_two_controls('SynLethDB CRISPR/CRISPRi [pre-specified primary]', _sl_real, _sl_ctrl_orig, _ctrl_p[_ctrl_p['tag'] == 'sl'])
if len(_hc_p):
    _ms += evaluate_two_controls('SynLethDB high-confidence [POST HOC: added after the CRISPR result was seen]', _hc_p, _sl_ctrl_orig, _ctrl_p[_ctrl_p['tag'] == 'hc'])
    _hc_p['in_bg'] = _hc_p['pair_id'].map(_in_bg_map)
    _ctrl_hc = _ctrl_p[_ctrl_p['tag'] == 'hc'].merge(_hc_real[['target', 'readout', 'pair_id']].rename(columns={'readout': 'for_readout'}), on=['target', 'for_readout'], how='left')
    _ctrl_hc['in_bg'] = _ctrl_hc['pair_id'].map(_in_bg_map)
    _ms += evaluate_two_controls('SynLethDB high-confidence, NOT in BioGRID [POST HOC; the independent subset]', _hc_p[~_hc_p['in_bg']], _sl_ctrl_orig, _ctrl_hc[~_ctrl_hc['in_bg'].fillna(True).astype(bool)])
matched_summary = pd.DataFrame(_ms)
print("\n(Chance AUC = 0.5. Matched controls remove the readout-property difference between real and random tests; an AUC that moves up toward 0.5 or above means the original controls were biased.)")

# %% [CELL 11f of 17 -- optional] What else could the BioGRID protein-layer signal be? Alternative-mechanism diagnostics
# #####################################################################
# CELL 11f (optional): quantify alternative explanations for the BioGRID SL-LT protein-cohort enrichment
#   FIGURE OUTPUT: (no figure; supplementary table / Discussion text)
#   DATA / DEPENDS ON: needs Cell 11e run first (re-uses _bg_real, _ctrl_p, make_cluster_sampler, _auc_low_p_is_real)
# #####################################################################
# DESIGN (fixed before running; every covariate and analysis is listed here, none is added after seeing the output):
#   Question: is the enrichment in the protein cohort test (primary slice, readout-matched controls) carried by something other than
#   "low target -> readout more essential because of synthetic lethality"? For every real and matched-control test (target A, readout B)
#   compute, over the shared cell lines (Spearman):
#     codep      |rho(dependency of A, dependency of B)|      co-essentiality: same complex / pathway (DEPENDENCY data)
#     abund_cov  |rho(abundance of A, abundance of B)|        protein co-variation: co-complex / co-regulated proteins (ABUNDANCE data)
#     rna_cov    |rho(expression of A, expression of B)|      co-expression / shared regulation
#     t_fitness  |rho(abundance of A, line mean dependency)|  target tracks a global fitness / proliferation axis (DEPENDENCY data)
#     t_pc1      |rho(abundance of A, PC1 of dependency)|     target tracks the dominant dependency axis (DEPENDENCY data)
#     r_pc1      |rho(dependency of B, PC1 of dependency)|    readout loads on the dominant dependency axis (DEPENDENCY data)
#     r_tissue   eta-squared of readout dependency by tissue  readout is lineage-driven (DEPENDENCY data)
#     same_stem  1 if the gene symbols share a family stem (e.g. SMARCA2 / SMARCA4): a crude paralog / family flag
#   Analyses: (1) real vs matched controls on each covariate (AUC, 0.5 = no difference); (2) logistic regression of real-vs-control on the
#   SimKO score (-log10 p) alone and with all covariates, coefficient on the score with a pair-clustered bootstrap CI; (3) the hit-rate
#   excess over controls after removing the top third of codep, abund_cov, rna_cov and every same_stem test.
import re
from scipy.stats import rankdata
from sklearn.linear_model import LogisticRegression
for _need_name in ('_bg_real', '_ctrl_p', 'make_cluster_sampler', '_auc_low_p_is_real'):
    if _need_name not in globals():
        raise RuntimeError(f"Cell 11f needs '{_need_name}' -- run Cells 11b to 11e first.")

_c11 = [c for c in abundance.columns if c in expression.columns and c in dependency.columns]
def _clean(df):
    return df.loc[~df.index.duplicated(), _c11]
_A11, _E11, _D11 = _clean(abundance), _clean(expression), _clean(dependency)
_Dv = _D11.to_numpy(float)
_line_fit = np.nanmean(_Dv, axis=0)
_Dfill = np.where(np.isnan(_Dv), np.nanmean(_Dv, axis=1, keepdims=True), _Dv)
_Dc = _Dfill - _Dfill.mean(axis=1, keepdims=True)
_pc1 = np.linalg.svd(_Dc, full_matrices=False)[2][0]              # line scores on the dominant dependency axis
_tissue = np.array([c.split('_', 1)[1] if '_' in c else 'NA' for c in _c11])
_tcodes = pd.factorize(_tissue)[0]

def _sp(x, y):
    m = ~(np.isnan(x) | np.isnan(y))
    if m.sum() < 20:
        return np.nan
    rx, ry = rankdata(x[m]), rankdata(y[m])
    if rx.std() == 0 or ry.std() == 0:
        return np.nan
    return float(np.corrcoef(rx, ry)[0, 1])
def _eta2(y):
    m = ~np.isnan(y)
    if m.sum() < 20:
        return np.nan
    yy, tc = y[m], _tcodes[m]
    tot = ((yy - yy.mean()) ** 2).sum()
    if tot == 0:
        return np.nan
    return float(sum(((yy[tc == k] - yy.mean()) ** 2).sum() * 0 + (tc == k).sum() * (yy[tc == k].mean() - yy.mean()) ** 2 for k in np.unique(tc)) / tot)
def _stem(g):
    return re.sub(r'[0-9]+[A-Z]?$', '', str(g))

_apos, _epos, _dpos = ({g: i for i, g in enumerate(df.index)} for df in (_A11, _E11, _D11))
_COVS = ['codep', 'abund_cov', 'rna_cov', 't_fitness', 't_pc1', 'r_pc1', 'r_tissue', 'same_stem']
_target_cache = {}
def _target_covs(t):
    if t not in _target_cache:
        a = _A11.to_numpy(float)[_apos[t]] if t in _apos else None
        _target_cache[t] = (a, abs(_sp(a, _line_fit)) if a is not None else np.nan, abs(_sp(a, _pc1)) if a is not None else np.nan)
    return _target_cache[t]
_readout_cache = {}
def _readout_covs(r):
    if r not in _readout_cache:
        y = _Dv[_dpos[r]] if r in _dpos else None
        _readout_cache[r] = (y, abs(_sp(y, _pc1)) if y is not None else np.nan, _eta2(y) if y is not None else np.nan)
    return _readout_cache[r]
def _covariates(df):
    rows = []
    for t, r in zip(df['target'], df['readout']):
        a, tf, tp = _target_covs(t)
        y, rp, et = _readout_covs(r)
        yt = _Dv[_dpos[t]] if t in _dpos else None
        rows.append({'codep': abs(_sp(yt, y)) if (yt is not None and y is not None) else np.nan,
                     'abund_cov': abs(_sp(a, _A11.to_numpy(float)[_apos[r]])) if (a is not None and r in _apos) else np.nan,
                     'rna_cov': abs(_sp(_E11.to_numpy(float)[_epos[t]], _E11.to_numpy(float)[_epos[r]])) if (t in _epos and r in _epos) else np.nan,
                     't_fitness': tf, 't_pc1': tp, 'r_pc1': rp, 'r_tissue': et,
                     'same_stem': float(len(_stem(t)) >= 3 and _stem(t) == _stem(r))})
    return pd.DataFrame(rows, index=df.index)

_ctrl_bg = _ctrl_p[_ctrl_p['tag'] == 'bg'].copy()
_real11 = pd.concat([_bg_real.reset_index(drop=True), _covariates(_bg_real.reset_index(drop=True))], axis=1)
_ctrl11 = pd.concat([_ctrl_bg.reset_index(drop=True), _covariates(_ctrl_bg.reset_index(drop=True))], axis=1)
_ctrl11['pair_id'] = ['ctrl%d' % i for i in range(len(_ctrl11))]
_all11 = pd.concat([_real11.assign(y=1), _ctrl11.assign(y=0)], ignore_index=True)
_all11['z'] = -np.log10(_all11['p_protein_cohort'].clip(lower=1 / (N_PERM_METHODS + 1)))
print(f"\nPrimary slice: {len(_real11)} real tests ({_real11['pair_id'].nunique()} pairs) vs {len(_ctrl11)} readout-matched controls.")

# ---- (1) real vs matched controls on each covariate ---------------------------------------------------------------------------
print("\n(1) Do real tests differ from matched controls on alternative-mechanism covariates?  AUC 0.5 = no difference; >0.5 = real tests higher.")
print(f"  {'covariate':<12}{'median real':>13}{'median ctrl':>13}{'AUC':>7}   {'95% CI (pairs resampled)':>26}")
_rng11 = np.random.default_rng(11)
_samp11 = make_cluster_sampler(_real11['pair_id'].to_numpy())
_cov_rows = []
for c in _COVS:
    xr, xc = _real11[c].to_numpy(float), _ctrl11[c].to_numpy(float)
    def _auc(a, b):
        a, b = a[~np.isnan(a)], b[~np.isnan(b)]
        if len(a) < 10 or len(b) < 10:
            return np.nan
        r = rankdata(np.concatenate([a, b]))
        return (r[:len(a)].sum() - len(a) * (len(a) + 1) / 2) / (len(a) * len(b))
    au = _auc(xr, xc)
    bs = [_auc(xr[_samp11(_rng11)], xc[_rng11.integers(0, len(xc), len(xc))]) for _ in range(300)]
    lo, hi = np.nanpercentile(bs, [2.5, 97.5])
    print(f"  {c:<12}{np.nanmedian(xr):>13.3f}{np.nanmedian(xc):>13.3f}{au:>7.3f}   [{lo:.3f}, {hi:.3f}]")
    _cov_rows.append({'covariate': c, 'median_real': np.nanmedian(xr), 'median_ctrl': np.nanmedian(xc), 'auc': au, 'ci_lo': lo, 'ci_hi': hi})
alt_mechanism_covariates = pd.DataFrame(_cov_rows)

# ---- (2) does the SimKO score still separate real from control once the covariates are in the model? ---------------------------
def _fit_coef(df, use_cov):
    X = df[['z'] + (_COVS if use_cov else [])].copy()
    for c in X.columns:
        X[c] = X[c].fillna(X[c].median())
        sd = X[c].std()
        X[c] = (X[c] - X[c].mean()) / (sd if sd > 0 else 1)
    m = LogisticRegression(C=1e3, class_weight='balanced', max_iter=2000).fit(X.to_numpy(), df['y'].to_numpy())
    return float(m.coef_[0][0])
print("\n(2) Logistic regression, real vs matched control, coefficient on the SimKO score (-log10 p, standardised; per SD).")
_b0, _b1 = _fit_coef(_all11, False), _fit_coef(_all11, True)
_boot0, _boot1 = [], []
_ctrl_part = _all11[_all11['y'] == 0].reset_index(drop=True)
_real_part = _all11[_all11['y'] == 1].reset_index(drop=True)
for _ in range(200):
    d_ = pd.concat([_real_part.iloc[_samp11(_rng11)], _ctrl_part.iloc[_rng11.integers(0, len(_ctrl_part), len(_ctrl_part))]], ignore_index=True)
    _boot0.append(_fit_coef(d_, False)); _boot1.append(_fit_coef(d_, True))
l0, h0 = np.percentile(_boot0, [2.5, 97.5]); l1, h1 = np.percentile(_boot1, [2.5, 97.5])
print(f"  score alone                 coefficient {_b0:+.3f}  95% CI [{l0:+.3f}, {h0:+.3f}]")
print(f"  score + all 8 covariates    coefficient {_b1:+.3f}  95% CI [{l1:+.3f}, {h1:+.3f}]")
print(f"  share of the score's coefficient that remains after adjustment: {_b1 / _b0:.0%}" if _b0 else "")
print("  (A coefficient that stays clearly above 0 means the protein-cohort score carries information the covariates do not.)")
_cov_only = _all11.copy()
_cov_only['z'] = 0.0
# which single covariate removes the most signal?
print("  leave-one-in: coefficient on the score with ONE covariate added at a time")
_single = []
for c in _COVS:
    X = _all11[['z', c]].copy()
    for k in X.columns:
        X[k] = X[k].fillna(X[k].median()); sd = X[k].std(); X[k] = (X[k] - X[k].mean()) / (sd if sd > 0 else 1)
    m = LogisticRegression(C=1e3, class_weight='balanced', max_iter=2000).fit(X.to_numpy(), _all11['y'].to_numpy())
    _single.append({'covariate': c, 'score_coef': float(m.coef_[0][0]), 'cov_coef': float(m.coef_[0][1])})
    print(f"    + {c:<11} score coef {m.coef_[0][0]:+.3f}   covariate coef {m.coef_[0][1]:+.3f}")
alt_mechanism_single = pd.DataFrame(_single)

# ---- (3) hit-rate excess after removing tests with obvious shared structure --------------------------------------------------
print("\n(3) Hit-rate excess (p<0.05, protein cohort) over matched controls after removing tests with obvious shared structure.")
_thr = {c: np.nanpercentile(_ctrl11[c].to_numpy(float), 66.7) for c in ('codep', 'abund_cov', 'rna_cov')}
def _clean_mask(df):
    k = np.ones(len(df), bool)
    for c, t in _thr.items():
        v = df[c].to_numpy(float)
        k &= ~(v > t)                     # NaN counts as kept
    k &= df['same_stem'].to_numpy(float) != 1
    return k
_rk, _ck = _clean_mask(_real11), _clean_mask(_ctrl11)
_rc, _cc = _real11[_rk].reset_index(drop=True), _ctrl11[_ck].reset_index(drop=True)
print(f"  kept {len(_rc)} of {len(_real11)} real tests ({_rc['pair_id'].nunique()} pairs) and {len(_cc)} of {len(_ctrl11)} controls "
      f"(thresholds = 66.7th percentile of controls: " + ", ".join(f"{c} {t:.2f}" for c, t in _thr.items()) + ")")
if len(_rc) >= 30 and len(_cc) >= 30:
    sm = make_cluster_sampler(_rc['pair_id'].to_numpy())
    bs = np.array([(_rc['p_protein_cohort'].to_numpy()[sm(_rng11)] < 0.05).mean() - (_cc['p_protein_cohort'].to_numpy()[_rng11.integers(0, len(_cc), len(_cc))] < 0.05).mean() for _ in range(500)])
    hr, hc = (_rc['p_protein_cohort'] < 0.05).mean(), (_cc['p_protein_cohort'] < 0.05).mean()
    print(f"  hit rate {hr:.3f} vs matched controls {hc:.3f}; excess {hr - hc:+.3f}  95% CI [{np.percentile(bs, 2.5):+.3f}, {np.percentile(bs, 97.5):+.3f}]  "
          f"one-sided pair-clustered p={(np.sum(bs <= 0) + 1) / (len(bs) + 1):.4f}   AUC={_auc_low_p_is_real(_rc['p_protein_cohort'].to_numpy(), _cc['p_protein_cohort'].to_numpy()):.3f}")
else:
    print("  (too few tests left to evaluate)")
print("\n(Interpretation guide: the primary slice excess before this filter was +0.053 [+0.018, +0.092]. If it stays near that value, shared complex / "
      "co-variation / family structure does not explain it. If it collapses, the signal is mostly that structure -- still useful, but not 'synthetic lethality' as such.)")

# %% [CELL 12 of 17] Exploratory: HSP110-family trigenic redundancy
# #####################################################################
# CELL 12 of 17: Exploratory: HSP110-family trigenic redundancy
#   FIGURE OUTPUT: (no figure; Discussion text only)
#   DATA / DEPENDS ON: none external
# #####################################################################
# --- Exploratory: HSPA4-family (HSPA4, HSPA4L, HSPH1) trigenic redundancy ---
# Real, published finding (in4mer platform, Nat Commun 2024, Gonatopoulos-
# Pournatzis et al.): combined CRISPR knockout of all three HSP110-family
# chaperones gives a severe fitness phenotype in a pooled screen, while any
# single or pairwise knockout gives only a moderate one. SimKO can't directly
# replicate an active combinatorial knockout screen -- there is no single
# "combined fitness" readout for a trio in DepMap's per-gene dependency data
# the way a real triple-knockout screen provides. This instead asks a
# related, weaker question SimKO's own mechanism can actually test: do cell
# lines NATURALLY low in all three simultaneously show a more distinctive
# Dependency-layer signature than cell lines low in just one or two? Cohorts
# are formed on Abundance (natural protein level); the downstream readout is
# Dependency -- a different layer combination from get_simko_differential's
# single-layer design, so cohort formation and the differential/permutation
# steps are wired separately here rather than reusing that function.
HSP_GENES = ['HSPA4', 'HSPA4L', 'HSPH1']
covered = [g for g in HSP_GENES if g in abundance.index]
print(f"{len(covered)} of {len(HSP_GENES)} HSP110-family genes covered in the Abundance matrix: {covered}")

if len(covered) < 2:
    print("Fewer than 2 covered -- cannot test pairwise/triple cohorts.")
else:
    from itertools import combinations

    # Abundance (cohort-forming) and Dependency (readout) do not necessarily
    # cover the exact same cell lines (DepMap CRISPR screening and CCLE
    # proteomics panels overlap only partially) -- restrict both to their
    # shared cell lines BEFORE forming cohorts, so every cell line the
    # cohort split produces is guaranteed to also exist in the readout
    # layer. Confirmed necessary directly: skipping this raised an
    # IndexError deep inside permutation shuffling, from pandas silently
    # returning fewer columns than requested rather than erroring loudly.
    shared_cls = [c for c in abundance.columns if c in set(dependency.columns)]
    n_dropped = abundance.shape[1] - len(shared_cls)
    if n_dropped:
        print(f"{n_dropped} of {abundance.shape[1]} Abundance cell lines are absent from Dependency "
              f"-- restricting this analysis to the {len(shared_cls)} cell lines common to both.")
    abundance_hsp = abundance[shared_cls]
    dependency_hsp = dependency[shared_cls]
    n_hsp = compute_cohort_n(len(shared_cls))
    hsp_results = {}

    def _cohort_dependency_fdr(target_list, label):
        print(f"\n--- {label} ---")
        class_df = get_classes_by_joint_value(target_list, abundance_hsp, n_hsp, mode='max')
        low_class = list(class_df.loc[class_df['class'] == 'low'].index)
        median_class = list(class_df.loc[class_df['class'] == 'median'].index)
        observed = get_differentials(class_df, dependency_hsp, len(low_class), len(median_class))
        observed = observed.loc[~observed.index.isin(target_list)]  # exclude targets' own trivial rows
        null_p = permutation_null_pvalues(dependency_hsp, low_class, median_class, n_permutations=200, seed=42)
        fdr_curve = empirical_fdr_curve(observed['p'].values, null_p)
        print(fdr_curve.to_string(index=False))
        return fdr_curve

    for g in covered:
        hsp_results[(g,)] = _cohort_dependency_fdr([g], f"Single: {g}")
    if len(covered) >= 2:
        for pair in combinations(covered, 2):
            hsp_results[pair] = _cohort_dependency_fdr(list(pair), f"Pair: {'+'.join(pair)}")
    if len(covered) == 3:
        hsp_results[tuple(covered)] = _cohort_dependency_fdr(covered, f"Triple: {'+'.join(covered)}")

    print(f"\n{'=' * 70}\nSummary: empirical FDR at p<0.01, by cohort\n{'=' * 70}")
    print("(Lower = more real downstream Dependency-layer signal in that cohort. The trigenic")
    print(" hypothesis predicts the triple cohort should stand out from every single/pairwise one.)")
    for combo, fdr_curve in hsp_results.items():
        row = fdr_curve[fdr_curve['nominal_p_threshold'] == 0.01]
        if len(row) and pd.notna(row.iloc[0]['empirical_fdr']):
            print(f"{'+'.join(combo):<25} empirical FDR={row.iloc[0]['empirical_fdr']:.4f}  ({row.iloc[0]['observed_hits']:.0f} hits)")
        else:
            print(f"{'+'.join(combo):<25} (no hits at this threshold)")

# %% [CELL 13 of 17] Targeted paralog-compensation battery: HSP110 family
# #####################################################################
# CELL 13 of 17: Targeted paralog-compensation battery: HSP110 family
#   FIGURE OUTPUT: (no figure; Discussion text only)
#   DATA / DEPENDS ON: needs Cell 12
# #####################################################################
# --- Targeted paralog-compensation test: HSP110 family (DepMine-style) ---
# DepMine's own most directly comparable finding isn't a triple-family
# fitness test -- it's a targeted, single-readout paralog-compensation
# test: does loss of one SWI/SNF paralog (ARID1A, or SMARCA4) increase
# dependency SPECIFICALLY on its partner paralog (ARID1B, or SMARCA2)?
# They found a strong, and notably ASYMMETRIC, signal there (ARID1A-loss
# -> ARID1B-dependency was strong; the reverse was not). This runs the
# direct analog for the three HSP110-family genes: every directional
# single-target pair, plus each pair-low cohort against the one
# remaining family member -- 9 tests total, each far more targeted (and
# better powered) than the exploratory genome-wide scan in the cell
# above. DepMine's own signal came from actual genetic LoF/deletion, a
# stronger, more complete perturbation than SimKO's continuous natural-
# abundance proxy -- a weaker or null result here isn't necessarily a
# problem with SimKO's method, just the expected cost of a gentler proxy.
HSP_GENES = ['HSPA4', 'HSPA4L', 'HSPH1']
covered = [g for g in HSP_GENES if g in abundance.index]

if len(covered) < 2:
    print("Fewer than 2 covered -- cannot run any directional test.")
else:
    from itertools import permutations, combinations

    shared_cls = [c for c in abundance.columns if c in set(dependency.columns)]
    n_dropped = abundance.shape[1] - len(shared_cls)
    if n_dropped:
        print(f"{n_dropped} of {abundance.shape[1]} Abundance cell lines are absent from Dependency "
              f"-- restricting this analysis to the {len(shared_cls)} cell lines common to both.")
    abundance_hsp = abundance[shared_cls]
    dependency_hsp = dependency[shared_cls]

    compensation_results = {}

    # Single-target directional tests: low-X -> dependency on Y, every ordered pair
    for target, readout in permutations(covered, 2):
        compensation_results[(target, readout)] = test_single_gene_dependency_shift(
            target, readout, abundance_hsp, dependency_hsp, n_permutations=200
        )

    # Pair-low -> dependency on the one remaining family member (only
    # meaningful with all 3 covered)
    if len(covered) == 3:
        for pair in combinations(covered, 2):
            remaining = [g for g in covered if g not in pair][0]
            compensation_results[(tuple(pair), remaining)] = test_single_gene_dependency_shift(
                list(pair), remaining, abundance_hsp, dependency_hsp, n_permutations=200
            )

    print(f"\n{'=' * 70}\nSummary: HSP110-family paralog-compensation tests\n{'=' * 70}")
    print("(One-sided empirical p: does the low cohort show HIGHER dependency on")
    print(" the readout gene, as DepMine's ARID1A->ARID1B pattern predicts? Note")
    print(" whether any real signal is similarly asymmetric by direction.)")
    for (target, readout), r in compensation_results.items():
        target_label = target if isinstance(target, str) else '+'.join(target)
        print(f"{target_label:<20} -> {readout:<10} diff={r['observed_diff']:+.4f}  "
              f"one-sided p={r['empirical_p_one_sided']:.4f}")

# %% [CELL 14 of 17] Family-wise joint permutation across the 9 HSP110 tests
# #####################################################################
# CELL 14 of 17: Family-wise joint permutation across the 9 HSP110 tests
#   FIGURE OUTPUT: (no figure; Discussion text only)
#   DATA / DEPENDS ON: needs Cell 13
# #####################################################################
# --- Family-wise significance check: joint permutation across all 9 tests ---
# The battery above computed 9 marginal p-values independently -- treating
# them as independent (e.g. via Benjamini-Hochberg) ignores that several
# share a cohort definition or a readout gene, and are therefore
# correlated. This reruns the same 9 tests through ONE joint permutation
# procedure that preserves that correlation, then asks the real question
# directly: how often would a null with this exact correlation structure
# produce as many "hits" as we actually saw, purely by chance?
if len(covered) == 3:
    from itertools import permutations as _permutations, combinations as _combinations
    joint_tests = [((t,), r) for t, r in _permutations(covered, 2)]
    for pair in _combinations(covered, 2):
        remaining = [g for g in covered if g not in pair][0]
        joint_tests.append((tuple(pair), remaining))

    joint_result = joint_permutation_dependency_battery(
        joint_tests, abundance, dependency, n_permutations=1000, seed=42
    )

# %% [CELL 15 of 17] ME2 -> ME3 collateral lethality replication
# #####################################################################
# CELL 15 of 17: ME2 -> ME3 collateral lethality replication
#   FIGURE OUTPUT: (no figure; Discussion text only)
#   DATA / DEPENDS ON: none external
# #####################################################################
# --- Further single-target replication: ME2 -> ME3 collateral lethality ---
# Real, established relationship (Dey et al. 2017, Nature): pancreatic
# cancer cells with genomic loss of ME2 (often collateral to SMAD4
# deletion) become specifically dependent on the paralogous ME3 for
# survival, via impaired NADPH/ROS homeostasis. ME1, the third malic
# enzyme paralog, was explicitly shown NOT to participate the same way --
# a built-in negative control for the same test. Unlike the HSPA4-family
# test above, this is a single-target relationship, so it reuses the
# exact mechanism already validated for BRCA1/BRCA2 -> PARP1 directly,
# with no new code.
ME_TARGETS = ['ME1', 'ME2']
me_covered = [g for g in ME_TARGETS if g in abundance.index]
if 'ME3' not in dependency.index or not me_covered:
    print("ME3 (readout) or neither of ME1/ME2 (targets) covered -- cannot run this test.")
else:
    me_results = {}
    for target in me_covered:
        me_results[target] = test_single_gene_dependency_shift(target, 'ME3', abundance, dependency, n_permutations=200)

    pancreas_cls = get_tissue_cell_lines(abundance.columns, 'PANCREAS')
    print(f"\n{len(pancreas_cls)} PANCREAS cell lines available for tissue-restricted comparison")
    me_tissue_results = {}
    if len(pancreas_cls) >= 2:
        for target in me_covered:
            me_tissue_results[target] = test_single_gene_dependency_shift(
                target, 'ME3', abundance, dependency, active_cell_lines=pancreas_cls, n_permutations=200
            )
    else:
        print("Too few PANCREAS cell lines -- skipping tissue-restricted comparison.")

    print(f"\n{'=' * 70}\nSummary: ME2/ME1 -> ME3 collateral lethality\n{'=' * 70}")
    print("(Real, established relationship: ME2-low predicts ME3 dependency;")
    print(" ME1-low should NOT, per Dey et al. 2017 -- a built-in negative control.)")
    for target, r in me_results.items():
        print(f"Pan-cancer:  {target} -> ME3  diff={r['observed_diff']:+.4f}  one-sided p={r['empirical_p_one_sided']:.4f}")
    for target, r in me_tissue_results.items():
        print(f"PANCREAS:    {target} -> ME3  diff={r['observed_diff']:+.4f}  one-sided p={r['empirical_p_one_sided']:.4f}")

# %% [CELL 16 of 17] Composite deficiency score vs Abundance alone
# #####################################################################
# CELL 16 of 17: Composite deficiency score vs Abundance alone
#   FIGURE OUTPUT: (no figure; Discussion text only)
#   DATA / DEPENDS ON: needs Cell 3
# #####################################################################
# --- Composite 'deficiency score' cohort: does combining Abundance + Expression
# + Mutation beat Abundance alone? ---
# Quick, deliberately lean first look -- not a fully-weighted, tuned composite,
# just equal-weighted z-scores of Abundance, Expression, and Mutation burden for
# ARID1A itself (excluding Dependency: that's a consequence of deficiency, not
# evidence for it, so including it here would be circular against the test
# below). Reuses the exact same real RPE1 ARID1A-KO proteomics ground truth and
# concordance methodology already established above -- this cell must be run
# after that one (needs real_diff_arid1a, result_arid1a, perm_arid1a).
TARGET = 'ARID1A'
composite_df, composite_cls = get_deficiency_composite(TARGET, abundance, expression, mutation)
n_composite = compute_cohort_n(len(composite_cls))

class_df_composite = get_classes_by_joint_value([TARGET], composite_df, n_composite, mode='max')
n_low = (class_df_composite['class'] == 'low').sum()
n_median = (class_df_composite['class'] == 'median').sum()
print(f"\nComposite-score cohort: {n_low} low (deficient), {n_median} median")

abundance_restricted = abundance[composite_cls]
simko_diff_composite = get_differentials(class_df_composite, abundance_restricted, n_low, n_median)
simko_diff_composite = simko_diff_composite.loc[~simko_diff_composite.index.isin([TARGET])]

print(f"\n{'=' * 70}\nComposite-score cohort vs real RPE1 ARID1A-KO proteomics\n{'=' * 70}")
result_composite = run_concordance_analysis(simko_diff_composite, real_diff_arid1a, top_n=200,
                                             output_prefix='concordance_ARID1A_composite')

# permutation_null_for_concordance recomputes its OWN 'observed' from raw Abundance
# ranking internally -- not what we want here, since our real cohort came from the
# composite ranking instead. Its NULL arrays are still valid and reusable exactly as
# generated (the null is a completely random cohort split of the same size, applied
# to the same genome-wide Abundance readout -- it doesn't depend on how the REAL
# cohort was formed), so we take those and compute our own empirical p against our
# own composite-based observed statistics instead of trusting this call's returned
# 'observed_*'/'empirical_p_*' fields directly.
perm_composite_raw = permutation_null_for_concordance(
    [TARGET], abundance_restricted, real_diff_arid1a, n=n_composite, n_permutations=200, seed=42
)
emp_p_rho_composite = float(np.mean(np.abs(perm_composite_raw['null_rhos']) >= abs(result_composite['spearman_rho'])))
emp_p_up_composite = float(np.mean(perm_composite_raw['null_ups'] >= result_composite['up_overlap']))
emp_p_down_composite = float(np.mean(perm_composite_raw['null_downs'] >= result_composite['down_overlap']))

print(f"\n{'=' * 70}\nHead-to-head: Abundance-only cohort vs Composite-score cohort\n{'=' * 70}")
print(f"{'Metric':<20}{'Abundance-only':<20}{'Composite':<20}")
print(f"{'UP overlap':<20}{result_arid1a['up_overlap']:<20}{result_composite['up_overlap']:<20}")
print(f"{'  empirical p':<20}{perm_arid1a['empirical_p_up']:<20.4f}{emp_p_up_composite:<20.4f}")
print(f"{'DOWN overlap':<20}{result_arid1a['down_overlap']:<20}{result_composite['down_overlap']:<20}")
print(f"{'  empirical p':<20}{perm_arid1a['empirical_p_down']:<20.4f}{emp_p_down_composite:<20.4f}")
print(f"{'rho':<20}{result_arid1a['spearman_rho']:<20.4f}{result_composite['spearman_rho']:<20.4f}")
print(f"{'  empirical p':<20}{perm_arid1a['empirical_p_rho']:<20.4f}{emp_p_rho_composite:<20.4f}")


# %% [CELL 17 of 17] ALL MANUSCRIPT FIGURES -- the ONLY cell that writes figure files
# #####################################################################
# CELL 17 of 17: ALL MANUSCRIPT FIGURES -- the ONLY cell that writes figure files
#   FIGURE OUTPUT: writes FIGURES 1-8, 8b and 9 to FIG_DIR
#   DATA / DEPENDS ON: reads results already in memory; skips any figure whose cell has not been run
# #####################################################################
# =====================================================================
# ALL MANUSCRIPT FIGURES -- run this cell last (or any time after the
# cells it draws from). Every figure is regenerated from results already
# in memory, so changing FIG_STYLE / FIG_DPI / the titles below and
# re-running THIS cell is all that's needed; nothing is recomputed.
# A figure whose inputs haven't been computed yet is skipped with a note
# naming the missing variable, rather than stopping the others.
# =====================================================================
INCLUDE_PROTEIN_PANEL_IN_FIG4 = True    # panel A: Abundance layer vs real RPE1 proteomics (Cells 3 and 5)
INCLUDE_ABUNDANCE_VS_RNA_PANEL = True   # panel B: Abundance layer vs real HCT116 RNA-seq (Cell 7)
INCLUDE_EXPRESSION_PANEL = False        # optional extra panel: Expression layer vs real HCT116 RNA-seq (Cell 7). Off by default --
                                        # the Abundance layer is the headline; Expression-vs-Abundance is already Figure 5.

def _scope_labels():
    labels = {'pancancer': 'Pan-cancer (%d lines)' % abundance.shape[1]}
    for tissue in phase0_results:
        if tissue != 'pancancer':
            n_lines = len(get_tissue_cell_lines(abundance.columns, tissue))
            labels[tissue.lower()] = '%s (%d lines)' % (tissue.replace('_', ' ').title(), n_lines)
    return labels

def _fig4():
    candidates = [
        (INCLUDE_PROTEIN_PANEL_IN_FIG4, 'SimKO Abundance (protein) vs real RPE1 knockout proteomics\n[protein vs protein]',
         'perm_arid1a', 'perm_corr_arid1a', 'combined_arid1a'),
        (INCLUDE_ABUNDANCE_VS_RNA_PANEL, 'SimKO Abundance (protein) vs real HCT116 knockout RNA-seq\n[protein vs RNA]',
         'perm_simko_abund_rna', 'perm_corr_abund_rna', 'combined_abund_rna'),
        (INCLUDE_EXPRESSION_PANEL, 'SimKO Expression (RNA) vs real HCT116 knockout RNA-seq\n[RNA vs RNA]',
         'perm_simko_expr', 'perm_corr_expr', 'combined_expr'),
    ]
    chosen = [c for c in candidates if c[0]]
    panels = []
    for k, (_, title, simko_name, corr_name, comb_name) in enumerate(chosen):
        panels.append({'title': '%s. %s' % ('ABC'[k], title),
                       'simko': globals()[simko_name], 'correlation': globals()[corr_name],
                       'combined': globals()[comb_name]})
    plot_benchmark_figure(panels)

FIGURE_JOBS = [
    ('Figure 1', ['phase0_results'],
     lambda: plot_real_target_fdr_curves(
         phase0_results['pancancer']['raw']['real'], FIG_DIR / 'fig1_real_target_fdr_curves.png',
         title='Empirical FDR for four biologically-motivated targets\n(Abundance layer, pan-cancer, %d cell lines)' % abundance.shape[1])),
    ('Figure 2', ['comparison_table', 'phase0_results'],
     lambda: plot_multipanel_multitissue_comparison(
         comparison_table, FIG_DIR / 'fig2_multipanel_multitissue_comparison.png',
         scope_labels=_scope_labels(), title='Real targets vs. non-specific panels, across tissue scopes')),
    ('Figure 3', ['tissue_bias_df', 'target_categories'],
     lambda: plot_tissue_bias_ratio(
         tissue_bias_df, target_categories, FIG_DIR / 'fig3_tissue_bias_ratio_per_gene.png',
         title='Per-gene tissue-bias ratio, pan-cancer (Abundance layer)')),
    ('Figure 4', ['perm_arid1a', 'perm_corr_arid1a', 'combined_arid1a',
                  'perm_simko_abund_rna', 'perm_corr_abund_rna', 'combined_abund_rna'], _fig4),
    ('Figure 5', ['perm_arid1a', 'perm_expr_vs_protein'],
     lambda: plot_cross_layer_grid(perm_arid1a, perm_expr_vs_protein)),
    ('Figure 6', ['abundance', 'dependency'],
     lambda: plot_brca_parp1_boxplot(get_cohort_dependency_values('BRCA1', 'PARP1', abundance, dependency),
                                     get_cohort_dependency_values('BRCA2', 'PARP1', abundance, dependency))),
    ('Figure 7', ['biogrid_combined'], lambda: plot_category_throughput_breakdown(biogrid_combined)),
    ('Figure 8', ['biogrid_combined'], lambda: plot_enrichment_curve(biogrid_combined)),
    ('Figure 9 (BioGRID: protein vs RNA vs correlation)', ['biogrid_method_summary'],
     lambda: plot_biogrid_method_comparison(biogrid_method_summary)),
    ('Figure 8b (supplementary ROC)', ['biogrid_combined'],
     lambda: print("  ROC-AUCs (quoted in the Results text):", {k: round(v, 3) for k, v in plot_roc_curves_secondary(biogrid_combined).items()})),
]

print("Writing figures to %s\n" % FIG_DIR)
for label, needs, draw in FIGURE_JOBS:
    missing = [v for v in needs if v not in globals()]
    if missing:
        print("%s: SKIPPED -- not computed yet, run the cell that creates: %s" % (label, ', '.join(missing)))
        continue
    try:
        draw()
    except Exception as e:
        print("%s: FAILED -- %s: %s" % (label, type(e).__name__, e))