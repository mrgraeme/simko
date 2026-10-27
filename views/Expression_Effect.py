import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from simko_utils import (
    get_abundance_data, get_expression_data, get_mutation_data, get_dependency_data,
    get_classes_by_joint_value, get_differentials, get_mutation_differentials,
    get_diff_summary, render_raw_table, get_differentials_boxplot, select_outliers,
    compute_tissue_variance_explained, style_results_table, get_tissue_variance_column_config,
    get_correlation_scores,
)

abundance = get_abundance_data()
expression = get_expression_data()
mutation = get_mutation_data()
dependency = get_dependency_data()
cmap = plt.cm.get_cmap('RdYlBu_r')

st.write("### Explore expression change effects 🧬")

ko_targets = st.session_state.get('saved_ko_proteins', [])
focus_proteins = st.session_state.get('saved_protein_list', [])
active_tissues = st.session_state.get('saved_tissues', [])
active_cell_lines = st.session_state.get('saved_cell_lines', [])

strict_multi_target = st.session_state.get('saved_strict_multi_target', True)
combo_mode = 'max' if strict_multi_target else 'mean'

correction_scope_choice = st.session_state.get('saved_correction_scope_choice', '')
if correction_scope_choice.startswith('All genome-wide'):
    scope_mode = None
elif correction_scope_choice.startswith('Features with overall'):
    scope_mode = 'overall_range'
else:
    scope_mode = 'diff'
prefilter_ab = st.session_state.get('saved_correction_threshold_ab', 0.25)
prefilter_ex = st.session_state.get('saved_correction_threshold_ex', 0.25)
prefilter_dep = st.session_state.get('saved_correction_threshold_dep', 0.25)

diff_floor = st.session_state.get('saved_diff_floor', 0.25)
sig_filter_enabled = st.session_state.get('saved_sig_filter_enabled', False)
sig_col = st.session_state.get('saved_sig_col', 'p') if sig_filter_enabled else None
sig_threshold = st.session_state.get('saved_sig_threshold', 0.01)
sort_by = st.session_state.get('saved_sort_by', 'diff')

# Cell-line coverage is not guaranteed identical across matrices (see Table
# 1 in the manuscript -- "cell line count to confirm"), so filter against
# the EXPRESSION matrix's own columns rather than assuming they match
# Abundance.
all_cols = list(expression.columns)
if active_cell_lines:
    active_cls = [c for c in active_cell_lines if c in all_cols]
elif active_tissues:
    active_cls = [s for s in all_cols if any(xs in s for xs in active_tissues)]
else:
    active_cls = all_cols

abundance = abundance.filter(active_cls)
expression = expression.filter(active_cls)
mutation = mutation.filter(active_cls)
dependency = dependency.filter(active_cls)

eta_abund = compute_tissue_variance_explained(abundance)
eta_exp = compute_tissue_variance_explained(expression)
eta_dep = compute_tissue_variance_explained(dependency)
eta_mut = compute_tissue_variance_explained(mutation)

# Correlation view alongside the cohort-based differential -- only defined
# for a single target (see get_correlation_scores' docstring for why: a
# "composite correlation" across multiple targets hasn't been validated
# the way the single-target case has, so this is left unset for
# multi-target selections rather than guessing at a definition).
if ko_targets and len(ko_targets) == 1:
    corr_abund = get_correlation_scores(ko_targets[0], abundance)
    corr_exp = get_correlation_scores(ko_targets[0], expression)
    corr_dep = get_correlation_scores(ko_targets[0], dependency)
else:
    corr_abund = corr_exp = corr_dep = None

if ko_targets:
    st.caption(f"Analyzing simulated KO effects of {ko_targets} across {len(active_cls)} cell lines (Expression-based classification).")

    n = ((expression.shape[1] - 1) // 3) if expression.shape[1] < 60 else 20
    class_df = get_classes_by_joint_value(ko_targets, expression, n, mode=combo_mode)

    tab1, tab2 = st.tabs(["Classification Map", "Cohort Breakdowns"])
    with tab1:
        fig = plt.figure(figsize=(10, 3))
        sns.heatmap(
            class_df[['joint_low_score']].sort_values('joint_low_score').T.astype(float).round(2),
            square=True, cmap="vlag", annot=True, annot_kws={'size': 6}, cbar=False
        )
        st.pyplot(fig)
    with tab2:
        st.table(class_df.style.background_gradient(cmap=cmap, vmin=-6, vmax=6, axis=None))

    st.markdown("**Mean abundance, expression, mutation count and dependency change across median and low classes**")

    local_additions = st.multiselect(
        'Select additional downstream proteins to inspect',
        options=expression.index,
        placeholder='Supplement the global protein tracking pool'
    )
    show_proteins = list(set(ko_targets + focus_proteins + local_additions))

    if show_proteins:
        diff_abund_track = get_differentials(class_df, abundance, n, n, prefilter_mode='protein_list', prefilter_protein_list=show_proteins)
        diff_exp_track = get_differentials(class_df, expression, n, n, prefilter_mode='protein_list', prefilter_protein_list=show_proteins)
        diff_mut_track = get_mutation_differentials(class_df, mutation, prefilter_protein_list=show_proteins)
        diff_dep_track = get_differentials(class_df, dependency, n, n, prefilter_mode='protein_list', prefilter_protein_list=show_proteins)

        diff_summary = get_diff_summary(diff_abund_track, diff_exp_track, diff_mut_track, show_proteins, diff_dep_df=diff_dep_track, eta_sq=eta_abund, corr_abund_df=corr_abund, corr_exp_df=corr_exp, corr_dep_df=corr_dep)

        tab_dep, tab_ab, tab_ex, tab_mut = st.tabs(["Gene Dependency", "Abundance", "Expression", "Mutation"])
        with tab_dep:
            st.pyplot(get_differentials_boxplot(class_df, dependency, show_proteins, diff_df=diff_dep_track, title='Dependency', ylabel='dependency score'))
        with tab_ab:
            st.pyplot(get_differentials_boxplot(class_df, abundance, show_proteins, diff_df=diff_abund_track, title='Abundance', ylabel='abundance'))
        with tab_ex:
            st.pyplot(get_differentials_boxplot(class_df, expression, show_proteins, diff_df=diff_exp_track, title='Expression', ylabel='expression'))
        with tab_mut:
            st.pyplot(get_differentials_boxplot(class_df, mutation, show_proteins, diff_df=diff_mut_track, title='Mutation', ylabel='mutation count'))

        st.caption(
            f"q-values above are corrected across only your {len(show_proteins)} chosen proteins "
            "(m = list size), not genome-wide — see Global Settings for details."
        )

        tbl_summary, tbl_raw_ab, tbl_raw_ex, tbl_raw_mut, tbl_raw_dep = st.tabs(
            ["Summary", "Raw Abundance", "Raw Expression", "Raw Mutation", "Raw Dependency"]
        )
        with tbl_summary:
            st.dataframe(style_results_table(diff_summary, cmap), column_config=get_tissue_variance_column_config())
        with tbl_raw_ab:
            render_raw_table(abundance, show_proteins, class_df, cmap, 'Abundance', f"abundance_raw_tracked_{'_'.join(ko_targets)}", eta_sq=eta_abund)
        with tbl_raw_ex:
            render_raw_table(expression, show_proteins, class_df, cmap, 'Expression', f"expression_raw_tracked_{'_'.join(ko_targets)}", eta_sq=eta_exp)
        with tbl_raw_mut:
            render_raw_table(mutation, show_proteins, class_df, cmap, 'Mutation', f"mutation_raw_tracked_{'_'.join(ko_targets)}", eta_sq=eta_mut)
        with tbl_raw_dep:
            render_raw_table(dependency, show_proteins, class_df, cmap, 'Dependency', f"dependency_raw_tracked_{'_'.join(ko_targets)}", eta_sq=eta_dep)

    st.markdown("**Proteins most significantly affected by simulated KO of targets**")

    diff_abund_df = get_differentials(class_df, abundance, n, n, prefilter_mode=scope_mode, prefilter_min=prefilter_ab)
    diff_exp_df = get_differentials(class_df, expression, n, n, prefilter_mode=scope_mode, prefilter_min=prefilter_ex)
    diff_mut_df = get_mutation_differentials(class_df, mutation)
    diff_dep_df = get_differentials(class_df, dependency, n, n, prefilter_mode=scope_mode, prefilter_min=prefilter_dep)

    filter_factor = st.selectbox('Rank & filter most-changed proteins by', options=['Dependency', 'Abundance', 'Expression', 'Mutation'])
    n_outliers = st.slider('Number of outliers by differential', value=5, min_value=0, max_value=50)

    factor_diff_map = {'Dependency': diff_dep_df, 'Abundance': diff_abund_df, 'Expression': diff_exp_df, 'Mutation': diff_mut_df}
    selected_diff_df = factor_diff_map[filter_factor]

    diff_proteins, n_pass, n_total = select_outliers(
        selected_diff_df,
        filter_col=sig_col if sig_col else 'p', threshold=sig_threshold if sig_col else 1.01,
        n_outliers=n_outliers, min_abs_diff=diff_floor, sort_by=sort_by,
    )
    st.caption(
        f"{n_pass} of {n_total} {filter_factor.lower()} features pass |diff| ≥ {diff_floor:g}"
        + (f" and {'p' if sig_col=='p' else 'q'} < {sig_threshold:g}" if sig_col else "") + ". "
        f"Table below shows up to {2 * n_outliers} of those. Adjust these thresholds in "
        "Global Settings → Statistical Analysis Settings."
    )

    # See 1_Abundance.py for why this recomputation is needed: without it,
    # genes selected via one layer's ranking can be missing from another
    # layer's independently diff-prefiltered table, showing up as
    # unlabelled boxplots and NaN ("blacked out") summary cells.
    diff_abund_top = get_differentials(class_df, abundance, n, n, prefilter_mode='protein_list', prefilter_protein_list=diff_proteins)
    diff_exp_top = get_differentials(class_df, expression, n, n, prefilter_mode='protein_list', prefilter_protein_list=diff_proteins)
    diff_mut_top = get_mutation_differentials(class_df, mutation, prefilter_protein_list=diff_proteins)
    diff_dep_top = get_differentials(class_df, dependency, n, n, prefilter_mode='protein_list', prefilter_protein_list=diff_proteins)

    tab_dep, tab_ab, tab_ex, tab_mut = st.tabs(["Gene Dependency", "Abundance", "Expression", "Mutation"])
    with tab_dep:
        st.pyplot(get_differentials_boxplot(class_df, dependency, diff_proteins, diff_df=diff_dep_top, title='Dependency', ylabel='dependency score', use_global_order=False))
    with tab_ab:
        st.pyplot(get_differentials_boxplot(class_df, abundance, diff_proteins, diff_df=diff_abund_top, title='Abundance', ylabel='abundance', use_global_order=False))
    with tab_ex:
        st.pyplot(get_differentials_boxplot(class_df, expression, diff_proteins, diff_df=diff_exp_top, title='Expression', ylabel='expression', use_global_order=False))
    with tab_mut:
        st.pyplot(get_differentials_boxplot(class_df, mutation, diff_proteins, diff_df=diff_mut_top, title='Mutation', ylabel='mutation count', use_global_order=False))

    diff_top = get_diff_summary(diff_abund_top, diff_exp_top, diff_mut_top, diff_proteins, diff_dep_df=diff_dep_top, eta_sq=eta_abund, corr_abund_df=corr_abund, corr_exp_df=corr_exp, corr_dep_df=corr_dep)

    tbl_summary, tbl_raw_ab, tbl_raw_ex, tbl_raw_mut, tbl_raw_dep = st.tabs(
        ["Summary", "Raw Abundance", "Raw Expression", "Raw Mutation", "Raw Dependency"]
    )
    with tbl_summary:
        st.dataframe(style_results_table(diff_top, cmap), column_config=get_tissue_variance_column_config())
        st.download_button(
            f'Download All {filter_factor} Foldchanges',
            selected_diff_df.to_csv(index=True).encode('utf-8'),
            f"{filter_factor.lower()}_foldchange_{'_'.join(ko_targets)}.csv",
            "text/csv",
            key="download_outlier_foldchange"
        )
    with tbl_raw_ab:
        render_raw_table(abundance, diff_proteins, class_df, cmap, 'Abundance', f"abundance_raw_{filter_factor.lower()}_{'_'.join(ko_targets)}", eta_sq=eta_abund)
    with tbl_raw_ex:
        render_raw_table(expression, diff_proteins, class_df, cmap, 'Expression', f"expression_raw_{filter_factor.lower()}_{'_'.join(ko_targets)}", eta_sq=eta_exp)
    with tbl_raw_mut:
        render_raw_table(mutation, diff_proteins, class_df, cmap, 'Mutation', f"mutation_raw_{filter_factor.lower()}_{'_'.join(ko_targets)}", eta_sq=eta_mut)
    with tbl_raw_dep:
        render_raw_table(dependency, diff_proteins, class_df, cmap, 'Dependency', f"dependency_raw_{filter_factor.lower()}_{'_'.join(ko_targets)}", eta_sq=eta_dep)

else:
    st.warning("⚠️ No Knockout targets set. Go to **Global Settings** to pick your target proteins for KO simulation.")