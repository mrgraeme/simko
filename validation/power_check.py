"""
Power check for the ARID1B tissue-matched concordance result.

The LARGE_INTESTINE-restricted cohort (27 cell lines) gave SimKO's dynamic
cohort-sizing formula n=8 per side -- well under half the pan-cancer n=20.
Before reading the tissue-matched result (Spearman rho=0.019, p=0.12,
essentially null) as evidence against real concordance, it's worth asking a
simpler question first: does a gene we ALREADY KNOW has a strong, real
effect still show that effect clearly at n=8? If yes, n=8 alone isn't the
explanation for ARID1B's null-ish tissue-matched result. If PCNA's own
known-strong signal also collapses at n=8, that's a much more mundane
explanation than "no real concordance" -- it just means the tissue-matched
test never had enough power to detect anything, real or not.

Uses PCNA specifically because Phase 0 already established it has the
strongest, cleanest real signal of any target tested (94% estimated
real-hit rate at p<0.01, pan-cancer, Abundance layer). Run here on the
EXPRESSION layer instead, matching the layer actually used in the ARID1B
concordance work, so the comparison is apples-to-apples with that result
rather than mixing layers.
"""
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
APP_ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(APP_ROOT))

from validation import load_matrices, run_permutation_validation


if __name__ == '__main__':
    abundance, expression, mutation, dependency = load_matrices()

    print(f"{'=' * 70}\nPCNA, Expression layer, pan-cancer, standard n=20\n{'=' * 70}")
    result_n20 = run_permutation_validation(['PCNA'], expression, n=20, n_permutations=200, seed=42)
    print(result_n20['fdr_curve'].to_string(index=False))

    print(f"\n{'=' * 70}\nPCNA, Expression layer, pan-cancer, forced n=8 (matching the "
          f"LARGE_INTESTINE tissue-matched cohort size)\n{'=' * 70}")
    result_n8 = run_permutation_validation(['PCNA'], expression, n=8, n_permutations=200, seed=42)
    print(result_n8['fdr_curve'].to_string(index=False))

    print(f"\n{'=' * 70}\nComparison at p<0.01\n{'=' * 70}")
    fdr20 = result_n20['fdr_curve']
    fdr8 = result_n8['fdr_curve']
    row20 = fdr20[fdr20['nominal_p_threshold'] == 0.01].iloc[0]
    row8 = fdr8[fdr8['nominal_p_threshold'] == 0.01].iloc[0]
    print(f"n=20: empirical FDR = {row20['empirical_fdr']:.4f} ({row20['observed_hits']} hits)")
    print(f"n=8:  empirical FDR = {row8['empirical_fdr']:.4f} ({row8['observed_hits']} hits)")
    print(
        "\nIf n=8's empirical FDR is still clearly well below the ~0.85-0.96 pure-noise "
        "ceiling (established in Phase 0's synthetic testing), PCNA's known-real signal "
        "survives the reduced cohort size reasonably well -- meaning n=8 alone probably "
        "isn't why ARID1B's tissue-matched concordance came back null, and that result is "
        "more likely telling us something real about the comparison itself. If n=8's "
        "empirical FDR has risen sharply toward that noise ceiling, low power is a "
        "genuine, mundane explanation for ARID1B's null tissue-matched result, and the "
        "concordance question there remains genuinely open rather than answered."
    )