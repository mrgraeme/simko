"""
One-time preprocessing step for GSE101966_HCT116_FPKM.txt.gz (and similar
HOMER analyzeRepeats output): converts the RefSeq transcript IDs
('Transcript/RepeatID' column -- NM_/NR_/XM_/XR_ accessions) to gene
symbols, so the result can be joined against SimKO's own gene-symbol-indexed
expression matrix in concordance_analysis.py.

Run this LOCALLY, not in a sandboxed/network-restricted environment -- it
queries the mygene.info web service via the `mygene` package, which needs
outbound internet access.

    pip install mygene
    python convert_refseq_to_symbol.py GSE101966_HCT116_FPKM.txt.gz GSE101966_HCT116_FPKM_bysymbol.txt

Multiple RefSeq transcripts (isoforms) commonly map to the same gene
symbol -- FPKM values for transcripts sharing a symbol are SUMMED per
sample, a defensible (if not the only reasonable) convention for going
from isoform-level to gene-level expression. Rows that aren't RefSeq
transcript accessions at all (this file mixes transcripts and repeat
elements under one ID column, per its own column name) are dropped before
querying, since they were never going to have a gene symbol anyway.
"""
import sys
import re
import pandas as pd

try:
    import mygene
except ImportError:
    sys.exit("This script needs the mygene package: pip install mygene")


def convert_refseq_to_gene_symbol(fpkm_path, output_path, id_col='Transcript/RepeatID', sep='\t'):
    df = pd.read_csv(fpkm_path, sep=sep)
    if id_col not in df.columns:
        sys.exit(f"Column {id_col!r} not found. Actual columns: {df.columns.tolist()}")

    is_refseq = df[id_col].astype(str).str.match(r'^[NX][MR]_')
    print(f"{is_refseq.sum()} of {len(df)} rows look like RefSeq transcript IDs "
          f"(NM_/NR_/XM_/XR_ prefix); {(~is_refseq).sum()} other rows (likely repeat "
          f"elements, sharing this ID column per the file's own naming) will be dropped.")
    df = df[is_refseq].copy()

    mg = mygene.MyGeneInfo()
    unique_ids = df[id_col].unique().tolist()
    print(f"Querying mygene.info for {len(unique_ids)} unique RefSeq IDs (may take a minute)...")
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
          f"({n_no_hit} not found in mygene.info, {n_no_symbol} found but with no symbol field).")

    df['gene_symbol'] = df[id_col].map(id_to_symbol)
    n_before = len(df)
    df = df.dropna(subset=['gene_symbol'])
    print(f"Dropped {n_before - len(df)} rows with no resolved gene symbol; {len(df)} rows remain.")

    fpkm_cols = [c for c in df.columns if c.strip().endswith('FPKM')]
    if not fpkm_cols:
        sys.exit(f"No columns ending in 'FPKM' found. Actual columns: {df.columns.tolist()}")

    n_genes_before_collapse = df['gene_symbol'].nunique()
    gene_level = df.groupby('gene_symbol')[fpkm_cols].sum()
    print(f"Collapsed {len(df)} transcript rows into {n_genes_before_collapse} genes "
          f"(summing FPKM across isoforms sharing a symbol).")

    gene_level.to_csv(output_path, sep='\t')
    print(f"Saved {output_path} -- {gene_level.shape[0]} genes x {gene_level.shape[1]} sample columns.")
    return gene_level


if __name__ == '__main__':
    if len(sys.argv) != 3:
        sys.exit("Usage: python convert_refseq_to_symbol.py <input_fpkm_file> <output_path>")
    convert_refseq_to_gene_symbol(sys.argv[1], sys.argv[2])