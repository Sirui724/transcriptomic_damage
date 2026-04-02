import re
import pandas as pd
import numpy as np


# =========================================================
# Helper functions
# =========================================================

def compute_ratio(numerator_df, denominator_df):
    """
    Compute element-wise ratio with safe handling of zeros.

    Rules:
        - denominator == 0 → treated as NaN
        - NaN results → replaced with 0
    """
    ratio = numerator_df.div(denominator_df.replace(0, np.nan))
    return ratio.fillna(0)


def align_matrices(df1, df2):
    """
    Align two matrices by common genes (index) and samples (columns).
    """
    common_genes = df1.index.intersection(df2.index)
    common_samples = df1.columns.intersection(df2.columns)

    df1_common = df1.loc[common_genes, common_samples]
    df2_common = df2.loc[common_genes, common_samples]

    return df1_common, df2_common, common_genes, common_samples


# =========================================================
# Part 1. Junction / gene expression ratio
# =========================================================

def compute_junction_ratio(expr_file, junc_file, output_file):
    expr_df = pd.read_csv(expr_file, sep="\t", index_col="SYMBOL")
    junc_df = pd.read_csv(junc_file, sep="\t", index_col="name")

    # Normalize column names (remove "SJ")
    junc_df.columns = [col.replace("SJ", "") for col in junc_df.columns]

    expr_common, junc_common, genes, samples = align_matrices(expr_df, junc_df)

    result = compute_ratio(junc_common, expr_common)

    result.to_csv(output_file, sep="\t")

    print(f"[Junction] Done: {output_file}")
    print(f"  Genes: {len(genes)}")
    print(f"  Samples: {len(samples)}")


# =========================================================
# Part 2. Stopsite / gene expression ratio
# =========================================================

def compute_stopsite_ratio(expr_file, stopsite_file, output_file):
    expr_df = pd.read_csv(expr_file, sep="\t", index_col="SYMBOL")
    stop_df = pd.read_csv(stopsite_file, sep="\t", index_col="gene")

    expr_common, stop_common, genes, samples = align_matrices(expr_df, stop_df)

    result = compute_ratio(stop_common, expr_common)

    result.to_csv(output_file, sep="\t")

    print(f"[Stopsite] Done: {output_file}")
    print(f"  Genes: {len(genes)}")
    print(f"  Samples: {len(samples)}")


# =========================================================
# Part 3. Normalize merge_re by mapped reads
# =========================================================

def normalize_merge_re(merge_re_file, damage_info_file, output_file, scale=1e8):
    """
    Normalize merge_re counts by mapped reads.

    Formula:
        normalized = count * scale / mapped_reads
    """
    damage_df = pd.read_csv(damage_info_file, sep="\t")
    mapped_dict = dict(zip(damage_df["sample"], damage_df["mapped_reads"]))

    merge_re = pd.read_csv(merge_re_file, sep="\t", index_col="name")

    # Normalize column names: remove trailing letters (e.g., "626-1a" → "626-1")
    merge_re.columns = [
        re.sub(r'([0-9]+-[0-9]+)[a-zA-Z]*', r'\1', col)
        for col in merge_re.columns
    ]

    for col in merge_re.columns:
        if col in mapped_dict:
            divisor = mapped_dict[col]

            if divisor == 0:
                merge_re[col] = 0
            else:
                merge_re[col] = merge_re[col] * scale / divisor
        else:
            merge_re[col] = 0

    result = merge_re.round(6).fillna(0)
    result.to_csv(output_file, sep="\t")

    print(f"[merge_re normalization] Done: {output_file}")


# =========================================================
# Main
# =========================================================

if __name__ == "__main__":

    # File paths
    expr_file = "merged_expected_counts_SYMBOL.txt"

    # 1. Junction ratio
    compute_junction_ratio(
        expr_file=expr_file,
        junc_file="junc_gene.tsv",
        output_file="junc_gene_ratio_result.txt"
    )

    # 2. Stopsite ratio
    compute_stopsite_ratio(
        expr_file=expr_file,
        stopsite_file="stopsite_merge_gene.txt",
        output_file="stopsite_ratio_result.txt"
    )

    # 3. merge_re normalization
    normalize_merge_re(
        merge_re_file="merge_re.txt",
        damage_info_file="damage_info.txt",
        output_file="normalized_merge_re.txt"
    )
