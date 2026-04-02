#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Human damage pipeline helper script

This script performs several tasks:
1. Merge exon-domain BED summaries across samples
2. Annotate merged regions with overlapping gene names from a GTF
3. Keep only protein-coding genes
4. Merge splice-site (ss) summary files
5. Merge retained-event (re) summary files
6. Merge RSEM expected counts across samples
7. Convert ENSEMBL gene IDs to gene symbols
"""

import bisect
import csv
import glob
import gzip
import os
import re
from collections import defaultdict

import pandas as pd


# =========================================================
# Configuration
# =========================================================

HUMAN_DAMAGE_DIR = "/IGF1/sirui/human_damage_out"
PIPELINE_OUT_DIR = "/Entropy/siruizhang/damage_pipeline/human"

GTF_PATH = "/Entropy/siruizhang/reference/hg38/Homo_sapiens.GRCh38.115.gtf"
SYMBOL_MAP_PATH = "/Entropy/siruizhang/damage_pipeline/human/ENSEMBL_to_SYMBOL_hg38.csv"

MERGED_BED_OUTPUT = "merged_result.tsv"
ANNOTATED_BED_OUTPUT = "annotated_result.tsv"
CODING_BED_OUTPUT = "coding_genes_result.tsv"

MERGED_SS_OUTPUT = os.path.join(PIPELINE_OUT_DIR, "merge_ss.txt")
MERGED_RE_OUTPUT = os.path.join(PIPELINE_OUT_DIR, "merge_re.txt")
MERGED_EXPECTED_COUNTS_OUTPUT = os.path.join(PIPELINE_OUT_DIR, "merged_expected_counts.txt")
MERGED_EXPECTED_COUNTS_SYMBOL_OUTPUT = os.path.join(
    PIPELINE_OUT_DIR, "merged_expected_counts_SYMBOL.txt"
)


# =========================================================
# Part 1. Merge exon-domain BED files
# =========================================================

def merge_bed_files(input_pattern: str, output_path: str, na_threshold_fraction: float = 1.0) -> pd.DataFrame:
    """
    Merge all matching BED-like files using column 4 as region name
    and column 5 as value.

    Parameters
    ----------
    input_pattern : str
        Glob pattern for input files.
    output_path : str
        Output TSV path.
    na_threshold_fraction : float
        Keep rows with NA count < n_columns * fraction.
        Original script effectively used all columns as threshold.
    """
    files = glob.glob(input_pattern)
    merged_df = pd.DataFrame()

    for file_path in files:
        sample_name = os.path.basename(file_path.split("_exon")[0])
        df = pd.read_csv(
            file_path,
            sep="\t",
            header=None,
            usecols=[3, 4],
            names=["name", sample_name]
        )
        df = df.drop_duplicates(subset=["name"], keep="first")
        df = df.set_index("name")

        if merged_df.empty:
            merged_df = df
        else:
            merged_df = merged_df.merge(df, how="outer", left_index=True, right_index=True)

    if not merged_df.empty:
        max_na = merged_df.shape[1] * na_threshold_fraction
        merged_df = merged_df.loc[merged_df.isna().sum(axis=1) < max_na]

    merged_df.fillna("NA").to_csv(output_path, sep="\t")
    return merged_df


# =========================================================
# Part 2. Build gene index from GTF and annotate merged table
# =========================================================

def build_gene_index(gtf_path: str):
    """
    Parse GTF and build a chromosome-wise gene interval index.

    Returns
    -------
    dict
        {chrom: (starts, [(start, end, gene_name), ...])}
    """
    gene_indices = defaultdict(lambda: ([], []))

    with open(gtf_path, "r") as handle:
        for line in handle:
            if line.startswith("#"):
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue

            chrom = parts[0]
            start = int(parts[3])
            end = int(parts[4])
            attrs = parts[8]

            match = re.search(r'gene_name "([^"]+)"', attrs)
            if not match:
                continue

            gene_name = match.group(1)

            starts, genes = gene_indices[chrom]
            idx = bisect.bisect_left(starts, start)
            starts.insert(idx, start)
            genes.insert(idx, (start, end, gene_name))

    return gene_indices


def annotate_genes(input_tsv: str, output_tsv: str, gene_indices: dict) -> None:
    """
    Annotate TSV row names of form chrom:start-end with overlapping gene names.
    """
    with open(input_tsv, "r") as f_in, open(output_tsv, "w") as f_out:
        header = f_in.readline().rstrip("\n")
        f_out.write(f"{header}\n")

        for line in f_in:
            line = line.rstrip("\n")
            if not line:
                continue

            parts = line.split("\t")
            original_name = parts[0]

            try:
                chrom, pos = original_name.split(":", 1)
                start_str, end_str = pos.split("-", 1)
                start = int(start_str)
                end = int(end_str)
            except ValueError:
                f_out.write(f"{line}\n")
                continue

            candidate_chroms = [chrom]
            if chrom.startswith("chr"):
                candidate_chroms.append(chrom[3:])
            else:
                candidate_chroms.append(f"chr{chrom}")

            gene_names = []

            for candidate_chrom in candidate_chroms:
                if candidate_chrom not in gene_indices:
                    continue

                starts, genes = gene_indices[candidate_chrom]
                idx = bisect.bisect_right(starts, end) - 1

                seen = set()
                while idx >= 0:
                    g_start, g_end, gene_name = genes[idx]
                    if g_end < start:
                        break
                    if gene_name not in seen:
                        gene_names.append(gene_name)
                        seen.add(gene_name)
                    idx -= 1

                if gene_names:
                    break

            if gene_names:
                parts[0] = f"{','.join(gene_names)}@chr{original_name}"

            f_out.write("\t".join(parts) + "\n")


# =========================================================
# Part 3. Extract protein-coding genes and filter
# =========================================================

def extract_protein_coding_genes(gtf_path: str) -> set:
    """
    Extract protein-coding gene names from a GTF/GTF.GZ file.
    """
    protein_coding_genes = set()
    open_func = gzip.open if gtf_path.endswith(".gz") else open

    with open_func(gtf_path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue

            attributes = parts[8]
            attr_dict = {
                m.group(1): m.group(2)
                for m in re.finditer(r'(\w+)\s+"([^"]+)"', attributes)
            }

            gene_type = attr_dict.get("gene_biotype", attr_dict.get("gene_type", ""))
            if gene_type == "protein_coding":
                gene_name = attr_dict.get("gene_name")
                if gene_name:
                    protein_coding_genes.add(gene_name)

    return protein_coding_genes


def filter_protein_coding_genes(input_tsv: str, output_tsv: str, protein_coding_genes: set) -> None:
    """
    Keep only rows whose annotated gene name is protein-coding.
    """
    with open(input_tsv, "r") as f_in, open(output_tsv, "w", newline="") as f_out:
        reader = csv.reader(f_in, delimiter="\t")
        writer = csv.writer(f_out, delimiter="\t")

        header = next(reader)
        writer.writerow(header)

        for row in reader:
            gene_name = row[0].split("@")[0]
            if gene_name in protein_coding_genes:
                writer.writerow(row)


# =========================================================
# Part 4. Merge ss / re summary files
# =========================================================

def merge_event_files(
    input_dir: str,
    suffix: str,
    output_path: str,
    include_gene_column: bool = True,
    add_event_type_to_gene: str = None,
    na_threshold_fraction: float = 1.0,
    max_row_len: int = None
) -> None:
    """
    Merge per-sample event summary files such as *_ss.txt or *_re.txt.

    Expected file format:
        col0 = event name
        col1 = gene annotation
        col2 = value
    """
    files = sorted(os.listdir(input_dir))
    event_values = {}
    annotations = {}
    sample_names = []
    sample_idx = 0

    for file_name in files:
        if not file_name.endswith(suffix):
            continue

        sample_name = file_name.split(suffix)[0]
        sample_names.append(sample_name)

        temp = {}
        file_path = os.path.join(input_dir, file_name)

        with open(file_path, "r") as handle:
            for line in handle:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 3:
                    continue

                event_name = parts[0]
                gene_anno = parts[1]
                value = parts[2]

                if add_event_type_to_gene is not None:
                    gene_anno = f"{gene_anno}\t{add_event_type_to_gene}"

                annotations[event_name] = gene_anno
                temp[event_name] = value

        if not event_values:
            for event_name, value in temp.items():
                event_values[event_name] = [value]
        else:
            for event_name, value in temp.items():
                if event_name in event_values:
                    event_values[event_name].append(value)
                else:
                    event_values[event_name] = ["NA"] * sample_idx + [value]

            for event_name in event_values:
                if event_name not in temp:
                    event_values[event_name].append("NA")

        sample_idx += 1

    with open(output_path, "w") as out:
        if include_gene_column:
            out.write("name\tgene\t" + "\t".join(sample_names) + "\n")
        else:
            out.write("name\t" + "\t".join(sample_names) + "\n")

        for event_name, values in event_values.items():
            na_count = values.count("NA")
            if na_count >= sample_idx * na_threshold_fraction:
                continue

            if max_row_len is not None and len(values) >= max_row_len:
                continue

            if include_gene_column:
                out.write(f"{event_name}\t{annotations.get(event_name, 'NA')}\t" + "\t".join(values) + "\n")
            else:
                out.write(f"{event_name}\t" + "\t".join(values) + "\n")


# =========================================================
# Part 5. Merge RSEM expected counts
# =========================================================

def merge_expected_counts(input_pattern: str, output_path: str) -> pd.DataFrame:
    """
    Merge expected_count from *.genes.results files by gene_id.
    """
    files = glob.glob(input_pattern)
    merged_df = pd.DataFrame()

    for file_path in files:
        sample_name = file_path.split("out/")[1].split(".gene")[0]

        df = pd.read_csv(
            file_path,
            sep="\t",
            usecols=["gene_id", "expected_count"]
        ).set_index("gene_id")

        df.rename(columns={"expected_count": sample_name}, inplace=True)

        if merged_df.empty:
            merged_df = df
        else:
            merged_df = merged_df.join(df, how="outer")

    merged_df = merged_df.fillna(0)
    merged_df.to_csv(output_path, sep="\t", index=True)
    return merged_df


# =========================================================
# Part 6. Convert ENSEMBL IDs to gene symbols
# =========================================================

def convert_ensembl_to_symbol(merged_path: str, symbol_map_path: str, output_path: str) -> pd.DataFrame:
    """
    Convert merged expected counts from ENSEMBL gene IDs to SYMBOL.
    Duplicate symbols are summed.
    """
    merged_df = pd.read_csv(merged_path, sep="\t", index_col="gene_id")
    symbol_df = pd.read_csv(symbol_map_path)

    merged_with_symbol = merged_df.merge(
        symbol_df,
        left_index=True,
        right_on="gene_id",
        how="left"
    )

    merged_with_symbol["SYMBOL"] = merged_with_symbol["SYMBOL"].fillna("Unknown")

    numeric_cols = merged_df.columns.tolist()
    merged_with_symbol = merged_with_symbol.groupby("SYMBOL")[numeric_cols].sum()

    merged_with_symbol.to_csv(output_path, sep="\t")
    return merged_with_symbol


# =========================================================
# Main
# =========================================================

def main():
    os.makedirs(PIPELINE_OUT_DIR, exist_ok=True)

    # 1. Merge exon-domain BED summaries
    merge_bed_files(
        input_pattern=os.path.join(HUMAN_DAMAGE_DIR, "*SJ_exon_domain_uniq.bed"),
        output_path=MERGED_BED_OUTPUT,
        na_threshold_fraction=1.0
    )

    # 2. Annotate merged regions with genes
    gene_indices = build_gene_index(GTF_PATH)
    annotate_genes(MERGED_BED_OUTPUT, ANNOTATED_BED_OUTPUT, gene_indices)

    # 3. Keep protein-coding genes only
    protein_coding_genes = extract_protein_coding_genes(GTF_PATH)
    print(f"Found {len(protein_coding_genes)} protein-coding genes.")
    filter_protein_coding_genes(ANNOTATED_BED_OUTPUT, CODING_BED_OUTPUT, protein_coding_genes)
    print("Protein-coding gene filtering completed.")

    # 4. Merge ss files
    merge_event_files(
        input_dir=HUMAN_DAMAGE_DIR,
        suffix="_ss.txt",
        output_path=MERGED_SS_OUTPUT,
        include_gene_column=True,
        add_event_type_to_gene=None,
        na_threshold_fraction=1.0
    )
    print("ss merge done.")

    # 5. Merge re files
    merge_event_files(
        input_dir=HUMAN_DAMAGE_DIR,
        suffix="_re.txt",
        output_path=MERGED_RE_OUTPUT,
        include_gene_column=False,
        add_event_type_to_gene="re",
        na_threshold_fraction=0.5,
        max_row_len=200
    )
    print("re merge done.")

    # 6. Merge expected counts
    merge_expected_counts(
        input_pattern=os.path.join(HUMAN_DAMAGE_DIR, "*.genes.results"),
        output_path=MERGED_EXPECTED_COUNTS_OUTPUT
    )
    print("Expected count merge done.")

    # 7. Convert ENSEMBL gene IDs to SYMBOL
    convert_ensembl_to_symbol(
        merged_path=MERGED_EXPECTED_COUNTS_OUTPUT,
        symbol_map_path=SYMBOL_MAP_PATH,
        output_path=MERGED_EXPECTED_COUNTS_SYMBOL_OUTPUT
    )
    print("ENSEMBL to SYMBOL conversion done.")


if __name__ == "__main__":
    main()
