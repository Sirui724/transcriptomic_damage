import csv


def safe_float(value, default=0.0):
    """Convert a string to float; return default if conversion fails."""
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def safe_int(value, default=0):
    """Convert a string to int; return default if conversion fails."""
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def merge_genes(input_file, output_file):
    """
    Merge junction-level values to gene-level values.

    Input format:
        First column: gene annotation with genomic position, e.g. "GENE1,GENE2@chr1:100-200"
        Remaining columns: sample values

    Rules:
        - 'NA' is treated as 0
        - if one junction maps to multiple genes, its value is added to each gene
    """
    gene_dict = {}

    with open(input_file, "r", newline="") as f_in:
        reader = csv.reader(f_in, delimiter="\t")
        headers = next(reader)
        sample_columns = headers[1:]

        for row in reader:
            if not row:
                continue

            gene_location = row[0]
            genes_part = gene_location.split("@")[0]
            genes = [g.strip() for g in genes_part.split(",") if g.strip()]

            values = []
            for val in row[1:]:
                if val.strip().upper() == "NA":
                    values.append(0.0)
                else:
                    values.append(safe_float(val, default=0.0))

            for gene in genes:
                if gene not in gene_dict:
                    gene_dict[gene] = [0.0] * len(sample_columns)

                for i, value in enumerate(values):
                    gene_dict[gene][i] += value

    with open(output_file, "w", newline="") as f_out:
        writer = csv.writer(f_out, delimiter="\t")
        writer.writerow(["name"] + sample_columns)

        for gene in sorted(gene_dict):
            counts = [f"{val:.1f}" for val in gene_dict[gene]]
            writer.writerow([gene] + counts)


def merge_stopsite_counts(input_file, output_file, max_count=10):
    """
    Merge stop-site counts at the gene level.

    Input format:
        Column 1: event name
        Column 2: gene name
        Column 3+: sample counts

    Rules:
        - header line starts with 'name'
        - 'NA' is treated as 0
        - only counts <= max_count are kept
        - counts > max_count are treated as 0
    """
    gene_counts = {}
    sample_header = None

    with open(input_file, "r", newline="") as f_in:
        reader = csv.reader(f_in, delimiter="\t")

        for row in reader:
            if not row:
                continue

            if row[0] == "name":
                sample_header = row[1:]
                continue

            if len(row) < 3:
                continue

            gene_name = row[1]

            if gene_name not in gene_counts:
                gene_counts[gene_name] = [0] * (len(row) - 2)

            for i, value in enumerate(row[2:]):
                if value == "NA":
                    continue

                count = safe_int(value, default=0)
                if count <= max_count:
                    gene_counts[gene_name][i] += count

    with open(output_file, "w", newline="") as f_out:
        writer = csv.writer(f_out, delimiter="\t")

        if sample_header is None:
            raise ValueError("Header row starting with 'name' was not found in the input file.")

        writer.writerow(sample_header)

        for gene in sorted(gene_counts):
            writer.writerow([gene] + gene_counts[gene])


if __name__ == "__main__":
    merge_genes("coding_genes_result.tsv", "junc_gene.tsv")
    merge_stopsite_counts("merge_ss.txt", "stopsite_merge_gene.txt", max_count=10)
