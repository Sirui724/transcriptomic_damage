# =========================================================
# Gene fusion analysis in mouse bulk RNA-seq
# Author: Sirui
# Description:
#   1. Load gene fusion counts
#   2. Normalize by library size / unique reads
#   3. Test correlation with age
#   4. Perform GO enrichment for age-associated genes
#   5. Normalize fusion counts by gene expression
#   6. Repeat association analysis on expression-normalized values
# =========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(reshape2)
  library(ggplot2)
  library(clusterProfiler)
  library(org.Mm.eg.db)
})

# -----------------------------
# File paths
# -----------------------------
fusion_file <- "/Lasagna/sirui/damage/gene_fusion/mouse/merge_gf_gene.txt"
metadata_file <- "/IGF1/sirui/GSE_GSE132040_sc/GSE132040_MACA_Bulk_metadata.csv"
gene_count_file <- "/Kimchi/sirui/GSE132040/GSE132040_genecount.txt"

out_dir <- "/Lasagna/sirui/damage/gene_fusion/mouse"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Helper functions
# -----------------------------
run_feature_age_cor <- function(df, feature_cols, age_col) {
  res_list <- lapply(feature_cols, function(feat) {
    x <- suppressWarnings(as.numeric(df[[feat]]))
    y <- suppressWarnings(as.numeric(df[[age_col]]))

    keep <- complete.cases(x, y)

    if (sum(keep) < 3) {
      return(data.frame(
        gene = feat,
        P = NA_real_,
        cor = NA_real_,
        stringsAsFactors = FALSE
      ))
    }

    test <- cor.test(x[keep], y[keep], use = "complete.obs")

    data.frame(
      gene = feat,
      P = test$p.value,
      cor = unname(test$estimate),
      stringsAsFactors = FALSE
    )
  })

  bind_rows(res_list)
}

add_significance_columns <- function(df, p_cutoff = 0.05, cor_cutoff = 0.1) {
  df %>%
    mutate(
      neg_log10_p = -log10(P),
      significance = case_when(
        !is.na(P) & P < p_cutoff & abs(cor) > cor_cutoff ~ "Significant",
        TRUE ~ "Not Significant"
      )
    )
}

plot_volcano_like_cor <- function(df, label_sig = FALSE) {
  p <- ggplot(df, aes(x = cor, y = neg_log10_p, color = significance)) +
    geom_point(alpha = 0.7, size = 2.5) +
    scale_color_manual(values = c("Significant" = "red", "Not Significant" = "grey70")) +
    labs(
      title = NULL,
      x = "Correlation with age",
      y = expression(-log[10](p-value)),
      color = "Significance"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "top")

  if (label_sig) {
    p <- p +
      geom_text(
        data = df %>% filter(significance == "Significant"),
        aes(label = gene),
        hjust = 0,
        vjust = 0,
        size = 3,
        color = "black"
      )
  }

  p
}

run_go_enrichment <- function(gene_vec) {
  gene_vec <- unique(na.omit(gene_vec))

  if (length(gene_vec) == 0) {
    return(NULL)
  }

  enrichGO(
    gene = gene_vec,
    OrgDb = org.Mm.eg.db,
    ont = "ALL",
    keyType = "SYMBOL"
  )
}

# -----------------------------
# Step 1. Load fusion matrix
# -----------------------------
data_gf <- read.table(fusion_file, header = TRUE, check.names = FALSE)

# TODO:
# The original script uses:
#   data_gf <- data_gf[, c("name", data$raw.file)]
#   mapply("/", data_gf[, 2:858], data$uniq_reads)
# but the object `data` is not defined in the provided script.
# Replace `sample_info` below with the correct object containing:
#   - raw.file
#   - uniq_reads

# sample_info <- ...
# data_gf <- data_gf[, c("name", sample_info$raw.file)]
# temp <- as.data.frame(mapply("/", data_gf[, -1], sample_info$uniq_reads))
# data_gf[, -1] <- temp

# -----------------------------
# Step 2. Reshape fusion matrix
# -----------------------------
data_gf_t <- as.data.frame(t(data_gf), stringsAsFactors = FALSE)
colnames(data_gf_t) <- data_gf_t[1, ]
data_gf_t <- data_gf_t[-1, , drop = FALSE]
data_gf_t$raw.file <- rownames(data_gf_t)

# -----------------------------
# Step 3. Load and clean metadata
# -----------------------------
d_sample <- read.csv(metadata_file, check.names = FALSE)

d_sample <- separate(
  d_sample,
  source.name,
  into = c("tissue"),
  sep = "_",
  remove = FALSE,
  convert = FALSE
)

d_sample <- filter(d_sample, characteristics..sex != "missing")
d_sample <- d_sample[, c(1, 3, 4, 6, 8, 9, 12)]

# Original script used: merge(..., bt = "raw.file")
# That is almost certainly a typo and should be `by = "raw.file"`
data_gf_ <- merge(data_gf_t, d_sample, by = "raw.file")

# -----------------------------
# Step 4. Correlation analysis:
# raw fusion burden vs age
# -----------------------------
feature_cols_raw <- colnames(data_gf_)[2:17497]

B <- run_feature_age_cor(
  df = data_gf_,
  feature_cols = feature_cols_raw,
  age_col = "characteristics..age"
)

B <- add_significance_columns(B)

p1 <- plot_volcano_like_cor(B, label_sig = FALSE)
print(p1)

write.table(
  B,
  file = file.path(out_dir, "gene_correlation.txt"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

B_pos_raw <- filter(B, cor > 0.1, P < 0.05)
B_neg_raw <- filter(B, cor < -0.1, P < 0.05)

# -----------------------------
# Step 5. GO enrichment:
# positively age-associated genes
# -----------------------------
go_raw <- run_go_enrichment(B_pos_raw$gene)

if (!is.null(go_raw)) {
  print(dotplot(go_raw, split = "ONTOLOGY") + facet_grid(ONTOLOGY ~ ., scale = "free"))
  print(barplot(go_raw) +
    scale_fill_gradient2(low = "darkred", mid = "orange", high = "yellow", midpoint = 0.02))
  print(cnetplot(go_raw))
}

# -----------------------------
# Step 6. Normalize fusion counts
# by matched gene expression
# -----------------------------
d_gene <- read.table(gene_count_file, header = TRUE, check.names = FALSE)

d_gene <- d_gene[which(d_gene$gene %in% data_gf$name), ]
d_gene_long <- melt(
  d_gene,
  id.vars = "gene",
  variable.name = "raw.file",
  value.name = "gene_count"
)
colnames(d_gene_long)[1] <- "name"

fusion_long <- melt(
  data_gf,
  id.vars = "name",
  variable.name = "raw.file",
  value.name = "count"
)

temp <- merge(fusion_long, d_gene_long, by = c("name", "raw.file"))
temp <- unique(temp)

temp[is.na(temp)] <- 0
temp$count_norm <- temp$count / temp$gene_count
temp[is.na(temp)] <- 0
temp[temp == "Inf"] <- NA
temp$count_norm[temp$count_norm >= 1] <- NA

data_gf_norm <- dcast(
  unique(temp[, c(1, 2, 5)]),
  name ~ raw.file,
  value.var = "count_norm",
  fun.aggregate = mean
)

data_gf_t_norm <- as.data.frame(t(data_gf_norm), stringsAsFactors = FALSE)
colnames(data_gf_t_norm) <- data_gf_t_norm[1, ]
data_gf_t_norm <- data_gf_t_norm[-1, , drop = FALSE]
data_gf_t_norm$raw.file <- rownames(data_gf_t_norm)

data_gf_n <- merge(data_gf_t_norm, d_sample, by = "raw.file")
data_gf_n[data_gf_n == "         Inf"] <- NA
data_gf_n <- data_gf_n[, colSums(!is.na(data_gf_n)) > 5]

write.table(
  data_gf_n,
  file = file.path(out_dir, "gf_merge_gene_norm.txt"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# -----------------------------
# Step 7. Correlation analysis:
# expression-normalized fusion burden vs age
# -----------------------------
feature_cols_norm <- colnames(data_gf_n)[2:11985]

B_norm <- run_feature_age_cor(
  df = data_gf_n,
  feature_cols = feature_cols_norm,
  age_col = "characteristics..age"
)

B_norm <- add_significance_columns(B_norm)

p2 <- plot_volcano_like_cor(B_norm, label_sig = FALSE)
print(p2)

p3 <- plot_volcano_like_cor(B_norm, label_sig = TRUE)
print(p3)

write.table(
  B_norm,
  file = file.path(out_dir, "gene_norm_correlation.txt"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

B_pos_norm <- filter(B_norm, cor > 0.1, P < 0.05)
B_neg_norm <- filter(B_norm, cor < -0.1, P < 0.05)

# -----------------------------
# Step 8. GO enrichment:
# normalized positively associated genes
# -----------------------------
go_norm <- run_go_enrichment(B_pos_norm$gene)

if (!is.null(go_norm)) {
  print(
    barplot(go_norm, label_format = 50) +
      scale_fill_gradient2(low = "darkred", mid = "orange", high = "yellow") +
      facet_grid(ONTOLOGY ~ ., scale = "free", space = "free")
  )
  print(cnetplot(go_norm))
}

# -----------------------------
# Step 9. Keep positively associated features
# -----------------------------
B_pos <- filter(B_norm, cor > 0)

data_gf_n_pos <- data_gf_n[, c("raw.file", B_pos$gene), drop = FALSE]
colnames(data_gf_n_pos)[-1] <- paste0(colnames(data_gf_n_pos)[-1], "_gf")
