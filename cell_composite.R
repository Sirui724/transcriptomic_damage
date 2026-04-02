devtools::install_github('AlmogAngel/xCell2')
library(xCell2)

data("MouseRNAseqData.xCell2Ref", package = "xCell2")

d_gene <- read.table('/Kimchi/sirui/GSE132040/GSE132040_genecount.txt', head = TRUE, row.names = 1)
d_info <- read.csv('/Lasagna/sirui/damage/GSE132040_MACA_Bulk_metadata.csv')

d_info <- data_2[, c(1, 8, 10, 12, 22)]
colnames(d_info)[5] <- 'all'

d_gene[is.na(d_gene)] <- 0
d_gene <- d_gene[, colSums(d_gene != 0) > 10000]

xcell2_results <- xCell2::xCell2Analysis(
  mix = log2(d_gene + 1),
  xcell2object = MouseRNAseqData.xCell2Ref,
  spillover = TRUE
)

data("TabulaMurisBlood.xCell2Ref", package = "xCell2")

d_gene_blood <- d_gene[, which(colnames(d_gene) %in%
  filter(d_info, tissue == 'Spleen' | tissue == 'Marrow' |
                   tissue == 'Bone' | tissue == 'WBC')$raw.file)]

xcell2_results <- xCell2::xCell2Analysis(
  mix = log2(d_gene_blood + 1),
  xcell2object = TabulaMurisBlood.xCell2Ref,
  spillover = TRUE
)

blood_xcell <- as.data.frame(t(xcell2_results))
blood_xcell$raw.file <- rownames(blood_xcell)
blood_xcell <- merge(blood_xcell, d_info, by = 'raw.file')

write.table(
  blood_xcell,
  file = '/Lasagna/sirui/damage/xcell/blood_xcell2.txt',
  quote = FALSE,
  sep = '\t',
  row.names = FALSE
)

AA2 <- {}

for (i in 2:7) {
  d_cor <- as.data.frame(ddply(blood_xcell, .(tissue), function(xx) {
    data.frame(
      COR = cor(xx$AGE, xx[, i]),
      P = cor.test(xx$AGE, xx[, i])$p.value
    )
  }))
  d_cor$cell <- colnames(blood_xcell)[i]
  AA2 <- rbind(AA2, d_cor)
}

AA2$COR <- as.numeric(AA2$COR)
AA2 <- AA2 %>% mutate(significance = ifelse(P < 0.05, "*", ""))

ggplot(AA2, aes(x = cell, y = tissue, fill = COR)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "blue", high = "red", mid = "white",
    midpoint = 0, name = "Cor", na.value = "white"
  ) +
  geom_text(aes(label = significance), color = "black", size = 5) +
  labs(title = "cor_with_damage", x = "", y = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

tissues <- c('Spleen', 'Marrow', 'WBC', 'Bone')

A <- {}

for (j_ in 1:4) {
  tissue_ <- tissues[j_]
  print(tissue_)

  cell_temp <- subset(blood_xcell, tissue == tissue_)
  cell_temp[, 2:7] <- scale(cell_temp[, 2:7])

  temp <- c()

  for (i in 2:7) {
    mod2 <- lm(cell_temp$all ~ cell_temp[, i] + cell_temp$AGE, data = cell_temp)

    p_age <- summary(mod2)$coefficients[2, 4]
    b_age <- summary(mod2)$coefficients[2, 1]

    temp1 <- c(colnames(cell_temp)[i], b_age, p_age)
    temp <- rbind(temp, temp1)
  }

  temp <- as.data.frame(temp)
  colnames(temp) <- c('cell', 'b', 'P')
  temp$tissue <- tissue_

  A <- rbind(A, temp)
}

A <- as.data.frame(A)
A$b <- as.numeric(A$b)

A <- A %>% mutate(significance = ifelse(as.numeric(P) < 0.05, "*", ""))

ggplot(A, aes(x = cell, y = tissue, fill = b)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "blue", high = "red", mid = "white",
    midpoint = 0, name = "Cor"
  ) +
  geom_text(aes(label = significance), size = 5) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

library(zellkonverter)

sce <- readRDS("/Wafer/siruizhang/damage/mouse_singlecell/facs.normalized.Liver.rds")
count_matrix <- as.data.frame(assay(sce, "X"))

dice_labels <- as.data.frame(sce@colData)

dice_labels$ont <- dice_labels$cell_ontology_id
dice_labels$label <- dice_labels$cell_ontology_class
dice_labels$sample <- rownames(dice_labels)
dice_labels$dataset <- "TMS_liver"

dice_labels <- dice_labels[, 15:18]

dice_labels$ont <- as.character(dice_labels$ont)
dice_labels$label <- as.character(dice_labels$label)

dice_labels$ont[dice_labels$ont == 'nan'] <- NA

set.seed(123)

count_matrix <- as.matrix(count_matrix)

Liver.xCell2Ref <- xCell2::xCell2Train(
  ref = count_matrix,
  labels = dice_labels,
  refType = "rnaseq"
)
