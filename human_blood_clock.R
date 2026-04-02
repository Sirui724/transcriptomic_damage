# =========================================================
# Aging clock pipeline (multi-omics integration)
# =========================================================

library(data.table)
library(dplyr)
library(purrr)
library(glmnet)
library(lightgbm)
library(ggplot2)

# =========================================================
# 0. Helper functions
# =========================================================

rmse <- function(a, b) sqrt(mean((a - b)^2, na.rm = TRUE))
mae  <- function(a, b) mean(abs(a - b), na.rm = TRUE)

make_folds <- function(n, K = 5, seed = 1) {
  set.seed(seed)
  sample(rep(1:K, length.out = n))
}

logit_transform <- function(M, eps = 1e-6) {
  M <- pmin(pmax(M, eps), 1 - eps)
  log(M / (1 - M))
}

# =========================================================
# 1. Data loading (统一格式)
# =========================================================

load_dataset <- function(file, info) {
  df <- read.table(file, header = TRUE, sep = "\t")
  df <- as.data.frame(t(df))
  colnames(df) <- df[1,]
  df <- df[-1,]
  df$sample <- rownames(df)
  merge(df, info, by = "sample")
}

# metadata
d_info <- read.table(
  "/Entropy/siruizhang/damage_pipeline/human/out_from_eristwo/data_info.txt",
  header = TRUE, sep = "\t"
)[, c(2,7)]

# load all datasets
data_stop <- load_dataset("stopsite_ratio_result.txt", d_info)
data_junc <- load_dataset("junc_gene_ratio_result.txt", d_info)
data_gf   <- load_dataset("gf_by_sample_CDS_fusion_normalized_by_expression.filtered.tsv", d_info)
data_re   <- load_dataset("normalized_merge_re.txt", d_info)

# =========================================================
# 2. Align samples
# =========================================================

prep_df <- function(df){
  df$age <- as.numeric(df$age)
  df <- df[!is.na(df$age), ]
  df
}

data_stop <- prep_df(data_stop)
data_junc <- prep_df(data_junc)
data_gf   <- prep_df(data_gf)
data_re   <- prep_df(data_re)

common_samples <- Reduce(intersect, list(
  data_stop$sample,
  data_junc$sample,
  data_gf$sample,
  data_re$sample
))

subset_common <- function(df){
  df <- df[df$sample %in% common_samples, ]
  df[match(common_samples, df$sample), ]
}

data_stop <- subset_common(data_stop)
data_junc <- subset_common(data_junc)
data_gf   <- subset_common(data_gf)
data_re   <- subset_common(data_re)

age_vec <- data_stop$age

# =========================================================
# 3. Sub-clock (LightGBM)
# =========================================================

fit_subclock <- function(df, train_samples, test_samples,
                         transform = "log1p", prev_min = 0.05) {

  sample_col <- "sample"
  age_col    <- "age"

  feature_cols <- setdiff(colnames(df), c(sample_col, age_col))

  df_tr <- df[df$sample %in% train_samples, ]
  df_te <- df[df$sample %in% test_samples, ]

  Xtr <- as.matrix(df_tr[, feature_cols])
  Xte <- as.matrix(df_te[, feature_cols])

  Xtr[is.na(Xtr)] <- 0
  Xte[is.na(Xte)] <- 0

  # prevalence filter
  keep <- colMeans(Xtr != 0) >= prev_min
  Xtr <- Xtr[, keep]
  Xte <- Xte[, keep]

  # transform
  if (transform == "log1p") {
    Xtr <- log1p(Xtr)
    Xte <- log1p(Xte)
  } else if (transform == "logit") {
    Xtr <- logit_transform(Xtr)
    Xte <- logit_transform(Xte)
  }

  ytr <- df_tr$age
  yte <- df_te$age

  dtrain <- lgb.Dataset(Xtr, label = ytr)

  model <- lgb.train(
    params = list(objective = "regression"),
    data = dtrain,
    nrounds = 500
  )

  pred_tr <- predict(model, Xtr)
  pred_te <- predict(model, Xte)

  list(
    train = data.frame(sample = df_tr$sample, age = ytr, pred = pred_tr),
    test  = data.frame(sample = df_te$sample, age = yte, pred = pred_te)
  )
}

# =========================================================
# 4. Train/test split
# =========================================================

set.seed(1)
train_samples <- sample(common_samples, 0.8 * length(common_samples))
test_samples  <- setdiff(common_samples, train_samples)

# =========================================================
# 5. Train 4 sub-clocks
# =========================================================

fit_re   <- fit_subclock(data_re,   train_samples, test_samples, "log1p")
fit_gf   <- fit_subclock(data_gf,   train_samples, test_samples, "logit")
fit_stop <- fit_subclock(data_stop, train_samples, test_samples, "logit")
fit_junc <- fit_subclock(data_junc, train_samples, test_samples, "logit")

# =========================================================
# 6. Merge sub-clock outputs
# =========================================================

merge_scores <- function(lst){
  reduce(lst, left_join, by = c("sample","age"))
}

sub_train <- merge_scores(list(
  rename(fit_re$train, score_re = pred),
  rename(fit_gf$train, score_gf = pred),
  rename(fit_stop$train, score_stop = pred),
  rename(fit_junc$train, score_junc = pred)
))

sub_test <- merge_scores(list(
  rename(fit_re$test, score_re = pred),
  rename(fit_gf$test, score_gf = pred),
  rename(fit_stop$test, score_stop = pred),
  rename(fit_junc$test, score_junc = pred)
))

# =========================================================
# 7. Final clock (Ridge)
# =========================================================

Xtr <- as.matrix(sub_train[, -c(1,2)])
ytr <- sub_train$age

Xte <- as.matrix(sub_test[, -c(1,2)])
yte <- sub_test$age

fit_final <- cv.glmnet(Xtr, ytr, alpha = 0)

pred_test <- predict(fit_final, Xte, s = "lambda.min")

cat(sprintf("Test RMSE=%.3f, cor=%.3f\n",
            rmse(yte, pred_test),
            cor(yte, pred_test)))

# =========================================================
# 8. Visualization
# =========================================================

plot_df <- data.frame(
  age = yte,
  pred = as.numeric(pred_test)
)

ggplot(plot_df, aes(age, pred)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm") +
  theme_classic()
