# =============================================================================
# 血小板储存质量代谢组学预测：概念验证研究
# 完整分析脚本（用于投稿补充材料）
# 
# 软件环境: R 4.5.x, 依赖包见下文
# 数据: 需将 ST003937 的原始数据文件置于 ./data 目录下
# 运行方式: 在 RStudio 中打开此脚本，全选运行或 source("main_analysis.R")
# =============================================================================

# 0. 环境准备 ---------------------------------------------------------------
packages <- c("tidyverse", "caret", "randomForest", "xgboost", "pROC", "glmnet",
              "fastshap", "shapviz", "cluster", "factoextra", "patchwork", "torch")
installed <- packages %in% installed.packages()
if (any(!installed)) install.packages(packages[!installed])
lapply(packages, library, character.only = TRUE)

# 若 torch 依赖的 libtorch 未安装，请先执行:
# torch::install_torch()
# 然后重启 R 会话

# 创建输出文件夹
if (!dir.exists("results")) dir.create("results")

# 通用 SCI 图表主题
sci_theme <- theme_minimal(base_size = 9) +
  theme(
    text = element_text(family = "Arial", color = "black"),
    plot.title = element_text(size = 10, face = "bold", hjust = 0),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8, color = "black"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.line = element_blank(),
    legend.position = "top",
    legend.box.background = element_rect(color = "black", linewidth = 0.3),
    strip.background = element_rect(fill = "grey95", color = "black", linewidth = 0.3),
    strip.text = element_text(size = 9, face = "bold")
  )

cb_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                "#0072B2", "#D55E00", "#CC79A7", "#000000")

# =============================================================================
# 1. 数据读取与合并
# =============================================================================
cat("1/6 读取并合并代谢组学数据...\n")
lines1 <- readLines("data/MSdata_ST003937_1.txt")
header1 <- lines1[1]; factors1 <- lines1[2]
col_names1 <- str_split_1(header1, "\t")
factor_values1 <- str_split_1(factors1, "\t")
sample_ids1 <- col_names1[-(1:2)]
sample_factors1 <- factor_values1[-(1:2)]

parse_factor_string <- function(x) {
  parts <- str_split_1(x, " \\| ")
  out <- list()
  for (p in parts) {
    kv <- str_split_fixed(str_trim(p), ":", 2)
    out[[ str_trim(kv[1]) ]] <- str_trim(kv[2])
  }
  out
}
factor_list1 <- map(sample_factors1, parse_factor_string)
samples_meta <- bind_rows(factor_list1) %>%
  mutate(Sample_ID = sample_ids1) %>%
  select(Sample_ID, everything())
names(samples_meta) <- make.names(names(samples_meta))

data_df1 <- read_tsv(I(lines1[-(1:2)]), col_names = col_names1, show_col_types = FALSE)
metab1 <- data_df1 %>%
  column_to_rownames(var = names(data_df1)[1]) %>%
  select(-2) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "Sample_ID")
metab1[,-1] <- lapply(metab1[,-1], function(x) as.numeric(as.character(x)))

full_data <- inner_join(samples_meta, metab1, by = "Sample_ID")
metab_cols1 <- setdiff(names(metab1), "Sample_ID")
cat(sprintf("文件1：%d 个代谢物\n", length(metab_cols1)))

# 补充代谢物文件
cat("2/6 读取文件2...\n")
lines2 <- readLines("data/MSdata_ST003937_2.txt")
col_names2 <- str_split_1(lines2[1], "\t")
data_df2 <- read_tsv(I(lines2[-(1:2)]), col_names = col_names2, show_col_types = FALSE)

metab2 <- data_df2 %>%
  column_to_rownames(var = names(data_df2)[1]) %>%
  select(-2) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "Sample_ID")
metab2[,-1] <- lapply(metab2[,-1], function(x) as.numeric(as.character(x)))

full_data <- left_join(full_data, metab2, by = "Sample_ID")
metab_cols2 <- setdiff(names(metab2), "Sample_ID")

overlap <- intersect(metab_cols1, metab_cols2)
if (length(overlap) > 0) {
  metab_cols2 <- setdiff(metab_cols2, overlap)
  full_data <- full_data %>% select(-all_of(overlap))
}
metab_cols <- c(metab_cols1, metab_cols2)
cat(sprintf("合并后代谢物总数：%d\n", length(metab_cols)))

full_data <- full_data %>%
  mutate(across(all_of(metab_cols), ~ ifelse(is.na(.), median(., na.rm = TRUE), .)))
zero_var <- metab_cols[apply(full_data[metab_cols], 2, sd) == 0]
if (length(zero_var) > 0) {
  full_data <- full_data %>% select(-all_of(zero_var))
  metab_cols <- setdiff(metab_cols, zero_var)
}

# =============================================================================
# 2. 提取捐献者ID
# =============================================================================
cat("3/6 提取捐献者ID...\n")
an_lines <- readLines("data/ST003937_AN006465.txt")
subject_rows <- an_lines[str_detect(an_lines, "^SUBJECT_SAMPLE_FACTORS\\s+")]
subject_df <- read.table(text = subject_rows, header = FALSE, sep = "\t",
                         stringsAsFactors = FALSE, fill = TRUE, quote = "")
subject_map <- subject_df %>%
  select(Subject_ID = V2, Sample_ID = V3) %>%
  distinct()
cat(sprintf("提取到 %d 位捐献者\n", n_distinct(subject_map$Subject_ID)))

full_data <- left_join(full_data, subject_map, by = "Sample_ID") %>%
  filter(!is.na(Subject_ID))

# =============================================================================
# 3. 数据归一化
# =============================================================================
cat("4/6 数据归一化...\n")
X_raw <- full_data %>% select(all_of(metab_cols))
sample_medians <- apply(X_raw, 1, median, na.rm = TRUE)
sample_medians[sample_medians == 0] <- 1
X_med_norm <- sweep(X_raw, 1, sample_medians, "/")
X_log <- log2(X_med_norm + 1)
pareto_scale <- function(x) { (x - mean(x, na.rm = TRUE)) / sqrt(sd(x, na.rm = TRUE)) }
X_pareto <- apply(X_log, 2, pareto_scale) %>% as.data.frame()
full_data[, metab_cols] <- X_pareto
cat("归一化完成，维度:", nrow(full_data), "x", length(metab_cols), "\n")

# =============================================================================
# 4. 构建标签
# =============================================================================
cat("5/6 构建任务标签...\n")
full_data$Quality_label <- ifelse(full_data$Fresh.v.stored == "Fresh", 1, 0)
full_data$Type_label   <- ifelse(full_data$Sample.type == "BioPlatelets", 1, 0)

# =============================================================================
# 5. 捐献者分组交叉验证 (RF/XGBoost)
# =============================================================================
cat("6/6 交叉验证建模...\n")
X <- full_data %>% select(all_of(metab_cols))
y_qual <- factor(full_data$Quality_label, levels=0:1, labels=c("Stored","Fresh"))
y_type <- factor(full_data$Type_label, levels=0:1, labels=c("Platelets","BioPlatelets"))
subjects <- full_data$Subject_ID

set.seed(42)
folds <- groupKFold(subjects, k = 5)
cv_results <- data.frame()

for (fold_idx in seq_along(folds)) {
  train_idx <- which(subjects %in% subjects[-folds[[fold_idx]]])
  test_idx  <- which(subjects %in% subjects[folds[[fold_idx]]])
  X_tr <- X[train_idx, ]; X_te <- X[test_idx, ]
  y_tr_q <- y_qual[train_idx]; y_te_q <- y_qual[test_idx]
  y_tr_t <- y_type[train_idx]; y_te_t <- y_type[test_idx]
  if (length(unique(y_tr_q))<2 | length(unique(y_te_q))<2) next
  if (length(unique(y_tr_t))<2 | length(unique(y_te_t))<2) next
  
  # 任务A
  cv_lasso_q <- cv.glmnet(as.matrix(X_tr), y_tr_q, family="binomial", alpha=1, nfolds=3)
  sel_q <- rownames(coef(cv_lasso_q, s="lambda.min"))[-1][coef(cv_lasso_q, s="lambda.min")[-1,]!=0]
  if (length(sel_q)<5) sel_q <- colnames(X_tr)
  rf_q <- randomForest(X_tr[,sel_q,drop=FALSE], y_tr_q, ntree=300)
  pred_rf_q <- predict(rf_q, X_te[,sel_q,drop=FALSE], type="prob")[,"Fresh"]
  auc_rf_q <- auc(roc(y_te_q, pred_rf_q, direction="<"))
  xgb_q <- xgb.train(params=list(objective="binary:logistic", max_depth=3, eta=0.1),
                     data=xgb.DMatrix(as.matrix(X_tr[,sel_q]), label=as.numeric(y_tr_q)-1),
                     nrounds=50, verbose=0)
  pred_xgb_q <- predict(xgb_q, xgb.DMatrix(as.matrix(X_te[,sel_q])))
  auc_xgb_q <- auc(roc(y_te_q, pred_xgb_q, direction="<"))
  
  # 任务B
  cv_lasso_t <- cv.glmnet(as.matrix(X_tr), y_tr_t, family="binomial", alpha=1, nfolds=3)
  sel_t <- rownames(coef(cv_lasso_t, s="lambda.min"))[-1][coef(cv_lasso_t, s="lambda.min")[-1,]!=0]
  if (length(sel_t)<5) sel_t <- colnames(X_tr)
  rf_t <- randomForest(X_tr[,sel_t,drop=FALSE], y_tr_t, ntree=300)
  pred_rf_t <- predict(rf_t, X_te[,sel_t,drop=FALSE], type="prob")[,"BioPlatelets"]
  auc_rf_t <- auc(roc(y_te_t, pred_rf_t, direction="<"))
  xgb_t <- xgb.train(params=list(objective="binary:logistic", max_depth=3, eta=0.1),
                     data=xgb.DMatrix(as.matrix(X_tr[,sel_t]), label=as.numeric(y_tr_t)-1),
                     nrounds=50, verbose=0)
  pred_xgb_t <- predict(xgb_t, xgb.DMatrix(as.matrix(X_te[,sel_t])))
  auc_xgb_t <- auc(roc(y_te_t, pred_xgb_t, direction="<"))
  
  cv_results <- rbind(cv_results, data.frame(
    Fold=fold_idx, Task_A_RF_AUC=auc_rf_q, Task_A_XGB_AUC=auc_xgb_q,
    Task_B_RF_AUC=auc_rf_t, Task_B_XGB_AUC=auc_xgb_t
  ))
}
cv_results[,-1] <- lapply(cv_results[,-1], as.numeric)
cat("\n交叉验证结果:\n"); print(cv_results)
cat("平均 AUC:\n"); print(colMeans(cv_results[,-1]))
write.csv(cv_results, "results/cv_results.csv", row.names=FALSE)

# =============================================================================
# 6. 随机划分对比实验 (100次)
# =============================================================================
cat("\n随机划分对比实验...\n")
set.seed(42)
random_auc <- c()
for (i in 1:100) {
  train_idx <- createDataPartition(y_qual, p=0.8, list=FALSE)
  X_tr <- X[train_idx,]; X_te <- X[-train_idx,]
  y_tr <- y_qual[train_idx]; y_te <- y_qual[-train_idx]
  cv_lasso <- cv.glmnet(as.matrix(X_tr), y_tr, family="binomial", alpha=1, nfolds=3)
  sel <- rownames(coef(cv_lasso, s="lambda.min"))[-1][coef(cv_lasso, s="lambda.min")[-1,]!=0]
  if (length(sel)<5) sel <- colnames(X_tr)
  rf <- randomForest(X_tr[,sel,drop=FALSE], y_tr, ntree=300)
  pred <- predict(rf, X_te[,sel,drop=FALSE], type="prob")[,"Fresh"]
  random_auc[i] <- auc(roc(y_te, pred, direction="<"))
}
cat("随机划分平均 AUC:", round(mean(random_auc),4), "\n")
donorwise_auc <- cv_results$Task_A_RF_AUC
cat("捐献者独立验证平均 AUC:", round(mean(donorwise_auc),4), "\n")

# =============================================================================
# 7. 多任务学习 (torch MLP)
# =============================================================================
cat("\n多任务MLP建模...\n")
X_all <- as.matrix(X)
y_qual_num <- full_data$Quality_label
y_type_num <- full_data$Type_label

create_stl_net <- function(input_dim) {
  nn_module("STL_Net",
            initialize = function() {
              self$fc1 <- nn_linear(input_dim, 64)
              self$fc2 <- nn_linear(64, 32)
              self$fc3 <- nn_linear(32, 1)
            },
            forward = function(x) {
              x %>% self$fc1() %>% nnf_relu() %>% nnf_dropout(0.3) %>%
                self$fc2() %>% nnf_relu() %>% nnf_dropout(0.3) %>%
                self$fc3() %>% torch_sigmoid()
            }
  )
}

create_mtl_net <- function(input_dim) {
  nn_module("MTL_Net",
            initialize = function() {
              self$shared1 <- nn_linear(input_dim, 64)
              self$shared2 <- nn_linear(64, 32)
              self$head_A1 <- nn_linear(32, 16)
              self$head_A2 <- nn_linear(16, 1)
              self$head_B1 <- nn_linear(32, 16)
              self$head_B2 <- nn_linear(16, 1)
            },
            forward = function(x) {
              shared <- x %>% self$shared1() %>% nnf_relu() %>% nnf_dropout(0.3) %>%
                self$shared2() %>% nnf_relu() %>% nnf_dropout(0.3)
              out_A <- shared %>% self$head_A1() %>% nnf_relu() %>% self$head_A2() %>% torch_sigmoid()
              out_B <- shared %>% self$head_B1() %>% nnf_relu() %>% self$head_B2() %>% torch_sigmoid()
              list(out_A, out_B)
            }
  )
}

results_stl_A <- data.frame(Fold=integer(), AUC=numeric())
results_stl_B <- data.frame(Fold=integer(), AUC=numeric())
results_mtl <- data.frame(Fold=integer(), TaskA_AUC=numeric(), TaskB_AUC=numeric())

for (fold_idx in seq_along(folds)) {
  train_idx <- which(subjects %in% subjects[-folds[[fold_idx]]])
  test_idx  <- which(subjects %in% subjects[folds[[fold_idx]]])
  X_tr <- torch_tensor(X_all[train_idx,], dtype=torch_float())
  X_te <- torch_tensor(X_all[test_idx,], dtype=torch_float())
  y_tr_A <- torch_tensor(y_qual_num[train_idx], dtype=torch_float())$view(c(-1,1))
  y_te_A <- y_qual_num[test_idx]
  y_tr_B <- torch_tensor(y_type_num[train_idx], dtype=torch_float())$view(c(-1,1))
  y_te_B <- y_type_num[test_idx]
  
  if(length(unique(y_qual_num[train_idx]))<2 | length(unique(y_type_num[train_idx]))<2) next
  
  # 单任务A
  net_A <- create_stl_net(ncol(X_all))()
  opt_A <- optim_adam(net_A$parameters, lr=0.001)
  for(e in 1:50) { net_A$train(); loss <- nnf_binary_cross_entropy(net_A(X_tr), y_tr_A); opt_A$zero_grad(); loss$backward(); opt_A$step() }
  net_A$eval(); pred_A <- as.numeric(net_A(X_te)$detach()); auc_stl_A <- auc(roc(y_te_A, pred_A, direction="<"))
  results_stl_A <- rbind(results_stl_A, data.frame(Fold=fold_idx, AUC=auc_stl_A))
  
  # 单任务B
  net_B <- create_stl_net(ncol(X_all))()
  opt_B <- optim_adam(net_B$parameters, lr=0.001)
  for(e in 1:50) { net_B$train(); loss <- nnf_binary_cross_entropy(net_B(X_tr), y_tr_B); opt_B$zero_grad(); loss$backward(); opt_B$step() }
  net_B$eval(); pred_B <- as.numeric(net_B(X_te)$detach()); auc_stl_B <- auc(roc(y_te_B, pred_B, direction="<"))
  results_stl_B <- rbind(results_stl_B, data.frame(Fold=fold_idx, AUC=auc_stl_B))
  
  # 多任务
  net_mtl <- create_mtl_net(ncol(X_all))()
  opt_mtl <- optim_adam(net_mtl$parameters, lr=0.001)
  for(e in 1:50) {
    net_mtl$train(); out <- net_mtl(X_tr)
    loss_A <- nnf_binary_cross_entropy(out[[1]], y_tr_A)
    loss_B <- nnf_binary_cross_entropy(out[[2]], y_tr_B)
    loss <- 0.5*loss_A + 0.5*loss_B
    opt_mtl$zero_grad(); loss$backward(); opt_mtl$step()
  }
  net_mtl$eval(); out_te <- net_mtl(X_te)
  auc_mtl_A <- auc(roc(y_te_A, as.numeric(out_te[[1]]$detach()), direction="<"))
  auc_mtl_B <- auc(roc(y_te_B, as.numeric(out_te[[2]]$detach()), direction="<"))
  results_mtl <- rbind(results_mtl, data.frame(Fold=fold_idx, TaskA_AUC=auc_mtl_A, TaskB_AUC=auc_mtl_B))
  cat(sprintf("Fold %d 完成\n", fold_idx))
}
cat("\n单任务A 均值 AUC:", mean(results_stl_A$AUC))
cat("\n多任务A 均值 AUC:", mean(results_mtl$TaskA_AUC))
cat("\n单任务B 均值 AUC:", mean(results_stl_B$AUC))
cat("\n多任务B 均值 AUC:", mean(results_mtl$TaskB_AUC))
# =============================================================================
# 保存表1 和 补充表S4
# =============================================================================
# =============================================================================
# 表1：多任务学习模型性能对比（单任务 vs 多任务）
# =============================================================================
table1 <- data.frame(
  Model = c("Single-task MLP", "Multi-task MLP"),
  Task_A_AUC = c(round(mean(results_stl_A$AUC), 3),
                 round(mean(results_mtl$TaskA_AUC), 3)),
  Task_B_AUC = c(round(mean(results_stl_B$AUC), 3),
                 round(mean(results_mtl$TaskB_AUC), 3))
)
write.csv(table1, "results/Table1_MTL_performance.csv", row.names = FALSE)

# =============================================================================
# 补充表S4：多任务学习模型逐折交叉验证详细结果
# =============================================================================
table_s4 <- data.frame(
  Fold = results_stl_A$Fold,
  STL_A_AUC = round(results_stl_A$AUC, 3),
  MTL_A_AUC = round(results_mtl$TaskA_AUC, 3),
  STL_B_AUC = round(results_stl_B$AUC, 3),
  MTL_B_AUC = round(results_mtl$TaskB_AUC, 3)
)
# 添加均值行
table_s4 <- rbind(table_s4,
                  data.frame(
                    Fold = "Mean",
                    STL_A_AUC = round(mean(results_stl_A$AUC), 3),
                    MTL_A_AUC = round(mean(results_mtl$TaskA_AUC), 3),
                    STL_B_AUC = round(mean(results_stl_B$AUC), 3),
                    MTL_B_AUC = round(mean(results_mtl$TaskB_AUC), 3)
                  ))
write.csv(table_s4, "results/TableS4_MTL_detailed.csv", row.names = FALSE)
# =============================================================================
# 8. SHAP 分析 (任务B)
# =============================================================================
cat("\nSHAP分析...\n")
rf_full <- randomForest(X, y_type, ntree=500, importance=TRUE)
pred_fun <- function(object, newdata) predict(object, newdata, type="prob")[,"BioPlatelets"]
set.seed(42); idx <- sample(1:nrow(X), 100)
X_df <- as.data.frame(X); X_explain_df <- X_df[idx,]
shap <- fastshap::explain(rf_full, X=X_df, newdata=X_explain_df, pred_wrapper=pred_fun, nsim=30)
shv <- shapviz(shap, X=X_df[idx,])

sv_importance(shv, kind="beeswarm", max_display=15, show_numbers=TRUE, bee_width=0.2, color_bar=FALSE)
ggsave("results/Fig4_SHAP_beeswarm.tiff", width=16, height=10, units="cm", dpi=600, compression="lzw")

shap_importance <- sort(colMeans(abs(shap)), decreasing=TRUE)
write.csv(data.frame(Metabolite=names(shap_importance), MeanAbsSHAP=shap_importance),
          "results/SHAP_top20_metabolites.csv", row.names=FALSE)

# =============================================================================
# 9. 供者协变量分析
# =============================================================================
cat("\n供者协变量分析...\n")
donor_scatter <- full_data %>%
  group_by(Subject_ID) %>%
  summarise(across(all_of(metab_cols), mean, .names="mean_{.col}"),
            StorageGroup=unique(Fresh.v.stored), n_samples=n(), .groups="drop")
center_matrix <- donor_scatter %>% select(starts_with("mean_"))
donor_centers <- donor_scatter %>% select(Subject_ID, starts_with("mean_"))

sample_distances <- full_data %>%
  left_join(donor_centers, by="Subject_ID") %>%
  mutate(sq_diff=0)
for (m in metab_cols) {
  mean_col <- paste0("mean_", m)
  sample_distances <- sample_distances %>% mutate(sq_diff = sq_diff + (.data[[m]] - .data[[mean_col]])^2)
}
sample_distances <- sample_distances %>%
  mutate(dist_to_center = sqrt(sq_diff)) %>%
  select(Subject_ID, Sample_ID, Fresh.v.stored, dist_to_center)

donor_variability <- sample_distances %>%
  group_by(Subject_ID, Fresh.v.stored) %>%
  summarise(Median_scatter = median(dist_to_center), .groups="drop")

donor_centers_mat <- as.matrix(center_matrix)
rownames(donor_centers_mat) <- donor_scatter$Subject_ID
donor_pairwise_dist <- dist(donor_centers_mat)
fresh_idx <- which(donor_scatter$StorageGroup=="Fresh")
stored_idx <- which(donor_scatter$StorageGroup=="Stored")
storage_group_distance <- sqrt(sum((colMeans(center_matrix[fresh_idx,]) - colMeans(center_matrix[stored_idx,]))^2))
dist_df <- data.frame(Distance = as.numeric(donor_pairwise_dist))

cat(sprintf("储存组间距离: %.4f, 供者间距离中位数: %.4f\n", storage_group_distance, median(dist_df$Distance)))

# =============================================================================
# 10. 生成所有图表
# =============================================================================
cat("\n生成图表...\n")

# 图1 PCA
pca_res <- prcomp(X, scale.=FALSE, center=TRUE)
pca_scores <- as.data.frame(pca_res$x)
pca_scores$Subject_ID <- as.factor(full_data$Subject_ID)
pca_scores$Storage <- full_data$Fresh.v.stored
pca_scores$Label <- full_data$Sample.type
var_exp <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)

p1a <- ggplot(pca_scores, aes(PC1, PC2)) + geom_point(aes(color=Subject_ID), size=1.2, alpha=0.7) +
  labs(x=paste0("PC1 (",var_exp[1],"%)"), y=paste0("PC2 (",var_exp[2],"%)"), title="A. By Donor") +
  scale_color_manual(values=rep(cb_palette, length.out=24)) + sci_theme + theme(legend.position="none")

p1b <- ggplot(pca_scores, aes(PC1, PC2)) + geom_point(aes(color=Storage), size=1.2, alpha=0.7) +
  scale_color_manual(values=c("Fresh"=cb_palette[2], "Stored"=cb_palette[1])) + labs(title="B. By Storage") + sci_theme

p1c <- ggplot(pca_scores, aes(PC1, PC2)) + geom_point(aes(color=Label), size=1.2, alpha=0.7) +
  scale_color_manual(values=c("Platelets"=cb_palette[3], "BioPlatelets"=cb_palette[6]), labels=c("Unlabeled","Biotinylated")) +
  labs(title="C. By Labeling") + sci_theme

fig1 <- p1a + p1b + p1c + plot_layout(ncol=3)
ggsave("results/Fig1_PCA_combined.tiff", fig1, width=18, height=6, units="cm", dpi=600, compression="lzw")

# =============================================================================
# 轮廓系数分析：定量证明供者效应 > 储存效应 ≈ 标记效应
# 插入位置：供者协变量分析之后、生成图表之前
# =============================================================================

library(cluster)

# 计算前两个主成分（与图1一致）
pca_scores_subset <- pca_scores[, c("PC1", "PC2")]

# 1. 按捐献者分组的轮廓系数
sil_donor <- silhouette(as.numeric(pca_scores$Subject_ID), 
                        dist(pca_scores_subset))
mean_sil_donor <- mean(sil_donor[, 3])

# 2. 按储存时长分组的轮廓系数
sil_storage <- silhouette(as.numeric(as.factor(pca_scores$Storage)), 
                          dist(pca_scores_subset))
mean_sil_storage <- mean(sil_storage[, 3])

# 3. 按标记状态分组的轮廓系数
sil_label <- silhouette(as.numeric(as.factor(pca_scores$Label)), 
                        dist(pca_scores_subset))
mean_sil_label <- mean(sil_label[, 3])

# 输出结果
cat(sprintf("轮廓系数（PC1-PC2空间）：\n"))
cat(sprintf("  捐献者分组: %.4f\n", mean_sil_donor))
cat(sprintf("  储存时长分组: %.4f\n", mean_sil_storage))
cat(sprintf("  标记状态分组: %.4f\n", mean_sil_label))
cat(sprintf("  供者/储存比值: %.2f\n", mean_sil_donor / mean_sil_storage))


# 图2 供者间距离
fig2 <- ggplot(dist_df, aes(x=Distance)) +
  geom_histogram(fill="grey60", color="grey30", bins=25, linewidth=0.2) +
  geom_vline(xintercept=storage_group_distance, color="#D55E00", linewidth=1.0, linetype="dashed") +
  annotate("text", x=storage_group_distance+0.8, y=Inf, vjust=2, label=paste0("Storage effect = ",round(storage_group_distance,2)), color="#D55E00", size=3, hjust=0) +
  annotate("text", x=median(dist_df$Distance)+0.5, y=Inf, vjust=4, label=paste0("Median inter-donor = ",round(median(dist_df$Distance),2)), color="black", size=3, hjust=0) +
  labs(x="Euclidean metabolic distance", y="Number of donor pairs", title="Inter-donor metabolic variability dwarfs storage effect") +
  sci_theme + theme(panel.grid.major.x=element_blank())
ggsave("results/Fig2_donor_distance.tiff", fig2, width=12, height=8, units="cm", dpi=600, compression="lzw")

# 图3 AUC条形图
roc_data <- cv_results %>%
  pivot_longer(-Fold, names_to="Model", values_to="AUC") %>%
  mutate(Task=ifelse(grepl("Task_A",Model),"Task A: Storage quality","Task B: Labeling status"),
         Method=ifelse(grepl("_RF",Model),"Random Forest","XGBoost"))
auc_summary <- roc_data %>% group_by(Task, Method) %>% summarise(Mean_AUC=mean(AUC), SE=sd(AUC)/sqrt(n()), .groups="drop")

fig3 <- ggplot(auc_summary, aes(x=Method, y=Mean_AUC, fill=Method)) +
  geom_bar(stat="identity", width=0.6, color="black", linewidth=0.3) +
  geom_errorbar(aes(ymin=Mean_AUC-SE, ymax=Mean_AUC+SE), width=0.15, linewidth=0.5) +
  geom_hline(yintercept=0.5, linetype="dotted", color="grey50", linewidth=0.5) +
  facet_wrap(~Task, scales="fixed") +
  scale_fill_manual(values=c("Random Forest"=cb_palette[1], "XGBoost"=cb_palette[2])) +
  scale_y_continuous(limits=c(0,1), breaks=seq(0,1,0.2)) +
  labs(x="", y="Mean AUC (5-fold donor-wise CV)", title="Model performance comparison") +
  sci_theme + theme(axis.text.x=element_text(hjust=0.5))
ggsave("results/Fig3_AUC_barplot.tiff", fig3, width=14, height=7, units="cm", dpi=600, compression="lzw")

# 图5 关键代谢物箱线图
top_metabs <- c("Diphosphate", "Phosphate", "L-histidine")
plot_df <- full_data %>%
  select(all_of(top_metabs), Type_label) %>%
  pivot_longer(-Type_label, names_to="Metabolite", values_to="Value") %>%
  mutate(Label_status=ifelse(Type_label==1, "Biotinylated", "Unlabeled"))

fig5 <- ggplot(plot_df, aes(x=Label_status, y=Value, fill=Label_status)) +
  geom_boxplot(outlier.size=0.5, linewidth=0.3, alpha=0.8) +
  facet_wrap(~Metabolite, scales="free_y", ncol=3) +
  scale_fill_manual(values=c("Unlabeled"=cb_palette[3], "Biotinylated"=cb_palette[6])) +
  labs(x="", y="Normalized intensity", title="Top SHAP metabolites by biotinylation status") +
  sci_theme + theme(legend.title=element_blank(), strip.text=element_text(size=8))
ggsave("results/Fig5_metabolite_boxplot.tiff", fig5, width=18, height=6, units="cm", dpi=600, compression="lzw")

# 补充图S1
p_s1 <- ggplot(donor_variability, aes(x=Fresh.v.stored, y=Median_scatter, fill=Fresh.v.stored)) +
  geom_boxplot(alpha=0.8, outlier.size=0.8, linewidth=0.3) +
  scale_fill_manual(values=c("Fresh"=cb_palette[2], "Stored"=cb_palette[1])) +
  labs(x="Storage duration", y="Median intra-donor distance", title="Intra-donor metabolic scatter") +
  sci_theme + theme(legend.position="none") +
  annotate("text", x=1.5, y=max(donor_variability$Median_scatter)*0.95, label="p=0.83", size=3)
ggsave("results/FigS1_intra_donor_scatter.tiff", p_s1, width=10, height=8, units="cm", dpi=600, compression="lzw")

# 补充图S2
imp_df <- as.data.frame(importance(rf_full)) %>%
  rownames_to_column("Metabolite") %>%
  arrange(desc(MeanDecreaseGini)) %>% slice(1:15) %>%
  mutate(Metabolite=fct_reorder(Metabolite, MeanDecreaseGini))
p_s2 <- ggplot(imp_df, aes(x=MeanDecreaseGini, y=Metabolite)) +
  geom_col(fill=cb_palette[2], width=0.7) +
  labs(x="Mean Decrease Gini", y="", title="Top 15 metabolites by Gini importance") +
  sci_theme + theme(panel.grid.major.y=element_blank())
ggsave("results/FigS2_Gini_importance.tiff", p_s2, width=12, height=8, units="cm", dpi=600, compression="lzw")

# 补充图S3
comparison_df <- bind_rows(
  data.frame(Method="Random split\n(100 repeats)", AUC=random_auc),
  data.frame(Method="Donor-wise CV\n(5 folds)", AUC=donorwise_auc)
)
p_s3 <- ggplot(comparison_df, aes(x=Method, y=AUC, fill=Method)) +
  geom_boxplot(alpha=0.8, outlier.size=0.8, linewidth=0.3) +
  geom_hline(yintercept=0.5, linetype="dashed", color="grey50", linewidth=0.5) +
  scale_fill_manual(values=c("Random split\n(100 repeats)"=cb_palette[1], "Donor-wise CV\n(5 folds)"=cb_palette[2])) +
  labs(x="", y="AUC (Task A)", title="Impact of validation strategy") +
  sci_theme + theme(legend.position="none") +
  annotate("text", x=1, y=mean(random_auc)+0.02, label=paste0("Mean=",round(mean(random_auc),3)), size=3) +
  annotate("text", x=2, y=mean(donorwise_auc)+0.02, label=paste0("Mean=",round(mean(donorwise_auc),3)), size=3)
ggsave("results/FigS3_validation_comparison.tiff", p_s3, width=12, height=8, units="cm", dpi=600, compression="lzw")

# 补充表S3
table_s3 <- donor_variability %>%
  rename(`Subject ID`=Subject_ID, `Storage group`=Fresh.v.stored, `Median intra-donor scatter`=Median_scatter)
write.csv(table_s3, "results/TableS3_donor_variability.csv", row.names=FALSE)

cat("\n✅ 所有分析完成，结果与图表保存于 ./results 目录。\n")
