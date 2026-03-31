# ============================================================
# ПРАКТИКА 3. Обогащение GO/KEGG + классификатор подтипов BC
# Датасет: GSE45827
# ============================================================

# ── 1. Загрузка библиотек ────────────────────────────────────
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("enrichR"), ask = FALSE)
install.packages(c("randomForest", "caret", "pROC",
                   "ggplot2", "dplyr", "gridExtra"))

library(enrichR)
library(randomForest)
library(caret)
library(pROC)
library(ggplot2)
library(dplyr)
library(gridExtra)

# ── 2. Загрузка объектов из Практики 2 ───────────────────────
load("GSE45827_modules.RData")
# Доступны: datExpr, pheno, subtype, MEs,
#           module_colors, kME, hub_genes, target_module

cat("Hub-гены:", unique(hub_genes$symbol), "\n")

# ── ЧАСТЬ A: ФУНКЦИОНАЛЬНОЕ ОБОГАЩЕНИЕ ──────────────────────

# ── 3. Подготовка списка hub-генов ───────────────────────────
hub_symbols <- unique(hub_genes$symbol)
hub_symbols <- hub_symbols[!is.na(hub_symbols) & hub_symbols != ""]
cat("Генов для обогащения:", length(hub_symbols), "\n")

# ── 4. Подключение к Enrichr ─────────────────────────────────
setEnrichrSite("Enrichr")
dbs_available <- listEnrichrDbs()

# Используемые базы
dbs <- c("GO_Biological_Process_2025", "KEGG_2026")

# ── 5. Обогащение GO Biological Process ──────────────────────
enrich_res <- enrichr(hub_symbols, dbs)

go_res <- enrich_res[["GO_Biological_Process_2025"]] |>
  filter(Adjusted.P.value < 0.05) |>
  arrange(Adjusted.P.value) |>
  head(15)

cat("Топ-15 GO BP терминов:\n")
print(go_res[, c("Term", "Overlap", "Adjusted.P.value")])

# ── 6. Обогащение KEGG ───────────────────────────────────────
kegg_res <- enrich_res[["KEGG_2021"]] |>
  filter(Adjusted.P.value < 0.05) |>
  arrange(Adjusted.P.value) |>
  head(15)

cat("Топ-15 KEGG путей:\n")
print(kegg_res[, c("Term", "Overlap", "Adjusted.P.value")])

# ── 7. Визуализация GO BP ────────────────────────────────────
# Укорачиваем длинные названия терминов
go_res$Term_short <- substr(go_res$Term, 1, 50)
go_res$log10p     <- -log10(go_res$Adjusted.P.value)

p_go <- ggplot(go_res,
               aes(x = reorder(Term_short, log10p),
                   y = log10p)) +
  geom_col(fill = "#2563eb", alpha = 0.85) +
  coord_flip() +
  labs(
    title = paste("GO Biological Process — модуль", target_module),
    x     = NULL,
    y     = "-log10(adj. p-value)"
  ) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

ggsave("go_bp_barplot.png", p_go,
       width = 10, height = 7, dpi = 150)
cat("График GO BP сохранён: go_bp_barplot.png\n")

# ── 8. Визуализация KEGG ─────────────────────────────────────
kegg_res$Term_short <- substr(kegg_res$Term, 1, 50)
kegg_res$log10p     <- -log10(kegg_res$Adjusted.P.value)

p_kegg <- ggplot(kegg_res,
                 aes(x = reorder(Term_short, log10p),
                     y = log10p)) +
  geom_col(fill = "#16a34a", alpha = 0.85) +
  coord_flip() +
  labs(
    title = paste("KEGG Pathways — модуль", target_module),
    x     = NULL,
    y     = "-log10(adj. p-value)"
  ) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

ggsave("kegg_barplot.png", p_kegg,
       width = 10, height = 7, dpi = 150)
cat("График KEGG сохранён: kegg_barplot.png\n")

# ── 9. Сохранение таблиц обогащения ─────────────────────────
write.csv(go_res,   "enrichment_GO_BP.csv",   row.names = FALSE)
write.csv(kegg_res, "enrichment_KEGG.csv",    row.names = FALSE)
cat("Таблицы обогащения сохранены\n")

# ── ЧАСТЬ B: КЛАССИФИКАТОР ПОДТИПОВ ─────────────────────────

# ── 10. Матрица признаков: экспрессия hub-генов ──────────────
# Извлекаем столбцы с hub-генами из datExpr
# datExpr содержит зонды — используем probe ID из hub_genes
hub_probes <- hub_genes$gene   # исходные probe ID

# Проверяем наличие зондов в матрице
available_probes <- hub_probes[hub_probes %in% colnames(datExpr)]
cat("Доступно hub-зондов в datExpr:", length(available_probes), "\n")

X <- datExpr[, available_probes]
colnames(X) <- hub_genes$symbol[match(available_probes, hub_genes$gene)]

# Целевая переменная
y <- factor(subtype)
cat("Распределение подтипов:\n")
print(table(y))

# Объединяем в датафрейм
df_model <- data.frame(X, subtype = y)

# ── 11. Разбивка train/test ──────────────────────────────────
set.seed(123)
train_idx <- createDataPartition(
  df_model$subtype,
  p     = 0.75,
  list  = FALSE
)

train_df <- df_model[train_idx, ]
test_df  <- df_model[-train_idx, ]

cat("Train:", nrow(train_df), "| Test:", nrow(test_df), "\n")

# ── 12. Обучение Random Forest ───────────────────────────────
set.seed(123)
rf_model <- randomForest(
  subtype ~ .,
  data     = train_df,
  ntree    = 500,
  mtry     = floor(sqrt(ncol(train_df) - 1)),
  importance = TRUE
)

print(rf_model)

# ── 13. Предсказание и оценка качества ───────────────────────
pred_class <- predict(rf_model, test_df)
pred_prob  <- predict(rf_model, test_df, type = "prob")

# Confusion matrix
cm <- confusionMatrix(pred_class, test_df$subtype)
print(cm)

# Макро-F1
macro_f1 <- mean(cm$byClass[, "F1"], na.rm = TRUE)
cat("Macro F1:", round(macro_f1, 3), "\n")

# Сохранение метрик
metrics_df <- data.frame(
  Metric   = c("Overall Accuracy", "Macro F1"),
  Value    = c(round(cm$overall["Accuracy"], 3),
               round(macro_f1, 3))
)
write.csv(metrics_df, "classifier_metrics.csv", row.names = FALSE)
cat("Метрики сохранены: classifier_metrics.csv\n")

# ── 14. Визуализация confusion matrix ───────────────────────
cm_table <- as.data.frame(cm$table)

p_cm <- ggplot(cm_table,
               aes(x = Reference, y = Prediction, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Freq), size = 5, fontface = "bold") +
  scale_fill_gradient(low = "white", high = "#2563eb") +
  labs(
    title = "Confusion Matrix — Random Forest (hub-гены)",
    x     = "Истинный подтип",
    y     = "Предсказанный подтип"
  ) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        plot.title  = element_text(face = "bold"))

ggsave("confusion_matrix.png", p_cm,
       width = 7, height = 6, dpi = 150)
cat("Confusion matrix сохранена: confusion_matrix.png\n")

# ── 15. Feature importance ───────────────────────────────────
imp_df <- importance(rf_model) |>
  as.data.frame() |>
  tibble::rownames_to_column("gene") |>
  arrange(desc(MeanDecreaseGini)) |>
  head(20)

p_imp <- ggplot(imp_df,
                aes(x = reorder(gene, MeanDecreaseGini),
                    y = MeanDecreaseGini)) +
  geom_col(fill = "#dc2626", alpha = 0.85) +
  coord_flip() +
  labs(
    title = "Feature Importance — топ hub-гены",
    x     = "Ген (HGNC)",
    y     = "Mean Decrease Gini"
  ) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggsave("feature_importance.png", p_imp,
       width = 8, height = 6, dpi = 150)
cat("Feature importance сохранена: feature_importance.png\n")

# ── 16. Итоговый вывод ───────────────────────────────────────
cat("\n=== Итоги Практики 3 ===\n")
cat("GO BP термины (топ-15):  enrichment_GO_BP.csv\n")
cat("KEGG пути (топ-15):      enrichment_KEGG.csv\n")
cat("Accuracy классификатора: ",
    round(cm$overall["Accuracy"], 3), "\n")
cat("Macro F1:                ", round(macro_f1, 3), "\n")
cat("Файлы графиков:\n")
cat("  go_bp_barplot.png, kegg_barplot.png\n")
cat("  confusion_matrix.png, feature_importance.png\n")