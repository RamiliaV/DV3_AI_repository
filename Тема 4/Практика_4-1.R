# ============================================================
# ПРАКТИКА 1. Подготовка данных и выбор soft-threshold power
# Датасет: GSE45827 (подтипы рака молочной железы)
# ============================================================

# ── 0. Установка пакетов (один раз) ─────────────────────────
# options(repos = "https://mirror.truenetwork.ru/CRAN/")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("GEOquery", "Biobase"), ask = FALSE)
BiocManager::install("WGCNA")
# install.packages(c("ggplot2", "dplyr"))

# ── 1. Загрузка библиотек ────────────────────────────────────
library(GEOquery)
library(Biobase)
library(WGCNA)
library(ggplot2)
library(dplyr)

options(stringsAsFactors = FALSE)
allowWGCNAThreads()   # многопоточность

# ── 2. Загрузка данных с GEO ─────────────────────────────────
gse <- getGEO("GSE45827", GSEMatrix = TRUE, AnnotGPL = TRUE)
gse <- gse[[1]]   # объект ExpressionSet

# Матрица экспрессии: гены (строки) × образцы (столбцы)
expr_raw <- exprs(gse)
pheno    <- pData(gse)

cat("Размерность исходной матрицы:", dim(expr_raw), "\n")
# Ожидается ~30000 генов × 151 образец

# ── 3. Разметка подтипов ─────────────────────────────────────
# Переменная с подтипом находится в столбце "characteristics_ch1.1"
# или аналогичном — проверяем
head(pheno[, grep("characteristics|subtype|type", colnames(pheno),
                  ignore.case = TRUE)])

# Извлечение подтипа (адаптируйте имя столбца по выводу выше)
subtype_raw <- pheno$`characteristics_ch1.1`
subtype <- gsub("tissue: ", "", subtype_raw)
subtype <- trimws(subtype)
table(subtype)

# Удаление клеточных линий: оставляем только тканевые образцы
keep_samples <- subtype != "cell line"
expr_raw <- expr_raw[, keep_samples]
pheno    <- pheno[keep_samples, ]
subtype  <- subtype[keep_samples]

cat("После удаления клеточных линий:", ncol(expr_raw), "образцов\n")
# Ожидается ~150 образцов

# ── 4. Контроль качества образцов ───────────────────────────
# 4.1 Box-plot: распределение экспрессии по образцам
boxplot(expr_raw[, 1:30],
        las = 2, cex.axis = 0.5,
        main = "Распределение экспрессии (первые 30 образцов)",
        ylab = "log2 intensity")

# 4.2 Кластеризация образцов для выявления выбросов
sampleTree <- hclust(dist(t(expr_raw)), method = "average")
plot(sampleTree,
     main = "Кластеризация образцов (выявление выбросов)",
     xlab = "", sub = "",
     cex = 0.6, cex.main = 1.2)

# 4.3 Удаление выбросов по порогу высоты дерева (если нужно)
# Визуально определите порог; пример: abline(h = 250, col = "red")
abline(h = 250, col = "red", lty = 2)
cut_height <- 250
keep <- sampleTree$height < cut_height

expr_raw <- expr_raw[, keep]
pheno    <- pheno[keep, ]
subtype  <- subtype[keep]

# ── 5. Транспонирование: WGCNA требует образцы × гены ────────
datExpr_full <- t(expr_raw)
cat("datExpr_full:", dim(datExpr_full), "\n")

# ── 6. Проверка goodSamplesGenes ─────────────────────────────
gsg <- goodSamplesGenes(datExpr_full, verbose = 3)
gsg$allOK

if (!gsg$allOK) {
  # Удаляем проблемные гены и образцы
  if (sum(!gsg$goodGenes) > 0)
    cat("Удалено генов:", sum(!gsg$goodGenes), "\n")
  if (sum(!gsg$goodSamples) > 0)
    cat("Удалено образцов:", sum(!gsg$goodSamples), "\n")
  datExpr_full <- datExpr_full[gsg$goodSamples, gsg$goodGenes]
  subtype      <- subtype[gsg$goodSamples]
}

# ── 7. Фильтрация по дисперсии — топ-5000 генов ──────────────
gene_var  <- apply(datExpr_full, 2, var)
top_genes <- order(gene_var, decreasing = TRUE)[1:5000]
datExpr   <- datExpr_full[, top_genes]

cat("Финальная матрица для WGCNA:", dim(datExpr), "\n")  # 137 × 5000

# Визуализация распределения дисперсий
var_df <- data.frame(variance = sort(gene_var, decreasing = TRUE),
                     rank = seq_along(gene_var))

ggplot(var_df[1:5000, ], aes(x = rank, y = variance)) +
  geom_line(color = "steelblue") +
  geom_vline(xintercept = 5000, linetype = "dashed", color = "red") +
  labs(title = "Дисперсия генов (топ-5000)",
       x = "Ранг гена", y = "Дисперсия") +
  theme_bw()

# ── 8. Выбор soft-threshold power ────────────────────────────
powers <- 1:20

sft <- pickSoftThreshold(
  datExpr,
  powerVector   = powers,
  networkType   = "signed",
  verbose       = 5
)

# 8.1 График 1: Scale-free fit (R²)
sft_df <- data.frame(
  Power      = sft$fitIndices$Power,
  R2         = -sign(sft$fitIndices$slope) * sft$fitIndices$SFT.R.sq,
  MeanK      = sft$fitIndices$mean.k.
)

p1 <- ggplot(sft_df, aes(x = Power, y = R2)) +
  geom_point(size = 3, color = "steelblue") +
  geom_line(color = "steelblue") +
  geom_text(aes(label = Power), vjust = -0.8, size = 3) +
  geom_hline(yintercept = 0.85, linetype = "dashed", color = "red") +
  annotate("text", x = 18, y = 0.83, label = "R² = 0.85",
           color = "red", size = 3.5) +
  labs(title = "Scale-free fit",
       x = "Soft-threshold power", y = "R² (scale-free fit)") +
  theme_bw()

# 8.2 График 2: Mean connectivity
p2 <- ggplot(sft_df, aes(x = Power, y = MeanK)) +
  geom_point(size = 3, color = "darkorange") +
  geom_line(color = "darkorange") +
  geom_text(aes(label = Power), vjust = -0.8, size = 3) +
  labs(title = "Mean connectivity",
       x = "Soft-threshold power", y = "Средняя связность") +
  theme_bw()

# Вывод двух графиков рядом
library(gridExtra)
grid.arrange(p1, p2, ncol = 2)

# ── 9. Фиксация выбранного power ─────────────────────────────
# Выберите минимальный power, при котором R² ≥ 0.85
softPower <- sft_df$Power[which(sft_df$R2 >= 0.85)[1]]
cat("Выбранный soft-threshold power:", softPower, "\n")

# Если автоматический выбор не сработал — задайте вручную:
# softPower <- 12

# ── 10. Сохранение объектов для Практики 2 ──────────────────
save(datExpr, pheno, subtype, softPower,
     file = "GSE45827_prepared.RData")

cat("Объекты сохранены в GSE45827_prepared.RData\n")
cat("  datExpr:   ", dim(datExpr),    "\n")
cat("  subtype:   ", table(subtype),  "\n")
cat("  softPower: ", softPower,       "\n")