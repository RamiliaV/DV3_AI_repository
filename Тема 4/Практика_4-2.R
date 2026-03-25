# ============================================================
# ПРАКТИКА 2. Построение модулей и hub-генов (WGCNA)
# Датасет: GSE45827 (подтипы рака молочной железы)
# ============================================================

# ── 1. Загрузка библиотек ────────────────────────────────────
# options(repos = "https://mirror.truenetwork.ru/CRAN/")
library(WGCNA)
library(ggplot2)
library(dplyr)
library(gridExtra)

options(stringsAsFactors = FALSE)
allowWGCNAThreads()

# ── 2. Загрузка объектов из Практики 1 ───────────────────────
load("GSE45827_prepared.RData")
# Доступны: datExpr, pheno, subtype, softPower

cat("Матрица экспрессии:", dim(datExpr), "\n")
cat("Подтипы:", table(subtype), "\n")
cat("softPower:", softPower, "\n")

# ── 3. Запуск blockwiseModules ───────────────────────────────
# Занимает 5–15 минут в зависимости от машины
bwnet <- blockwiseModules(
  datExpr,
  power             = softPower,
  networkType       = "signed",
  TOMType           = "signed",
  minModuleSize     = 30,
  mergeCutHeight    = 0.25,
  numericLabels     = FALSE,   # цвета вместо цифр
  pamRespectsDendro = FALSE,
  saveTOMs          = FALSE,
  verbose           = 3
)

# Сохраняем сразу — пересчёт дорогой
save(bwnet, file = "GSE45827_bwnet.RData")

# ── 4. Таблица модулей ───────────────────────────────────────
module_colors <- bwnet$colors
mod_table     <- sort(table(module_colors), decreasing = TRUE)
print(mod_table)

# Крупнейший модуль (кроме grey — «мусорный» модуль)
largest_module <- names(mod_table)[names(mod_table) != "grey"][1]
cat("Крупнейший модуль:", largest_module,
    "(", mod_table[largest_module], "генов )\n")

# ── 5. Дендрограмма с раскрашенными ветвями ─────────────────
png("dendro_modules.png", width = 1200, height = 600, res = 120)
plotDendroAndColors(
  bwnet$dendrograms[[1]],
  module_colors[bwnet$blockGenes[[1]]],
  "Модули",
  dendroLabels = FALSE,
  hang         = 0.03,
  addGuide     = TRUE,
  guideHang    = 0.05,
  main         = "Дендрограмма генов и цвета модулей (GSE45827)"
)
dev.off()
cat("Дендрограмма сохранена: dendro_modules.png\n")

# ── 6. Module Eigengenes (ME) ────────────────────────────────
MEs_obj <- moduleEigengenes(datExpr, module_colors)
MEs     <- MEs_obj$eigengenes
MEs     <- orderMEs(MEs)

cat("Число модулей (включая grey):", ncol(MEs), "\n")

# ── 7. Связь ME с подтипами — ANOVA ─────────────────────────
# Кодируем подтипы как числа для корреляции
subtype_num <- as.numeric(factor(subtype))

# Корреляция ME с подтипом
cor_ME  <- cor(MEs, subtype_num, use = "pairwise.complete.obs")
pval_ME <- corPvalueStudent(cor_ME, nrow(MEs))

ME_subtype <- data.frame(
  Module  = rownames(cor_ME),
  r       = round(cor_ME[, 1], 3),
  p.value = round(pval_ME[, 1], 4)
) |> arrange(r)

print(head(ME_subtype, 10))

# ── 8. Heatmap ME × подтип ───────────────────────────────────
# Матрица: образцы × ME
ME_plot <- data.frame(MEs, subtype = subtype)

# Средние значения ME по подтипам
ME_means <- ME_plot |>
  group_by(subtype) |>
  summarise(across(starts_with("ME"), mean)) |>
  as.data.frame()

rownames(ME_means) <- ME_means$subtype
ME_means$subtype   <- NULL
ME_matrix          <- as.matrix(ME_means)

png("heatmap_ME_subtype.png", width = 1000, height = 600, res = 120)
heatmap(
  ME_matrix,
  scale   = "column",
  col     = colorRampPalette(c("#2563eb", "white", "#dc2626"))(100),
  margins = c(10, 8),
  main    = "Средние Module Eigengenes по подтипам BC"
)
dev.off()
cat("Heatmap сохранён: heatmap_ME_subtype.png\n")

# ── 9. Выбор целевого модуля ─────────────────────────────────
# Берём модуль с наименьшим p-value по связи с подтипом
# (исключая grey)
target_module <- ME_subtype |>
  filter(Module != "MEgrey") |>
  slice(1) |>
  pull(Module) |>
  gsub("ME", "", x = _)

cat("Целевой модуль для hub-генов:", target_module, "\n")

# ── 10. Расчёт kME (module membership) ──────────────────────
kME <- signedKME(datExpr, MEs, outputColumnName = "kME")

# kME для целевого модуля
kme_col    <- paste0("kME", target_module)
kME_target <- data.frame(
  gene = colnames(datExpr),
  kME  = kME[, kme_col]
) |>
  arrange(desc(kME)) |>
  filter(module_colors[match(gene, colnames(datExpr))] == target_module)

cat("Топ-10 hub-генов по kME в модуле", target_module, ":\n")
print(head(kME_target, 10))

# ── 11. Топ-20 hub-генов → сохранение ───────────────────────
hub_genes <- head(kME_target, 20)

BiocManager::install("hgu133plus2.db", ask = FALSE)
library(hgu133plus2.db)

hub_symbols <- mapIds(
  hgu133plus2.db,
  keys     = hub_genes$gene,
  column   = "SYMBOL",
  keytype  = "PROBEID",
  multiVals = "first"
)

hub_genes$symbol <- hub_symbols

# Убираем зонды без аннотации
hub_genes <- hub_genes[!is.na(hub_genes$symbol), ]

cat("Hub-генов с аннотацией:", nrow(hub_genes), "\n")
print(hub_genes[, c("gene", "symbol", "kME")])

# Обновляем CSV
write.csv(hub_genes, "hub_genes_GSE45827.csv", row.names = FALSE)

# ── 12. Визуализация kME топ-20 ─────────────────────────────
p_hub <- ggplot(hub_genes, aes(x = reorder(symbol, kME), y = kME)) +
  geom_col(fill = target_module, alpha = 0.85) +
  coord_flip() +
  labs(
    title = paste("Топ-20 hub-генов модуля", target_module),
    x     = "Ген",
    y     = "kME (module membership)"
  ) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggsave("hub_genes_barplot.png", p_hub,
       width = 8, height = 6, dpi = 150)
cat("График hub-генов сохранён: hub_genes_barplot.png\n")

# ── 13. Сохранение объектов для Практики 3 ──────────────────
save(datExpr, pheno, subtype, MEs, module_colors,
     kME, hub_genes, target_module,
     file = "GSE45827_modules.RData")

cat("\n=== Итоги Практики 2 ===\n")
cat("Модулей построено:   ", ncol(MEs) - 1, "(без grey)\n")
cat("Целевой модуль:      ", target_module, "\n")
cat("Hub-генов отобрано:  ", nrow(hub_genes), "\n")
cat("Файлы: GSE45827_modules.RData, hub_genes_GSE45827.csv\n")
cat("       dendro_modules.png, heatmap_ME_subtype.png\n")
cat("       hub_genes_barplot.png\n")
