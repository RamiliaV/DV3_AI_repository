############################################################
# ПРАКТИКА 3. КЛАСТЕРИЗАЦИЯ (K-MEANS И ИЕРАРХИЧЕСКАЯ) НА BCW
# Датасет: Breast Cancer Wisconsin (Diagnostic), mlbench::BreastCancer
# Цели:
#  - стандартизация признаков
#  - выбор числа кластеров (Elbow, Silhouette)
#  - K-means и иерархическая кластеризация
#  - визуализация через PCA
#  - сравнение с истинными классами (ARI)
############################################################

## Блок 0. Подготовка окружения
## Пакеты:
## - mlbench: датасет BreastCancer
## - tidyverse: dplyr, ggplot2 и др.
## - cluster: silhouette
## - mclust: Adjusted Rand Index (ARI)
## - factoextra: удобная визуализация PCA/кластеров
## - pheatmap: heatmap с дендрограммой

# Установить зеркало по умолчанию (рекомендуется)
# options(repos = "https://mirror.truenetwork.ru/CRAN/")

# install.packages("mlbench")
# install.packages("tidyverse")
# install.packages("cluster")
# install.packages("mclust")
# install.packages("factoextra")
# install.packages("pheatmap")

library(mlbench)
library(tidyverse)
library(cluster)
library(mclust)
library(factoextra)
library(pheatmap)

# Фиксируем seed для воспроизводимости кластеризации (K-means, выбор подвыборки и т.п.)
set.seed(123)

############################################################
## Блок 1. Загрузка и первичный осмотр данных
############################################################

# Загружаем датасет BreastCancer
data(BreastCancer)

# Копируем в bcw для дальнейшей работы
bcw <- BreastCancer

cat("Размер датасета (строки x столбцы):\n")
print(dim(bcw))      # число образцов x число колонок

cat("\nСтруктура датасета:\n")
str(bcw)             # типы переменных

cat("\nРаспределение классов (benign / malignant):\n")
print(table(bcw$Class))   # сколько доброкачественных/злокачественных

############################################################
## Блок 2. Препроцессинг: удаление Id, NA, перевод признаков в numeric
############################################################

# Цель: подготовить числовую матрицу признаков для кластеризации.
# Шаги:
#  - удалить идентификатор Id (не информативен)
#  - убрать строки с пропусками
#  - перевести все диагностические признаки (кроме Class) в числовой формат

bcw_clean <- bcw %>%
  select(-Id) %>%   # удаляем колонку Id
  drop_na() %>%     # убираем любые строки с NA
  mutate(
    across(-Class, as.numeric)  # все колонки кроме Class -> numeric
  )

cat("\nПосле очистки наблюдений:\n")
print(nrow(bcw_clean))

cat("\nРаспределение классов (benign / malignant) после очистки:\n")
print(table(bcw_clean$Class))

############################################################
## Блок 3. Выделение признаков и истинных меток
############################################################

# Истинные классы (диагноз) будем использовать только для оценки качества кластеризации
true_labels <- bcw_clean$Class

# Матрица признаков: все колонки кроме Class
bcw_features <- bcw_clean %>%
  select(-Class)

cat("\nЧисло признаков для кластеризации:\n")
print(ncol(bcw_features))

############################################################
## Блок 4. Стандартизация признаков
############################################################

# Почему это важно:
# - K-means использует евклидово расстояние
# - признаки из разных шкал (1–10, 1–1000) должны быть приведены к сопоставимому масштабу
# scale() -> Z-score: (x - mean) / sd

bcw_scaled <- bcw_features %>%
  scale() %>%
  as.matrix()   # K-means и dist работают с матрицами/матрицами-like

cat("\nПроверка: средние (~0):\n")
print(round(colMeans(bcw_scaled), 3))

cat("\nПроверка: sd (~1):\n")
print(round(apply(bcw_scaled, 2, sd), 3))

############################################################
## Блок 5. Выбор числа кластеров: Elbow Method
############################################################

# Elbow Method:
#  - для k = 1..K считаем tot.withinss (WCSS — суммарная внутрикластерная сумма квадратов)
#  - ищем "локоть" — точку, после которой уменьшение WCSS замедляется
#  - часто даёт разумную оценку k

max_k <- 10

wcss <- sapply(1:max_k, function(k) {
  kmeans(bcw_scaled, centers = k, nstart = 25)$tot.withinss
})

elbow_df <- tibble(
  k    = 1:max_k,
  WCSS = wcss
)

cat("\nWCSS для разных k (Elbow Method):\n")
print(elbow_df)

ggplot(elbow_df, aes(x = k, y = WCSS)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(title = "Elbow Method для выбора числа кластеров",
       x = "Число кластеров (k)",
       y = "WCSS (tot.withinss)") +
  scale_x_continuous(breaks = 1:10) +
  theme_minimal()

############################################################
## Блок 6. Выбор числа кластеров: Silhouette Score
############################################################

# Silhouette Score:
#  - для каждой точки i: s(i) = (b(i) - a(i)) / max(a(i), b(i))
#    a(i) — среднее расстояние до точек в своём кластере
#    b(i) — среднее расстояние до closest других кластеров
#  - средний silhouette по датасету > 0.5 — хорошие кластеры
#  - мы считаем средний silhouette для k = 2..K и выбираем максимум

max_k_sil <- 6

sil_scores <- sapply(2:max_k_sil, function(k) {
  km  <- kmeans(bcw_scaled, centers = k, nstart = 25)
  sil <- silhouette(km$cluster, dist(bcw_scaled))
  mean(sil[, 3])  # индекс 3 — значение silhouette для каждой точки
})

sil_df <- tibble(
  k          = 2:max_k_sil,
  Silhouette = sil_scores
)

cat("\nСредние Silhouette для разных k:\n")
print(sil_df)

ggplot(sil_df, aes(x = k, y = Silhouette)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(title = "Silhouette Score для выбора числа кластеров",
       x = "Число кластеров (k)",
       y = "Средний Silhouette") +
  theme_minimal()

############################################################
## Блок 7. K-means кластеризация (k = 2)
############################################################

# Для BCW чаще всего разумное k = 2 (benign vs malignant).
# Здесь явно фиксируем k_opt = 2.
k_opt <- 2

set.seed(123)
km_res <- kmeans(bcw_scaled, centers = k_opt, nstart = 25)

# Номера кластеров для каждого наблюдения (1..k)
clusters_kmeans <- km_res$cluster

cat("\nРазмеры кластеров (K-means, k = 2):\n")
print(table(clusters_kmeans))

cat("\nЦентры кластеров (в стандартизованных координатах):\n")
print(round(km_res$centers, 2))

############################################################
## Блок 8. PCA и визуализация кластеров K-means
############################################################

# PCA применяется для визуализации высокомерных данных (много признаков)
#  - prcomp на стандартизованной матрице
#  - используем первые две главные компоненты (PC1, PC2)

pca_res <- prcomp(bcw_scaled, scale. = FALSE)

# Собираем tibble с координатами PC и информацией о кластерах и истинных классах
pca_df <- tibble(
  PC1       = pca_res$x[, 1],
  PC2       = pca_res$x[, 2],
  Cluster   = factor(clusters_kmeans),
  TrueClass = true_labels
)

# Визуализация K-means кластеров на плоскости PC1–PC2
ggplot(pca_df, aes(x = PC1, y = PC2, color = Cluster)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(title = "K-means (k = 2), визуализация через PCA",
       x = "Главная компонента 1",
       y = "Главная компонента 2",
       color = "Кластер") +
  theme_minimal()

############################################################
## Блок 9. Иерархическая кластеризация и дендрограмма
############################################################

# Иерархическая кластеризация:
#  - начинаем с каждого объекта в отдельном кластере
#  - последовательно объединяем ближайшие кластеры
#  - получаем дендрограмму (дерево)
# Здесь используем:
#  - расстояние: евклидово
#  - метод связи: average (среднее расстояние между кластерами)

dist_mat <- dist(bcw_scaled, method = "euclidean")

hc_res <- hclust(dist_mat, method = "average")

# Рисуем дендрограмму
plot(hc_res, labels = FALSE,
     main = "Дендрограмма (average linkage)",
     xlab = "", ylab = "Высота", sub = "")

# Разрезаем дерево на 2 кластера (k = 2)
hc_clusters <- cutree(hc_res, k = 2)

# Добавляем прямоугольники вокруг кластеров на дендрограмме
rect.hclust(hc_res, k = 2, border = "red")

cat("\nРазмеры кластеров (иерархическая, k = 2):\n")
print(table(hc_clusters))

############################################################
## Блок 10. Сравнение K-means и иерархической кластеризации
############################################################

# Добавляем результаты иерархической кластеризации в pca_df
pca_df <- pca_df %>%
  mutate(HierCluster = factor(hc_clusters))

# Сравниваем визуально на PCA: слева K-means, справа Hierarchical
install.packages("gridExtra")
library(gridExtra)

p1 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Cluster)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(title = "K-means (k = 2)") +
  theme_minimal()

p2 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = HierCluster)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(title = "Иерархическая (k = 2)") +
  theme_minimal()

grid.arrange(p1, p2, ncol = 2)

############################################################
## Блок 11. Сравнение кластеров с истинными классами
############################################################

# Таблицы сопряжённости: кластеры vs истинный диагноз
cat("\nK-means (k = 2) vs истинные классы:\n")
print(table(KMeans = clusters_kmeans, Truth = true_labels))

cat("\nHierarchical (k = 2) vs истинные классы:\n")
print(table(Hier = hc_clusters, Truth = true_labels))

# Визуализация истинных классов на PCA
ggplot(pca_df, aes(x = PC1, y = PC2, color = TrueClass)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(title = "Истинные диагнозы (benign / malignant)",
       color = "Диагноз") +
  theme_minimal()

############################################################
## Блок 12. Оценка качества кластеризации: Adjusted Rand Index (ARI)
############################################################

# ARI измеряет согласованность двух разбиений:
#  - 1 = идеальное совпадение
#  - 0 = случайное совпадение
#  - < 0 = хуже случайного

true_factor <- factor(true_labels)

ari_kmeans <- adjustedRandIndex(clusters_kmeans, true_factor)
ari_hc     <- adjustedRandIndex(hc_clusters,       true_factor)

cat(sprintf("\nARI (K-means vs истина): %.3f\n", ari_kmeans))
cat(sprintf("ARI (Hierarchical vs истина): %.3f\n", ari_hc))

############################################################
## Блок 13. Heatmap с дендрограммой
############################################################

# Для наглядности возьмём первые 100 образцов
n_show <- min(100, nrow(bcw_scaled))
heat_data <- bcw_scaled[1:n_show, ]

# Аннотация строк: кластеры K-means, иерархические и истинный класс
annotation_row <- data.frame(
  KMeans = factor(clusters_kmeans[1:n_show]),
  Hier   = factor(hc_clusters[1:n_show]),
  Class  = true_labels[1:n_show]
)
rownames(annotation_row) <- rownames(heat_data)

pheatmap(heat_data,
         clustering_distance_rows = "euclidean",
         clustering_method        = "average",
         annotation_row           = annotation_row,
         show_rownames            = FALSE,
         main = "Heatmap: первые 100 образцов, BCW")
