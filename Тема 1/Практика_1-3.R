## Блок 0. Подготовка окружения

install.packages("mlbench")
install.packages("ggplot2")
install.packages("cluster")
install.packages("mclust")
install.packages("factoextra")
install.packages("pheatmap")
install.packages("tidyverse")

library(mlbench)
library(ggplot2)
library(cluster)
library(mclust)
library(factoextra)
library(pheatmap)
library(tidyverse)

set.seed(123)

## Блок 1. Загрузка и первичный осмотр

data(BreastCancer)
bcw <- BreastCancer

cat("Размер датасета (строки x столбцы):\n")
print(dim(bcw))

cat("\nСтруктура датасета:\n")
str(bcw)

cat("\nРаспределение классов (benign / malignant):\n")
print(table(bcw$Class))

## Блок 2. Препроцессинг

bcw_clean <- bcw %>%
  select(-Id) %>%
  drop_na() %>%
  mutate(
    across(-Class, as.numeric)
  )

cat("\nПосле очистки наблюдений:\n")
print(nrow(bcw_clean))

cat("\nРаспределение классов (benign / malignant):\n")
print(table(bcw_clean$Class))

## Блок 3. Признаки и метки

true_labels <- bcw_clean$Class

bcw_features <- bcw_clean %>%
  select(-Class)

cat("\nЧисло признаков:\n")
print(ncol(bcw_features))

## Блок 4. Стандартизация признаков

bcw_scaled <- bcw_features %>%
  scale() %>%
  as.matrix()

cat("\nПроверка: средние (~0):\n")
print(round(colMeans(bcw_scaled), 3))

cat("\nПроверка: sd (~1):\n")
print(round(apply(bcw_scaled, 2, sd), 3))

## Блок 5. Elbow Method

max_k <- 10

wcss <- sapply(1:max_k, function(k) {
  kmeans(bcw_scaled, centers = k, nstart = 25)$tot.withinss
})

elbow_df <- tibble(
  k = 1:max_k,
  WCSS = wcss
)

print(elbow_df)

ggplot(elbow_df, aes(x = k, y = WCSS)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(title = "Elbow Method для выбора числа кластеров",
       x = "Число кластеров (k)",
       y = "WCSS") +
  theme_minimal()

## Блок 6. Silhouette Score

max_k_sil <- 6

sil_scores <- sapply(2:max_k_sil, function(k) {
  km  <- kmeans(bcw_scaled, centers = k, nstart = 25)
  sil <- silhouette(km$cluster, dist(bcw_scaled))
  mean(sil[, 3])
})

sil_df <- tibble(
  k = 2:max_k_sil,
  Silhouette = sil_scores
)

print(sil_df)

ggplot(sil_df, aes(x = k, y = Silhouette)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(title = "Silhouette Score для выбора числа кластеров",
       x = "Число кластеров (k)",
       y = "Средний Silhouette") +
  theme_minimal()

## Блок 7. K-means (k = 2)

k_opt <- 2

set.seed(123)
km_res <- kmeans(bcw_scaled, centers = k_opt, nstart = 25)

clusters_kmeans <- km_res$cluster

cat("\nРазмеры кластеров (K-means, k = 2):\n")
print(table(clusters_kmeans))

cat("\nЦентры кластеров:\n")
print(round(km_res$centers, 2))

## Блок 8. PCA и визуализация K-means

pca_res <- prcomp(bcw_scaled, scale. = FALSE)

pca_df <- tibble(
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  Cluster   = factor(clusters_kmeans),
  TrueClass = true_labels
)

ggplot(pca_df, aes(x = PC1, y = PC2, color = Cluster)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(title = "K-means (k = 2), PCA",
       x = "PC1", y = "PC2",
       color = "Кластер") +
  theme_minimal()

## Блок 9. Иерархическая кластеризация

dist_mat <- dist(bcw_scaled, method = "euclidean")

hc_res <- hclust(dist_mat, method = "average")

plot(hc_res, labels = FALSE,
     main = "Дендрограмма (average linkage)",
     xlab = "", ylab = "Высота", sub = "")

hc_clusters <- cutree(hc_res, k = 2)
rect.hclust(hc_res, k = 2, border = "red")

cat("\nРазмеры кластеров (иерархическая, k = 2):\n")
print(table(hc_clusters))

## Блок 10. Сравнение K-means и иерархической

pca_df <- pca_df %>%
  mutate(HierCluster = factor(hc_clusters))

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

## Блок 11. Сравнение с истинными диагнозами

cat("\nK-means vs истинные классы:\n")
print(table(KMeans = clusters_kmeans, Truth = true_labels))

cat("\nHierarchical vs истинные классы:\n")
print(table(Hier = hc_clusters, Truth = true_labels))

ggplot(pca_df, aes(x = PC1, y = PC2, color = TrueClass)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(title = "Истинные диагнозы (benign / malignant)",
       color = "Диагноз") +
  theme_minimal()

## Блок 12. ARI

true_factor <- factor(true_labels)

ari_kmeans <- adjustedRandIndex(clusters_kmeans, true_factor)
ari_hc     <- adjustedRandIndex(hc_clusters,       true_factor)

cat(sprintf("\nARI (K-means vs истина): %.3f\n", ari_kmeans))
cat(sprintf("ARI (Hierarchical vs истина): %.3f\n", ari_hc))

## Блок 13. Heatmap (бонус)

n_show <- min(100, nrow(bcw_scaled))
heat_data <- bcw_scaled[1:n_show, ]

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

## Блок 14. Упражнение: другие k

set.seed(123)
km3 <- kmeans(bcw_scaled, centers = 3, nstart = 25)
ari_km3 <- adjustedRandIndex(km3$cluster, true_factor)
cat(sprintf("\nARI (K-means, k = 3 vs истина): %.3f\n", ari_km3))