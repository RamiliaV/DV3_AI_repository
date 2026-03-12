## Практика 2.2. Grid Search для Random Forest
## Данные: PimaIndiansDiabetes (mlbench)

## 1. Подготовка окружения ---------------------------------------------

# Установить зеркало по умолчанию (рекомендуется)
# options(repos = "https://mirror.truenetwork.ru/CRAN/")

# Установите пакеты при необходимости:
# install.packages(c("mlbench", "dplyr", "randomForest", "caret", "pROC"))

library(mlbench)
library(dplyr)
library(randomForest)
library(caret)
library(pROC)

## 2. Загрузка и подготовка данных -------------------------------------

data(PimaIndiansDiabetes)
pid <- PimaIndiansDiabetes

# Удаляем строки с пропусками, кодируем исход как 0/1
pid_clean <- pid %>%
  filter(complete.cases(.)) %>%
  mutate(diabetes = ifelse(diabetes == "pos", "pos", "neg"))

# Преобразуем в фактор с понятными уровнями для caret
pid_clean$diabetes <- factor(pid_clean$diabetes, levels = c("neg", "pos"))

table(pid_clean$diabetes)

## 3. Разбиение на train / test (70 / 30) -------------------------------

set.seed(123)
n <- nrow(pid_clean)
train_idx <- sample(seq_len(n), size = floor(0.7 * n))

train <- pid_clean[train_idx, ]
test  <- pid_clean[-train_idx, ]

table(train$diabetes)
table(test$diabetes)

## 4. Базовая модель Random Forest (для сравнения) ---------------------

rf_baseline <- randomForest(
  diabetes ~ .,
  data       = train,
  ntree      = 500,
  mtry       = 2,
  importance = TRUE
)

print(rf_baseline)

## Функция для F1 (положительный класс = "pos")
f1_score <- function(true, pred_class, positive = "pos") {
  true <- factor(true, levels = c("neg", "pos"))
  pred_class <- factor(pred_class, levels = c("neg", "pos"))
  cm <- table(pred_class, true)
  tp <- cm[positive, positive]
  fp <- sum(cm[positive, ]) - tp
  fn <- sum(cm[, positive]) - tp
  precision <- tp / (tp + fp)
  recall    <- tp / (tp + fn)
  2 * precision * recall / (precision + recall)
}

# Оценка на тесте
pred_baseline_prob <- predict(rf_baseline, newdata = test, type = "prob")[, "pos"]
roc_baseline <- roc(test$diabetes, pred_baseline_prob)
auc_baseline <- auc(roc_baseline)
auc_baseline

# Классы для baseline
pred_baseline_class <- ifelse(pred_baseline_prob > 0.5, "pos", "neg")
f1_baseline <- f1_score(test$diabetes, pred_baseline_class)
f1_baseline

## 5. Grid Search: определение сетки гиперпараметров -------------------

# Сетка для перебора
tune_grid <- expand.grid(
  mtry = c(2, 4, 6)
)

# Примечание: caret не поддерживает ntree в tuneGrid для метода "rf"
# ntree задаётся отдельно через параметр ntree в train()
# Мы будем варьировать mtry, а ntree зафиксируем

## 6. Настройка cross-validation ---------------------------------------

ctrl <- trainControl(
  method          = "cv",
  number          = 5,
  summaryFunction = twoClassSummary,
  classProbs      = TRUE,
  savePredictions = "final",
  verboseIter     = TRUE
)

## 7. Grid Search с caret ----------------------------------------------

# Обучение с перебором гиперпараметров
grid_search <- train(
  diabetes ~ .,
  data      = train,
  method    = "rf",
  trControl = ctrl,
  tuneGrid  = tune_grid,
  metric    = "ROC",
  ntree     = 500  # фиксируем число деревьев
)

print(grid_search)
plot(grid_search)

# Лучшая комбинация параметров
best_params <- grid_search$bestTune
best_params

## 8. Оценка финальной модели на тесте ---------------------------------

pred_best_prob <- predict(grid_search, newdata = test, type = "prob")[, "pos"]
roc_best <- roc(test$diabetes, pred_best_prob)
auc_best <- auc(roc_best)
auc_best

# Классы для лучшей модели (тот же порог 0.5)
pred_best_class <- ifelse(pred_best_prob > 0.5, "pos", "neg")
f1_best <- f1_score(test$diabetes, pred_best_class)
f1_best

# Сравнение базовой и настроенной модели
comparison <- data.frame(
  model = c("Baseline (mtry=2, ntree=500)", "Grid Search Best"),
  AUC   = c(as.numeric(auc_baseline), as.numeric(auc_best))
)
comparison

## 9. Random Search (опциональный блок) --------------------------------

# Random Search: случайный выбор комбинаций из более широкого пространства

ctrl_random <- trainControl(
  method          = "cv",
  number          = 5,
  search          = "random",  # включаем random search
  summaryFunction = twoClassSummary,
  classProbs      = TRUE,
  savePredictions = "final",
  verboseIter     = TRUE
)

set.seed(123)
random_search <- train(
  diabetes ~ .,
  data      = train,
  method    = "rf",
  trControl = ctrl_random,
  tuneLength = 10,  # количество случайных комбинаций
  metric    = "ROC",
  ntree     = 500
)

print(random_search)
plot(random_search)

# Оценка на тесте
pred_random_prob <- predict(random_search, newdata = test, type = "prob")[, "pos"]
roc_random <- roc(test$diabetes, pred_random_prob)
auc_random <- auc(roc_random)
auc_random

pred_random_class <- ifelse(pred_random_prob > 0.5, "pos", "neg")
f1_random <- f1_score(test$diabetes, pred_random_class)
f1_random

## 10. Итоговое сравнение всех моделей ---------------------------------

final_comparison <- data.frame(
  model = c("Baseline", "Grid Search", "Random Search"),
  AUC   = c(as.numeric(auc_baseline), 
            as.numeric(auc_best), 
            as.numeric(auc_random))
)
final_comparison

## 11. Визуализация результатов Grid Search ----------------------------

# График ROC-кривых для сравнения
plot(roc_baseline, col = "gray", lwd = 2, main = "ROC Curves Comparison")
lines(roc_best, col = "blue", lwd = 2)
lines(roc_random, col = "red", lwd = 2)
legend("bottomright", 
       legend = c("Baseline", "Grid Search", "Random Search"),
       col = c("gray", "blue", "red"),
       lwd = 2)

## 12. Детальные результаты Grid Search --------------------------------

# Посмотреть все результаты по fold'ам
grid_search$results

# Визуализация зависимости AUC от mtry
ggplot(grid_search$results, aes(x = mtry, y = ROC)) +
  geom_line(color = "blue", size = 1.2) +
  geom_point(color = "blue", size = 3) +
  labs(title = "Grid Search: AUC-ROC vs mtry",
       x = "mtry (число признаков на разбиение)",
       y = "Mean AUC-ROC (5-fold CV)") +
  theme_minimal()

## 13. Финальная модель с лучшими параметрами --------------------------

# Переобучим финальную модель на всём train с лучшими параметрами
final_rf <- randomForest(
  diabetes ~ .,
  data  = train,
  mtry  = best_params$mtry,
  ntree = 500,
  importance = TRUE
)

print(final_rf)

# Важность признаков в финальной модели
varImpPlot(final_rf, main = "Feature Importance (Final Model)")

## 14. Confusion Matrix на тесте ---------------------------------------

pred_final_class <- predict(final_rf, newdata = test, type = "class")
cm_final <- confusionMatrix(pred_final_class, test$diabetes)
cm_final

## 15. Выводы ----------------------------------------------------------

cat("\n=== ИТОГИ ПРАКТИКИ 2.2 ===\n")
cat("Базовая модель (mtry=2): AUC =", round(auc_baseline, 3), "\n")
cat("Grid Search (лучший mtry):", best_params$mtry, ", AUC =", round(auc_best, 3), "\n")
cat("Random Search: AUC =", round(auc_random, 3), "\n")
cat("\nУлучшение качества:", 
    ifelse(auc_best > auc_baseline, "ДА", "НЕТ"), "\n")
cat("Разница AUC:", round(auc_best - auc_baseline, 4), "\n")
