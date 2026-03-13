## Практика 2.2. Grid Search для Random Forest
## Данные: PimaIndiansDiabetes (mlbench)

## 1. Подготовка окружения ---------------------------------------------

# Установить зеркало по умолчанию (рекомендуется)
# options(repos = "https://mirror.truenetwork.ru/CRAN/")

# Установите пакеты при необходимости:
install.packages(c("mlbench", "dplyr", "randomForest", "caret", "pROC"))

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

## 12. Финальная модель с лучшими параметрами --------------------------

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

## 13. SVM с RBF-ядром и Grid Search ----------------------------------
## Требуется пакет kernlab (используется внутри caret)

# install.packages("kernlab")  # при необходимости
library(kernlab)

# 13.1. Сетка гиперпараметров для svmRadial
# C — штраф за ошибки; sigma — параметр ядра RBF
grid_svm <- expand.grid(
  C     = c(0.25, 1, 4),
  sigma = c(0.01, 0.05, 0.1)
)

# 13.2. Обучение SVM с 5-fold CV по ROC
set.seed(123)
svm_grid <- train(
  diabetes ~ .,
  data      = train,
  method    = "svmRadial",
  trControl = ctrl,       # тот же trainControl, что и для RF
  tuneGrid  = grid_svm,
  metric    = "ROC"
)

print(svm_grid)
plot(svm_grid)

# Лучшая комбинация C и sigma
svm_best <- svm_grid$bestTune
svm_best

# 13.3. Оценка лучшей SVM на тесте ------------------------------------

# Предсказанные вероятности и классы
prob_svm <- predict(svm_grid, newdata = test, type = "prob")[, "pos"]
pred_svm <- predict(svm_grid, newdata = test, type = "raw")

# AUC
roc_svm <- roc(test$diabetes, prob_svm)
auc_svm <- auc(roc_svm)
auc_svm

# F1 (используем ту же функцию f1_score, что и для RF)
f1_svm <- f1_score(test$diabetes, pred_svm)
f1_svm

# 13.4. Сравнение RF vs SVM -------------------------------------------

comparison_svm_rf <- data.frame(
  model = c("Random Forest (best mtry)", "SVM RBF (best C, sigma)"),
  AUC   = c(as.numeric(auc_best), as.numeric(auc_svm)),  # auc_best из RF
  F1    = c(f1_best, f1_svm)                             # f1_best из RF
)

comparison_svm_rf

# 13.5. Совместный ROC-график для RF и SVM ----------------------------

plot(roc_best, col = "steelblue", lwd = 2,
     main = "ROC-кривые: Random Forest vs SVM (Pima)")
lines(roc_svm, col = "darkorange", lwd = 2)
abline(a = 0, b = 1, lty = 2, col = "gray")

legend("bottomright",
       legend = c("Random Forest", "SVM RBF"),
       col    = c("steelblue", "darkorange"),
       lwd    = 2)
