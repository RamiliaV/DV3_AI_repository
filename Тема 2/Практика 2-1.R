## Практика 2.1. Random Forest и важность признаков
## Данные: PimaIndiansDiabetes (mlbench)

# github!

## 1. Подготовка окружения ---------------------------------------------

# Установить зеркало по умолчанию (рекомендуется)
# options(repos = "https://mirror.truenetwork.ru/CRAN/")

# Установите пакеты при необходимости:
install.packages(c("mlbench", "dplyr", "randomForest", "pROC"))

library(mlbench)
library(dplyr)
library(randomForest)
library(pROC)

set.seed(123)

## 2. Загрузка и подготовка данных -------------------------------------

data(PimaIndiansDiabetes)
pid <- PimaIndiansDiabetes

# Удаляем строки с пропусками, кодируем исход как 0/1
pid_clean <- pid %>%
  filter(complete.cases(.)) %>%
  mutate(diabetes = ifelse(diabetes == "pos", 1, 0))

table(pid_clean$diabetes)

## 3. Разбиение на train / test (70 / 30) -------------------------------

n <- nrow(pid_clean)
train_idx <- sample(seq_len(n), size = floor(0.7 * n))

train <- pid_clean[train_idx, ]
test  <- pid_clean[-train_idx, ]

table(train$diabetes)
table(test$diabetes)

## 4. Базовый Random Forest на всех признаках --------------------------

# Целевая переменная как фактор
train$diabetes_factor <- factor(train$diabetes, levels = c(0, 1))
test$diabetes_factor  <- factor(test$diabetes,  levels = c(0, 1))

rf_full <- randomForest(
  diabetes_factor ~ . - diabetes,   # исключаем числовой дубль исхода
  data       = train,
  ntree      = 500,
  mtry       = floor(sqrt(ncol(train) - 2)), # без diabetes и diabetes_factor
  importance = TRUE
)

print(rf_full)

## 5. Качество базовой модели на тесте ---------------------------------

# Предсказанные вероятности класса "1"
test_prob_full <- predict(rf_full, newdata = test, type = "prob")[, "1"]
test_pred_full <- ifelse(test_prob_full > 0.5, 1, 0)

cm_full <- table(pred = test_pred_full, true = test$diabetes)
cm_full

accuracy_full <- sum(diag(cm_full)) / sum(cm_full)
accuracy_full

roc_full <- roc(test$diabetes, test_prob_full)
auc_full <- auc(roc_full)
auc_full

## 6. Важность признаков: Gini -----------------------------------------

# type = 2 — MeanDecreaseGini
imp_gini_mat <- importance(rf_full, type = 2)
imp_gini <- imp_gini_mat[, "MeanDecreaseGini"]
imp_gini <- sort(imp_gini, decreasing = TRUE)

imp_gini

## 7. Важность признаков: permutation (по AUC) -------------------------

metric_auc <- function(y, p) {
  as.numeric(auc(y, p))
}

perm_importance <- function(model, data, target, metric_fun) {
  # базовое качество
  base_prob <- predict(model, newdata = data, type = "prob")[, "1"]
  base_metric <- metric_fun(data[[target]], base_prob)

  res <- c()
  feature_names <- setdiff(names(data), c(target, paste0(target, "_factor")))

  for (feat in feature_names) {
    data_perm <- data
    data_perm[[feat]] <- sample(data_perm[[feat]])
    prob_perm <- predict(model, newdata = data_perm, type = "prob")[, "1"]
    res[feat] <- base_metric - metric_fun(data_perm[[target]], prob_perm)
  }

  sort(res, decreasing = TRUE)
}

imp_perm <- perm_importance(
  model      = rf_full,
  data       = test,
  target     = "diabetes",
  metric_fun = metric_auc
)

imp_perm

## 8. Сводная таблица важности признаков -------------------------------

imp_df <- data.frame(
  feature          = names(imp_gini),
  gini_importance  = as.numeric(imp_gini),
  row.names        = NULL
)

imp_df$perm_importance <- imp_perm[match(imp_df$feature, names(imp_perm))]

imp_df <- imp_df %>%
  mutate(
    rank_gini = rank(-gini_importance, ties.method = "min"),
    rank_perm = rank(-perm_importance, ties.method = "min")
  ) %>%
  arrange(rank_perm)

imp_df

## 9. Отбор топ-5 признаков по permutation importance -----------------

top5 <- head(imp_df$feature[!is.na(imp_df$perm_importance)], 5)
top5

# Формируем обучающие и тестовые наборы с топ-10 признаками
train_top5 <- train[, c(top5, "diabetes_factor", "diabetes")]
test_top5  <- test[,  c(top5, "diabetes_factor", "diabetes")]

rf_top5 <- randomForest(
  diabetes_factor ~ . - diabetes,
  data       = train_top5,
  ntree      = 500,
  mtry       = floor(sqrt(length(top5))),
  importance = TRUE
)

print(rf_top5)

## 10. Качество модели с топ-10 признаками -----------------------------

test_prob_top5 <- predict(rf_top5, newdata = test_top5, type = "prob")[, "1"]
test_pred_top5 <- ifelse(test_prob_top5 > 0.5, 1, 0)

cm_top5 <- table(pred = test_pred_top5, true = test_top5$diabetes)
cm_top5

accuracy_top5 <- sum(diag(cm_top5)) / sum(cm_top5)
accuracy_top5

roc_top5 <- roc(test_top5$diabetes, test_prob_top5)
auc_top5 <- auc(roc_top5)
auc_top5

## 11. Краткое сравнение моделей --------------------------------------

plot(roc_full, col = "steelblue", lwd = 2,
     main = "ROC-кривые: полный набор vs топ-5")
lines(roc_top5, col = "darkorange", lwd = 2)
abline(a = 0, b = 1, lty = 2, col = "gray")
legend("bottomright",
       legend = c("Все признаки", "Топ-5 признаков"),
       col    = c("steelblue", "darkorange"),
       lwd    = 2)

results <- data.frame(
  model      = c("RF_full", "RF_top5"),
  accuracy   = c(accuracy_full, accuracy_top5),
  auc        = c(as.numeric(auc_full), as.numeric(auc_top5))
)

results

## 12. Дополнительный блок: SHAP для Random Forest --------------------
## (опционально)

# install.packages("iml")   # выполнить один раз при необходимости
library(iml)

## 12.1. Переобучаем RF без формулы -----------------------------------

# Матрица признаков и цель
X_train <- train[, setdiff(names(train), c("diabetes", "diabetes_factor"))]
y_train <- factor(train$diabetes, levels = c(0, 1))

# Обучение RF в виде x/y, без формулы
rf_xy <- randomForest(
  x         = X_train,
  y         = y_train,
  ntree     = 500,
  mtry      = floor(sqrt(ncol(X_train))),
  importance = TRUE
)

## 12.2. Обёртка-предиктор для iml -------------------------------------

predict_fun <- function(model, newdata) {
  # newdata — data.frame только с признаками
  p <- predict(model, newdata = newdata, type = "prob")
  as.numeric(p[, "1"])
}

predictor <- Predictor$new(
  model            = rf_xy,
  data             = X_train,
  y                = y_train,
  predict.function = predict_fun
)

## 12.3. Локальные SHAP-значения --------------------------------------

x_interest <- X_train[1, , drop = FALSE]

shap_local <- Shapley$new(
  predictor   = predictor,
  x.interest  = x_interest
)

plot(shap_local)