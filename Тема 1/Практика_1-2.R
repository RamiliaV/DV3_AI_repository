############################################################
# ПРАКТИКА 2. СРАВНЕНИЕ МОДЕЛЕЙ КЛАССИФИКАЦИИ НА PIMA
# Модели: логистическая регрессия, дерево решений, Random Forest
# Цель: сравнить качество по Accuracy, Sensitivity, Specificity, AUC
############################################################

## Блок 0. Подготовка окружения
## Устанавливаем и подключаем пакеты:
## - mlbench: датасет PimaIndiansDiabetes
## - pROC: построение ROC-кривых и AUC
## - rpart, rpart.plot: дерево решений и его визуализация
## - randomForest: Random Forest
## - tidyverse: удобный синтаксис (dplyr, pipes)

# Установить зеркало по умолчанию (рекомендуется)
# options(repos = "https://mirror.truenetwork.ru/CRAN/")

install.packages("mlbench")
install.packages("pROC")
install.packages("rpart")
install.packages("rpart.plot")
install.packages("randomForest")
install.packages("tidyverse")

library(mlbench)
library(pROC)
library(rpart)
library(rpart.plot)
library(randomForest)
library(tidyverse)

# Фиксируем seed для воспроизводимости (одно и то же разбиение train/test)
set.seed(123)

############################################################
## Блок 1. Загрузка и первичный осмотр данных
############################################################

# Загружаем датасет PimaIndiansDiabetes из пакета mlbench
data(PimaIndiansDiabetes)

# Копируем в объект pid, чтобы не работать с оригиналом напрямую
pid <- PimaIndiansDiabetes

cat("Размер датасета (строки x столбцы):\n")
print(dim(pid))  # число пациентов x число признаков

cat("\nСтруктура датасета:\n")
str(pid)         # типы переменных (numeric, factor и т.д.)

cat("\nРаспределение классов (neg/pos):\n")
print(table(pid$diabetes))  # сколько neg (нет диабета) и pos (есть диабет)

############################################################
## Блок 2. Препроцессинг (очистка, кодирование класса 0/1)
############################################################

# Цель: убрать возможные пропуски и перевести целевой признак в 0/1
# Используем цепочку dplyr: filter + mutate

pid_clean <- pid %>%
  # оставляем только полные строки (без NA во всех столбцах)
  filter(complete.cases(.)) %>%
  # переводим diabetes: pos -> 1, neg -> 0
  mutate(
    diabetes = ifelse(diabetes == "pos", 1, 0)
  )

cat("\nПосле очистки наблюдений:\n")
print(nrow(pid_clean))

cat("\nРаспределение классов (0/1):\n")
print(table(pid_clean$diabetes))

############################################################
## Блок 3. Разбиение на train / test
############################################################

# Стандартное разбиение: 70% на обучение, 30% на тест
n <- nrow(pid_clean)
train_size <- floor(0.7 * n)

# Случайным образом выбираем индексы для обучающей выборки
train_idx <- sample(seq_len(n), size = train_size)

# Формируем обучающую и тестовую выборки
train <- pid_clean[train_idx, ]
test  <- pid_clean[-train_idx, ]

cat("\nTrain:", nrow(train), "наблюдений\n")
cat("Test :", nrow(test),  "наблюдений\n\n")

cat("Классы в train (0/1):\n")
print(table(train$diabetes))

cat("\nКлассы в test (0/1):\n")
print(table(test$diabetes))

############################################################
## Блок 4. Модель 1 — логистическая регрессия
############################################################

# Строим логистическую регрессию:
# diabetes (0/1) ~ . (все остальные признаки)
model_logit <- glm(
  diabetes ~ .,
  data   = train,
  family = binomial(link = "logit")
)

# summary показывает:
# - коэффициенты β
# - стандартные ошибки
# - p-значения
summary(model_logit)

# Считаем вероятности класса 1 (диабет) для наблюдений из test
test <- test %>%
  mutate(
    prob_logit = predict(model_logit, newdata = ., type = "response")
  )

############################################################
## Блок 5. Модель 2 — дерево решений (rpart)
############################################################

# Строим дерево решений:
# factor(diabetes) ~ . — целевая переменная как фактор
# method = "class" — задача классификации
# control — параметры сложности дерева (cp, minsplit и т.д.)
model_tree <- rpart(
  factor(diabetes) ~ .,
  data   = train,
  method = "class",
  control = rpart.control(cp = 0.01, minsplit = 20)
)

# Визуализация дерева: структура правил "если-то"
rpart.plot(model_tree, cex = 0.5)

# Предсказанные вероятности для класса 1 (диабет) на test
prob_tree <- predict(model_tree, newdata = test, type = "prob")[, "1"]

test <- test %>%
  mutate(
    prob_tree = prob_tree
  )

############################################################
## Блок 6. Модель 3 — Random Forest
############################################################

# Строим Random Forest:
# - ntree: число деревьев в ансамбле
# - mtry: число признаков, случайно выбираемых в каждом сплите
model_rf <- randomForest(
  factor(diabetes) ~ .,
  data  = train,
  ntree = 100,   # 100 деревьев — разумный старт
  mtry  = 3      # пример значения mtry (можно настраивать)
)

# Предсказанные вероятности для класса 1 (диабет) на test
prob_rf <- predict(model_rf, newdata = test, type = "prob")[, "1"]

test <- test %>%
  mutate(
    prob_rf = prob_rf
  )

############################################################
## Блок 7. Функция для расчёта метрик по одной модели
############################################################

# Функция принимает:
# - probs: предсказанные вероятности класса 1
# - true: истинные метки (0/1)
# - threshold: порог, по которому вероятности превращаются в классы
# Возвращает вектор из 4 метрик: Accuracy, Sensitivity, Specificity, AUC

metrics_for_model <- function(probs, true, threshold = 0.5) {
  # Переводим вероятности в классы 0/1
  pred <- ifelse(probs >= threshold, 1, 0)
  
  # Строим матрицу ошибок
  cm <- table(pred = pred, true = true)
  
  # Accuracy = доля правильных предсказаний
  accuracy <- sum(diag(cm)) / sum(cm)
  
  # Sensitivity = TP / (TP + FN) — доля правильно найденных больных
  sensitivity <- cm["1", "1"] / sum(cm[, "1"])
  
  # Specificity = TN / (TN + FP) — доля правильно найденных здоровых
  specificity <- cm["0", "0"] / sum(cm[, "0"])
  
  # ROC и AUC по всем порогам
  roc_obj <- roc(true, probs)
  auc_val <- as.numeric(auc(roc_obj))
  
  # Возвращаем named vector
  c(Accuracy   = accuracy,
    Sensitivity = sensitivity,
    Specificity = specificity,
    AUC        = auc_val)
}

############################################################
## Блок 8. Сравнение трёх моделей в одной таблице
############################################################

true <- test$diabetes  # истинные классы для удобства

# Считаем метрики для каждой модели
res_logit <- metrics_for_model(test$prob_logit, true)
res_tree  <- metrics_for_model(test$prob_tree,  true)
res_rf    <- metrics_for_model(test$prob_rf,    true)

# Объединяем в одну таблицу (строки — модели, столбцы — метрики)
metrics_table <- rbind(
  LogReg = res_logit,
  DesTree   = res_tree,
  RF     = res_rf
)

# Округляем до 3 знаков после запятой для компактного вывода
round(metrics_table, 3)

############################################################
## Блок 9. Общий ROC‑график для всех трёх моделей
############################################################

# ROC‑объекты для каждой модели
roc_logit <- roc(true, test$prob_logit)
roc_tree  <- roc(true, test$prob_tree)
roc_rf    <- roc(true, test$prob_rf)

# Рисуем ROC логистической регрессии
plot(roc_logit, col = "blue", lwd = 2,
     main = "ROC-кривые: LogReg, Tree, RF")

# Добавляем ROC для дерева и RF на тот же график
plot(roc_tree, col = "red",      lwd = 2, add = TRUE)
plot(roc_rf,   col = "darkgreen", lwd = 2, add = TRUE)

# Легенда
legend("bottomright",
       legend = c("LogReg", "Tree", "RF"),
       col    = c("blue", "red", "darkgreen"),
       lwd    = 2)

############################################################
## Блок 10. Упражнение: влияние порога для всех моделей
############################################################

# Функция: считает метрики для одной модели при одном пороге
compute_metrics_threshold <- function(model_name, thr, probs, true_class) {
  pred_thr <- ifelse(probs >= thr, 1, 0)
  cm_thr   <- table(pred = pred_thr, true = true_class)
  
  accuracy_thr    <- sum(diag(cm_thr)) / sum(cm_thr)
  sensitivity_thr <- cm_thr["1", "1"] / sum(cm_thr[, "1"])
  specificity_thr <- cm_thr["0", "0"] / sum(cm_thr[, "0"])
  
  data.frame(
    Model       = model_name,
    Threshold   = thr,
    Accuracy    = accuracy_thr,
    Sensitivity = sensitivity_thr,
    Specificity = specificity_thr
  )
}

thresholds <- c(0.3, 0.5, 0.7)
true <- test$diabetes

results_list <- list()

for (thr in thresholds) {
  res_logit <- compute_metrics_threshold("LogReg", thr, test$prob_logit, true)
  res_tree  <- compute_metrics_threshold("Tree",   thr, test$prob_tree,  true)
  res_rf    <- compute_metrics_threshold("RF",     thr, test$prob_rf,    true)
  
  results_list[[length(results_list) + 1]] <- res_logit
  results_list[[length(results_list) + 1]] <- res_tree
  results_list[[length(results_list) + 1]] <- res_rf
}

threshold_metrics <- bind_rows(results_list)

# Округлим для компактности
threshold_metrics_round <- threshold_metrics %>%
  mutate(
    Accuracy    = round(Accuracy, 3),
    Sensitivity = round(Sensitivity, 3),
    Specificity = round(Specificity, 3)
  )

threshold_metrics_round
