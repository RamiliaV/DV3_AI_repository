############################################################
# ПРАКТИКА 1. ЛОГИСТИЧЕСКАЯ РЕГРЕССИЯ: ПРИМЕР НА ДИАБЕТЕ
############################################################

## Блок 0. Подготовка окружения
## Устанавливаем и подключаем нужные пакеты.
## tidyverse включает dplyr, ggplot2, tibble и др.

# Установить зеркало по умолчанию (рекомендуется)
# options(repos = "https://mirror.truenetwork.ru/CRAN/")

install.packages("mlbench")    # датасеты, в т.ч. PimaIndiansDiabetes
install.packages("pROC")       # ROC-кривая и AUC
install.packages("tidyverse")  # dplyr, ggplot2 и др.

library(mlbench)
library(pROC)
library(tidyverse)

# Фиксируем seed, чтобы разбиение на train/test было воспроизводимым.
set.seed(123)

## Блок 1. Загрузка и первичный осмотр данных
## Загружаем датасет PimaIndiansDiabetes из mlbench и смотрим его структуру.

data(PimaIndiansDiabetes)
pid <- PimaIndiansDiabetes

cat("Размер датасета (строки x столбцы):\n")
print(dim(pid))   # число наблюдений и число признаков

cat("\nСтруктура датасета:\n")
str(pid)          # типы столбцов, пример значений

cat("\nРаспределение классов (neg/pos):\n")
print(table(pid$diabetes))  # баланс классов исходно (neg/pos)

## Блок 2. Препроцессинг (очистка, кодирование класса 0/1)
## Цель: убрать строки с пропусками и перевести целевой признак в 0/1.
## Используем dplyr-цепочку: filter + mutate.

pid_clean <- pid %>%
  # Оставляем только полные случаи (без NA в любых столбцах)
  filter(complete.cases(.)) %>%
  # Переводим diabetes в числовой формат: pos -> 1, neg -> 0
  mutate(
    diabetes = ifelse(diabetes == "pos", 1, 0)
  )

cat("\nПосле очистки наблюдений:\n")
print(nrow(pid_clean))

cat("\nРаспределение классов (0/1):\n")
print(table(pid_clean$diabetes))

## Блок 3. Разбиение на train / test
## Цель: разделить данные на обучающую (70%) и тестовую (30%) выборки.

n <- nrow(pid_clean)               # всего наблюдений
train_size <- floor(0.7 * n)       # размер обучающей выборки

# Выбираем случайные индексы для train
train_idx <- sample(seq_len(n), size = train_size)

# Формируем train и test как обычные data.frame
train <- pid_clean[train_idx, ]
test  <- pid_clean[-train_idx, ]

cat("\nTrain:", nrow(train), "наблюдений\n")
cat("Test :", nrow(test),  "наблюдений\n\n")

cat("Классы в train (0/1):\n")
print(table(train$diabetes))

cat("\nКлассы в test (0/1):\n")
print(table(test$diabetes))

## Блок 4. Обучение логистической регрессии
## Строим модель: diabetes ~ . (все остальные переменные как предикторы).
## family = binomial(link = "logit") — классическая логистическая регрессия.

model_logit <- glm(
  diabetes ~ .,              # слева — целевой признак, справа — все X
  data   = train,            # обучающие данные
  family = binomial(link = "logit")
)

# Сводка по модели: коэффициенты, значимость предикторов и т.д.
summary(model_logit)

## Сигмоидная функция (логистическая)
# Извлекаем линейный предиктор η = β0 + βX и вероятность p
test$eta <- predict(model_logit, newdata = test, type = "link")      # линейный предиктор
test$prob <- predict(model_logit, newdata = test, type = "response") # вероятность (у вас уже есть)

# Строим гладкую сигмоиду
sigmoid <- function(x) 1 / (1 + exp(-x))

eta_range <- range(test$eta)
x_grid <- seq(eta_range[1] - 1, eta_range[2] + 1, length.out = 400)
y_grid <- sigmoid(x_grid)

plot(x_grid, y_grid,
     type = "l", lwd = 2, col = "blue",
     xlab = "Линейный предиктор η",
     ylab = "Вероятность p",
     main = "Сигмоидная функция и реальные предсказания")

# Добавляем реальные точки (η, p) для тестовых наблюдений
points(test$eta, test$prob,
       pch = 16, col = rgb(1, 0, 0, 0.5))

# Линия порога p = 0.5 (η = 0)
abline(h = 0.5, lty = 2, col = "gray")
abline(v = 0,   lty = 2, col = "gray")

cols <- ifelse(test$diabetes == 1, "red", "darkgreen")
points(test$eta, test$prob, pch = 16, col = adjustcolor(cols, alpha.f = 0.6))
legend("topleft",
       legend = c("diabetes = 1", "diabetes = 0"),
       col    = c("red", "darkgreen"),
       pch    = 16)


## Блок 5. Прогноз на тесте и классификация при пороге 0.5
## Сначала получаем предсказанные вероятности диабета,
## затем превращаем их в классы 0/1 по порогу.

# Добавляем в test столбец prob — вероятность диабета по модели
test <- test %>%
  mutate(
    prob = predict(model_logit, newdata = ., type = "response")
  )

# Фиксируем порог классификации
threshold <- 0.5

# pred = 1, если prob >= 0.5; иначе pred = 0
test <- test %>%
  mutate(
    pred = ifelse(prob >= threshold, 1, 0)
  )

cat("\nПервые 10 наблюдений (истина, p, предсказание):\n")
print(
  test %>%
    select(diabetes, prob, pred) %>%
    head(10)
)

## Блок 6. Матрица ошибок и базовые метрики
## Строим confusion matrix и считаем Accuracy, Sensitivity, Specificity.

cm <- table(pred = test$pred, true = test$diabetes)

cat("\nМатрица ошибок (threshold =", threshold, "):\n")
print(cm)

# Accuracy = (TP + TN) / (всего)
accuracy <- sum(diag(cm)) / sum(cm)

# Sensitivity (Recall) = доля правильно найденных больных:
# TP / (TP + FN) = cm["1","1"] / (cm["1","1"] + cm["0","1"])
sensitivity <- cm["1", "1"] / sum(cm[, "1"])

# Specificity = доля правильно найденных здоровых:
# TN / (TN + FP) = cm["0","0"] / (cm["0","0"] + cm["1","0"])
specificity <- cm["0", "0"] / sum(cm[, "0"])

cat(sprintf("\nAccuracy : %.3f\n", accuracy))
cat(sprintf("Sensitivity: %.3f\n", sensitivity))
cat(sprintf("Specificity: %.3f\n", specificity))

## Блок 7. ROC-кривая и AUC
## ROC показывает качество модели при всех возможных порогах.

roc_obj <- roc(
  response  = test$diabetes,  # истинные значения (0/1)
  predictor = test$prob       # предсказанные вероятности
)

auc_value <- auc(roc_obj)
cat(sprintf("\nAUC: %.3f\n", auc_value))

plot(roc_obj,
     col  = "blue",
     lwd  = 2,
     main = paste("ROC (Pima, AUC =", round(auc_value, 3), ")"))

## Блок 8. Упражнение: влияние порога (0.3 и 0.7)
## Показываем trade-off между Sensitivity и Specificity при разных порогах.

check_threshold <- function(thr, probs, true_class) {
  # Перевод вероятностей в классы по порогу thr
  pred_thr <- ifelse(probs >= thr, 1, 0)
  cm_thr   <- table(pred = pred_thr, true = true_class)

  accuracy_thr    <- sum(diag(cm_thr)) / sum(cm_thr)
  sensitivity_thr <- cm_thr["1", "1"] / sum(cm_thr[, "1"])
  specificity_thr <- cm_thr["0", "0"] / sum(cm_thr[, "0"])

  cat("\n=== Порог:", thr, "===\n")
  print(cm_thr)
  cat(sprintf("Accuracy : %.3f\n",  accuracy_thr))
  cat(sprintf("Sensitivity: %.3f\n", sensitivity_thr))
  cat(sprintf("Specificity: %.3f\n", specificity_thr))
}

check_threshold(0.3, test$prob, test$diabetes)
check_threshold(0.5, test$prob, test$diabetes)
check_threshold(0.7, test$prob, test$diabetes)

plot(x_grid, y_grid,
     type = "l", lwd = 2, col = "blue",
     xlab = "Линейный предиктор η",
     ylab = "Вероятность p",
     main = "Сигмоидная функция и реальные предсказания")

# Добавляем реальные точки (η, p) для тестовых наблюдений
points(test$eta, test$prob,
       pch = 16, col = rgb(1, 0, 0, 0.5))

# Линия порога p = 0.5 (η = 0)
abline(h = 0.5, lty = 2, col = "gray")
abline(v = 0,   lty = 2, col = "gray")
abline(h = 0.3, lty = 2, col = "red")
abline(h = 0.7, lty = 2, col = "blue")

cols <- ifelse(test$diabetes == 1, "red", "darkgreen")
points(test$eta, test$prob, pch = 16, col = adjustcolor(cols, alpha.f = 0.6))
legend("topleft",
       legend = c("diabetes = 1", "diabetes = 0"),
       col    = c("red", "darkgreen"),
       pch    = 16)
