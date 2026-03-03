## Блок 0. Подготовка окружения

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

set.seed(123)

## Блок 1. Загрузка и первичный осмотр

data(PimaIndiansDiabetes)
pid <- PimaIndiansDiabetes

cat("Размер датасета (строки x столбцы):\n")
print(dim(pid))

cat("\nСтруктура:\n")
str(pid)

cat("\nРаспределение классов (neg/pos):\n")
print(table(pid$diabetes))

## Блок 2. Препроцессинг

pid_clean <- pid %>%
  filter(complete.cases(.)) %>%
  mutate(
    diabetes = ifelse(diabetes == "pos", 1, 0)
  )

cat("\nПосле очистки наблюдений:\n")
print(nrow(pid_clean))

cat("\nРаспределение классов (0/1):\n")
print(table(pid_clean$diabetes))

## Блок 3. Train / Test split

n <- nrow(pid_clean)
train_size <- floor(0.7 * n)

train_idx <- sample(seq_len(n), size = train_size)

train <- pid_clean[train_idx, ]
test  <- pid_clean[-train_idx, ]

cat("\nTrain:", nrow(train), "наблюдений\n")
cat("Test :", nrow(test),  "наблюдений\n\n")

cat("Классы в train (0/1):\n")
print(table(train$diabetes))

cat("\nКлассы в test (0/1):\n")
print(table(test$diabetes))

## Блок 4. Модель 1 — логистическая регрессия

model_logit <- glm(
  diabetes ~ .,
  data   = train,
  family = binomial(link = "logit")
)

summary(model_logit)

test <- test %>%
  mutate(
    prob_logit = predict(model_logit, newdata = ., type = "response")
  )

## Блок 5. Модель 2 — дерево решений

model_tree <- rpart(
  factor(diabetes) ~ .,
  data   = train,
  method = "class",
  control = rpart.control(cp = 0.01, minsplit = 20)
)

rpart.plot(model_tree)

prob_tree <- predict(model_tree, newdata = test, type = "prob")[, "1"]
test <- test %>%
  mutate(prob_tree = prob_tree)

## Блок 6. Модель 3 — Random Forest

model_rf <- randomForest(
  factor(diabetes) ~ .,
  data  = train,
  ntree = 100,
  mtry  = 3
)

prob_rf <- predict(model_rf, newdata = test, type = "prob")[, "1"]
test <- test %>%
  mutate(prob_rf = prob_rf)

## Блок 7. Функция для метрик

metrics_for_model <- function(probs, true, threshold = 0.5) {
  pred <- ifelse(probs >= threshold, 1, 0)
  cm   <- table(pred = pred, true = true)

  accuracy    <- sum(diag(cm)) / sum(cm)
  sensitivity <- cm["1", "1"] / sum(cm[, "1"])
  specificity <- cm["0", "0"] / sum(cm[, "0"])

  roc_obj <- roc(true, probs)
  auc_val <- as.numeric(auc(roc_obj))

  c(Accuracy = accuracy,
    Sensitivity = sensitivity,
    Specificity = specificity,
    AUC = auc_val)
}

## Блок 8. Сводная таблица по трём моделям

true <- test$diabetes

res_logit <- metrics_for_model(test$prob_logit, true)
res_tree  <- metrics_for_model(test$prob_tree,  true)
res_rf    <- metrics_for_model(test$prob_rf,    true)

metrics_table <- rbind(
  LogReg = res_logit,
  Tree   = res_tree,
  RF     = res_rf
)

round(metrics_table, 3)

## Блок 9. Общий ROC‑график

roc_logit <- roc(true, test$prob_logit)
roc_tree  <- roc(true, test$prob_tree)
roc_rf    <- roc(true, test$prob_rf)

plot(roc_logit, col = "blue", lwd = 2,
     main = "ROC-кривые: LogReg (синий), Tree (красный), RF (зелёный)")
plot(roc_tree, col = "red", lwd = 2, add = TRUE)
plot(roc_rf, col = "darkgreen", lwd = 2, add = TRUE)
abline(a = 0, b = 1, lty = 2, col = "gray")

legend("bottomright",
       legend = c("LogReg", "Tree", "RF"),
       col    = c("blue", "red", "darkgreen"),
       lwd    = 2)

## Блок 10. Упражнение: влияние порога

check_threshold <- function(thr, probs, true_class, model_name) {
  pred_thr <- ifelse(probs >= thr, 1, 0)
  cm_thr   <- table(pred = pred_thr, true = true_class)

  accuracy_thr    <- sum(diag(cm_thr)) / sum(cm_thr)
  sensitivity_thr <- cm_thr["1", "1"] / sum(cm_thr[, "1"])
  specificity_thr <- cm_thr["0", "0"] / sum(cm_thr[, "0"])

  cat("\n=== Модель:", model_name, "| Порог:", thr, "===\n")
  print(cm_thr)
  cat(sprintf("Accuracy : %.3f\n",  accuracy_thr))
  cat(sprintf("Sensitivity: %.3f\n", sensitivity_thr))
  cat(sprintf("Specificity: %.3f\n", specificity_thr))
}

thresholds <- c(0.3, 0.5, 0.7)

for (thr in thresholds) {
  check_threshold(thr, test$prob_logit, true, "LogReg")
  check_threshold(thr, test$prob_tree,  true, "Tree")
  check_threshold(thr, test$prob_rf,    true, "RF")
}
