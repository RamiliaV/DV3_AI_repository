## Практика 3. Сравнение моделей и Grid Search на Cleveland Heart Disease
## Модели: Logistic Regression, Decision Tree, Random Forest
## Метрики: AUC-ROC, F1

## 1. Подготовка окружения ---------------------------------------------

# options(repos = "https://mirror.truenetwork.ru/CRAN/")
# install.packages(c("kmed", "dplyr", "caret", "randomForest", "pROC", "rpart", "corrplot"))

library(kmed)
library(dplyr)
library(caret)
library(randomForest)
library(pROC)
library(rpart)
library(corrplot)

set.seed(123)

## 2. Загрузка и подготовка данных (Cleveland Heart Disease) -----------

data("heart")
hd <- heart

# Проверяем структуру
str(hd)
# Целевая переменная обычно называется 'class' в kmed::heart
# 0 = нет болезни, 1-4 = есть болезнь

hd_clean <- hd %>%
  filter(complete.cases(.)) %>%
  mutate(
    disease = ifelse(class > 0, "disease", "healthy"),
    disease = factor(disease, levels = c("healthy", "disease"))
  ) %>%
  select(-class)

# Целевая переменная и признаки
target_col   <- "disease"
feature_cols <- setdiff(names(hd_clean), target_col)

table(hd_clean$disease)
prop.table(table(hd_clean$disease))

## 2b. EDA: базовая статистика -----------------------------------------

summary(hd_clean[, feature_cols])

## 3. Train / Test split (70 / 30) -------------------------------------

train_idx <- createDataPartition(hd_clean$disease, p = 0.7, list = FALSE)

train_raw <- hd_clean[train_idx, ]
test_raw  <- hd_clean[-train_idx, ]

table(train_raw$disease)
table(test_raw$disease)
prop.table(table(train_raw$disease))
prop.table(table(test_raw$disease))

## 3b. Препроцессинг: center/scale для числовых признаков -------------

num_cols <- names(train_raw)[sapply(train_raw, is.numeric)]
num_cols <- setdiff(num_cols, target_col)

pre_proc <- preProcess(train_raw[, num_cols],
                       method = c("center", "scale"))

train <- train_raw
train[, num_cols] <- predict(pre_proc, train_raw[, num_cols])

test <- test_raw
test[, num_cols] <- predict(pre_proc, test_raw[, num_cols])

## 3c. EDA: матрица корреляций по числовым признакам -------------------

if (length(num_cols) > 1) {
  corr_mat <- cor(train[, num_cols])
  corrplot(corr_mat, method = "color", tl.cex = 0.7)
}

## 4. Общая настройка caret: k-fold CV и метрики -----------------------

ctrl <- trainControl(
  method          = "cv",
  number          = 5,
  summaryFunction = twoClassSummary,
  classProbs      = TRUE,
  savePredictions = "final"
)

## Функция для F1 ------------------------------------------------------

f1_score <- function(true, pred_class, positive = "disease") {
  true <- factor(true, levels = c("healthy", "disease"))
  pred_class <- factor(pred_class, levels = c("healthy", "disease"))
  cm <- table(pred_class, true)
  tp <- cm[positive, positive]
  fp <- sum(cm[positive, ]) - tp
  fn <- sum(cm[, positive]) - tp
  precision <- tp / (tp + fp)
  recall    <- tp / (tp + fn)
  2 * precision * recall / (precision + recall)
}

## 5. Модель 1: Логистическая регрессия -------------------------------

model_glm <- train(
  disease ~ .,
  data      = train,
  method    = "glm",
  family    = binomial,
  trControl = ctrl,
  metric    = "ROC"
)

model_glm

prob_glm <- predict(model_glm, newdata = test, type = "prob")[, "disease"]
pred_glm <- predict(model_glm, newdata = test, type = "raw")

roc_glm <- roc(test$disease, prob_glm)
auc_glm <- auc(roc_glm)
f1_glm  <- f1_score(test$disease, pred_glm)

auc_glm
f1_glm

## 6. Модель 2: Decision Tree (rpart) ----------------------------------

grid_rpart <- expand.grid(cp = c(0.001, 0.01, 0.05))

model_rpart <- train(
  disease ~ .,
  data      = train,
  method    = "rpart",
  trControl = ctrl,
  tuneGrid  = grid_rpart,
  metric    = "ROC"
)

model_rpart
plot(model_rpart)

prob_rpart <- predict(model_rpart, newdata = test, type = "prob")[, "disease"]
pred_rpart <- predict(model_rpart, newdata = test, type = "raw")

roc_rpart <- roc(test$disease, prob_rpart)
auc_rpart <- auc(roc_rpart)
f1_rpart  <- f1_score(test$disease, pred_rpart)

auc_rpart
f1_rpart

## 7. Модель 3: Random Forest + Grid Search по mtry --------------------

p <- ncol(train) - 1
grid_rf <- expand.grid(mtry = c(2, 4, 6, 8))

model_rf <- train(
  disease ~ .,
  data      = train,
  method    = "rf",
  trControl = ctrl,
  tuneGrid  = grid_rf,
  metric    = "ROC",
  ntree     = 200
)

model_rf
plot(model_rf)

best_mtry <- model_rf$bestTune$mtry
best_mtry

prob_rf <- predict(model_rf, newdata = test, type = "prob")[, "disease"]
pred_rf <- predict(model_rf, newdata = test, type = "raw")

roc_rf <- roc(test$disease, prob_rf)
auc_rf <- auc(roc_rf)
f1_rf  <- f1_score(test$disease, pred_rf)

auc_rf
f1_rf

## 8. Итоговое сравнение моделей --------------------------------------

results <- data.frame(
  model = c("Logistic Regression", "Decision Tree (rpart)", "Random Forest"),
  AUC   = c(as.numeric(auc_glm),
            as.numeric(auc_rpart),
            as.numeric(auc_rf)),
  F1    = c(f1_glm, f1_rpart, f1_rf)
)

results

## 9. Совместный график ROC-кривых ------------------------------------

plot(roc_glm,   col = "steelblue", lwd = 2,
     main = "ROC-кривые: сравнение моделей (Heart Disease)")
lines(roc_rpart, col = "darkorange",  lwd = 2)
lines(roc_rf,    col = "forestgreen", lwd = 2)
abline(a = 0, b = 1, lty = 2, col = "gray")

legend("bottomright",
       legend = c("LogReg", "Decision Tree", "Random Forest"),
       col    = c("steelblue", "darkorange", "forestgreen"),
       lwd    = 2)

## 10. Feature importance для Random Forest ----------------------------

imp_rf_gini <- varImp(model_rf, scale = FALSE)
imp_rf_gini
plot(imp_rf_gini, top = 10,
     main = "Random Forest: Gini feature importance")

## 11. Permutation importance для RF (по AUC) --------------------------

perm_importance_auc <- function(model, data, target = "disease") {
  base_prob <- predict(model, newdata = data, type = "prob")[, "disease"]
  base_auc  <- as.numeric(auc(roc(data[[target]], base_prob)))
  
  res <- c()
  feature_names <- setdiff(names(data), target)
  
  for (feat in feature_names) {
    data_perm <- data
    data_perm[[feat]] <- sample(data_perm[[feat]])
    prob_perm <- predict(model, newdata = data_perm, type = "prob")[, "disease"]
    auc_perm  <- as.numeric(auc(roc(data_perm[[target]], prob_perm)))
    res[feat] <- base_auc - auc_perm
  }
  
  sort(res, decreasing = TRUE)
}

imp_rf_perm <- perm_importance_auc(model_rf, test, target = "disease")
imp_rf_perm

## 12. Сводная таблица важности признаков ------------------------------

gini_df <- imp_rf_gini$importance
gini_df$feature <- rownames(gini_df)
colnames(gini_df)[1] <- "gini_importance"

imp_df <- gini_df %>%
  mutate(
    perm_importance = imp_rf_perm[feature],
    rank_gini       = rank(-gini_importance, ties.method = "min"),
    rank_perm       = rank(-perm_importance, ties.method = "min")
  ) %>%
  arrange(rank_perm)

imp_df

## 12b. EDA: boxplot для самого важного признака -----------------------

top_feature <- imp_df$feature[1]
top_feature

if (is.numeric(train[[top_feature]])) {
  boxplot(
    train[[top_feature]] ~ train$disease,
    main = paste("Распределение", top_feature, "по классам"),
    xlab = "disease",
    ylab = top_feature,
    col  = c("lightblue", "lightpink")
  )
}

## 13. Топ-k признаков и RF на топ-k ----------------------------------

k <- min(5, length(feature_cols))

topk <- head(imp_df$feature[!is.na(imp_df$perm_importance)], k)
topk

train_topk <- train[, c(topk, "disease")]
test_topk  <- test[,  c(topk, "disease")]

model_rf_topk <- train(
  disease ~ .,
  data      = train_topk,
  method    = "rf",
  trControl = ctrl,
  tuneGrid  = grid_rf,
  metric    = "ROC",
  ntree     = 200
)

model_rf_topk

prob_rf_topk <- predict(model_rf_topk, newdata = test_topk, type = "prob")[, "disease"]
pred_rf_topk <- predict(model_rf_topk, newdata = test_topk, type = "raw")

roc_rf_topk <- roc(test_topk$disease, prob_rf_topk)
auc_rf_topk <- auc(roc_rf_topk)
f1_rf_topk  <- f1_score(test_topk$disease, pred_rf_topk)

data.frame(
  model = c("RF_all_features", paste0("RF_top", k)),
  AUC   = c(as.numeric(auc_rf), as.numeric(auc_rf_topk)),
  F1    = c(f1_rf, f1_rf_topk)
)

## 14. Эксперимент: различное число деревьев ---------------------------

rf_50 <- train(
  disease ~ .,
  data      = train,
  method    = "rf",
  trControl = ctrl,
  tuneGrid  = data.frame(mtry = best_mtry),
  metric    = "ROC",
  ntree     = 50
)

rf_200 <- model_rf

rf_1000 <- train(
  disease ~ .,
  data      = train,
  method    = "rf",
  trControl = ctrl,
  tuneGrid  = data.frame(mtry = best_mtry),
  metric    = "ROC",
  ntree     = 1000
)

models_ntree <- list(
  rf_50   = rf_50,
  rf_200  = rf_200,
  rf_1000 = rf_1000
)

ntree_results <- lapply(names(models_ntree), function(name) {
  m <- models_ntree[[name]]
  prob <- predict(m, newdata = test, type = "prob")[, "disease"]
  pred <- predict(m, newdata = test, type = "raw")
  data.frame(
    model = name,
    AUC   = as.numeric(auc(roc(test$disease, prob))),
    F1    = f1_score(test$disease, pred)
  )
})

ntree_results <- do.call(rbind, ntree_results)
ntree_results
