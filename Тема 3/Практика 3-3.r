## Практика 3. Многомерная Cox, проверка PH, C-index,
##            интеграция с классификацией
## Данные: TCGA-LUAD (клинические данные из RTCGA)

## 1. Подготовка окружения ------------------------------------------------

# options(repos = "https://mirror.truenetwork.ru/CRAN/")
# install.packages(c("survival", "survminer", "dplyr", "ggplot2",
#                     "broom", "caret", "pROC", "survcomp",
#                     "timeROC", "randomForest"))
# if (!requireNamespace("BiocManager")) install.packages("BiocManager")
# BiocManager::install("RTCGA")
# BiocManager::install("RTCGA.clinical")
# BiocManager::install("survcomp")

library(survival)
library(survminer)
library(dplyr)
library(ggplot2)
library(broom)

set.seed(123)

## 2. Загрузка и препроцессинг TCGA-LUAD ----------------------------------

library(RTCGA)
library(RTCGA.clinical)

data(LUAD.clinical)
df <- LUAD.clinical

df_surv <- df %>%
  mutate(
    time = dplyr::coalesce(
      as.numeric(patient.days_to_death),
      as.numeric(patient.days_to_last_followup)
    ),
    vital_clean = tolower(trimws(patient.vital_status)),
    status = ifelse(vital_clean == "dead", 1L, 0L),
    age = as.numeric(patient.age_at_initial_pathologic_diagnosis),
    stage_raw = tolower(trimws(patient.stage_event.pathologic_stage))
  ) %>%
  filter(!is.na(time), !is.na(status)) %>%
  mutate(
    stage_clean = gsub("^stage\\s+", "", stage_raw),
    stage_core  = sub("([ivx]+).*", "\\1", stage_clean),
    stage = case_when(
      stage_core == "i"   ~ "I",
      stage_core == "ii"  ~ "II",
      stage_core == "iii" ~ "III",
      stage_core == "iv"  ~ "IV",
      TRUE ~ NA_character_
    ),
    stage = factor(stage, levels = c("I", "II", "III", "IV")),
    smoking_raw = tolower(trimws(patient.tobacco_smoking_history)),
    smoking = case_when(
      grepl("lifelong non-smoker", smoking_raw) ~ "Never",
      grepl("current smoker", smoking_raw) ~ "Smoker",
      grepl("current reformed smoker", smoking_raw) ~ "Former",
      TRUE ~ NA_character_
    ),
    gender = case_when(
      tolower(patient.gender) == "male" ~ "Male",
      tolower(patient.gender) == "female" ~ "Female",
      TRUE ~ NA_character_
    ),
    time_years = time / 365.25
  )

df_surv$smoking <- factor(df_surv$smoking, levels = c("Never", "Former", "Smoker"))
df_surv$gender <- factor(df_surv$gender, levels = c("Female", "Male"))
# Полный набор (без NA в ключевых предикторах)
df_complete <- df_surv %>%
  filter(!is.na(age), !is.na(stage), !is.na(smoking), !is.na(gender))

cat("Полный набор для многомерной модели:", nrow(df_complete), "пациентов\n")
cat("Событий:", sum(df_complete$status), "\n\n")

## =====================================================================
## ЧАСТЬ 1. МНОГОМЕРНАЯ COX-МОДЕЛЬ
## =====================================================================

## 3. Построение многомерной Cox ------------------------------------------

cox_multi <- coxph(
  Surv(time_years, status) ~ age + stage + smoking + gender,
  data = df_complete
)

cat("--- Результаты многомерной Cox ---\n")
print(summary(cox_multi))

## 4. Forest plot многомерной модели --------------------------------------

ggforest(cox_multi, data = df_complete,
         main = "Многомерная Cox: TCGA-LUAD")

## =====================================================================
## ЧАСТЬ 2. C-INDEX
## =====================================================================

## 5. C-index из модели ---------------------------------------------------

# Concordance index — встроенная метрика из coxph
concordance_result <- concordance(cox_multi)
c_index <- concordance_result$concordance
c_se    <- sqrt(concordance_result$var)

cat(sprintf("C-index:          %.4f\n", c_index))
cat(sprintf("95%% CI:           [%.4f; %.4f]\n",
            c_index - 1.96 * c_se,
            c_index + 1.96 * c_se))
cat(sprintf("SE:               %.4f\n\n", c_se))

cat("Интерпретация C-index:\n")
cat("  0.5 — модель не лучше случайной\n")
cat("  0.6-0.7 — слабая модель\n")
cat("  0.7-0.8 — приемлемая модель\n")
cat("  >= 0.8 — отличная модель\n\n")

## 5b. C-index для одномерных моделей (сравнение) -------------------------

cox_age_uni   <- coxph(Surv(time_years, status) ~ age,   data = df_complete)
cox_stage_uni <- coxph(Surv(time_years, status) ~ stage, data = df_complete)

c_age   <- concordance(cox_age_uni)$concordance
c_stage <- concordance(cox_stage_uni)$concordance
c_multi <- c_index

cat("Сравнение C-index:\n")
cat(sprintf("  Univariate (age):   %.4f\n", c_age))
cat(sprintf("  Univariate (stage): %.4f\n", c_stage))
cat(sprintf("  Multivariate:       %.4f\n\n", c_multi))

## =====================================================================
## ЧАСТЬ 3. ПРОВЕРКА PH-ПРЕДПОЛОЖЕНИЯ
## =====================================================================

## 6. cox.zph — тест пропорциональности рисков ----------------------------

ph_test <- cox.zph(cox_multi)
cat("--- Тест Шоенфельда ---\n")
print(ph_test)

cat("\nИнтерпретация:\n")
cat("  p > 0.05 — PH-предположение НЕ нарушено (OK)\n")
cat("  p < 0.05 — PH-предположение нарушено (проблема)\n\n")

## 7. Графики остатков Шоенфельда -----------------------------------------

par(mfrow = c(2, 2))
plot(ph_test)
par(mfrow = c(1, 1))

# Отдельные графики для каждой переменной
for (i in seq_len(ncol(ph_test$y))) {
  plot(ph_test, var = i,
       main = paste("Schoenfeld residuals:", colnames(ph_test$y)[i]))
}

## 8. Обсуждение нарушений PH ---------------------------------------------

cat("--- Нарушения PH ---\n")
ph_pvals <- ph_test$table[, "p"]
violated <- names(ph_pvals)[ph_pvals < 0.05]
if (length(violated) > 0) {
  cat("Переменные с нарушением PH (p < 0.05):\n")
  for (v in violated) {
    cat(sprintf("  %s: p = %.4f\n", v, ph_pvals[v]))
  }
  cat("\nВозможные решения:\n")
  cat("  1. Стратификация по этой переменной\n")
  cat("  2. Добавление time-varying coefficient\n")
  cat("  3. Разбиение анализа по временным интервалам\n")
} else {
  cat("PH-предположение не нарушено ни для одной переменной\n")
}

## 9. Проблема отделения и редких уровней ----------------------------------

cat("\n--- Проверка редких уровней ---\n")
cat("Стадии:\n")
print(table(df_complete$stage))
cat("Курение:\n")
print(table(df_complete$smoking))

# Если стадия IV очень мала, можно объединить III и IV
stage_counts <- table(df_complete$stage)
if (any(stage_counts < 30)) {
  cat("Рекомендация: объединить (collapsing) редкие категории\n")
  
  # Пример collapsing
  df_complete <- df_complete %>%
    mutate(
      stage_collapsed = case_when(
        stage %in% c("I")    ~ "I",
        stage %in% c("II")    ~ "II",
        stage %in% c("III", "IV")  ~ "Late (III-IV)"
      ),
      stage_collapsed = factor(stage_collapsed,
                               levels = c("I", "II", "Late (III-IV)"))
    )
  
  cox_collapsed <- coxph(
    Surv(time_years, status) ~ age + stage_collapsed + smoking + gender,
    data = df_complete
  )
  cat("\nМодель с объединёнными стадиями:\n")
  print(summary(cox_collapsed))
}

## =====================================================================
## ЧАСТЬ 4. ИНТЕГРАЦИЯ С КЛАССИФИКАЦИЕЙ
## =====================================================================

## 10. Подготовка данных для классификации ---------------------------------

# Бинарный исход: умер в течение 5 лет vs жив/цензурирован > 5 лет
horizon <- 5  # лет

df_class <- df_complete %>%
  filter(time_years >= horizon | status == 1) %>%
  mutate(
    outcome = ifelse(status == 1 & time_years <= horizon, 1L, 0L),
    outcome = factor(outcome, levels = c(0, 1), labels = c("alive", "dead"))
  )

cat("Классификация (горизонт =", horizon, "лет):\n")
cat("N =", nrow(df_class), "\n")
print(table(df_class$outcome))

## 11. Логистическая регрессия (классификатор) -----------------------------

library(caret)
library(pROC)

# Train/Test split
train_idx <- createDataPartition(df_class$outcome, p = 0.7, list = FALSE)
train_data <- df_class[train_idx, ]
test_data  <- df_class[-train_idx, ]

# Логистическая регрессия
logit_model <- glm(
  outcome ~ age + stage + smoking + gender,
  data = train_data,
  family = binomial
)
cat("\n--- Логистическая регрессия ---\n")
print(summary(logit_model))

# Предсказания
prob_logit <- predict(logit_model, newdata = test_data, type = "response")
pred_logit <- ifelse(prob_logit > 0.5, "dead", "alive")
pred_logit <- factor(pred_logit, levels = c("alive", "dead"))

# Метрики классификатора
cm <- confusionMatrix(pred_logit, test_data$outcome, positive = "dead")
cat("\nConfusion Matrix:\n")
print(cm)

roc_logit <- roc(test_data$outcome, prob_logit)
auc_logit <- auc(roc_logit)
cat(sprintf("\nAUC (логистическая регрессия): %.4f\n", auc_logit))

## 12. Random Forest (классификатор) --------------------------------------

library(randomForest)

rf_data_train <- train_data %>%
  select(age, stage, smoking, gender, outcome) %>%
  na.omit()
rf_data_test <- test_data %>%
  select(age, stage, smoking, gender, outcome) %>%
  na.omit()

rf_model <- randomForest(
  outcome ~ age + stage + smoking + gender,
  data  = rf_data_train,
  ntree = 200,
  importance = TRUE
)

prob_rf <- predict(rf_model, newdata = rf_data_test, type = "prob")[, "dead"]
pred_rf <- predict(rf_model, newdata = rf_data_test)

roc_rf  <- roc(rf_data_test$outcome, prob_rf)
auc_rf  <- auc(roc_rf)
cat(sprintf("AUC (random forest):          %.4f\n", auc_rf))

## 13. Cox-модель на тех же данных ----------------------------------------

cox_model_cls <- coxph(
  Surv(time_years, status) ~ age + stage+ smoking + gender,
  data = df_class
)

c_index_cls <- concordance(cox_model_cls)$concordance
cat(sprintf("C-index (Cox на тех же данных): %.4f\n", c_index_cls))

## 14. Сравнение метрик: Cox vs Классификация ------------------------------

comparison <- data.frame(
  Model = c("Cox regression", "Logistic Regression", "Random Forest"),
  Metric = c("C-index", "AUC", "AUC"),
  Value = round(c(c_index_cls, as.numeric(auc_logit), as.numeric(auc_rf)), 4)
)
print(comparison)

cat("\nОбсуждение различий:\n")
cat("  Cox: моделирует ВРЕМЯ до события (hazard)\n")
cat("  Классификация: моделирует СТАТУС на фиксированном горизонте\n")
cat("  C-index ~ AUC, но считается по разному\n")
cat("  Cox учитывает цензурирование, классификатор — нет\n\n")

## 15. ROC-кривые сравнения -----------------------------------------------

plot(roc_logit, col = "steelblue", lwd = 2,
     main = "ROC: классификация vs Cox (2-летний горизонт)")
lines(roc_rf, col = "forestgreen", lwd = 2)
abline(a = 0, b = 1, lty = 2, col = "gray")

legend("bottomright",
       legend = c(
         sprintf("LogReg (AUC=%.3f)", auc_logit),
         sprintf("RF (AUC=%.3f)", auc_rf)
       ),
       col = c("steelblue", "forestgreen"),
       lwd = 2)

## =====================================================================
## ЧАСТЬ 5. МИНИ-ПРОЕКТ: РИСК-СКОР
## =====================================================================

## 16. Экспорт риск-скора из классификатора --------------------------------

# Риск-скор = вероятность из логистической регрессии
df_class$risk_score_logit <- predict(logit_model,
                                     newdata = df_class,
                                     type = "response")

## 17. Стратификация KM по риск-группам -----------------------------------

# Разделение на high/low risk по медиане
median_risk <- median(df_class$risk_score_logit, na.rm = TRUE)

df_class <- df_class %>%
  mutate(
    risk_group = ifelse(risk_score_logit > median_risk, "High risk", "Low risk"),
    risk_group = factor(risk_group, levels = c("Low risk", "High risk"))
  )

cat("Группы риска:\n")
print(table(df_class$risk_group))

# KM по группам риска
km_risk <- survfit(Surv(time_years, status) ~ risk_group, data = df_class)

ggsurvplot(
  km_risk,
  data         = df_class,
  pval         = TRUE,
  pval.method  = TRUE,
  conf.int     = TRUE,
  risk.table   = TRUE,
  risk.table.height = 0.2,
  palette      = c("#2E9FDF", "#FC4E07"),
  legend.labs  = c("Low risk", "High risk"),
  legend.title = "Группа риска",
  xlab = "Время (годы)",
  ylab = "Вероятность выживаемости",
  title = "KM по группам риска (на основе логистической регрессии)"
)

## 18. Риск-скор как ковариата в Cox ---------------------------------------

cat("\n--- Cox с риск-скором как ковариатой ---\n")
cox_with_risk <- coxph(
  Surv(time_years, status) ~ risk_score_logit + age + stage + smoking + gender,
  data = df_class
)
print(summary(cox_with_risk))

c_risk <- concordance(cox_with_risk)$concordance
cat(sprintf("C-index (Cox + risk score): %.4f\n", c_risk))

library(timeROC)

## 1. Линейный предиктор Cox (risk score)
df_complete$risk_score <- predict(cox_multi, type = "lp")

## 2. Выбираем временные точки (в годах)
times_years <- c(2, 5, 10, 12, 15)

## Переводим в те же единицы, что и time (если time в годах — можно так и оставить)
times <- times_years

## 3. Расчёт time-dependent ROC и AUC
roc_cox <- timeROC(
  T      = df_complete$time_years,
  delta  = df_complete$status,
  marker = df_complete$risk_score,
  cause  = 1,                  # событие = смерть
  times  = times,
  iid    = TRUE
)

roc_cox$AUC    # AUC(t) для выбранных времён
roc_cox$times  # какие t реально использованы

cols <- c("steelblue", "darkorange", "forestgreen", "red", "yellow")

# первая кривая
plot(
  roc_cox$FP[, 1],  roc_cox$TP[, 1],
  type = "l", lwd = 2, col = cols[1],
  xlab = "Специфичность",
  ylab = "Чувствительность",
  main = "Time-dependent ROC для Cox-модели (1, 3, 5 лет)"
)

# остальные
lines(roc_cox$FP[, 2], roc_cox$TP[, 2], col = cols[2], lwd = 2)
lines(roc_cox$FP[, 3], roc_cox$TP[, 3], col = cols[3], lwd = 2)
lines(roc_cox$FP[, 4], roc_cox$TP[, 4], col = cols[4], lwd = 2)
lines(roc_cox$FP[, 5], roc_cox$TP[, 5], col = cols[5], lwd = 2)

abline(0, 1, lty = 2, col = "gray")

legend(
  "bottomright",
  legend = paste0("t = ", times, ", AUC = ",
                  sprintf("%.2f", roc_cox$AUC)),
  col = cols, lwd = 2, bty = "n"
)

## 19. Итоговое сравнение всех подходов ------------------------------------

cat("\n============================\n")
cat("  ИТОГОВОЕ СРАВНЕНИЕ\n")
cat("============================\n\n")

final_table <- data.frame(
  Approach = c(
    "Multivariate Cox (clinical)",
    "Cox + risk score from LogReg",
    "Logistic Regression (5y)",
    "Random Forest (5y)"
  ),
  Metric = c("C-index", "C-index", "AUC", "AUC"),
  Value = round(c(c_index_cls, c_risk,
                   as.numeric(auc_logit), as.numeric(auc_rf)), 4)
)
print(final_table)

cat("\nВыводы:\n")
cat("1. Cox-модель моделирует время, классификатор — статус\n")
cat("2. C-index сопоставим с AUC, но учитывает цензурирование\n")
cat("3. Риск-скор из классификатора можно использовать в Cox\n")
cat("4. Стратификация KM по risk score визуализирует различия\n")
cat("5. PH-предположение нужно проверять (cox.zph)\n")
cat("6. Редкие категории вызывают проблему separation\n")

