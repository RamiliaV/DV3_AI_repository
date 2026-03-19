## Практика 2. Одномерные модели Кокса
## Данные: TCGA-LUAD (клинические данные из RTCGA)
## Цель: Surv-объект, univariate Cox для age/stage/ER/PR/HER2,
##       оценка HR, 95% CI, p-value, таблица значимых предикторов

## 1. Подготовка окружения ------------------------------------------------

# options(repos = "https://mirror.truenetwork.ru/CRAN/")
# install.packages(c("survival", "survminer", "dplyr", "ggplot2", "broom"))
# if (!requireNamespace("BiocManager")) install.packages("BiocManager")
# BiocManager::install("RTCGA")
# BiocManager::install("RTCGA.clinical")

library(survival)
library(survminer)
library(dplyr)
library(ggplot2)
library(broom)       # tidy() для красивых таблиц

set.seed(123)

## 2. Загрузка данных TCGA-LUAD -----------------------------------------

library(RTCGA)
library(RTCGA.clinical)

data(LUAD.clinical)
df <- LUAD.clinical

## 3. Препроцессинг -------------------------------------------------------

df_surv <- df %>%
  mutate(
    # время: дни до смерти или последнего визита
    time = dplyr::coalesce(
      as.numeric(patient.days_to_death),
      as.numeric(patient.days_to_last_followup)
    ),
    # статус: 1 = dead, 0 = цензурирован
    vital_clean = tolower(trimws(patient.vital_status)),
    status = ifelse(vital_clean == "dead", 1L, 0L),
    # стадия из pathologic_stage
    stage_raw = patient.stage_event.pathologic_stage
  ) %>%
  filter(!is.na(time),
         !is.na(status),
         !is.na(stage_raw))

smoking_cols <- names(df)[grepl("smoking|tobacco", names(df), ignore.case = TRUE)]
gender_cols <- names(df)[grepl("gender|sex", names(df), ignore.case = TRUE)]

cat("Smoking переменные:\n")
cat(smoking_cols, sep = "\n")
cat("\nGender переменные:\n")
cat(gender_cols, sep = "\n")

## Препроцессинг smoking и gender
df_surv <- df_surv %>%
  mutate(
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
    )
  ) %>%
  filter(!is.na(time),
         !is.na(status),
         !is.na(stage_raw)) %>%
  mutate(age = patient.age_at_initial_pathologic_diagnosis)

df_surv$smoking <- factor(df_surv$smoking, levels = c("Never", "Former", "Smoker"))
df_surv$gender <- factor(df_surv$gender, levels = c("Female", "Male"))

cat("Пациентов после базовой фильтрации:", nrow(df_surv), "\n")

## 3b. Упрощение стадий --------------------------------------------------

df_surv <- df_surv %>%
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
    stage = factor(stage, levels = c("I", "II", "III", "IV"))
  )

# Время в годах
df_surv$time_years <- df_surv$time / 365.25

df_surv$age <- as.numeric(df_surv$age)

cat("\nСтруктура данных для анализа:\n")
cat("N =", nrow(df_surv), "\n")
cat("Событий:", sum(df_surv$status), "\n")
cat("Возраст: mean =", round(mean(df_surv$age, na.rm = TRUE), 1), "\n")
cat("Стадии:\n"); print(table(df_surv$stage, useNA = "ifany"))
cat("Smoking:\n"); print(table(df_surv$smoking, useNA = "ifany"))
cat("Gender:\n"); print(table(df_surv$gender, useNA = "ifany"))

## 4. Объект Surv ---------------------------------------------------------

surv_obj <- Surv(time = df_surv$time_years, event = df_surv$status)

## 5. Одномерные модели Кокса ---------------------------------------------

## 5a. Age ---------------------------------------------------------------

cat("\n--- Cox: age ---\n")
cox_age <- coxph(Surv(time_years, status) ~ age, data = df_surv)
print(summary(cox_age))

## 5b. Stage --------------------------------------------------------------

cat("\n--- Cox: stage ---\n")
cox_stage <- coxph(Surv(time_years, status) ~ stage, data = df_surv)
print(summary(cox_stage))

## 5c. Cox по smoking ----------------------------------------------------

cat("--- Cox smoking ---\n")
cox_smoking <- coxph(Surv(time_years, status) ~ smoking, data = df_surv)
print(summary(cox_smoking))

## 5d. Cox по gender ----------------------------------------------------

cat("--- Cox gender ---\n")
cox_gender <- coxph(Surv(time_years, status) ~ gender, data = df_surv)
print(summary(cox_gender))

## 6. Сводная таблица: HR, 95% CI, p-value --------------------------------

extract_cox <- function(model, name) {
  s <- summary(model)
  coef_df <- as.data.frame(s$conf.int)
  p_df    <- as.data.frame(s$coefficients)
  
  data.frame(
    predictor = rownames(coef_df),
    model_name = name,
    HR        = round(coef_df[, "exp(coef)"], 3),
    lower_CI  = round(coef_df[, "lower .95"], 3),
    upper_CI  = round(coef_df[, "upper .95"], 3),
    p_value   = signif(p_df[, "Pr(>|z|)"], 3),
    stringsAsFactors = FALSE
  )
}

cox_results <- bind_rows(
  extract_cox(cox_age,   "age"),
  extract_cox(cox_stage, "stage"),
  extract_cox(cox_smoking, "smoking"),
  extract_cox(cox_gender, "gender")
)

cat("\n============================\n")
cat("  СВОДНАЯ ТАБЛИЦА РЕЗУЛЬТАТОВ\n")
cat("============================\n\n")
print(cox_results, row.names = FALSE)

## 7. Отбор значимых предикторов ------------------------------------------

significant <- cox_results %>%
  filter(p_value < 0.05) %>%
  arrange(p_value)

cat("\n=== ЗНАЧИМЫЕ ПРЕДИКТОРЫ (p < 0.05) ===\n")
print(significant, row.names = FALSE)

## 8. Визуализация: forest plot для одномерных Cox -------------------------

# Forest plot для stage
cat("\n--- Forest plot: stage ---\n")
ggforest(cox_stage, data = df_surv,
         main = "Cox: стадия (TCGA-LUAD)")

cat("--- Forest plot smoking ---\n")
ggforest(cox_smoking, data = df_surv, main = "Cox: статус курения (TCGA-LUAD)")

cat("--- Forest plot gender ---\n")
ggforest(cox_gender, data = df_surv, main = "Cox: пол (TCGA-LUAD)")

## 9. Интерпретация HR ----------------------------------------------------

cat("\n============================\n")
cat("  ИНТЕРПРЕТАЦИЯ HR\n")
cat("============================\n\n")

cat("Что означает HR:\n")
cat("  HR = 1.0  — нет эффекта (нулевая гипотеза)\n")
cat("  HR > 1.0  — увеличение риска (фактор риска)\n")
cat("  HR < 1.0  — снижение риска (защитный фактор)\n")
cat("  HR = 2.0  — риск события в 2 раза выше\n")
cat("  HR = 0.5  — риск события в 2 раза ниже\n\n")

cat("Пример интерпретации для возраста:\n")
age_hr <- cox_results$HR[cox_results$model_name == "age"]
cat(sprintf("  HR(age) = %.3f\n", age_hr))
cat(sprintf("  Каждый дополнительный год возраста увеличивает\n"))
cat(sprintf("  риск смерти на %.1f%%\n", (age_hr - 1) * 100))

cat("Курение (Smoker vs Never):\n")
smoking_row <- cox_results[cox_results$predictor == "smokingSmoker", ]

if (nrow(smoking_row) > 0) {
  hr  <- smoking_row$HR[1]
  lcl <- smoking_row$lowerCI[1]
  ucl <- smoking_row$upperCI[1]
  p   <- smoking_row$p_value[1]

  if (!is.na(p)) {
    cat(sprintf(
      "HR(smoker vs never) = %.2f (95%% CI %.2f–%.2f), p = %.3f\n",
      hr, lcl, ucl, p
    ))

    if (p >= 0.05) {
      cat("Различия статистически незначимы (p >= 0.05),\n",
          "поэтому интерпретировать эффект курения как реальный нельзя.\n",
          sep = "")
    } else if (hr > 1) {
      cat(sprintf(
        "Оценка модели: курильщики имеют риск на %.0f%% выше.\n",
        (hr - 1) * 100
      ))
    } else if (hr < 1) {
      cat(sprintf(
        "Оценка модели: курильщики имеют риск на %.0f%% ниже.\n",
        (1 - hr) * 100
      ))
    } else {
      cat("Оценка модели: риск у курильщиков и некурящих примерно одинаков.\n")
    }
  } else {
    cat("p-value для эффекта курения не вычислен (NA).\n")
  }
} else {
  cat("Строка для smokingSmoker не найдена в cox_results.\n")
}



## 10. Связь с Заданием 6 ------------------------------------------------

cat("\n============================\n")
cat("  ПОДГОТОВКА К ЗАДАНИЮ 6\n")
cat("============================\n\n")
cat("Для Задания 6 (одномерный анализ):\n")
cat("1. Постройте одномерные Cox для age, stage, ER, PR, HER2\n")
cat("2. Извлеките HR, 95% CI, p-value\n")
cat("3. Определите значимые предикторы (p < 0.05)\n")
cat("4. Интерпретируйте HR клинически\n")
cat("5. Далее: многомерная модель на Практике 3\n")
