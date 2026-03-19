## Практика 1. Kaplan–Meier и Log-rank тест
## Данные: TCGA-LUAD (клинические данные из RTCGA)
## Цель: KM-кривые по стадиям, log-rank тест, медианы выживаемости

## 1. Подготовка окружения ------------------------------------------------

# options(repos = "https://mirror.truenetwork.ru/CRAN/")
# install.packages(c("survival", "survminer", "dplyr", "ggplot2"))
# if (!requireNamespace("BiocManager")) install.packages("BiocManager")
# BiocManager::install("RTCGA")
# BiocManager::install("RTCGA.clinical")

library(survival)    # Surv, survfit, survdiff
library(survminer)   # ggsurvplot — KM-кривые
library(dplyr)       # препроцессинг
library(ggplot2)

set.seed(123)

## 2. Загрузка данных TCGA-LUAD -----------------------------------------

library(RTCGA)
library(RTCGA.clinical)

data(LUAD.clinical)
df <- LUAD.clinical

cat("Размер датасета:", nrow(df), "x", ncol(df), "\n")

# Посмотрим на доступные переменные survival
cat("Ключевые переменные:\n")
cat(names(df)[grep("days_to|vital|stage", names(df))], sep = "\n")

## 3. Препроцессинг: формируем time, status, stage -----------------------

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
         !is.na(stage_raw)) %>%
  mutate(
    time = as.numeric(time),
    status = as.integer(status),
    stage_raw = tolower(trimws(stage_raw))
  )

cat("После базовой фильтрации:", nrow(df_surv), "пациентов\n")

## 4. Упрощение стадий: I, II, III, IV -----------------------------------

df_surv <- df_surv %>%
  mutate(
    stage_clean = gsub("^stage\\s+", "", stage_raw),
    stage_clean = trimws(stage_clean),
    stage_core  = sub("([ivx]+).*", "\\1", stage_clean),
    stage_simple = case_when(
      stage_core == "i"   ~ "I",
      stage_core == "ii"  ~ "II",
      stage_core == "iii" ~ "III",
      stage_core == "iv"  ~ "IV",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(stage_simple)) %>%
  mutate(stage_simple = factor(stage_simple,
                               levels = c("I", "II", "III", "IV")))

cat("После упрощения стадий:", nrow(df_surv), "пациентов\n")
cat("Стадии:\n")
print(table(df_surv$stage_simple))
cat("Событий (смертей):", sum(df_surv$status),
    "(", round(100 * mean(df_surv$status), 1), "%)\n\n")

## 5. Перевод времени в годы (для удобства графиков) ----------------------

df_surv <- df_surv %>%
  mutate(time_years = time / 365.25)

## 6. Объект Surv ---------------------------------------------------------

surv_obj <- Surv(time = df_surv$time_years, event = df_surv$status)

head(surv_obj, 20)

## 7. Kaplan–Meier: общая кривая (без стратификации) ----------------------

km_overall <- survfit(Surv(time_years, status) ~ 1, data = df_surv)

cat("=== ОБЩАЯ KM-КРИВАЯ ===\n")
print(km_overall)

# Медиана выживаемости (общая)
cat("\nМедиана выживаемости (общая):\n")
print(surv_median(km_overall))

# График
ggsurvplot(
  km_overall,
  data       = df_surv,
  conf.int   = TRUE,
  risk.table = TRUE,
  risk.table.height = 0.2,
  surv.median.line  = "hv",
  palette    = "#2E9FDF",
  xlab = "Время (годы)",
  ylab = "Вероятность выживаемости",
  title = "Kaplan–Meier: общая выживаемость TCGA-LUAD"
)

## 8. Kaplan–Meier по стадиям I–IV ----------------------------------------

km_stage <- survfit(surv_obj ~ stage_simple, data = df_surv)

cat("\n=== KAPLAN–MEIER ПО СТАДИЯМ ===\n")
print(km_stage)

# Медианы выживаемости по стадиям
cat("\n=== МЕДИАНЫ ВЫЖИВАЕМОСТИ ПО СТАДИЯМ ===\n")
medians <- surv_median(km_stage)
print(medians)

## 9. KM-график по стадиям (публикационное качество) ----------------------

km_plot <- ggsurvplot(
  km_stage,
  data         = df_surv,
  pval         = TRUE,              # p-value log-rank на графике
  pval.method  = TRUE,              # метка метода теста
  conf.int     = TRUE,              # доверительные интервалы
  risk.table   = TRUE,              # таблица числа пациентов под риском
  risk.table.height = 0.25,
  surv.median.line  = "hv",         # горизонтальные/вертикальные линии медиан
  palette      = c("#E7B800", "#2E9FDF", "#FC4E07", "#00BA38"),
  legend.labs  = levels(df_surv$stage_simple),
  legend.title = "Стадия",
  xlab = "Время (годы)",
  ylab = "Вероятность выживаемости",
  title = "Kaplan–Meier: выживаемость по стадиям (TCGA-LUAD)"
)

print(km_plot)

# Сохранение графика
png("km_stages_LUAD.png", width = 900, height = 750)
print(km_plot)
dev.off()
cat("График сохранён: km_stages_LUAD.png\n")

## 10. Log-rank тест ------------------------------------------------------

cat("\n=== LOG-RANK ТЕСТ ===\n")
logrank <- survdiff(surv_obj ~ stage_simple, data = df_surv)
print(logrank)

pval <- 1 - pchisq(logrank$chisq,
                    df = length(unique(df_surv$stage_simple)) - 1)
cat(sprintf("Chi-square: %.3f\n", logrank$chisq))
cat(sprintf("p-value:    %.2e\n", pval))

if (pval < 0.05) {
  cat("Значимые различия между стадиями (p < 0.05)\n\n")
} else {
  cat("Нет значимых различий (p >= 0.05)\n\n")
}

## 11. Попарные сравнения (pairwise log-rank) ------------------------------

cat("=== ПОПАРНЫЕ СРАВНЕНИЯ (log-rank) ===\n")
pairwise <- pairwise_survdiff(
  Surv(time_years, status) ~ stage_simple,
  data = df_surv,
  p.adjust.method = "BH"   # поправка Бенджамини–Хохберга
)
print(pairwise)

## 12. Сводная таблица результатов -----------------------------------------

cat("\n=== ИТОГОВАЯ ТАБЛИЦА ===\n")

summary_tbl <- df_surv %>%
  group_by(stage_simple) %>%
  summarise(
    n         = n(),
    events    = sum(status),
    pct_event = round(100 * mean(status), 1),
    .groups = "drop"
  )
print(summary_tbl)
cat(sprintf("\nОбщий log-rank p-value: %.2e\n", pval))

## 13. Дополнительно: KM для ранних vs поздних стадий ---------------------

df_surv <- df_surv %>%
  mutate(
    stage_group = ifelse(stage_simple %in% c("I", "II"), "Ранняя (I-II)",
                         "Поздняя (III-IV)"),
    stage_group = factor(stage_group,
                         levels = c("Ранняя (I-II)", "Поздняя (III-IV)"))
  )

km_group <- survfit(Surv(time_years, status) ~ stage_group, data = df_surv)

ggsurvplot(
  km_group,
  data         = df_surv,
  pval         = TRUE,
  pval.method  = TRUE,
  conf.int     = TRUE,
  risk.table   = TRUE,
  risk.table.height = 0.2,
  palette      = c("#2E9FDF", "#FC4E07"),
  legend.labs  = levels(df_surv$stage_group),
  legend.title = "Группа стадий",
  xlab = "Время (годы)",
  ylab = "Вероятность выживаемости",
  title = "KM: ранние vs поздние стадии (TCGA-LUAD)"
)
