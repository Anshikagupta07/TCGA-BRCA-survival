# ============================================================================
# Breast Cancer Survival Analysis — TCGA BRCA
# Author : Anshika Gupta
# Methods: Kaplan-Meier survival analysis + multivariate Cox regression
# Data   : TCGA-BRCA clinical cohort (public, via TCGAbiolinks)
# ============================================================================

# libraries
library(TCGAbiolinks)
library(dplyr)
library(tidyr)
library(tidyverse)
library(broom)
library(ggplot2)
library(survival)
library(survminer)
dir.create("figures", showWarnings = FALSE)

# ── Helper: ggsurvplot requires png() + print() — ggsave() only saves risk table
save_km <- function(plot, filename, width = 10, height = 8) {
  png(filename, width = width, height = height, units = "in", res = 150)
  print(plot)
  dev.off()
}


# 1. TCGA-BRCA clinical data
query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Clinical",
  data.type = "Clinical Supplement",
  data.format = "BCR Biotab"
)

GDCdownload(query)
clinical_data <- GDCprepare(query)
patient_data <- clinical_data$clinical_patient_brca

# 2. Data cleaning & Wrangling
df <- patient_data %>%
    
    filter(!bcr_patient_barcode %in% c("bcr_patient_barcode", "CDE_ID:")) %>%
    
    select(
      patient_id        = bcr_patient_barcode,
      age               = age_at_diagnosis,
      vital_status,
      days_last_contact = last_contact_days_to,
      days_to_death     = death_days_to,
      stage             = ajcc_pathologic_tumor_stage,
      er_status         = er_status_by_ihc,
      pr_status         = pr_status_by_ihc,
      her2_status       = her2_status_by_ihc,
      histological_type
    ) %>%
    
    mutate(across(c(age, days_last_contact, days_to_death),
                  ~ suppressWarnings(as.numeric(.)))) %>%
    
    mutate(
      OS_event    = ifelse(vital_status == "Dead", 1, 0),
      OS_time     = ifelse(OS_event == 1, days_to_death, days_last_contact),
      
      stage_clean = case_when(
        grepl("stage\\s*i(a|b)?$",   stage, ignore.case = TRUE) ~ "Stage I",
        grepl("stage\\s*ii(a|b|c)?$", stage, ignore.case = TRUE) ~ "Stage II",
        grepl("stage\\s*iii",         stage, ignore.case = TRUE) ~ "Stage III",
        grepl("stage\\s*iv",          stage, ignore.case = TRUE) ~ "Stage IV",
        TRUE                                                     ~ NA_character_
      )
    ) %>%
    
    mutate(across(c(er_status, pr_status, her2_status),
                  ~ ifelse(. %in% c("Positive", "Negative"), ., NA_character_))) %>%
    
    filter(!is.na(OS_time), OS_time > 0)
  
# Sanity check
cat("Dataset:", nrow(df), "patients |", sum(df$OS_event), "events\n")
table(df$stage_clean)
summary(df$OS_time)

# 3. EDA Plots
# PLOT 1: Age at diagnosis by stage
df %>%
  filter(!is.na(stage_clean)) %>%
  ggplot(aes(x = stage_clean, y = age, fill = stage_clean)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21) +
  scale_fill_manual(values = c("#B5D4F4", "#378ADD", "#185FA5", "#042C53")) +
  labs(title = "Age at diagnosis by cancer stage", x = "Stage", y = "Age (years)") +
  theme_minimal() +
  theme(legend.position = "none")

ggsave("figures/01_age_by_stage.png", width = 7, height = 5, dpi = 150)

# PLOT 2: Receptor status counts
df %>%
  select(er_status, pr_status, her2_status) %>%
  pivot_longer(everything(), names_to = "receptor", values_to = "status") %>%
  filter(!is.na(status)) %>%
  mutate(receptor = recode(receptor,
                           er_status = "ER", pr_status = "PR", her2_status = "HER2")) %>%
  ggplot(aes(x = receptor, fill = status)) +
  geom_bar(position = "fill", alpha = 0.85) +
  scale_fill_manual(values = c("Positive" = "#378ADD", "Negative" = "#B4B2A9")) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Receptor status distribution",
       x = "Receptor", y = "Proportion", fill = "Status") +
  theme_minimal()

ggsave("figures/02_receptor_status.png", width = 6, height = 5, dpi = 150)

# PLOT 3: Survival time distribution 
ggplot(df, aes(x = OS_time / 365.25)) +
  geom_histogram(bins = 40, fill = "#378ADD", alpha = 0.8, color = "white") +
  geom_vline(xintercept = median(df$OS_time, na.rm = TRUE) / 365.25,
             linetype = "dashed", color = "#A32D2D", linewidth = 0.8) +
  labs(title = "Distribution of follow-up time",
       subtitle = "Dashed line = median follow-up",
       x = "Years", y = "Number of patients") +
  theme_minimal()

ggsave("figures/03_followup_distribution.png", width = 7, height = 5, dpi = 150)

# 4. KM Curves
# PLOT 1: By cancer stage
df_stage <- df %>%
  filter(!is.na(stage_clean)) %>%
  mutate(stage_clean = factor(stage_clean,
                              levels = c("Stage I", "Stage II",
                                         "Stage III", "Stage IV")))

fit_stage <- survfit(Surv(OS_time, OS_event) ~ stage_clean, data = df_stage)

p_stage <- ggsurvplot(
  fit_stage,
  data               = df_stage,
  pval               = TRUE,
  pval.size          = 4,
  conf.int           = TRUE,
  conf.int.alpha     = 0.1,
  risk.table         = TRUE,
  risk.table.height  = 0.28,
  risk.table.title   = "Patients at risk",
  palette            = c("#B5D4F4", "#378ADD", "#185FA5", "#042C53"),
  legend.title       = "",
  legend.labs        = c("Stage I", "Stage II", "Stage III", "Stage IV"),
  title              = "Overall survival by cancer stage — TCGA BRCA (n=1,012)",
  xlab               = "Time (days)",
  ylab               = "Survival probability",
  xlim               = c(0, 7000),
  break.time.by      = 1000,
  surv.median.line   = "hv",        # draws median survival lines
  ggtheme            = theme_minimal(base_size = 13)
)

save_km(p_stage, "figures/04_km_by_stage.png")

# PLOT 2: By ER status 
df_er <- df %>%
  filter(!is.na(er_status)) %>%
  mutate(er_status = factor(er_status, levels = c("Negative", "Positive")))

fit_er <- survfit(Surv(OS_time, OS_event) ~ er_status, data = df_er)

p_er <- ggsurvplot(
  fit_er,
  data               = df_er,
  pval               = TRUE,
  pval.size          = 4,
  conf.int           = TRUE,
  conf.int.alpha     = 0.1,
  risk.table         = TRUE,
  risk.table.height  = 0.25,
  risk.table.title   = "Patients at risk",
  palette            = c("#E24B4A", "#042C53"),  # red = negative, navy = positive
  legend.title       = "",
  legend.labs        = c("Negative", "Positive"),
  title              = "Overall survival by ER status — TCGA BRCA",
  xlab               = "Time (days)",
  ylab               = "Survival probability",
  xlim               = c(0, 7000),
  break.time.by      = 1000,
  surv.median.line   = "hv",
  ggtheme            = theme_minimal(base_size = 13)
)

save_km(p_er, "figures/05_km_by_er_status.png")


# PLOT 3: By HER2 status 
df_her2 <- df %>%
  filter(!is.na(her2_status)) %>%
  mutate(her2_status = factor(her2_status, levels = c("Negative", "Positive")))

fit_her2 <- survfit(Surv(OS_time, OS_event) ~ her2_status, data = df_her2)

p_her2 <- ggsurvplot(
  fit_her2,
  data               = df_her2,
  pval               = TRUE,
  pval.size          = 4,
  conf.int           = TRUE,
  conf.int.alpha     = 0.1,
  risk.table         = TRUE,
  risk.table.height  = 0.25,
  risk.table.title   = "Patients at risk",
  palette            = c("#3B6D11", "#EF9F27"),  # green = negative, amber = positive
  legend.title       = "HER2 status",
  legend.labs        = c("Negative", "Positive"),
  title              = "Overall survival by HER2 status — TCGA BRCA",
  xlab               = "Time (days)",
  ylab               = "Survival probability",
  xlim               = c(0, 4000),
  break.time.by      = 1000,
  surv.median.line   = "hv",
  ggtheme            = theme_minimal(base_size = 13)
)

save_km(p_her2, "figures/06_km_by_her2_status.png")


# Quick check: print p-values to console
cat("\n── Stage log-rank p-value:\n")
print(survdiff(Surv(OS_time, OS_event) ~ stage_clean, data = df_stage))

cat("\n── ER status log-rank p-value:\n")
print(survdiff(Surv(OS_time, OS_event) ~ er_status, data = df_er))

cat("\n── HER2 status log-rank p-value:\n")
print(survdiff(Surv(OS_time, OS_event) ~ her2_status, data = df_her2))

# Cox model dataset
df_cox <- df %>%
  filter(!is.na(stage_clean), !is.na(er_status), !is.na(her2_status)) %>%
  mutate(
    stage_clean = factor(stage_clean,
                         levels = c("Stage I", "Stage II", "Stage III", "Stage IV")),
    er_status   = factor(er_status,   levels = c("Negative", "Positive")),
    her2_status = factor(her2_status, levels = c("Negative", "Positive")),
    age_scaled  = scale(age)[, 1]
  )

# Fit model 
cox_model <- coxph(
  Surv(OS_time, OS_event) ~ age_scaled + stage_clean + er_status + her2_status,
  data = df_cox
)
summary(cox_model)

# Extract results directly from model 
conc <- round(summary(cox_model)$concordance["C"], 3)

hr_data <- tidy(cox_model, exponentiate = TRUE, conf.int = TRUE) %>%
  mutate(
    variable = recode(term,
                      "age_scaled"           = "Age (per SD increase)",
                      "stage_cleanStage II"  = "Stage II vs Stage I",
                      "stage_cleanStage III" = "Stage III vs Stage I",
                      "stage_cleanStage IV"  = "Stage IV vs Stage I",
                      "er_statusPositive"    = "ER positive vs negative",
                      "her2_statusPositive"  = "HER2 positive vs negative"
    ),
    pval_label = ifelse(p.value < 0.001, "p<0.001",
                        paste0("p=", round(p.value, 3))),
    group    = case_when(
      grepl("Age",   variable) ~ "Clinical",
      grepl("Stage", variable) ~ "Stage",
      TRUE                     ~ "Receptor"
    ),
    sig      = ifelse(p.value < 0.05, "significant", "trend"),
    variable = factor(variable, levels = rev(variable))
  )

# Forest plot
library(grid)

p <- ggplot(hr_data, aes(estimate, variable, color = sig)) +
  
  geom_vline(xintercept = 1, linetype = "dashed", color = "#7a7a7a", linewidth = 0.6) +
  
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high),
                width = 0.2, orientation = "y", linewidth = 0.8) +
  
  geom_point(size = 3, shape = 18) +
  
  geom_text(aes(label = sprintf("%.2f (%.2f–%.2f)", estimate, conf.low, conf.high)),
            hjust = -0.05, size = 3, color = "black") +
  
  geom_text(aes(x = 4.4, label = pval_label),
            hjust = 1, size = 3, color = "black") +
  
  scale_x_log10(
    breaks = c(0.25, 0.5, 1, 2, 4),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  
  scale_color_manual(values = c("significant" = "#2C7BB6", "trend" = "#9E9E9E")) +
  
  facet_grid(group ~ ., scales = "free_y", space = "free_y", switch = "y") +
  
  labs(
    title = "Survival Analysis (Cox Model)",
    subtitle = "Hazard ratios with 95% confidence intervals",
    x = "Hazard Ratio (log scale)",
    y = NULL
  ) +
  
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    
    # ✅ compact facets
    panel.spacing.y = unit(0.2, "lines"),
    
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, face = "bold", size = 9),
    
    axis.text.y = element_text(size = 9),
    
    axis.text.x = element_text(margin = margin(t = -8)),
    axis.title.x = element_text(margin = margin(t = -4)),
    axis.ticks.length = unit(-0.15, "cm"),
    
    plot.margin = margin(t = 8, r = 15, b = 5, l = 8)
  ) +
  
  coord_cartesian(xlim = c(0.25, 4.5), clip = "off")
ggsave("forest_plot.png", p, width = 8, height = 5, dpi = 300)
