#install.packages("remotes")
#devtools::install_github("MRCIEU/TwoSampleMR")
#devtools::install_github("MRCIEU/MRInstruments")
#devtools::install_github("NightingaleHealth/ggforestplot")
#library(remotes)
#install_github("WSpiller/MVMR", build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE)
# Load Required Libraries

library(tidyverse)
library(TwoSampleMR)
library(glue)
library(data.table)
library(dplyr)

# Load exposure data: protein GWAS
exposure_data <- fread("pQTL_PCSK9_GWAS.txt")
dim(exposure_data)
exposure_data$Phenotype <- "PCSK9"
exposure_data <- exposure_data[!duplicated(exposure_data$SNP), ]
dim(exposure_data)
exposure_data_1 <- exposure_data %>% dplyr::filter(P < 5e-8)
dim(exposure_data_1)

# Save filtered exposure data
fwrite(exposure_data_1, "pQTL_PCSK9_GWAS_1.txt", sep = "\t", quote = FALSE)

# Load using TwoSampleMR format
exposure_data_1 <- read_exposure_data(
  filename = "pQTL_PCSK9_GWAS_1.txt",
  sep = "\t",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  eaf_col = "MAF",
  pval_col = "P",
  phenotype_col = "Phenotype"
)
exposure_data_1[1:3, ]

# Read and Format Outcome Data (Triglycerides GWAS)
outcome_data <- fread("Triglycerides_GWAS.txt")
dim(outcome_data)
outcome_data$Phenotype <- "Triglycerides"
outcome_data <- outcome_data[!duplicated(outcome_data$SNP), ]
dim(outcome_data)
outcome_data_1 <- outcome_data %>% dplyr::filter(P < 0.05)
dim(outcome_data_1)

# Ensure Z and SE are numeric
outcome_data_1$Z <- as.numeric(outcome_data_1$Z)
outcome_data_1$SE <- as.numeric(outcome_data_1$SE)

# Compute beta
outcome_data_1$beta <- outcome_data_1$Z * outcome_data_1$SE

# Save cleaned outcome data
fwrite(outcome_data_1, "Triglycerides_GWAS_processed_with_beta.txt", sep = "\t", quote = FALSE)

# Load using TwoSampleMR format
outcome_data_1 <- read_outcome_data(
  filename = "Triglycerides_GWAS_processed_with_beta.txt",
  sep = "\t",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "SE",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  eaf_col = "MAF",
  pval_col = "P",
  phenotype_col = "Phenotype"
)
outcome_data_1[1:3, ]

# Harmonise Exposure and Outcome Data
harmonised_data <- harmonise_data(
  exposure_dat = exposure_data_1,
  outcome_dat = outcome_data_1
)
dim(harmonised_data)
harmonised_data[1:3, ]

# Data frame with the three main methods only
mr_results <- data.frame(
  method = c("MR Egger", "Weighted median", "Inverse variance weighted"),
  b = c(-0.10166960, -0.04855068, 0.55805900),
  se = c(6.964044e-02, 5.808784e-03, 2.868964e-02)
)

# Calculate confidence intervals
mr_results$lower_ci <- mr_results$b - 1.96 * mr_results$se
mr_results$upper_ci <- mr_results$b + 1.96 * mr_results$se

# Plot
ggplot(mr_results, aes(x = method, y = b)) +
  geom_point(size = 3, color = "blue") +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.2, color = "black") +
  theme_minimal() +
  labs(
    title = "MR Estimates with 95% Confidence Intervals",
    x = "MR Method",
    y = "Beta (Effect Size)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))