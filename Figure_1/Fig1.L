############################################################
# CPTAC BRCA
# ENPP1 protein expression by CIN groups
# Median cutoff: global + subtype-specific
############################################################

rm(list = ls())

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggpubr)
})

############################################################
# 0) set working directory
############################################################

setwd("~/Desktop/brca_cptac_2020")

outdir <- "ENPP1_CIN_group_median_only"
dir.create(outdir, showWarnings = FALSE)

############################################################
# 1) helper: read cBioPortal clinical files
############################################################

read_cbio_clinical <- function(file) {
  x <- readLines(file)
  header_idx <- which(grepl("^(PATIENT_ID|SAMPLE_ID)\t", x))[1]
  if (is.na(header_idx)) stop(paste("Cannot find header row in", file))
  fread(file, skip = header_idx - 1, sep = "\t", header = TRUE)
}

standardize_subtype <- function(x) {
  y <- toupper(trimws(as.character(x)))
  y[y %in% c("", "NA", "N/A", "NULL")] <- NA
  
  out <- rep(NA_character_, length(y))
  out[grepl("LUMINAL A|LUMA", y)] <- "Luminal A"
  out[grepl("LUMINAL B|LUMB", y)] <- "Luminal B"
  out[grepl("HER2", y)] <- "HER2-enriched"
  out[grepl("BASAL|BASAL-LIKE|TNBC|TRIPLE", y)] <- "Basal/TNBC"
  out[grepl("NORMAL", y)] <- "Normal-like"
  out
}

############################################################
# 2) read data
############################################################

prot <- fread("data_protein_quantification.txt")
clin_pat <- read_cbio_clinical("data_clinical_patient.txt")
clin_sam <- read_cbio_clinical("data_clinical_sample.txt")

cat("===== patient columns =====\n")
print(colnames(clin_pat))

cat("===== sample columns =====\n")
print(colnames(clin_sam))

############################################################
# 3) extract ENPP1 protein
############################################################

enpp1 <- prot %>%
  filter(grepl("^ENPP1\\|", .data[[colnames(prot)[1]]])) %>%
  pivot_longer(
    cols = -1,
    names_to = "SampleID",
    values_to = "ENPP1_protein"
  ) %>%
  select(SampleID, ENPP1_protein)

cat("===== ENPP1 summary =====\n")
print(summary(enpp1$ENPP1_protein))

############################################################
# 4) build clinical tables
############################################################

clin_pat2 <- clin_pat %>%
  transmute(
    PatientID = PATIENT_ID,
    CIN = CHROMOSOME_INSTABILITY_INDEX_CIN_,
    TNBC_status = TNBC_UPDATED_CLINICAL_STATUS,
    ER_status = ER_UPDATED_CLINICAL_STATUS,
    PR_status = PR_CLINICAL_STATUS,
    HER2_status = ERBB2_UPDATED_CLINICAL_STATUS,
    Stemness = STEMNESS_SCORE,
    ESTIMATE_Immune = ESTIMATE_IMMUNE_SCORE,
    XCELL_Immune = XCELL_IMMUNE_SCORE,
    ESTIMATE_Stromal = ESTIMATE_STROMAL_SCORE,
    XCELL_Stromal = XCELL_STROMAL_SCORE
  )

clin_sam2 <- clin_sam %>%
  transmute(
    SampleID = SAMPLE_ID,
    PatientID = PATIENT_ID,
    PAM50_raw = PAM50,
    Subtype = standardize_subtype(PAM50)
  )

############################################################
# 5) merge
############################################################

dat <- enpp1 %>%
  left_join(clin_sam2, by = "SampleID") %>%
  left_join(clin_pat2, by = "PatientID")

cat("===== merged overview =====\n")
cat("n total:", nrow(dat), "\n")
print(table(dat$Subtype, useNA = "ifany"))
cat("n with CIN:", sum(!is.na(dat$CIN)), "\n")

fwrite(
  dat,
  file.path(outdir, "ENPP1_CPTAC_merged_table.tsv"),
  sep = "\t"
)

############################################################
# 6) keep CIN + ENPP1 analyzable samples
############################################################

df_cin <- dat %>%
  filter(!is.na(CIN), !is.na(ENPP1_protein))

cat("===== analyzable for CIN analysis =====\n")
cat("n =", nrow(df_cin), "\n")

############################################################
# 7) GLOBAL median CIN cutoff
############################################################

cin_median_cut <- median(df_cin$CIN, na.rm = TRUE)

df_cin <- df_cin %>%
  mutate(
    CIN_median_group = ifelse(
      CIN >= cin_median_cut,
      "CIN_high",
      "CIN_low"
    )
  )

df_cin$CIN_median_group <- factor(
  df_cin$CIN_median_group,
  levels = c("CIN_low", "CIN_high")
)

cat("===== global median cutoff =====\n")
cat("CIN median cutoff =", round(cin_median_cut, 4), "\n")
print(table(df_cin$CIN_median_group))

global_median_stats <- df_cin %>%
  group_by(CIN_median_group) %>%
  summarise(
    n = n(),
    median_CIN = median(CIN, na.rm = TRUE),
    mean_CIN = mean(CIN, na.rm = TRUE),
    median_ENPP1 = median(ENPP1_protein, na.rm = TRUE),
    mean_ENPP1 = mean(ENPP1_protein, na.rm = TRUE),
    .groups = "drop"
  )

fwrite(
  global_median_stats,
  file.path(outdir, "Global_CIN_median_group_summary.tsv"),
  sep = "\t"
)

global_wilcox <- wilcox.test(
  ENPP1_protein ~ CIN_median_group,
  data = df_cin
)

global_pval <- global_wilcox$p.value

global_test_res <- data.frame(
  comparison = "CIN_high vs CIN_low",
  cutoff_method = "Global median",
  CIN_median_cutoff = cin_median_cut,
  p_value = global_pval,
  n_low = sum(df_cin$CIN_median_group == "CIN_low"),
  n_high = sum(df_cin$CIN_median_group == "CIN_high")
)

fwrite(
  global_test_res,
  file.path(outdir, "Global_CIN_median_Wilcoxon_test.tsv"),
  sep = "\t"
)

p_global_median <- ggplot(
  df_cin,
  aes(
    x = CIN_median_group,
    y = ENPP1_protein,
    fill = CIN_median_group
  )
) +
  geom_violin(trim = FALSE, alpha = 0.5, color = NA) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.10, size = 1.5, alpha = 0.7) +
  stat_compare_means(method = "wilcox.test") +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "none"
  ) +
  labs(
    title = "ENPP1 expression by global CIN median group",
    subtitle = paste0(
      "Median CIN cutoff = ",
      round(cin_median_cut, 4)
    ),
    x = "",
    y = "Normalized ENPP1 protein expression"
  )

ggsave(
  file.path(outdir, "ENPP1_by_CIN_global_median_group.pdf"),
  p_global_median,
  width = 5,
  height = 5
)

############################################################
# 8) subtype-specific median CIN cutoff
############################################################

subtypes_to_run <- c(
  "Luminal A",
  "Luminal B",
  "HER2-enriched",
  "Basal/TNBC",
  "Normal-like"
)

sub_median_summary <- list()

for (st in subtypes_to_run) {
  
  df_sub <- df_cin %>%
    filter(Subtype == st)
  
  cat("===== subtype median:", st, "=====\n")
  cat("n =", nrow(df_sub), "\n")
  
  if (nrow(df_sub) < 10) {
    cat("Skip: too few samples\n")
    next
  }
  
  sub_median_cut <- median(df_sub$CIN, na.rm = TRUE)
  
  df_sub <- df_sub %>%
    mutate(
      CIN_median_group_sub = ifelse(
        CIN >= sub_median_cut,
        "CIN_high",
        "CIN_low"
      )
    )
  
  df_sub$CIN_median_group_sub <- factor(
    df_sub$CIN_median_group_sub,
    levels = c("CIN_low", "CIN_high")
  )
  
  if (length(unique(df_sub$CIN_median_group_sub)) < 2) {
    cat("Skip: only one CIN group\n")
    next
  }
  
  wt <- wilcox.test(
    ENPP1_protein ~ CIN_median_group_sub,
    data = df_sub
  )
  
  sub_median_summary[[st]] <- data.frame(
    Subtype = st,
    N = nrow(df_sub),
    Median_CIN_cutoff = sub_median_cut,
    P_value = wt$p.value,
    n_low = sum(df_sub$CIN_median_group_sub == "CIN_low"),
    n_high = sum(df_sub$CIN_median_group_sub == "CIN_high"),
    ENPP1_lowCIN_median = median(
      df_sub$ENPP1_protein[df_sub$CIN_median_group_sub == "CIN_low"],
      na.rm = TRUE
    ),
    ENPP1_highCIN_median = median(
      df_sub$ENPP1_protein[df_sub$CIN_median_group_sub == "CIN_high"],
      na.rm = TRUE
    ),
    ENPP1_lowCIN_mean = mean(
      df_sub$ENPP1_protein[df_sub$CIN_median_group_sub == "CIN_low"],
      na.rm = TRUE
    ),
    ENPP1_highCIN_mean = mean(
      df_sub$ENPP1_protein[df_sub$CIN_median_group_sub == "CIN_high"],
      na.rm = TRUE
    ),
    stringsAsFactors = FALSE
  )
  
  p_sub_median <- ggplot(
    df_sub,
    aes(
      x = CIN_median_group_sub,
      y = ENPP1_protein,
      fill = CIN_median_group_sub
    )
  ) +
    geom_violin(trim = FALSE, alpha = 0.5, color = NA) +
    geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.8) +
    geom_jitter(width = 0.10, size = 1.5, alpha = 0.7) +
    stat_compare_means(method = "wilcox.test") +
    theme_classic(base_size = 13) +
    theme(
      legend.position = "none"
    ) +
    labs(
      title = paste0(st, ": ENPP1 by subtype-specific CIN median"),
      subtitle = paste0(
        "Median CIN cutoff = ",
        round(sub_median_cut, 4)
      ),
      x = "",
      y = "Normalized ENPP1 protein expression"
    )
  
  ggsave(
    file.path(
      outdir,
      paste0(
        "ENPP1_by_CIN_median_",
        gsub("[/ ]", "_", st),
        ".pdf"
      )
    ),
    p_sub_median,
    width = 5,
    height = 5
  )
}

############################################################
# 9) save subtype summary
############################################################

if (length(sub_median_summary) > 0) {
  
  sub_median_summary_df <- bind_rows(sub_median_summary)
  
  fwrite(
    sub_median_summary_df,
    file.path(outdir, "Subtype_specific_CIN_median_summary.tsv"),
    sep = "\t"
  )
  
  print(sub_median_summary_df)
}

############################################################
# 10) save overall cutoff summary
############################################################

cutoff_summary <- data.frame(
  Global_median_CIN_cutoff = cin_median_cut
)

fwrite(
  cutoff_summary,
  file.path(outdir, "Global_CIN_median_cutoff_summary.tsv"),
  sep = "\t"
)

############################################################
# 11) done
############################################################

cat("===== DONE =====\n")
cat("Results saved to:\n")
cat(outdir, "\n")
print(list.files(outdir))
