############################################################
# CPTAC BRCA: ENPP1 protein by subtype and CIN
############################################################

rm(list = ls())

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggpubr)
})

setwd("~/Desktop/brca_cptac_2020")

outdir <- "ENPP1_CPTAC_subtype_CIN"
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

############################################################
# 6) ENPP1 by PAM50 subtype
############################################################
p_subtype <- ggplot(dat %>% filter(!is.na(Subtype)),
                    aes(x = Subtype, y = ENPP1_protein, fill = Subtype)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.08, size = 1.2, alpha = 0.6) +
  theme_classic(base_size = 13) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(
    title = "ENPP1 protein across PAM50 subtypes",
    x = "",
    y = "Normalized ENPP1 protein expression"
  )

ggsave(file.path(outdir, "ENPP1_by_PAM50_violin.pdf"), p_subtype, width = 7, height = 5)

############################################################
# 7) Basal/TNBC vs non-Basal
############################################################
dat2 <- dat %>%
  mutate(
    Basal_group = ifelse(Subtype == "Basal/TNBC", "Basal/TNBC", "Non-Basal")
  )

p_basal <- ggplot(dat2 %>% filter(!is.na(Basal_group)),
                  aes(x = Basal_group, y = ENPP1_protein, fill = Basal_group)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.10, size = 1.5, alpha = 0.7) +
  stat_compare_means(method = "wilcox.test") +
  theme_classic(base_size = 13) +
  labs(
    title = "ENPP1 protein: Basal/TNBC vs Non-Basal",
    x = "",
    y = "Normalized ENPP1 protein expression"
  )

ggsave(file.path(outdir, "ENPP1_Basal_vs_nonBasal.pdf"), p_basal, width = 5, height = 5)

############################################################
# 8) ENPP1 vs CIN
############################################################
df_cin <- dat %>%
  filter(!is.na(ENPP1_protein), !is.na(CIN))

cor_test <- cor.test(df_cin$ENPP1_protein, df_cin$CIN, method = "spearman")

capture.output(cor_test, file = file.path(outdir, "ENPP1_vs_CIN_spearman.txt"))

p_cin <- ggplot(df_cin, aes(x = ENPP1_protein, y = CIN)) +
  geom_point(alpha = 0.75, size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic(base_size = 13) +
  labs(
    title = paste0(
      "ENPP1 protein vs CIN (Spearman rho = ",
      round(cor_test$estimate, 3),
      ", p = ",
      signif(cor_test$p.value, 3),
      ")"
    ),
    x = "Normalized ENPP1 protein expression",
    y = "Chromosome instability index (CIN)"
  )

ggsave(file.path(outdir, "ENPP1_vs_CIN_scatter.pdf"), p_cin, width = 6, height = 5)

############################################################
# 9) ENPP1 high/low by median vs CIN
############################################################
median_cut <- median(df_cin$ENPP1_protein, na.rm = TRUE)

df_cin2 <- df_cin %>%
  mutate(
    ENPP1_group = ifelse(ENPP1_protein >= median_cut, "ENPP1_high", "ENPP1_low")
  )

writeLines(
  paste("Median cutoff for ENPP1 =", round(median_cut, 4)),
  file.path(outdir, "ENPP1_median_cutoff.txt")
)

p_cin_group <- ggplot(df_cin2, aes(x = ENPP1_group, y = CIN, fill = ENPP1_group)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.10, size = 1.5, alpha = 0.7) +
  stat_compare_means(method = "wilcox.test") +
  theme_classic(base_size = 13) +
  labs(
    title = "CIN by ENPP1 median group",
    x = "",
    y = "Chromosome instability index (CIN)"
  )

ggsave(file.path(outdir, "CIN_by_ENPP1_median_group.pdf"), p_cin_group, width = 5, height = 5)

############################################################
# 10) subtype-specific ENPP1 vs CIN
############################################################
subtypes_to_run <- c("Luminal A", "Luminal B", "HER2-enriched", "Basal/TNBC", "Normal-like")

sub_cor_summary <- list()

for (st in subtypes_to_run) {
  df_sub <- df_cin %>% filter(Subtype == st)
  
  if (nrow(df_sub) < 8) next
  
  ct <- cor.test(df_sub$ENPP1_protein, df_sub$CIN, method = "spearman")
  
  sub_cor_summary[[st]] <- data.frame(
    Subtype = st,
    N = nrow(df_sub),
    Spearman_rho = unname(ct$estimate),
    P_value = ct$p.value
  )
  
  p_sub <- ggplot(df_sub, aes(x = ENPP1_protein, y = CIN)) +
    geom_point(alpha = 0.75, size = 2) +
    geom_smooth(method = "lm", se = FALSE) +
    theme_classic(base_size = 13) +
    labs(
      title = paste0(
        st, ": ENPP1 vs CIN (rho = ",
        round(ct$estimate, 3),
        ", p = ",
        signif(ct$p.value, 3),
        ")"
      ),
      x = "Normalized ENPP1 protein expression",
      y = "Chromosome instability index (CIN)"
    )
  
  ggsave(
    file.path(outdir, paste0("ENPP1_vs_CIN_", gsub("[/ ]", "_", st), ".pdf")),
    p_sub, width = 6, height = 5
  )
}

if (length(sub_cor_summary) > 0) {
  sub_cor_summary_df <- bind_rows(sub_cor_summary)
  fwrite(sub_cor_summary_df,
         file.path(outdir, "Subtype_specific_ENPP1_vs_CIN.tsv"),
         sep = "\t")
}

############################################################
# 11) additional: TNBC clinical status
############################################################
p_tnbc <- ggplot(dat %>% filter(!is.na(TNBC_status)),
                 aes(x = TNBC_status, y = ENPP1_protein, fill = TNBC_status)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.10, size = 1.5, alpha = 0.7) +
  theme_classic(base_size = 13) +
  labs(
    title = "ENPP1 protein by TNBC clinical status",
    x = "",
    y = "Normalized ENPP1 protein expression"
  )

ggsave(file.path(outdir, "ENPP1_by_TNBC_status.pdf"), p_tnbc, width = 6, height = 5)

############################################################
# 12) save merged table
############################################################
fwrite(dat, file.path(outdir, "ENPP1_CPTAC_merged_table.tsv"), sep = "\t")

cat("===== analysis done =====\n")
print(list.files(outdir))
