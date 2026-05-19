############################################################
## SCAN-B Basal only
## ENPP1 expression vs BRCA altered / functional loss
############################################################

rm(list = ls())

############################################################
## 0. Packages
############################################################

need_pkgs <- c(
  "data.table",
  "readxl",
  "stringr",
  "GEOquery",
  "ggplot2"
)

for (p in need_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
}

library(data.table)
library(readxl)
library(stringr)
library(GEOquery)
library(ggplot2)

options(stringsAsFactors = FALSE)

############################################################
## 1. Choose files
############################################################

cat("Choose SupplementaryDataTable.xlsx ...\n")
cohort_xlsx <- file.choose()

cat("Choose WGS-RNA mapping txt/tsv file ...\n")
map_file <- file.choose()

cat("Choose GSE96058 gene_expression_transformed.csv ...\n")
expr_file <- file.choose()

cat("cohort_xlsx =", cohort_xlsx, "\n")
cat("map_file    =", map_file, "\n")
cat("expr_file   =", expr_file, "\n")

############################################################
## 2. Read PatientDataTable
############################################################

cs <- as.data.table(read_excel(cohort_xlsx, sheet = "PatientDataTable"))
cat("PatientDataTable dims:", dim(cs), "\n")

############################################################
## 3. Find PD_ID-like column
############################################################

find_pd_col <- function(dt) {
  for (nm in names(dt)) {
    v <- as.character(dt[[nm]])
    v <- v[!is.na(v)]
    if (length(v) == 0) next
    if (mean(grepl("^PD\\d+", v)) > 0.3) return(nm)
  }
  return(NA_character_)
}

pd_col_cs <- find_pd_col(cs)
cat("Detected PD_ID column:", pd_col_cs, "\n")
if (is.na(pd_col_cs)) stop("Cannot find PD_ID-like column in PatientDataTable")

############################################################
## 4. Check required columns
############################################################

required_cols <- c("PAM50_AIMS")
missing_cols <- setdiff(required_cols, names(cs))
if (length(missing_cols) > 0) {
  stop(paste("Missing columns in PatientDataTable:", paste(missing_cols, collapse = ", ")))
}

############################################################
## 5. Build clinical subtype table
############################################################

subtype_tab <- cs[, .(
  PD_ID = as.character(get(pd_col_cs)),
  PAM50_AIMS = as.character(PAM50_AIMS)
)]

subtype_tab <- unique(subtype_tab[!is.na(PD_ID) & PD_ID != ""])
cat("subtype_tab dims:", dim(subtype_tab), "\n")

############################################################
## 6. Read WGS-RNA mapping
############################################################

map <- fread(map_file)
cat("mapping raw dims:", dim(map), "\n")

if (!all(c("PD_ID", "External_ID_rba") %in% names(map))) {
  stop("Mapping file must contain columns: PD_ID and External_ID_rba")
}

mapping_table <- map[, .(
  PD_ID = as.character(PD_ID),
  RNAseq_external_id = as.character(External_ID_rba)
)]

mapping_table <- unique(mapping_table[
  !is.na(PD_ID) & PD_ID != "" &
    !is.na(RNAseq_external_id) & RNAseq_external_id != ""
])

cat("mapping_table dims:", dim(mapping_table), "\n")

############################################################
## 7. Merge subtype + mapping
############################################################

step1_merged <- merge(subtype_tab, mapping_table, by = "PD_ID")
cat("After subtype + mapping merge dims:", dim(step1_merged), "\n")
cat("Unique PD_ID after merge:", length(unique(step1_merged$PD_ID)), "\n")

############################################################
## 8. Read expression matrix and extract ENPP1
############################################################

expr_dt <- fread(expr_file)
cat("Expression matrix dims:", dim(expr_dt), "\n")

find_gene_col <- function(dt) {
  for (nm in names(dt)) {
    v <- as.character(dt[[nm]])
    v <- v[!is.na(v)]
    if (length(v) == 0) next
    if ("ENPP1" %in% v) return(nm)
  }
  return(NA_character_)
}

gene_col <- find_gene_col(expr_dt)
cat("Detected gene column:", gene_col, "\n")
if (is.na(gene_col)) stop("Cannot find a column containing ENPP1 in expression matrix")

enpp1_row <- expr_dt[get(gene_col) == "ENPP1"]
if (nrow(enpp1_row) == 0) stop("ENPP1 row not found in expression matrix")
if (nrow(enpp1_row) > 1) {
  cat("Multiple ENPP1 rows found; keeping first row.\n")
  enpp1_row <- enpp1_row[1]
}

sample_cols <- setdiff(names(enpp1_row), gene_col)

enpp1_long <- melt(
  enpp1_row,
  id.vars = gene_col,
  measure.vars = sample_cols,
  variable.name = "F_ID",
  value.name = "ENPP1"
)

enpp1_long[, ENPP1 := suppressWarnings(as.numeric(ENPP1))]
enpp1_long <- enpp1_long[!is.na(ENPP1)]
enpp1_long <- enpp1_long[grepl("^F\\d+", F_ID)]

# Collapse rep1/rep2 if present
enpp1_long[, F_core := str_replace(as.character(F_ID), "rep\\d+$", "")]
enpp1_F <- enpp1_long[, .(ENPP1 = mean(ENPP1, na.rm = TRUE)), by = F_core]
setnames(enpp1_F, "F_core", "F_ID")

cat("Collapsed ENPP1 sample count:", nrow(enpp1_F), "\n")

############################################################
## 9. GEO metadata: F_ID <-> RNAseq_external_id
############################################################

cat("Downloading GEO metadata for GSE96058 ...\n")
gse <- getGEO("GSE96058", GSEMatrix = TRUE)

pd_all <- rbindlist(lapply(gse, function(x) as.data.table(pData(x))), fill = TRUE)
cat("pd_all dims:", dim(pd_all), "\n")

allcols <- names(pd_all)
pd_all[, blob := do.call(paste, c(.SD, sep = " | ")), .SDcols = allcols]

pd_all[, F_ID := str_extract(blob, "F\\d+(?:rep\\d+)?")]
pd_all[, RNAseq_external_id := str_extract(blob, "S\\d{6}[^\\s\\|]*")]

map2 <- unique(pd_all[
  !is.na(F_ID) & !is.na(RNAseq_external_id),
  .(F_ID, RNAseq_external_id)
])

cat("map2 dims:", dim(map2), "\n")

############################################################
## 10. Final expression table
############################################################

tmp <- merge(step1_merged, map2, by = "RNAseq_external_id")
final_df <- merge(tmp, enpp1_F, by = "F_ID")
final_df <- unique(final_df, by = "PD_ID")

cat("final_df dims:", dim(final_df), "\n")

############################################################
## 11. Keep Basal only
############################################################

dat_basal <- final_df[PAM50_AIMS == "Basal"]
cat("Basal sample size:", nrow(dat_basal), "\n")

############################################################
## 12. Read BRCA sheet and merge
############################################################

brca_status <- as.data.table(
  read_excel(cohort_xlsx, sheet = "BRCA_PALB2_RAD51C_status")
)

cat("BRCA sheet dims:", dim(brca_status), "\n")

if (!"Assay" %in% names(brca_status)) {
  stop("BRCA_PALB2_RAD51C_status sheet must contain column 'Assay'")
}

setnames(brca_status, "Assay", "PD_ID")

brca_cols_needed <- c(
  "isBRCA1Biallelic",
  "isBRCA1Monoallelic",
  "isBRCA2Biallelic",
  "isBRCA2Monoallelic",
  "BRCA1promoterMet"
)

missing_brca_cols <- setdiff(brca_cols_needed, names(brca_status))
if (length(missing_brca_cols) > 0) {
  stop(paste("Missing BRCA columns:", paste(missing_brca_cols, collapse = ", ")))
}

dat_basal2 <- merge(dat_basal, brca_status, by = "PD_ID", all.x = TRUE)
cat("dat_basal2 dims:", dim(dat_basal2), "\n")

############################################################
## 13. Force BRCA columns numeric
############################################################

for (cn in brca_cols_needed) {
  dat_basal2[, (cn) := suppressWarnings(as.numeric(get(cn)))]
}

############################################################
## 14. Define altered groups
############################################################

dat_basal2[, BRCA1_altered := fifelse(
  isBRCA1Biallelic == 1 |
    isBRCA1Monoallelic == 1 |
    BRCA1promoterMet == 1,
  "BRCA1 altered", "BRCA1 WT"
)]

dat_basal2[, BRCA2_altered := fifelse(
  isBRCA2Biallelic == 1 |
    isBRCA2Monoallelic == 1,
  "BRCA2 altered", "BRCA2 WT"
)]

dat_basal2[, BRCA1or2_altered := fifelse(
  isBRCA1Biallelic == 1 |
    isBRCA1Monoallelic == 1 |
    BRCA1promoterMet == 1 |
    isBRCA2Biallelic == 1 |
    isBRCA2Monoallelic == 1,
  "BRCA1/2 altered", "BRCA1/2 WT"
)]

dat_basal2[, BRCA1_altered := factor(BRCA1_altered, levels = c("BRCA1 WT", "BRCA1 altered"))]
dat_basal2[, BRCA2_altered := factor(BRCA2_altered, levels = c("BRCA2 WT", "BRCA2 altered"))]
dat_basal2[, BRCA1or2_altered := factor(BRCA1or2_altered, levels = c("BRCA1/2 WT", "BRCA1/2 altered"))]

############################################################
## 15. Define functional groups (2-group version)
## Functional = WT + Monoallelic
## Nonfunctional = Biallelic loss / BRCA1 hypermethylation
############################################################

# Combined BRCA1/2 functional status
dat_basal2[, BRCA_function_group := fifelse(
  isBRCA1Biallelic == 1 |
    isBRCA2Biallelic == 1 |
    BRCA1promoterMet == 1,
  "Nonfunctional",
  "Functional"
)]
dat_basal2[, BRCA_function_group := factor(
  BRCA_function_group,
  levels = c("Functional", "Nonfunctional")
)]

# BRCA1 functional status
dat_basal2[, BRCA1_function_group := fifelse(
  isBRCA1Biallelic == 1 | BRCA1promoterMet == 1,
  "Nonfunctional",
  "Functional"
)]
dat_basal2[, BRCA1_function_group := factor(
  BRCA1_function_group,
  levels = c("Functional", "Nonfunctional")
)]

# BRCA2 functional status
dat_basal2[, BRCA2_function_group := fifelse(
  isBRCA2Biallelic == 1,
  "Nonfunctional",
  "Functional"
)]
dat_basal2[, BRCA2_function_group := factor(
  BRCA2_function_group,
  levels = c("Functional", "Nonfunctional")
)]

############################################################
## 16. Check group sizes
############################################################

cat("\n===== Altered groups =====\n")
print(table(dat_basal2$BRCA1_altered, useNA = "ifany"))
print(table(dat_basal2$BRCA2_altered, useNA = "ifany"))
print(table(dat_basal2$BRCA1or2_altered, useNA = "ifany"))

cat("\n===== Functional groups =====\n")
print(table(dat_basal2$BRCA_function_group, useNA = "ifany"))
print(table(dat_basal2$BRCA1_function_group, useNA = "ifany"))
print(table(dat_basal2$BRCA2_function_group, useNA = "ifany"))

############################################################
## 17. Helper functions
############################################################

plot_two_group <- function(dt, group_col, title_text, out_pdf, colors_named) {
  d <- copy(dt)[!is.na(get(group_col)) & !is.na(ENPP1)]
  d[[group_col]] <- factor(d[[group_col]], levels = names(colors_named))
  
  res <- wilcox.test(ENPP1 ~ get(group_col), data = d)
  
  n_tab <- table(d[[group_col]])
  label_txt <- paste(
    paste0(names(n_tab), " n=", as.integer(n_tab)),
    collapse = ", "
  )
  
  p <- ggplot(d, aes_string(x = group_col, y = "ENPP1", color = group_col)) +
    geom_boxplot(outlier.shape = NA, width = 0.6) +
    geom_jitter(width = 0.15, alpha = 0.65, size = 1.5) +
    scale_color_manual(values = colors_named) +
    annotate(
      "text",
      x = Inf, y = Inf,
      label = paste0("Wilcoxon p = ", signif(res$p.value, 2)),
      hjust = 1.05, vjust = 1.2, size = 4
    ) +
    theme_classic(base_size = 14) +
    theme(legend.position = "none") +
    labs(
      title = paste0(title_text, "\n", label_txt),
      x = "",
      y = "ENPP1 expression"
    )
  
  ggsave(out_pdf, p, width = 5.6, height = 4.5, device = "pdf")
  return(list(test = res, counts = n_tab, plot = p))
}

plot_three_group <- function(dt, group_col, title_text, out_pdf, colors_named) {
  d <- copy(dt)[!is.na(get(group_col)) & !is.na(ENPP1)]
  d[[group_col]] <- factor(d[[group_col]], levels = names(colors_named))
  
  res_kw <- kruskal.test(ENPP1 ~ get(group_col), data = d)
  res_pw <- pairwise.wilcox.test(
    d$ENPP1,
    d[[group_col]],
    p.adjust.method = "BH"
  )
  
  n_tab <- table(d[[group_col]])
  label_txt <- paste(
    paste0(names(n_tab), " n=", as.integer(n_tab)),
    collapse = ", "
  )
  
  p <- ggplot(d, aes_string(x = group_col, y = "ENPP1", color = group_col)) +
    geom_boxplot(outlier.shape = NA, width = 0.6) +
    geom_jitter(width = 0.15, alpha = 0.65, size = 1.5) +
    scale_color_manual(values = colors_named) +
    annotate(
      "text",
      x = Inf, y = Inf,
      label = paste0("Kruskal-Wallis p = ", signif(res_kw$p.value, 2)),
      hjust = 1.05, vjust = 1.2, size = 4
    ) +
    theme_classic(base_size = 14) +
    theme(legend.position = "none") +
    labs(
      title = paste0(title_text, "\n", label_txt),
      x = "",
      y = "ENPP1 expression"
    )
  
  ggsave(out_pdf, p, width = 6.2, height = 4.8, device = "pdf")
  return(list(test = res_kw, pairwise = res_pw, counts = n_tab, plot = p))
}

############################################################
## 18. Altered analyses
############################################################

res_brca1_alt <- plot_two_group(
  dat_basal2,
  group_col = "BRCA1_altered",
  title_text = "Basal only: ENPP1 by BRCA1 altered status",
  out_pdf = "Basal_ENPP1_by_BRCA1_altered.pdf",
  colors_named = c("BRCA1 WT" = "blue", "BRCA1 altered" = "red")
)

res_brca2_alt <- plot_two_group(
  dat_basal2,
  group_col = "BRCA2_altered",
  title_text = "Basal only: ENPP1 by BRCA2 altered status",
  out_pdf = "Basal_ENPP1_by_BRCA2_altered.pdf",
  colors_named = c("BRCA2 WT" = "blue", "BRCA2 altered" = "red")
)

res_brca12_alt <- plot_two_group(
  dat_basal2,
  group_col = "BRCA1or2_altered",
  title_text = "Basal only: ENPP1 by BRCA1/2 altered status",
  out_pdf = "Basal_ENPP1_by_BRCA1or2_altered.pdf",
  colors_named = c("BRCA1/2 WT" = "blue", "BRCA1/2 altered" = "red")
)

############################################################
## 19. Functional analyses
############################################################

res_brca_func <- plot_two_group(
  dat_basal2,
  group_col = "BRCA_function_group",
  title_text = "Basal only: ENPP1 by BRCA functional status",
  out_pdf = "Basal_ENPP1_by_Functional_BRCA_Status.pdf",
  colors_named = c(
    "Functional" = "blue",
    "Nonfunctional" = "red"
  )
)

res_brca1_func <- plot_two_group(
  dat_basal2,
  group_col = "BRCA1_function_group",
  title_text = "Basal only: ENPP1 by BRCA1 functional status",
  out_pdf = "Basal_ENPP1_by_Functional_BRCA1_Status.pdf",
  colors_named = c(
    "Functional" = "blue",
    "Nonfunctional" = "red"
  )
)

res_brca2_func <- plot_two_group(
  dat_basal2,
  group_col = "BRCA2_function_group",
  title_text = "Basal only: ENPP1 by BRCA2 functional status",
  out_pdf = "Basal_ENPP1_by_Functional_BRCA2_Status.pdf",
  colors_named = c(
    "Functional" = "blue",
    "Nonfunctional" = "red"
  )
)

############################################################
## 20. Save source data
############################################################

fwrite(dat_basal2, "dat_basal2_BRCA_ENPP1.csv")

############################################################
## 21. Save summary text
############################################################

sink("Basal_BRCA_ENPP1_summary.txt")

cat("===== Altered analyses =====\n\n")

cat("BRCA1 altered vs WT\n")
print(res_brca1_alt$counts)
print(res_brca1_alt$test)
cat("\n")

cat("BRCA2 altered vs WT\n")
print(res_brca2_alt$counts)
print(res_brca2_alt$test)
cat("\n")

cat("BRCA1/2 altered vs WT\n")
print(res_brca12_alt$counts)
print(res_brca12_alt$test)
cat("\n")

cat("===== Functional analyses =====\n\n")

cat("Combined BRCA functional groups\n")
print(res_brca_func$counts)
print(res_brca_func$test)
cat("\n")

cat("BRCA1 functional groups\n")
print(res_brca1_func$counts)
print(res_brca1_func$test)
cat("\n")

cat("BRCA2 functional groups\n")
print(res_brca2_func$counts)
print(res_brca2_func$test)
cat("\n")

sink()

cat("\nAll done. Files saved to:\n", getwd(), "\n")
