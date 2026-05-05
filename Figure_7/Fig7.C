############################################################
## SCAN-B Basal only
## ENPP1 vs OS within BRCA functional groups
## Functional = WT + Monoallelic
## Nonfunctional = Biallelic loss / BRCA1 hypermethylation
## ENPP1 grouping methods: median and quartile only
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
  "ggplot2",
  "survival",
  "survminer"
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
library(survival)
library(survminer)

options(stringsAsFactors = FALSE)

############################################################
## Output directory
############################################################

outdir <- path.expand("~/Desktop/ENPP1_BRCA_analysis_outputs")
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}
cat("Output folder:", outdir, "\n")

## Global KM axis settings: force all KM plots to share the same x/y axis.
KM_XLIM <- c(0, 3000)
KM_XBREAK <- 500
KM_YLIM <- c(0, 1)


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

required_cols <- c(
  "PAM50_AIMS",
  "TumorAssay",
  "OS", "OSbin",
  "Clinical.RFI", "Clinical.RFIbin",
  "Clinical.IDFS", "Clinical.IDFSbin"
)

missing_cols <- setdiff(required_cols, names(cs))
if (length(missing_cols) > 0) {
  stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
}

############################################################
## 5. Build clinical endpoint table
############################################################

clin_tab <- cs[, .(
  PD_ID = as.character(get(pd_col_cs)),
  TumorAssay = as.character(TumorAssay),
  PAM50_AIMS = as.character(PAM50_AIMS),

  OS_time_year   = suppressWarnings(as.numeric(OS)),
  OS_event       = suppressWarnings(as.integer(OSbin)),
)]

clin_tab <- unique(clin_tab[!is.na(PD_ID) & PD_ID != ""])
cat("clin_tab dims:", dim(clin_tab), "\n")

## year -> day
clin_tab[, OS_time_day   := OS_time_year * 365.25]

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
## 7. Merge clinical + mapping
############################################################

step1_merged <- merge(clin_tab, mapping_table, by = "PD_ID")
cat("After clinical + mapping merge dims:", dim(step1_merged), "\n")
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
## 10. Final merged table
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

for (cn in brca_cols_needed) {
  dat_basal2[, (cn) := suppressWarnings(as.numeric(get(cn)))]
}

############################################################
## 13. Define functional / nonfunctional groups
############################################################

# BRCA1: Functional = WT + Monoallelic; Nonfunctional = Biallelic + methylation
dat_basal2[, BRCA1_function_group := fifelse(
  isBRCA1Biallelic == 1 | BRCA1promoterMet == 1,
  "Nonfunctional",
  "Functional"
)]
dat_basal2[, BRCA1_function_group := factor(
  BRCA1_function_group,
  levels = c("Functional", "Nonfunctional")
)]

# BRCA2: Functional = WT + Monoallelic; Nonfunctional = Biallelic
dat_basal2[, BRCA2_function_group := fifelse(
  isBRCA2Biallelic == 1,
  "Nonfunctional",
  "Functional"
)]
dat_basal2[, BRCA2_function_group := factor(
  BRCA2_function_group,
  levels = c("Functional", "Nonfunctional")
)]

# Combined BRCA1/2
dat_basal2[, BRCA12_function_group := fifelse(
  isBRCA1Biallelic == 1 |
    isBRCA2Biallelic == 1 |
    BRCA1promoterMet == 1,
  "Nonfunctional",
  "Functional"
)]
dat_basal2[, BRCA12_function_group := factor(
  BRCA12_function_group,
  levels = c("Functional", "Nonfunctional")
)]

cat("\n===== Functional group counts =====\n")
cat("BRCA1 function group:\n")
print(table(dat_basal2$BRCA1_function_group, useNA = "ifany"))

cat("BRCA2 function group:\n")
print(table(dat_basal2$BRCA2_function_group, useNA = "ifany"))

cat("BRCA1/2 function group:\n")
print(table(dat_basal2$BRCA12_function_group, useNA = "ifany"))

############################################################
## 14. Helpers: ENPP1 grouping
############################################################

make_enpp1_group_median <- function(dt, expr_col = "ENPP1") {
  d <- copy(dt)[!is.na(get(expr_col))]

  if (nrow(d) == 0) return(NULL)

  cut_val <- median(d[[expr_col]], na.rm = TRUE)
  if (length(cut_val) == 0 || is.null(cut_val) || is.na(cut_val)) return(NULL)

  d[, ENPP1_group := fifelse(
    get(expr_col) >= cut_val,
    "ENPP1_high",
    "ENPP1_low"
  )]

  d[, ENPP1_group := factor(ENPP1_group, levels = c("ENPP1_low", "ENPP1_high"))]

  return(list(
    data = d,
    cutoff = cut_val,
    plot_type = "median"
  ))
}

## Q1 vs Q4 only: bottom 25% vs top 25%; middle 50% is excluded.
make_enpp1_group_quartile_extreme <- function(dt, expr_col = "ENPP1") {
  d <- copy(dt)[!is.na(get(expr_col))]

  if (nrow(d) == 0) return(NULL)
  if (nrow(d) < 12) {
    cat("Too few samples for Q1 vs Q4 grouping\n")
    return(NULL)
  }

  q1 <- as.numeric(quantile(d[[expr_col]], probs = 0.25, na.rm = TRUE, type = 7))
  q3 <- as.numeric(quantile(d[[expr_col]], probs = 0.75, na.rm = TRUE, type = 7))

  if (is.na(q1) || is.na(q3)) return(NULL)
  if (q1 >= q3) {
    cat("Quartile cutoff collapsed (q1 >= q3); skip Q1 vs Q4 grouping\n")
    return(NULL)
  }

  d <- d[get(expr_col) <= q1 | get(expr_col) >= q3]
  if (nrow(d) == 0) return(NULL)

  d[, ENPP1_group := fifelse(
    get(expr_col) >= q3,
    "ENPP1_high_Q4",
    "ENPP1_low_Q1"
  )]

  d[, ENPP1_group := factor(ENPP1_group, levels = c("ENPP1_low_Q1", "ENPP1_high_Q4"))]

  return(list(
    data = d,
    cutoff = c(Q1 = q1, Q3 = q3),
    plot_type = "quartile_extreme"
  ))
}

## Full quartile display: bottom 25%, middle 50%, top 25%.
make_enpp1_group_quartile_three <- function(dt, expr_col = "ENPP1") {
  d <- copy(dt)[!is.na(get(expr_col))]

  if (nrow(d) == 0) return(NULL)
  if (nrow(d) < 12) {
    cat("Too few samples for 3-level quartile grouping\n")
    return(NULL)
  }

  q1 <- as.numeric(quantile(d[[expr_col]], probs = 0.25, na.rm = TRUE, type = 7))
  q3 <- as.numeric(quantile(d[[expr_col]], probs = 0.75, na.rm = TRUE, type = 7))

  if (is.na(q1) || is.na(q3)) return(NULL)
  if (q1 >= q3) {
    cat("Quartile cutoff collapsed (q1 >= q3); skip 3-level quartile grouping\n")
    return(NULL)
  }

  d[, ENPP1_group := fifelse(
    get(expr_col) <= q1,
    "ENPP1_low_Q1",
    fifelse(get(expr_col) >= q3, "ENPP1_high_Q4", "ENPP1_mid_Q2_Q3")
  )]

  d[, ENPP1_group := factor(
    ENPP1_group,
    levels = c("ENPP1_low_Q1", "ENPP1_mid_Q2_Q3", "ENPP1_high_Q4")
  )]

  return(list(
    data = d,
    cutoff = c(Q1 = q1, Q3 = q3),
    plot_type = "quartile_three"
  ))
}

format_p <- function(p) {
  if (is.na(p)) return("p = NA")
  if (p < 0.001) return("p < 0.001")
  paste0("p = ", signif(p, 3))
}

pairwise_logrank_table <- function(d) {
  lv <- levels(d$ENPP1_group)
  lv <- lv[lv %in% as.character(unique(d$ENPP1_group))]

  if (length(lv) < 2) return(data.table())

  res <- rbindlist(lapply(seq_len(ncol(combn(lv, 2))), function(i) {
    pair <- combn(lv, 2)[, i]
    dd <- d[ENPP1_group %in% pair]
    dd[, ENPP1_group_pair := droplevels(factor(as.character(ENPP1_group), levels = pair))]

    if (length(unique(dd$tmp_event)) < 2) {
      pval <- NA_real_
      chisq <- NA_real_
    } else {
      lr <- survdiff(Surv(tmp_time, tmp_event) ~ ENPP1_group_pair, data = dd)
      chisq <- lr$chisq
      pval <- 1 - pchisq(lr$chisq, df = 1)
    }

    data.table(
      comparison = paste(pair, collapse = " vs "),
      group1 = pair[1],
      group2 = pair[2],
      n1 = sum(dd$ENPP1_group_pair == pair[1], na.rm = TRUE),
      n2 = sum(dd$ENPP1_group_pair == pair[2], na.rm = TRUE),
      chisq = chisq,
      p = pval,
      p_label = format_p(pval)
    )
  }))

  return(res)
}

pretty_group_name <- function(x) {
  x <- as.character(x)
  x <- gsub("ENPP1_low_Q1", "Low/Q1", x)
  x <- gsub("ENPP1_mid_Q2_Q3", "Mid/Q2-Q3", x)
  x <- gsub("ENPP1_high_Q4", "High/Q4", x)
  x <- gsub("ENPP1_low", "Low", x)
  x <- gsub("ENPP1_high", "High", x)
  x
}

############################################################
## 15. Generic endpoint runner
############################################################

run_enpp1_endpoint_in_subset <- function(
    dt,
    subset_name,
    time_col,
    event_col,
    endpoint_label,
    ylab_txt,
    method = c("median", "quartile_extreme", "quartile_three"),
    out_pdf
) {
  method <- match.arg(method)

  d0 <- copy(dt)[
    !is.na(ENPP1) &
      !is.na(get(time_col)) &
      !is.na(get(event_col))
  ]

  d0[, tmp_time := suppressWarnings(as.numeric(get(time_col)))]
  d0[, tmp_event := suppressWarnings(as.integer(get(event_col)))]

  d0 <- d0[!is.na(tmp_time) & !is.na(tmp_event) & tmp_time > 0]

  cat("\n============================\n")
  cat("Subset :", subset_name, "\n")
  cat("Endpoint:", endpoint_label, "\n")
  cat("Method :", method, "\n")
  cat("n before grouping:", nrow(d0), "\n")
  cat("Event table:\n")
  print(table(d0$tmp_event, useNA = "ifany"))

  if (nrow(d0) < 10) {
    cat("Skip: too few samples before grouping\n")
    return(NULL)
  }

  if (length(unique(d0$tmp_event)) < 2) {
    cat("Skip: endpoint event has only one class before grouping\n")
    return(NULL)
  }

  if (method == "median") {
    grp_res <- make_enpp1_group_median(d0, expr_col = "ENPP1")
  } else if (method == "quartile_extreme") {
    grp_res <- make_enpp1_group_quartile_extreme(d0, expr_col = "ENPP1")
  } else if (method == "quartile_three") {
    grp_res <- make_enpp1_group_quartile_three(d0, expr_col = "ENPP1")
  }

  if (is.null(grp_res)) {
    cat("Skip: grouping failed\n")
    return(NULL)
  }

  d <- copy(grp_res$data)
  cutoff_val <- grp_res$cutoff

  if (!"ENPP1_group" %in% names(d)) {
    cat("Skip: ENPP1_group not generated\n")
    return(NULL)
  }

  d <- d[
    !is.na(ENPP1_group) &
      !is.na(tmp_time) &
      !is.na(tmp_event) &
      tmp_time > 0
  ]

  cat("n after grouping:", nrow(d), "\n")
  cat("ENPP1_group table:\n")
  print(table(d$ENPP1_group, useNA = "ifany"))

  if (nrow(d) == 0) {
    cat("Skip: no observations after grouping\n")
    return(NULL)
  }

  n_tab <- table(d$ENPP1_group)

  if (sum(n_tab) == 0 || length(n_tab[n_tab > 0]) < 2) {
    cat("Skip: only one non-empty ENPP1 group after grouping\n")
    return(NULL)
  }

  if (length(unique(d$tmp_event)) < 2) {
    cat("Skip: event has only one class after grouping\n")
    return(NULL)
  }

  fit_km <- survfit(Surv(tmp_time, tmp_event) ~ ENPP1_group, data = d)
  pairwise_p <- pairwise_logrank_table(d)

  pairwise_label <- if (nrow(pairwise_p) == 0) {
    "Pairwise log-rank: NA"
  } else {
    paste0(
      "Pairwise log-rank (Mantel-Cox):\n",
      paste0(
        pretty_group_name(pairwise_p$group1), " vs ",
        pretty_group_name(pairwise_p$group2), ": ",
        pairwise_p$p_label,
        collapse = "\n"
      )
    )
  }

  ## Use fixed x-axis for every KM plot, regardless of the maximum follow-up time.
  xmax <- KM_XLIM[2]

  group_levels <- levels(d$ENPP1_group)
  group_levels <- group_levels[group_levels %in% names(n_tab)]
  group_levels <- group_levels[n_tab[group_levels] > 0]

  legend_labs <- paste0(pretty_group_name(group_levels), " (n=", as.integer(n_tab[group_levels]), ")")

  ## Color rule:
  ## low = #00BFC4, high = #F8766D; middle 50% = grey55
  if (length(group_levels) == 2) {
    palette_use <- c("#00BFC4", "#F8766D")
  } else {
    palette_use <- c("#00BFC4", "grey55", "#F8766D")
  }

  title_method <- switch(
    method,
    median = "median split",
    quartile_extreme = "Q1 vs Q4",
    quartile_three = "Q1 / middle 50% / Q4"
  )

  title_txt <- paste0(subset_name, " | ", endpoint_label, " | ENPP1 ", title_method)

  cutoff_label <- if (method == "median") {
    paste0("median = ", round(as.numeric(cutoff_val), 3))
  } else {
    paste0("Q1 = ", round(as.numeric(cutoff_val["Q1"]), 3),
           "; Q3 = ", round(as.numeric(cutoff_val["Q3"]), 3))
  }

  ## Auto n-number label shown directly on each KM plot
  n_label <- paste0(
    "n: ",
    paste0(pretty_group_name(group_levels), "=", as.integer(n_tab[group_levels]), collapse = "; ")
  )

  p <- ggsurvplot(
    fit_km,
    data = d,
    pval = FALSE,
    risk.table = FALSE,
    conf.int = FALSE,
    censor = TRUE,
    palette = palette_use,
    legend.title = "",
    legend.labs = legend_labs,
    ggtheme = theme_classic(base_size = 14),
    title = title_txt,
    xlab = "Time (days)",
    ylab = ylab_txt,
    break.time.by = KM_XBREAK,
    xlim = KM_XLIM,
    ylim = KM_YLIM
  )

  p$plot <- p$plot +
    labs(subtitle = n_label) +
    annotate(
      "text",
      x = xmax * 0.98, y = 0.18,
      hjust = 1,
      label = paste0(cutoff_label, "\n", pairwise_label),
      size = 3.5
    ) +
    theme(
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5),
      legend.position = "bottom",
      legend.box = "horizontal",
      plot.margin = margin(8, 28, 5, 10)
    )

  ggsave(
    filename = out_pdf,
    plot = p$plot,
    width = 7.8,
    height = 5.1,
    device = "pdf"
  )

  fwrite(d, sub("\\.pdf$", ".csv", out_pdf))
  fwrite(pairwise_p, sub("\\.pdf$", "__pairwise_logrank.csv", out_pdf))

  return(list(
    subset = subset_name,
    endpoint = endpoint_label,
    method = method,
    cutoff = cutoff_val,
    counts = n_tab,
    pairwise_logrank = pairwise_p
  ))
}

############################################################
## 16. Build subgroup list
############################################################

subset_list <- list(
  "BRCA1_Functional"    = dat_basal2[BRCA1_function_group == "Functional"],
  "BRCA1_Nonfunctional" = dat_basal2[BRCA1_function_group == "Nonfunctional"],

  "BRCA2_Functional"    = dat_basal2[BRCA2_function_group == "Functional"],
  "BRCA2_Nonfunctional" = dat_basal2[BRCA2_function_group == "Nonfunctional"],

  "BRCA1or2_Functional"    = dat_basal2[BRCA12_function_group == "Functional"],
  "BRCA1or2_Nonfunctional" = dat_basal2[BRCA12_function_group == "Nonfunctional"]
)

cat("\n===== Subset sizes =====\n")
for (nm in names(subset_list)) {
  cat("\n", nm, ":", nrow(subset_list[[nm]]), "\n")
  cat("OS event table:\n")
  print(table(subset_list[[nm]]$OS_event, useNA = "ifany"))
}

############################################################
## 17. Endpoint config
############################################################

endpoint_cfg <- list(
  OS = list(
    time_col = "OS_time_day",
    event_col = "OS_event",
    ylab = "Overall survival probability"
  )
)

############################################################
## 18. Run all analyses
############################################################

results_all <- list()

for (sub_nm in names(subset_list)) {
  dt_sub <- subset_list[[sub_nm]]

  for (ep in names(endpoint_cfg)) {
    time_col  <- endpoint_cfg[[ep]]$time_col
    event_col <- endpoint_cfg[[ep]]$event_col
    ylab_txt  <- endpoint_cfg[[ep]]$ylab

    ## median
    key1 <- paste(sub_nm, ep, "median", sep = "__")
    results_all[[key1]] <- run_enpp1_endpoint_in_subset(
      dt = dt_sub,
      subset_name = sub_nm,
      time_col = time_col,
      event_col = event_col,
      endpoint_label = ep,
      ylab_txt = ylab_txt,
      method = "median",
      out_pdf = file.path(outdir, paste0(sub_nm, "__", ep, "__ENPP1_median.pdf"))
    )

    ## quartile full display: Q1 / middle 50% / Q4
    key2 <- paste(sub_nm, ep, "quartile_three", sep = "__")
    results_all[[key2]] <- run_enpp1_endpoint_in_subset(
      dt = dt_sub,
      subset_name = sub_nm,
      time_col = time_col,
      event_col = event_col,
      endpoint_label = ep,
      ylab_txt = ylab_txt,
      method = "quartile_three",
      out_pdf = file.path(outdir, paste0(sub_nm, "__", ep, "__ENPP1_quartile_3groups.pdf"))
    )

    ## quartile extreme comparison: Q1 vs Q4 only
    key3 <- paste(sub_nm, ep, "quartile_extreme", sep = "__")
    results_all[[key3]] <- run_enpp1_endpoint_in_subset(
      dt = dt_sub,
      subset_name = sub_nm,
      time_col = time_col,
      event_col = event_col,
      endpoint_label = ep,
      ylab_txt = ylab_txt,
      method = "quartile_extreme",
      out_pdf = file.path(outdir, paste0(sub_nm, "__", ep, "__ENPP1_Q1_vs_Q4.pdf"))
    )
  }
}

############################################################
## 19. Save merged table
############################################################

fwrite(dat_basal2, file.path(outdir, "dat_basal2_with_BRCA_function_groups_and_ENPP1_all_endpoints.csv"))

############################################################
## 20. Summary report
############################################################

sink(file.path(outdir, "ENPP1_survival_within_BRCA_function_groups_summary.txt"))

cat("===== ENPP1 survival within BRCA functional groups =====\n\n")

for (nm in names(results_all)) {
  cat("---- ", nm, " ----\n", sep = "")

  res <- results_all[[nm]]

  if (is.null(res)) {
    cat("Skipped\n\n")
  } else {
    cat("Subset        :", res$subset, "\n")
    cat("Endpoint      :", res$endpoint, "\n")
    cat("Method        :", res$method, "\n")
    cat("Cutoff        :\n")
    print(res$cutoff)
    cat("Counts        :\n")
    print(res$counts)
    cat("Pairwise log-rank tests (Mantel-Cox; no overall p value shown):\n")
    print(res$pairwise_logrank)
    cat("\n")
  }
}

sink()

cat("\nAll done. Files saved to:\n", outdir, "\n")
