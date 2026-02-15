#### === TCGA-BRCA: ENPP1 vs HRD (Survival analysis)

## Install / load packages
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
bio_pkgs <- c("TCGAbiolinks","SummarizedExperiment","edgeR")
for (p in bio_pkgs) if (!requireNamespace(p, quietly=TRUE)) BiocManager::install(p, ask=FALSE, update=FALSE)
cran_pkgs <- c("dplyr","stringr","ggplot2","survival","survminer","readr","tibble")
to_install <- setdiff(cran_pkgs, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install, repos="https://cloud.r-project.org")

library(TCGAbiolinks)
library(SummarizedExperiment)
library(edgeR)
library(dplyr)
library(stringr)
library(ggplot2)
library(survival)
library(survminer)
library(readr)
library(tibble)

HRD_CUTOFF <- 42
outdir <- "results_tcga_enpp1_hrd"; dir.create(outdir, showWarnings = FALSE)

# helpers
.pick <- function(df, candidates) { x <- intersect(candidates, names(df)); if (length(x)) x[1] else NA_character_ }
quiet <- function(expr) suppressWarnings(suppressMessages(expr))

#1) Download the files
message("1) Downloading TCGA-BRCA STAR-Counts (Primary Tumor only)…")
query.exp <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = "Primary Tumor"
)
GDCdownload(query.exp, files.per.chunk = 50)
se.exp <- GDCprepare(query.exp)   # RangedSummarizedExperiment

## Clinical from colData(se.exp) → OS_days / OS_event (patient-level)
cd <- as.data.frame(colData(se.exp))
patient_col <- .pick(cd, c("patient","bcr_patient_barcode","submitter_id","case_submitter_id",
                           "cases.submitter_id","patient_submitter_id"))
death_col   <- .pick(cd, c("days_to_death","demographic.days_to_death","patient.days_to_death",
                           "days_to_death.1","vital_status_days_to_death"))
lfu_col     <- .pick(cd, c("days_to_last_follow_up","days_to_last_followup","last_contact_days_to",
                           "follow_up_days_to","days_to_last_known_alive","patient.days_to_last_followup"))
vital_col   <- .pick(cd, c("vital_status","patient.vital_status","vital_status_at_last_contact","vital_status.1"))

cd$OS_days <- suppressWarnings(as.numeric(dplyr::coalesce(
  if (!is.na(death_col)) cd[[death_col]] else NA_real_,
  if (!is.na(lfu_col))   cd[[lfu_col]]   else NA_real_
)))
cd$OS_event <- ifelse(
  grepl("deceased|dead", tolower(as.character(if (!is.na(vital_col)) cd[[vital_col]] else NA_character_))),
  1, 0
)
clin <- cd %>%
  transmute(patient12 = substr(.data[[patient_col]], 1, 12),
            OS_days = OS_days, OS_event = OS_event) %>%
  filter(!is.na(patient12)) %>%
  distinct(patient12, .keep_all = TRUE)

##2) ENPP1 log2-CPM (compute ONLY for ENPP1 to save RAM)
rdat <- as.data.frame(rowData(se.exp))
symcol <- .pick(rdat, c("external_gene_name","gene_name","GeneSymbol","symbol"))
if (is.na(symcol)) stop("No gene symbol column in rowData(se.exp)")

idx <- which(rdat[[symcol]] == "ENPP1")
stopifnot(length(idx) >= 1)

lib.size <- colSums(assay(se.exp))                               # vector
ENPP1_counts <- assay(se.exp)[idx[1], , drop = FALSE]             # 1 x n
ENPP1_lcpm_vec <- as.numeric(edgeR::cpm(ENPP1_counts, log = TRUE, prior.count = 1,
                                        lib.size = lib.size))
samps <- colnames(se.exp)

expr_df <- tibble(
  barcode     = samps,
  ENPP1_lcpm  = ENPP1_lcpm_vec,
  sample16    = substr(samps, 1, 16),
  patient12   = substr(samps, 1, 12)
)

# free bulky pieces
rm(ENPP1_counts, lib.size, rdat); gc()

##3) HRD scores (load local DDRscores.txt or download from mirrors)
message("2) Loading HRD (genomic scar) scores…")
hrd <- NULL
# Try local first (edit if you saved it elsewhere)
hrd_local_path <- "~/Downloads/TCGA_DDR_Data_Resources/DDRscores.txt"
hrd_local_path <- path.expand(hrd_local_path)

if (file.exists(hrd_local_path)) {
  message("Using local HRD file: ", hrd_local_path)
  hrd <- read_tsv(hrd_local_path, show_col_types = FALSE)
} else {
  # search in common places
  cand <- c(
    list.files(path.expand("~/Downloads"), pattern="^DDRscores\\.txt$", recursive=TRUE, full.names=TRUE),
    list.files(getwd(), pattern="^DDRscores\\.txt$", recursive=TRUE, full.names=TRUE)
  )
  if (length(cand)) {
    message("Found HRD file: ", cand[1])
    hrd <- read_tsv(cand[1], show_col_types = FALSE)
  } else {
    # mirrors
    message("No local HRD file; attempting download…")
    hrd_urls <- c(
      "https://raw.githubusercontent.com/GerkeLab/TCGAhrd/master/data/intermediateFIles/DDRscores.txt",
      "https://raw.githubusercontent.com/GerkeLab/TCGAhrd/master/data/DDRscores.txt",
      "https://raw.githubusercontent.com/GerkeLab/TCGAhrd/master/data-raw/DDR_Scores/DDRscores.txt"
    )
    for (u in hrd_urls) {
      message("Trying HRD URL: ", u)
      tmp <- tryCatch(readr::read_tsv(u, show_col_types = FALSE), error = function(e) NULL)
      if (!is.null(tmp) && nrow(tmp) > 0) { hrd <- tmp; break }
    }
    if (is.null(hrd)) stop("Could not load HRD scores. Provide DDRscores.txt locally or check network.")
  }
}

#4) Filter HRD to BRCA if a tumor/cohort column exists
pick_hrd_col <- function(df, cands) { x <- intersect(cands, names(df)); if (length(x)) x[1] else NA_character_ }
tumor_col <- pick_hrd_col(hrd, c("tumor","Tumor","cancer","Cancer","study","Study","project","Project","cohort"))
if (!is.na(tumor_col)) {
  hrd <- hrd[grepl("BRCA|TCGA[-_ ]?BRCA", hrd[[tumor_col]], ignore.case = TRUE), , drop = FALSE]
}

#5) Find barcode column (patient or sample)
find_barcode_col <- function(df){
  exact <- intersect(c("barcode","Barcode","bcr_patient_barcode","bcr_sample_barcode",
                       "Sample","sample","Tumor_Sample_Barcode","Aliquot_Barcode",
                       "submitter_id","case_submitter_id"), names(df))
  if (length(exact)) return(exact[1])
  for (nm in names(df)) {
    v <- df[[nm]]; if (is.factor(v)) v <- as.character(v)
    if (is.character(v) && any(grepl("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}", v))) return(nm)
  }
  NA_character_
}
barcode_col <- find_barcode_col(hrd)
if (is.na(barcode_col)) stop("No sample/patient barcode column detected in HRD table. Use names(hrd) to inspect.")

#6)Ensure HRDsum column
if (!"HRDsum" %in% names(hrd)) {
  if (all(c("LOH","LST","TAI") %in% names(hrd))) {
    hrd <- hrd %>% mutate(HRDsum = suppressWarnings(as.numeric(LOH)) +
                            suppressWarnings(as.numeric(LST)) +
                            suppressWarnings(as.numeric(TAI)))
  } else {
    alt <- .pick(hrd, c("HRD","HRD_score","HRDscore","HRD_Score"))
    if (is.na(alt)) stop("No HRDsum and no LOH/LST/TAI (or alternate HRD column) in HRD table.")
    hrd <- hrd %>% dplyr::rename(HRDsum = !!alt)
  }
}

#7) Normalize IDs and build HRD calls at patient/sample level
hrd_norm <- hrd %>%
  mutate(
    .raw_id   = as.character(.data[[barcode_col]]),
    patient12 = substr(.raw_id, 1, 12),
    sample16  = if_else(nchar(.raw_id) >= 16, substr(.raw_id, 1, 16), NA_character_),
    HRDsum    = suppressWarnings(as.numeric(HRDsum))
  ) %>%
  filter(!is.na(patient12)) %>%
  arrange(desc(HRDsum)) %>%                         # if duplicates, keep highest
  distinct(patient12, .keep_all = TRUE)

hrd_calls <- hrd_norm %>%
  mutate(HRD_status = case_when(
    is.na(HRDsum)        ~ NA_character_,
    HRDsum >= HRD_CUTOFF ~ "HRD-high",
    TRUE                 ~ "HRD-low"
  )) %>%
  select(patient12, sample16, HRDsum, HRD_status)

message("HRD table rows kept: ", nrow(hrd_calls))

##8) Merge (patient-level key), build expression and survival cohorts
dat_pre <- expr_df %>%
  left_join(hrd_calls, by = "patient12")
message("After HRD join: matched HRD for ", sum(!is.na(dat_pre$HRDsum)),
        " / ", nrow(dat_pre), " expression samples")

dat_all <- dat_pre %>%
  left_join(clin, by = "patient12")

dat_expr <- dat_all %>% filter(!is.na(ENPP1_lcpm), !is.na(HRD_status))
dat_surv <- dat_expr %>% filter(!is.na(OS_days), !is.na(OS_event)) %>%
  mutate(OS_months = OS_days / 30.44)

message("dat_expr N=", nrow(dat_expr),
        " | HRD-low=",  sum(dat_expr$HRD_status=="HRD-low", na.rm=TRUE),
        " | HRD-high=", sum(dat_expr$HRD_status=="HRD-high", na.rm=TRUE))
message("dat_surv N=", nrow(dat_surv), " (with OS)")


  
##9) --- ENPP1 ~ HRD violin (PDF, HRD-low first) ---
  # ensure order
  dat_expr <- dat_expr %>%
    dplyr::mutate(HRD_status = factor(HRD_status, levels = c("HRD-low","HRD-high")))
  
  if (length(na.omit(unique(dat_expr$HRD_status))) == 2) {
    wil <- wilcox.test(ENPP1_lcpm ~ HRD_status, data = dat_expr)
    p_lab <- paste0("Wilcoxon p = ", signif(wil$p.value, 3))
  } else {
    only <- as.character(unique(na.omit(dat_expr$HRD_status)))[1]
    n_only <- sum(dat_expr$HRD_status == only, na.rm = TRUE)
    p_lab <- paste0("Wilcoxon not computed: only ", only, " (n=", n_only, ")")
  }
  
  p_expr <- ggplot(dat_expr, aes(HRD_status, ENPP1_lcpm, fill = HRD_status)) +
    geom_violin(trim = FALSE, alpha = .35) +
    geom_boxplot(width = .15, outlier.shape = NA) +
    geom_jitter(width = .15, alpha = .6, size = 1) +
    labs(x = paste0("HRD status (HRDsum ≥ ", HRD_CUTOFF, " = high)"),
         y = "ENPP1 (log2-CPM)",
         title = "TCGA-BRCA: ENPP1 (log2-CPM) by HRD status") +
    annotate("text", x = 1.5, y = max(dat_expr$ENPP1_lcpm, na.rm = TRUE),
             label = p_lab, vjust = -0.5)
  
  pdf(file.path(outdir, "tcga_brca_ENPP1lcpm_by_HRD_status.pdf"),
      width = 6.4, height = 4.4)
  print(p_expr)
  dev.off()
  
##10)---------KM plot###
  run_surv <- function(df, label){
    cut <- median(df$ENPP1_lcpm, na.rm=TRUE)
    df$ENPP1_group <- ifelse(df$ENPP1_lcpm > cut, "High","Low")
    
    fit <- survfit(Surv(OS_months, OS_event) ~ ENPP1_group, data=df)
    cox <- coxph(Surv(OS_months, OS_event) ~ ENPP1_lcpm, data=df)
    
    g <- ggsurvplot(fit, data=df, pval=TRUE, risk.table=TRUE,
                    title=paste0("TCGA-BRCA ", label, ": OS by ENPP1 (median split)"),
                    xlab="Months", ylab="Overall survival probability")
    
    ## === only this line changed: write PDF instead of PNG ===
    ggsave(file.path(outdir, paste0("tcga_brca_OS_by_ENPP1_", gsub(" ","_",label), ".pdf")),
           g$plot, width = 6.4, height = 4.4, device = "pdf")
    
    tibble(group=label, n=nrow(df), median_cut=cut,
           HR=exp(coef(cox))[1],
           HR_low=exp(confint(cox))[1],
           HR_high=exp(confint(cox))[2],
           HR_p=summary(cox)$coefficients[,"Pr(>|z|)"][1])
  }
  surv_out <- list()
  if (nrow(filter(dat_surv, HRD_status=="HRD-low"))  >= 10) surv_out[["HRD-low"]]  <- run_surv(filter(dat_surv, HRD_status=="HRD-low"),  "HRD-low")
  if (nrow(filter(dat_surv, HRD_status=="HRD-high")) >= 10) surv_out[["HRD-high"]] <- run_surv(filter(dat_surv, HRD_status=="HRD-high"), "HRD-high")
  
