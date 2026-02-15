#### === TCGA-BRCA TNBC: ENPP1 expression in HRD low versus high 

# (A) Raise virtual memory cap (mac/Linux; harmless elsewhere)
if (.Platform$OS.type != "windows") Sys.setenv("R_MAX_VSIZE" = "50Gb")

# (B) Install / load packages
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
outdir_base <- "results_tcga_enpp1_hrd"; dir.create(outdir_base, showWarnings = FALSE)
outdir <- file.path(outdir_base, "TNBC"); dir.create(outdir, showWarnings = FALSE)

# --- helpers ---
.pick <- function(df, cands){ x <- intersect(cands, names(df)); if (length(x)) x[1] else NA_character_ }
quiet <- function(expr) suppressWarnings(suppressMessages(expr))

# 1) Identify TNBC patients via PAM50 Basal-like in subtype table
message("1) Identify TNBC patients (PAM50 Basal-like proxy)…")
sub <- TCGAbiolinks::TCGAquery_subtype("brca")

pat_col <- .pick(sub, c("patient","bcr_patient_barcode","submitter_id","case_submitter_id"))
pam_col <- .pick(sub, c("BRCA_Subtype_PAM50","PAM50","PAM50.mRNA"))
if (any(is.na(c(pat_col, pam_col)))) stop("Subtype table missing patient or PAM50 column. Run names(sub).")

tnbc_patients <- sub %>%
  transmute(patient12 = substr(.data[[pat_col]], 1, 12),
            PAM50 = as.character(.data[[pam_col]])) %>%
  filter(!is.na(PAM50), PAM50 %in% c("Basal-like","Basal")) %>%
  distinct(patient12) %>% pull(patient12)

message("TNBC (Basal-like) patients: ", length(tnbc_patients))

# 2) Build TNBC PRIMARY TUMOR file IDs straight from metadata
# (these column names match your meta screenshot)
pat_meta_col <- "cases.submitter_id"   # patient-level ID
stype_col    <- "sample_type"          # "Primary Tumor", etc.
aliquot_col  <- "submitter_id"         # full aliquot ID (not actually needed now)
fileid_col   <- "id"                   # GDC file UUID  <-- use this

stopifnot(all(c(pat_meta_col, stype_col, aliquot_col, fileid_col) %in% names(meta)))

sel <- meta[[stype_col]] == "Primary Tumor" &
  substr(meta[[pat_meta_col]], 1, 12) %in% tnbc_patients

file_ids_keep <- unique(na.omit(as.character(meta[[fileid_col]][sel])))
aliquots_keep <- unique(na.omit(as.character(meta[[aliquot_col]][sel])))

message("TNBC primary tumor files: ", length(file_ids_keep))
message("First 3 aliquots (for sanity): ", paste(head(aliquots_keep, 3), collapse = " | "))
stopifnot(length(file_ids_keep) > 0)

# 3) Query by FILE IDs (not barcode) to avoid matching issues
query.exp <- GDCquery(
  project       = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

# keep only our file IDs inside the query object
query.exp$results[[1]] <- dplyr::filter(query.exp$results[[1]], id %in% file_ids_keep)
GDCdownload(query.exp, files.per.chunk = 30)
se <- GDCprepare(query.exp)
# small TNBC-only object

# 4) Clinical from colData(se) → OS_days / OS_event (patient-level)
cd <- as.data.frame(colData(se))
patient_col_cd <- .pick(cd, c("patient","bcr_patient_barcode","submitter_id","case_submitter_id",
                              "cases.submitter_id","patient_submitter_id"))
death_col <- .pick(cd, c("days_to_death","demographic.days_to_death","patient.days_to_death",
                         "days_to_death.1","vital_status_days_to_death"))
lfu_col   <- .pick(cd, c("days_to_last_follow_up","days_to_last_followup","last_contact_days_to",
                         "follow_up_days_to","days_to_last_known_alive","patient.days_to_last_followup"))
vital_col <- .pick(cd, c("vital_status","patient.vital_status","vital_status_at_last_contact","vital_status.1"))

cd$OS_days <- suppressWarnings(as.numeric(dplyr::coalesce(
  if (!is.na(death_col)) cd[[death_col]] else NA_real_,
  if (!is.na(lfu_col))   cd[[lfu_col]]   else NA_real_
)))
cd$OS_event <- ifelse(
  grepl("deceased|dead", tolower(as.character(if (!is.na(vital_col)) cd[[vital_col]] else NA_character_))),
  1, 0
)
clin <- cd %>%
  transmute(patient12 = substr(.data[[patient_col_cd]], 1, 12),
            OS_days = OS_days, OS_event = OS_event) %>%
  filter(!is.na(patient12)) %>%
  distinct(patient12, .keep_all = TRUE)

# 5) ENPP1 log2-CPM (compute only that gene; memory-lean)
rdat <- as.data.frame(rowData(se))
symcol <- .pick(rdat, c("external_gene_name","gene_name","GeneSymbol","symbol"))
if (is.na(symcol)) stop("No gene symbol column in rowData(se). Run names(rowData(se)).")

idx <- which(rdat[[symcol]] == "ENPP1"); stopifnot(length(idx) >= 1)
lib.size <- colSums(assay(se))
ENPP1_counts <- assay(se)[idx[1], , drop = FALSE]
ENPP1_lcpm <- as.numeric(edgeR::cpm(ENPP1_counts, log = TRUE, prior.count = 1, lib.size = lib.size))
samps <- colnames(se)

expr_df <- tibble(
  barcode    = samps,                    # full aliquot barcodes
  ENPP1_lcpm = ENPP1_lcpm,
  sample16   = substr(samps, 1, 16),
  patient12  = substr(samps, 1, 12)
)
rm(ENPP1_counts, lib.size, rdat); gc()

# 6) HRD scores (HRDsum) → HRD_status (≥42 = high)
message("4) Loading HRDsum and building HRD_status…")
hrd <- NULL
hrd_local_path <- path.expand("~/Downloads/TCGA_DDR_Data_Resources/DDRscores.txt")
if (file.exists(hrd_local_path)) {
  hrd <- readr::read_tsv(hrd_local_path, show_col_types = FALSE)
} else {
  cand <- c(
    list.files(path.expand("~/Downloads"), pattern="^DDRscores\\.txt$", recursive=TRUE, full.names=TRUE),
    list.files(getwd(), pattern="^DDRscores\\.txt$", recursive=TRUE, full.names=TRUE)
  )
  if (length(cand)) {
    hrd <- readr::read_tsv(cand[1], show_col_types = FALSE)
  } else {
    urls <- c(
      "https://raw.githubusercontent.com/GerkeLab/TCGAhrd/master/data/intermediateFIles/DDRscores.txt",
      "https://raw.githubusercontent.com/GerkeLab/TCGAhrd/master/data/DDRscores.txt",
      "https://raw.githubusercontent.com/GerkeLab/TCGAhrd/master/data-raw/DDR_Scores/DDRscores.txt"
    )
    for (u in urls) {
      tmp <- tryCatch(readr::read_tsv(u, show_col_types = FALSE), error = function(e) NULL)
      if (!is.null(tmp) && nrow(tmp) > 0) { hrd <- tmp; break }
    }
    if (is.null(hrd)) stop("Could not load HRD scores. Provide DDRscores.txt locally or check network.")
  }
}

# keep BRCA rows if a cohort column exists
tumor_col <- .pick(hrd, c("tumor","Tumor","cancer","Cancer","study","Study","project","Project","cohort"))
if (!is.na(tumor_col)) hrd <- hrd[grepl("BRCA|TCGA[-_ ]?BRCA", hrd[[tumor_col]], ignore.case = TRUE), , drop = FALSE]

# detect barcode column and ensure HRDsum
barcode_col <- .pick(hrd, c("barcode","Barcode","bcr_patient_barcode","bcr_sample_barcode",
                            "Sample","sample","Tumor_Sample_Barcode","Aliquot_Barcode",
                            "submitter_id","case_submitter_id"))
if (is.na(barcode_col)) {
  # fallback: scan for any column with TCGA-like strings
  bc_cand <- names(hrd)[vapply(hrd, function(v) any(grepl("^TCGA-", as.character(v))), logical(1))]
  if (length(bc_cand)) barcode_col <- bc_cand[1]
}
if (is.na(barcode_col)) stop("No barcode column in HRD table. Run names(hrd).")

if (!"HRDsum" %in% names(hrd)) {
  if (all(c("LOH","LST","TAI") %in% names(hrd))) {
    hrd <- hrd %>% mutate(HRDsum = as.numeric(LOH) + as.numeric(LST) + as.numeric(TAI))
  } else {
    alt <- .pick(hrd, c("HRD","HRD_score","HRDscore","HRD_Score"))
    if (is.na(alt)) stop("No HRDsum and no LOH/LST/TAI or alternate HRD column.")
    hrd <- hrd %>% dplyr::rename(HRDsum = !!alt)
  }
}

hrd_calls <- hrd %>%
  mutate(
    .raw_id   = as.character(.data[[barcode_col]]),
    patient12 = substr(.raw_id, 1, 12),
    HRDsum    = suppressWarnings(as.numeric(HRDsum))
  ) %>%
  filter(patient12 %in% unique(expr_df$patient12)) %>%
  arrange(desc(HRDsum)) %>%
  distinct(patient12, .keep_all = TRUE) %>%
  mutate(HRD_status = case_when(
    is.na(HRDsum)        ~ NA_character_,
    HRDsum >= HRD_CUTOFF ~ "HRD-high",
    TRUE                 ~ "HRD-low"
  )) %>%
  select(patient12, HRDsum, HRD_status)

message("HRD matches in TNBC cohort: ", nrow(hrd_calls), " / ", length(unique(expr_df$patient12)))

# 7) Merge TNBC data
dat_pre <- expr_df %>% left_join(hrd_calls, by = "patient12")
dat_all <- dat_pre %>% left_join(clin, by = "patient12")

dat_expr <- dat_all %>% filter(!is.na(ENPP1_lcpm), !is.na(HRD_status))
dat_surv <- dat_expr %>% filter(!is.na(OS_days), !is.na(OS_event)) %>%
  mutate(OS_months = OS_days / 30.44)

message("TNBC expr N=", nrow(dat_expr),
        " | HRD-low=",  sum(dat_expr$HRD_status=="HRD-low", na.rm=TRUE),
        " | HRD-high=", sum(dat_expr$HRD_status=="HRD-high", na.rm=TRUE))
message("TNBC surv N=", nrow(dat_surv))

# 8) ENPP1 difference (HRD-high vs HRD-low)
if (length(na.omit(unique(dat_expr$HRD_status))) == 2) {
  wil <- wilcox.test(ENPP1_lcpm ~ HRD_status, data = dat_expr)
  tt  <- t.test(ENPP1_lcpm ~ HRD_status, data = dat_expr)
  out_diff <- tibble(
    test      = c("Wilcoxon","t-test"),
    pval      = c(wil$p.value, tt$p.value),
    LOW_mean  = mean(dat_expr$ENPP1_lcpm[dat_expr$HRD_status=="HRD-low"]),
    HIGH_mean = mean(dat_expr$ENPP1_lcpm[dat_expr$HRD_status=="HRD-high"]),
    LOW_n     = sum(dat_expr$HRD_status=="HRD-low"),
    HIGH_n    = sum(dat_expr$HRD_status=="HRD-high")
  )
  write_csv(out_diff, file.path(outdir, "TNBC_ENPP1lcpm_vs_HRD_stats.csv"))
  
  p_expr <- ggplot(dat_expr, aes(HRD_status, ENPP1_lcpm, fill=HRD_status)) +
    geom_violin(trim=FALSE, alpha=.35) +
    geom_boxplot(width=.15, outlier.shape=NA) +
    geom_jitter(width=.15, alpha=.6, size=1) +
    labs(x=paste0("HRD status (HRDsum ≥ ", HRD_CUTOFF, " = high)"),
         y="ENPP1 (log2-CPM)",
         title="TNBC (TCGA-BRCA): ENPP1 by HRD status") +
    annotate("text", x=1.5, y=max(dat_expr$ENPP1_lcpm, na.rm=TRUE),
             label=paste0("Wilcoxon p=", signif(wil$p.value,3)))
  ggsave(file.path(outdir, "TNBC_ENPP1lcpm_by_HRD_status.png"),
         p_expr, width=6.4, height=4.4, dpi=300)
} else {
  message("TNBC: skipping ENPP1~HRD test (need both HRD groups).")
}
## Export the data 
message("Done. Results in: ", normalizePath(outdir), "\n",
        "- TNBC_ENPP1lcpm_by_HRD_status.png\n",
        "- TNBC_ENPP1lcpm_vs_HRD_stats.csv)

