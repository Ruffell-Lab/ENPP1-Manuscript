######clean the data########
library(data.table)
library(stringr)
library(openxlsx)

setwd("~/Desktop/ISPY2")

## =========================================================
## 1. Read GSE194040 all-gene expression matrix
## =========================================================

expr_file <- list.files(
  pattern = "GSE194040.*geneLevel.*n988.*\\.txt$|GSE194040.*geneLevel.*n988.*\\.txt\\.gz$",
  full.names = TRUE
)[1]

expr_raw <- fread(expr_file)

dim(expr_raw)
expr_raw[1:5, 1:8]

## V1 = gene symbol
expr_raw <- expr_raw[!is.na(V1) & V1 != ""]
expr_raw[, gene := make.names(V1, unique = TRUE)]

sample_cols <- setdiff(names(expr_raw), c("V1", "gene"))

## transpose: genes x patients -> patients x genes
expr_mat <- as.data.table(t(as.matrix(expr_raw[, ..sample_cols])))

colnames(expr_mat) <- expr_raw$gene

expr_mat[, ResID := sample_cols]
setcolorder(expr_mat, "ResID")
expr_mat[, ResID := as.character(ResID)]

## check ENPP1
grep("^ENPP1$", names(expr_mat), value = TRUE)
summary(expr_mat$ENPP1)

## export expression-only file
fwrite(
  expr_mat,
  "ISPY2_ResID_ALL_gene_expression.csv"
)

## =========================================================
## 2. Read GSE196096 SOFT file and extract clinical annotation
## =========================================================

soft_file <- "GSE196096_family.soft"

x <- readLines(soft_file, warn = FALSE)

gsm_start <- grep("^\\^SAMPLE = GSM", x)
gsm_end <- c(gsm_start[-1] - 1, length(x))

sample_meta_list <- list()

for (i in seq_along(gsm_start)) {
  block <- x[gsm_start[i]:gsm_end[i]]
  
  gsm <- sub("^\\^SAMPLE = ", "", block[1])
  
  title_line <- block[grep("^!Sample_title = ", block)]
  title <- sub("^!Sample_title = ", "", title_line)
  
  chars <- block[grep("^!Sample_characteristics_ch", block)]
  chars <- sub("^!Sample_characteristics_ch[0-9]+ = ", "", chars)
  
  dt <- data.table(GSM = gsm, title = title)
  
  if (length(chars) > 0) {
    char_dt <- rbindlist(lapply(chars, function(z) {
      parts <- strsplit(z, ": ", fixed = TRUE)[[1]]
      
      if (length(parts) >= 2) {
        data.table(
          variable = parts[1],
          value = paste(parts[-1], collapse = ": ")
        )
      } else {
        data.table(variable = z, value = NA_character_)
      }
    }), fill = TRUE)
    
    char_wide <- dcast(char_dt, . ~ variable, value.var = "value")
    char_wide[, "." := NULL]
    
    dt <- cbind(dt, char_wide)
  }
  
  sample_meta_list[[i]] <- dt
}

sample_meta <- rbindlist(sample_meta_list, fill = TRUE)

dim(sample_meta)
names(sample_meta)
head(sample_meta)

## =========================================================
## 3. Clean clinical/sample annotation
## =========================================================

map196096 <- sample_meta[, .(
  ResID = as.character(`patient id`),
  GSM,
  title,
  arm,
  `arm (short name)`,
  her2,
  hr,
  mp,
  pcr,
  tissue,
  `adjustment method`,
  `rppa array (normalization key)`
)]

map196096 <- map196096[!is.na(ResID) & ResID != ""]

## one row per ResID
map196096 <- map196096[!duplicated(ResID)]

dim(map196096)
head(map196096)

## =========================================================
## 4. Merge clinical annotation with all-gene expression
## =========================================================

merged_clean <- merge(
  map196096,
  expr_mat,
  by = "ResID",
  all.x = TRUE
)

dim(merged_clean)
sum(is.na(merged_clean$ENPP1))
sum(duplicated(merged_clean$ResID))
summary(merged_clean$ENPP1)

## matched expression-only cohort
merged_expr_only <- merged_clean[!is.na(ENPP1)]

dim(merged_expr_only)

## =========================================================
## 5. Export final files
## =========================================================

fwrite(
  merged_clean,
  "ISPY2_GSE196096_annotation_ALL_gene_expression_all_annotation.csv"
)

fwrite(
  merged_expr_only,
  "ISPY2_GSE196096_annotation_ALL_gene_expression_matched_only.csv"
)

## gene list
gene_cols <- setdiff(names(expr_mat), "ResID")

fwrite(
  data.table(Gene = gene_cols),
  "ISPY2_all_gene_list.csv"
)

#####Differential Genes from non PCR vs PCR in GSE194040######
############################################################
## ISPY2 – VC volcano plot
############################################################

rm(list = ls())
gc()

library(data.table)
library(limma)
library(ggplot2)
library(ggrepel)

setwd("~/Desktop/ISPY2")

############################################################
## 1. Load data
############################################################

df <- fread("ISPY2_GSE196096_annotation_ALL_gene_expression_matched_only.csv")

df[, hr := as.integer(hr)]
df[, her2 := as.integer(her2)]
df[, pcr := as.integer(pcr)]

############################################################
## 2. Volcano function
############################################################

run_volcano <- function(df_sub, group_name, prefix) {
  
  df_sub <- df_sub[!is.na(pcr) & pcr %in% c(0, 1)]
  
  cat("\n====", group_name, "====\n")
  print(table(df_sub$pcr))
  
  clinical_cols <- c(
    "ResID", "GSM", "title", "arm", "arm (short name)",
    "her2", "hr", "mp", "pcr", "tissue",
    "adjustment method", "rppa array (normalization key)"
  )
  
  gene_cols <- setdiff(names(df_sub), clinical_cols)
  gene_cols <- gene_cols[sapply(df_sub[, ..gene_cols], is.numeric)]
  
  expr_patient_gene <- as.matrix(df_sub[, ..gene_cols])
  rownames(expr_patient_gene) <- df_sub$ResID
  
  expr <- t(expr_patient_gene)   # genes × samples
  
  group <- factor(df_sub$pcr, levels = c(0, 1))
  design <- model.matrix(~ group)
  
  fit <- lmFit(expr, design)
  fit <- eBayes(fit)
  
  res <- topTable(fit, coef = 2, number = Inf)
  res$logFC <- -res$logFC
  res$gene <- rownames(res)
  
  res$negLog10P <- -log10(res$P.Value)
  
  res$direction <- "NS"
  res$direction[res$P.Value < 0.05 & res$logFC > 0] <- "Higher in non-pCR"
  res$direction[res$P.Value < 0.05 & res$logFC < 0] <- "Higher in pCR"
  
  cat("\nENPP1:\n")
  print(res[res$gene == "ENPP1", ])
  
  top_n <- 15
  label_genes <- unique(c(
    "ENPP1",
    res[order(res$P.Value), "gene"][1:top_n]
  ))
  
  res$label <- ifelse(res$gene %in% label_genes, res$gene, NA)
  
  p <- ggplot(res, aes(x = logFC, y = negLog10P)) +
    geom_point(aes(color = direction), alpha = 0.75, size = 1.5) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_text_repel(
      aes(label = label),
      max.overlaps = Inf,
      size = 3
    ) +
    scale_color_manual(values = c(
      "Higher in non-pCR" = "#F8766D",
      "Higher in pCR" = "#00BFC4",
      "NS" = "grey80"
    )) +
    labs(
      title = paste0(group_name, ": non-pCR vs pCR"),
      subtitle = "logFC > 0 indicates higher expression in non-pCR responders",
      x = "logFC (non-pCR vs pCR)",
      y = "-log10(raw p value)",
      color = ""
    ) +
    theme_classic(base_size = 13)
  
  ggsave(
    paste0(prefix, "_volcano_nonpCR_vs_pCR.pdf"),
    p,
    width = 7,
    height = 6
  )
  
  ggsave(
    paste0(prefix, "_volcano_nonpCR_vs_pCR.png"),
    p,
    width = 7,
    height = 6,
    dpi = 300
  )
  
  write.csv(
    res,
    paste0(prefix, "_limma_results_nonpCR_vs_pCR.csv"),
    row.names = FALSE
  )
  
  return(res)
}

############################################################
## 3. VC arm
############################################################

df_vc <- df[
  grepl("veliparib|carboplatin|VC", arm, ignore.case = TRUE) |
    grepl("veliparib|carboplatin|VC", `arm (short name)`, ignore.case = TRUE)
]

res_vc <- run_volcano(
  df_sub = df_vc,
  group_name = "VC arm",
  prefix = "VC"
)
############################################################
## 4. ENPP1 comparison
############################################################

cat("\n==== ENPP1 comparison ====\n")

print(res_vc[res_vc$gene == "ENPP1", ])
enpp1_summary <- rbind(
  cbind(Arm = "VC", res_vc[res_vc$gene == "ENPP1", ]),
)

write.csv(
  enpp1_summary,
  "ENPP1_VC_volcano_summary.csv",
  row.names = FALSE
)

#####Differential Genes from non PCR vs PCR in GSE173839######
############################################################
## GSE173839 – durva + olaparib arm
## non-pCR vs pCR (same as VC)
############################################################

library(data.table)
library(dplyr)
library(limma)
library(ggplot2)
library(ggrepel)

setwd("~/Desktop/GSE173839")

############################################################
## 1.  read expression
############################################################

expr <- fread("GSE173839_ISPY2_AgilentGeneExp_durvaPlusCtr_FFPE_meanCol_geneLevel_n105.txt")

head(expr[,1:5])
colnames(expr)[1] <- "gene"

expr_mat <- as.data.frame(expr)
rownames(expr_mat) <- expr_mat$gene
expr_mat <- expr_mat[, -1]

############################################################
## 2. read clinical
############################################################

clin <- fread("GSE173839_ISPY2_DurvalumabOlaparibArm_biomarkers.csv.gz")

head(clin)
colnames(clin)

############################################################
## 3. pCR
############################################################
clin$pCR <- as.numeric(clin$pCR)
clin$group <- ifelse(clin$pCR == 1, "pCR", "non_pCR")
clin$group <- factor(clin$group, levels = c("pCR", "non_pCR"))

############################################################
## 4. Match samples（GSE173839)
############################################################
ncol(expr_mat)
nrow(clin)
min(ncol(expr_mat), nrow(clin))

## control
expr_use <- expr_control
clin_use <- clin_control

## or DOP
expr_use <- expr_dop
clin_use <- clin_dop
############################################################
## 4.5 Split arms
############################################################

## control arm
clin_control <- clin[clin$Arm == "control", ]
expr_control <- expr_mat[, clin$Arm == "control"]

## durvalumab/olaparib arm
clin_dop <- clin[clin$Arm == "durvalumab/olaparib", ]
expr_dop <- expr_mat[, clin$Arm == "durvalumab/olaparib"]

table(clin_control$group)
table(clin_dop$group)

############################################################
## FIX expression scale
############################################################
expr_use <- as.matrix(expr_use)
mode(expr_use) <- "numeric"
if (max(expr_use, na.rm = TRUE) > 50) {
  expr_use <- log2(expr_use + 1)
}
summary(expr_use)


############################################################
## 5. limma
############################################################

design <- model.matrix(~0 + clin$group)
colnames(design) <- levels(clin$group)

fit <- lmFit(expr_use, design)

contrast <- makeContrasts(
  nonpCR_vs_pCR = `non_pCR` - pCR,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)

res <- topTable(fit2,
                coef = "nonpCR_vs_pCR",
                number = Inf,
                sort.by = "P")

res$gene <- rownames(res)

############################################################
## 6. Volcano
############################################################

res$logP <- -log10(res$P.Value)

res$group <- "NS"
res$group[res$P.Value < 0.05 & res$logFC > 0] <- "Higher in non-pCR"
res$group[res$P.Value < 0.05 & res$logFC < 0] <- "Higher in pCR"

p <- ggplot(res, aes(x = logFC, y = logP, color = group)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c(
    "Higher in non-pCR" = "#F8766D",
    "Higher in pCR" = "#00BFC4",
    "NS" = "grey"
  )) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_classic() +
  labs(
    title = "Durvalumab + Olaparib arm",
    x = "logFC (non-pCR vs pCR)",
    y = "-log10(p value)"
  )
############################################################
## 7. Save
############################################################

write.csv(res, "DurvaOlaparib_limma_results_nonpCR_vs_pCR.csv", row.names = FALSE)

ggsave("DurvaOlaparib_volcano.pdf", width = 6, height = 5)
