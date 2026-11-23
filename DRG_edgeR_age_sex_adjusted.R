# UT Dallas DRG: chronic versus acute differential expression
#
# Input directory must contain:
#   subject_aggregated_raw_counts.tsv  - gene rows, subject columns
#   subject_metadata.tsv               - subject_id, condition, sex, age
#
# Run with:
#   Rscript 52_DRG_edgeR_age_sex_adjusted.R

suppressPackageStartupMessages(library(edgeR))

# Replace these two paths before running.
INPUT <- "path/to/input_directory"
OUT <- "path/to/output_directory"

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# 1. Load subject-level raw counts and metadata
counts <- read.delim(
  file.path(INPUT, "subject_aggregated_raw_counts.tsv"),
  row.names = 1,
  check.names = FALSE
)

counts <- as.matrix(counts)
storage.mode(counts) <- "numeric"

md <- read.delim(
  file.path(INPUT, "subject_metadata.tsv"),
  stringsAsFactors = FALSE
)

# Align metadata rows to the subject columns in the count matrix.
md <- md[match(colnames(counts), md$subject_id), ]

stopifnot(
  ncol(counts) == 13,
  !anyNA(md$subject_id),
  identical(colnames(counts), md$subject_id),
  !anyNA(counts),
  all(counts >= 0)
)

# 2. Prepare model variables
# Acute is the reference level for both comparisons.
md$condition <- factor(tolower(md$condition), levels = c("acute", "chronic"))
md$sex <- factor(tolower(md$sex))
md$age <- as.numeric(md$age)

stopifnot(
  !anyNA(md$condition),
  !anyNA(md$sex),
  !anyNA(md$age)
)

md$age_scaled <- as.numeric(scale(md$age))

cat("Subjects:", nrow(md), "\n")
print(table(md$condition))
cat("\nAge and sex by condition:\n")
print(md[, c("subject_id", "condition", "sex", "age")])

# 3. Construct and check both designs
design_unadjusted <- model.matrix(~ condition, data = md)
design_adjusted <- model.matrix(~ age_scaled + sex + condition, data = md)

cat(
  "\nUnadjusted design rank:", qr(design_unadjusted)$rank,
  "of", ncol(design_unadjusted), "\n"
)
cat(
  "Adjusted design rank:", qr(design_adjusted)$rank,
  "of", ncol(design_adjusted), "\n"
)

stopifnot(
  qr(design_unadjusted)$rank == ncol(design_unadjusted),
  qr(design_adjusted)$rank == ncol(design_adjusted)
)

# 4. Filter once using the unadjusted design and calculate TMM factors
dge <- DGEList(counts = counts, group = md$condition)

keep <- filterByExpr(dge, design = design_unadjusted)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")

cat("\nGenes retained:", nrow(dge), "\n")

logCPM <- cpm(dge, log = TRUE, prior.count = 2)

write.table(
  logCPM,
  file.path(OUT, "TMM_logCPM.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

# 5. Fit each model and test chronic versus acute
run_edger <- function(design, model_name) {
  y <- estimateDisp(dge, design, robust = TRUE)
  fit <- glmQLFit(y, design, robust = TRUE)
  test <- glmQLFTest(fit, coef = "conditionchronic")

  res <- topTags(test, n = Inf)$table
  res$gene_id <- rownames(res)
  res <- res[, c("gene_id", "logFC", "logCPM", "F", "PValue", "FDR")]

  write.table(
    res,
    file.path(OUT, paste0(model_name, "_edgeR.tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  cat("\nModel:", model_name, "\n")
  cat("Tested genes:", nrow(res), "\n")
  cat("FDR < 0.05:", sum(res$FDR < 0.05), "\n")
  cat("Higher in chronic:", sum(res$FDR < 0.05 & res$logFC > 0), "\n")
  cat("Higher in acute:", sum(res$FDR < 0.05 & res$logFC < 0), "\n")
}

run_edger(design_unadjusted, "unadjusted")
run_edger(design_adjusted, "age_sex_adjusted")

# 6. Save the metadata used in the models
write.table(
  md,
  file.path(OUT, "subject_metadata.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\nAnalysis completed.\nResults saved to:\n", OUT, "\n")
