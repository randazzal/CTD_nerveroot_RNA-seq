# Hallmark GSEA for 24 pooled patients versus 3 Mount Sinai nerve-root controls.
# Input: gene-level results from the final dataset-, age-, and sex-adjusted
# pooled edgeR model. All tested genes are eligible for ranking.

suppressPackageStartupMessages(library(fgsea))

# Replace these paths before running.
DE_FILE  <- "path/to/pooled27_strict_edgeR.tsv"
MAP_FILE <- "path/to/hg38.ensg_to_symbol_stripped.tsv"
GMT_FILE <- "path/to/h.all.v2024.1.Hs.symbols.gmt"
OUT      <- "path/to/Hallmark_GSEA_final_adjusted"

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# 1. Read the pooled edgeR results and map Ensembl IDs to gene symbols.
de <- read.delim(
  DE_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

mapping <- read.table(
  MAP_FILE,
  header = FALSE,
  sep = "",
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)

colnames(mapping) <- c("ENSG", "symbol")

de$ENSG <- sub("\\.[0-9]+$", "", de$gene_id)
de$symbol <- mapping$symbol[match(de$ENSG, mapping$ENSG)]

# 2. Rank genes using the direction of logFC and the edgeR F statistic.
de$rank_stat <- sign(de$logFC) * sqrt(pmax(de$F, 0))

valid <- (
  !is.na(de$symbol) &
  nzchar(de$symbol) &
  is.finite(de$rank_stat)
)

rank_df <- de[valid, ]

# For duplicate symbols, keep the strongest absolute ranking statistic.
# Ensembl ID breaks ties.
rank_df <- rank_df[order(-abs(rank_df$rank_stat), rank_df$ENSG), ]
rank_df <- rank_df[!duplicated(rank_df$symbol), ]

ranks <- rank_df$rank_stat
names(ranks) <- rank_df$symbol
ranks <- sort(ranks, decreasing = TRUE)

# 3. Run preranked GSEA using the human Hallmark v2024.1 gene sets.
pathways <- fgsea::gmtPathways(GMT_FILE)

set.seed(42)
gsea <- fgsea::fgseaMultilevel(
  pathways = pathways,
  stats = ranks,
  minSize = 15,
  maxSize = 500,
  eps = 1e-10
)

# 4. Save pathway results and the gene ranking used for GSEA.
gsea <- as.data.frame(gsea)
gsea$Enriched_group <- ifelse(
  gsea$NES > 0,
  "Pooled patients",
  "Mt Sinai NR controls"
)

gsea$leadingEdge <- vapply(
  gsea$leadingEdge,
  paste,
  collapse = ";",
  FUN.VALUE = character(1)
)

gsea <- gsea[order(gsea$padj, -abs(gsea$NES)), ]

write.table(
  gsea,
  file.path(OUT, "Hallmark_GSEA_all.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

gene_ranking <- data.frame(
  symbol = names(ranks),
  rank_stat = unname(ranks)
)

write.table(
  gene_ranking,
  file.path(OUT, "Hallmark_gene_ranking.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# 5. Summarize the pooled comparison.
sig <- !is.na(gsea$padj) & gsea$padj < 0.05

cat("Genes tested:", nrow(de), "\n")
cat("Genes with symbols:", sum(!is.na(de$symbol) & nzchar(de$symbol)), "\n")
cat("Unique ranked genes:", length(ranks), "\n")
cat("Hallmark gene sets loaded:", length(pathways), "\n")
cat("Pathways tested:", nrow(gsea), "\n")
cat("BH FDR < 0.05:", sum(sig), "\n")
cat("Enriched in pooled patients:", sum(sig & gsea$NES > 0), "\n")
cat("Enriched in Mt Sinai NR controls:", sum(sig & gsea$NES < 0), "\n")
