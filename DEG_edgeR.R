library(edgeR)
library(dplyr)
library(ggplot2)
library(stringr)
library(readr)
library(tidyr)
library(ggrepel)
library(statmod)
library(purrr)

BioCC_input = read.table(file = "/data2/lackey_lab/DownloadedSequenceData/randazza/nerve_root/rrna_removal/more_stringent/counts_NR_rRNAremove.txt", sep = "\t", header = T)
#change row names to gene id
rownames(BioCC_input) <- BioCC_input$Geneid
#delete unnecessary columns
extra <- BioCC_input[,1:6]
BioCC_input[, 1:6] <- NULL
colnames(BioCC_input) <- c("N6", "N7", "N8","P10", "P11", "C1", "C2", "C3", "C4",
                           "P4b", "P4a", "P5", "P6", "P8", "P9")

#assign samples to groups
meta <- data.frame(sample = c("N6", "N7", "N8", "NIH", "P10", "P11", "C1", "C2", "C3", "C4",
                              "P4b", "P4a", "P5", "P6", "P7", "P8", "P9"),
                   group = c("patient", "patient", "patient", "control", "patient", "patient", "control",
                             "control", "control", "control", "patient", "patient", "patient",
                             "patient", "patient", "patient", "patient"),
                   sex = c("F", "F", "F", "M", "M", "F", "M", "M", "F", "M",
                           "F", "F", "F", "M", "F", "F", "F"),
                   age = c("young", "young", "young", "young", "young", "young", "retired", "retired", "retired", "retired",
                           "young", "young", "young", "young", "young", "retired", "young")
)
meta <- meta[match(
  colnames(BioCC_input),
  meta$sample
), ]
meta$group <- factor(meta$group, levels = c("control","patient"))
design <- model.matrix(~ sex + age + group, data = meta)

rpk.norm.g <- DGEList(counts = BioCC_input)
cpm_data <- cpm(rpk.norm.g)

keep_f <- rowSums(cpm_data >= 4) >= 4
sum(keep_f) 

rpk.norm.g <- rpk.norm.g[keep, , keep.lib.sizes= FALSE]

###TPM calculations
counts <- rpk.norm.g$counts
gene_lengths <- extra[, c("Geneid", "Length")]
gene_lengths$KB <- gene_lengths$Length / 1000

gene_len_vec <- setNames(
  gene_lengths$KB,
  gene_lengths$Geneid
)
common_genes <- intersect(
  rownames(counts),
  names(gene_len_vec)
)

length(common_genes) 
counts <- counts[common_genes, , drop = FALSE]
gene_len_aligned <- gene_len_vec[rownames(counts)]
# Calculate RPK
rpk <- sweep(
  counts,
  1,
  gene_len_aligned,
  FUN = "/"
)

# Calculate TPM
tpm <- sweep(
  rpk,
  2,
  colSums(rpk, na.rm = TRUE) / 1e6,
  FUN = "/"
)

colSums(tpm) 
tpm_data <- as.data.frame(tpm)
tpm_data$Genes <- rownames(tpm_data)
write.table(tpm_data, "/data2/lackey_lab/randazza/EDS/short-read/rrna_removed_7to12/DGE/more_stringent/minus_outliers/all_TPM_keep_default.txt", sep = "\t")

rpk.norm.g <- calcNormFactors(rpk.norm.g, method = "TMM")
logCPM <- cpm(rpk.norm.g, log=TRUE, prior.count=1)
plotMDS(logCPM, labels=colnames(logCPM))

rpk.norm.g <- estimateDisp(rpk.norm.g,design, robust = TRUE)

mean(rpk.norm.g$tagwise.dispersion)
plotBCV(rpk.norm.g)
fit <- glmQLFit(rpk.norm.g, design, robust = TRUE)
plotQLDisp(fit)

test <- glmQLFTest(fit, coef = "grouppatient")
ddx_P1 <- topTags(test, n=Inf)
summary(dt_K700E_lrt<-decideTestsDGE(test,p.value = 0.05, adjust.method = "BH"))

#create ranks for GSEApreRanked analysis
tab <- topTags(test, n = Inf)$table
rnk <- sign(tab$logFC) * sqrt(pmax(tab$F, 0))
names(rnk) <- rownames(tab)
rnk <- sort(rnk, decreasing = TRUE)
rnk.df <- data.frame(
  Gene = names(rnk),
  Rank = as.numeric(rnk)
)
rnk.df$Gene <- sub("\\..*", "", rnk.df$Gene)
write.table(
  rnk.df,
  "minus_outliers/comboco_keepdefault_prerank_genename_Fstat.rnk",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

new <- as.data.frame(ddx_P1$table)
new$Gene.Symbol <- rownames(new)
rownames(new) <- NULL
new$negLogFDR <- -log10(new$FDR)
ggplot(new, aes(x = logFC, y = negLogFDR)) +
  geom_point(size = 0.8, alpha = 0.5) +
  theme_classic()+
  xlab("logFC") +
  ylab("-log10(FDR)")

new$Sig <- "Not Sig"
new$Sig[new$FDR < 0.05 & abs(new$logFC) < 1] <- "Sig, low FC"
new$Sig[new$FDR < 0.05 & abs(new$logFC) >= 1] <- "Sig, high FC"
gene_names <- read.delim("/data2/lackey_lab/GenomeReferences/hisat_index/gene_names.tab", sep = "\t", header = FALSE)
new <- merge(new, gene_names, by.x = "Gene.Symbol", by.y = "V1")

write.table(new, "/data2/lackey_lab/randazza/EDS/short-read/rrna_removed_7to12/DGE/more_stringent/minus_outliers/nerve_root_DGE_list_keepdef_comboco.txt", sep = "\t")

ggplot(new, aes(x = logFC, y = negLogFDR, color = Sig)) +
  geom_point(alpha = 0.4, size = 1) +
  theme_classic() +
  scale_color_manual(values = c("Not Sig" = "#cccccc", "Sig, low FC" = "#cccccc", "Sig, high FC" = "purple"))
