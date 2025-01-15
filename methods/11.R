library(tximport)
library(stringr)
library(dplyr)
library(edgeR)
library(ggplot2)
library(ggVennDiagram)

# create tx2gene file

gff <- read.delim("../material/reference.gff", header = FALSE, skip = 21)

gff[length(gff), ]

gff_genes <- gff[gff[[3]] == "gene", 9]
gff_mrna <- gff[gff[[3]] == "mRNA", 9]

get_id <- function(txt) {
  str_match(txt, "ID=(.*);Name.*")[, 2]
}

transcript_ids <- as.data.frame(get_id(gff_mrna))
gene_ids <- as.data.frame(get_id(gff_genes))

transcript_ids$gene_id <- str_match(transcript_ids[, 1], "(.*)_.*")[, 2]

names(gene_ids)[[1]] <- "gene_id"

tx2gene <- merge(transcript_ids, gene_ids, by = "gene_id", all = TRUE)
names(tx2gene)[[2]] <- "transcript_id"
tx2gene <- tx2gene[, c(2, 1)]

head(tx2gene)

# read abundance files

control_files <- list.files(path = "../results/control", pattern = "tsv$", recursive = TRUE, full.names = TRUE)
condition1_files <- list.files(path = "../results/condition1", pattern = "tsv$", recursive = TRUE, full.names = TRUE)
condition2_files <- list.files(path = "../results/condition2", pattern = "tsv$", recursive = TRUE, full.names = TRUE)

txi_control <- tximport(control_files,
  type = "kallisto", tx2gene = tx2gene,
  ignoreTxVersion = TRUE
)
txi_condition1 <- tximport(condition1_files,
  type = "kallisto", tx2gene = tx2gene,
  ignoreTxVersion = TRUE
)
txi_condition2 <- tximport(condition2_files,
  type = "kallisto", tx2gene = tx2gene,
  ignoreTxVersion = TRUE
)

txi_all <- list(
  control = txi_control,
  condition1 = txi_condition1,
  condition2 = txi_condition2
)

counts <- cbind(
  txi_control$counts,
  txi_condition1$counts,
  txi_condition2$counts
)

head(counts)

# perform differential gene expression analysis

group <- factor(c(rep("control", 3), rep("condition1", 3), rep("condition2", 3)))
group <- factor(group, levels = c("control", "condition1", "condition2"))

y <- DGEList(counts = counts, group = group)

keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes = FALSE]

y <- normLibSizes(y)

design <- model.matrix(~ 0 + group)
fit <- glmQLFit(y, design)

test_control_vs_cond1 <- glmQLFTest(fit, contrast = c(-1, 1, 0))
results_control_vs_cond1 <- topTags(test_control_vs_cond1, n = Inf)
write.csv(results_control_vs_cond1, "../results/control_vs_condition1_results.csv")

res_cond1 <- results_control_vs_cond1$table
res_cond1$sig <- ifelse(
  res_cond1$FDR < 0.05,
  "Significant",
  "Not Significant"
)

test_control_vs_cond2 <- glmQLFTest(fit, contrast = c(-1, 0, 1))
results_control_vs_cond2 <- topTags(test_control_vs_cond2, n = Inf)
write.csv(results_control_vs_cond2, "../results/control_vs_condition2_results.csv")

res_cond2 <- results_control_vs_cond2$table
res_cond2$sig <- ifelse(
  res_cond2$FDR < 0.05,
  "Significant",
  "Not Significant"
)

summ_table <- data.frame(
  gene_id = sort(row.names(res_cond1))
)

merge(tx2gene, read.delim(control_files[[1]]), by.x = "gene_id", by.y = "target_id", all = TRUE)

head(summ_table)

volcano_con1 <- ggplot(data = res_cond1, aes(x = logFC, y = -log10(PValue), color = sig)) +
  geom_point() +
  scale_color_manual(values = c("Not Significant" = "grey", "Significant" = "red")) +
  theme_minimal() +
  labs(
    x = "logFC",
    y = "-log10(PValue)",
    color = "Significance"
  )

volcano_con2 <- ggplot(data = res_cond2, aes(x = logFC, y = -log10(PValue), color = sig)) +
  geom_point() +
  scale_color_manual(values = c("Not Significant" = "grey", "Significant" = "red")) +
  theme_minimal() +
  labs(
    x = "logFC",
    y = expression(-log[10](PValue)),
    color = "Significance"
  )

ggsave("volcano_con1.pdf", plot = volcano_con1, device = "pdf", path = "../results")
ggsave("volcano_con2.pdf", plot = volcano_con2, device = "pdf", path = "../results")

# venn diagram

lfc_ge_2_cn1 <- row.names(res_cond1[abs(res_cond1$logFC) >= 2, ])
fdr_ad_p_le_0.05_cn1 <- row.names(res_cond1[p.adjust(res_cond1$FDR, "BH") <= 0.05, ])
fdr_le_0.01_cn1 <- row.names(res_cond1[res_cond1$FDR <= 0.01, ])

deg_list_cn1 <- list(
  LFC_greater_2 = lfc_ge_2_cn1,
  FDR_ad_p_value = fdr_ad_p_le_0.05_cn1,
  FDR_less_than_0.01 = fdr_le_0.01_cn1
)

ggVennDiagram(deg_list_cn1) +
  scale_fill_gradient(low = "white", high = "blue") +
  theme_minimal()

lfc_ge_2_cn2 <- row.names(res_cond1[abs(res_cond2$logFC) >= 2, ])
fdr_ad_p_le_0.05_cn2 <- row.names(res_cond2[p.adjust(res_cond1$FDR, "BH") <= 0.05, ])
fdr_le_0.01_cn2 <- row.names(res_cond2[res_cond1$FDR <= 0.01, ])

deg_list_cn2 <- list(
  LFC_greater_2 = lfc_ge_2_cn2,
  FDR_ad_p_value = fdr_ad_p_le_0.05_cn2,
  FDR_less_than_0.01 = fdr_le_0.01_cn2
)

ggVennDiagram(deg_list_cn2) +
  scale_fill_gradient(low = "white", high = "blue") +
  theme_minimal()

res_cond1_inc <- row.names(res_cond1[res_cond1$logFC > 0, ])
res_cond1_dec <- row.names(res_cond1[res_cond1$logFC < 0, ])

res_cond2_inc <- row.names(res_cond1[res_cond2$logFC > 0, ])
res_cond2_dec <- row.names(res_cond1[res_cond2$logFC < 0, ])

# top 20 differentially expressed geens
head(res_cond1[order(abs(res_cond1$logFC), decreasing = TRUE), ], 20)
head(res_cond2[order(abs(res_cond2$logFC), decreasing = TRUE), ], 20)

head(res_cond1)

# sdeg <- data.frame(
#   gene_id <-
# )
