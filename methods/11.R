library(tidyverse)
library(tximport)
library(stringr)
library(edgeR)
library(ggplot2)
library(ggVennDiagram)
library(limma)
library(wordcloud)

# map transcripts to genes

control_files <- list.files(path = "../material/control", pattern = "tsv$", recursive = TRUE, full.names = TRUE)
condition1_files <- list.files(path = "../material/condition1", pattern = "tsv$", recursive = TRUE, full.names = TRUE)
condition2_files <- list.files(path = "../material/condition2", pattern = "tsv$", recursive = TRUE, full.names = TRUE)

tx2gene <- read.table("../results/tx2gene", header = FALSE)

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
write.csv(results_control_vs_cond1, "../results/11_edgeR_control_vs_condition1_results.csv")

test_control_vs_cond2 <- glmQLFTest(fit, contrast = c(-1, 0, 1))
results_control_vs_cond2 <- topTags(test_control_vs_cond2, n = Inf)
write.csv(results_control_vs_cond2, "../results/11_edgeR_control_vs_condition2_results.csv")

# results

## volcano plots

res_cond1 <- results_control_vs_cond1$table
res_cond1$sig <- ifelse(
  res_cond1$FDR < 0.05,
  "Significant",
  "Not Significant"
)

res_cond2 <- results_control_vs_cond2$table
res_cond2$sig <- ifelse(
  res_cond2$FDR < 0.05,
  "Significant",
  "Not Significant"
)

edgeR_volcano_con1 <- ggplot(data = res_cond1, aes(x = logFC, y = -log10(PValue), color = sig)) +
  geom_point() +
  scale_color_manual(values = c("Not Significant" = "grey", "Significant" = "red")) +
  theme_minimal() +
  labs(
    x = "logFC",
    y = "-log(PValue)",
    color = "Significance"
  )
ggsave("11_edgeR_volcano_con1.pdf", plot = edgeR_volcano_con1, device = "pdf", path = "../results")

edgeR_volcano_con2 <- ggplot(data = res_cond2, aes(x = logFC, y = -log10(PValue), color = sig)) +
  geom_point() +
  scale_color_manual(values = c("Not Significant" = "grey", "Significant" = "red")) +
  theme_minimal() +
  labs(
    x = "logFC",
    y = "-log(PValue)",
    color = "Significance"
  )
ggsave("11_edgeR_volcano_con2.pdf", plot = edgeR_volcano_con2, device = "pdf", path = "../results")

# venn diagrams

lfc_ge_2_cn1 <- row.names(res_cond1[abs(res_cond1$logFC) >= 2, ])
fdr_le_0.05_cn1 <- row.names(res_cond1[res_cond1$FDR <= 0.05, ])
fdr_le_0.01_cn1 <- row.names(res_cond1[res_cond1$FDR <= 0.01, ])

deg_list_cn1 <- list(
  LFC_gr_2 = lfc_ge_2_cn1,
  FDR_le_0.05 = fdr_le_0.05_cn1,
  FDR_le_0.01 = fdr_le_0.01_cn1
)

edgeR_venn_diagram_cn1 <- ggVennDiagram(deg_list_cn1) +
  scale_fill_gradient(low = "white", high = "blue") +
  theme_minimal()
ggsave("../results/11_edgeR_venn_diagram_cn1.pdf", plot = edgeR_venn_diagram_cn1, device = "pdf")

lfc_ge_2_cn2 <- row.names(res_cond1[abs(res_cond2$logFC) >= 2, ])
fdr_le_0.05_cn2 <- row.names(res_cond2[res_cond1$FDR <= 0.05, ])
fdr_le_0.01_cn2 <- row.names(res_cond2[res_cond1$FDR <= 0.01, ])

deg_list_cn2 <- list(
  LFC_gr_2 = lfc_ge_2_cn2,
  FDR_le_0.05 = fdr_le_0.05_cn2,
  FDR_le_0.01 = fdr_le_0.01_cn2
)

edgeR_venn_diagram_cn2 <- ggVennDiagram(deg_list_cn2) +
  scale_fill_gradient(low = "white", high = "blue") +
  theme_minimal()
ggsave("../results/11_edgeR_venn_diagram_cn2.pdf", plot = edgeR_venn_diagram_cn2, device = "pdf")

# increased and decreased genes

edgeR_res_cond1_inc <- row.names(res_cond1[res_cond1$logFC > 0, ])
edgeR_res_cond1_dec <- row.names(res_cond1[res_cond1$logFC < 0, ])

edgeR_res_cond2_inc <- row.names(res_cond1[res_cond2$logFC > 0, ])
edgeR_res_cond2_dec <- row.names(res_cond1[res_cond2$logFC < 0, ])

# top 20 differentially expressed geens
top_up_cn1 <- head(res_cond1[order(res_cond1$logFC, decreasing = TRUE), ], 20)
top_up_cn2 <- head(res_cond2[order(res_cond2$logFC, decreasing = TRUE), ], 20)
top_down_cn1 <- head(res_cond1[order(res_cond1$logFC, decreasing = FALSE), ], 20)
top_down_cn2 <- head(res_cond2[order(res_cond2$logFC, decreasing = FALSE), ], 20)

top_up_cn1_edgeR <- row.names(top_up_cn1)
top_up_cn2_edgeR <- row.names(top_up_cn2)
top_down_cn1_edgeR <- row.names(top_down_cn1)
top_down_cn2_edgeR <- row.names(top_down_cn2)

go_annotations <- read_delim("../material/go_annotations", skip = 4, col_names = FALSE)

go_annotations$gene_id <- ifelse(
  is.na(str_match(go_annotations[[11]], ".*\\|(.*)")[, 2]),
  go_annotations[[11]],
  str_match(go_annotations[[11]], ".*\\|(.*)")[, 2]
)

names(go_annotations)[[3]] <- "protein_id"
names(go_annotations)[[10]] <- "protein_description"

top_up_cn1_desc <- subset(go_annotations, gene_id %in% rownames(top_up_cn1))
top_up_cn2_desc <- subset(go_annotations, gene_id %in% rownames(top_up_cn2))
top_down_cn1_desc <- subset(go_annotations, gene_id %in% rownames(top_down_cn1))
top_down_cn2_desc <- subset(go_annotations, gene_id %in% rownames(top_down_cn2))

top_up_cn1_list <- str_split(top_up_cn1_desc$protein_description, " ") %>%
  unlist() %>%
  tolower()
top_up_cn2_list <- str_split(top_up_cn2_desc$protein_description, " ") %>%
  unlist() %>%
  tolower()
top_down_cn1_list <- str_split(top_down_cn1_desc$protein_description, " ") %>%
  unlist() %>%
  tolower()
top_down_cn2_list <- str_split(top_down_cn2_desc$protein_description, " ") %>%
  unlist() %>%
  tolower()

bad_words_list <- c(
  "protein",
  "subunit",
  "complex",
  "factor",
  "and",
  "putative",
  "iii",
  "rna",
  "dna",
  "uncharacterized",
  "cell",
  "homolog",
  "with",
  "atp-dependent",
  "yel009c-a",
  "probable"
)

top_up_cn1_list <- top_up_cn1_list[!top_up_cn1_list %in% bad_words_list & nchar(top_up_cn1_list) > 2]
top_up_cn2_list <- top_up_cn2_list[!top_up_cn2_list %in% bad_words_list & nchar(top_up_cn2_list) > 2]
top_down_cn1_list <- top_down_cn1_list[!top_down_cn1_list %in% bad_words_list & nchar(top_down_cn1_list) > 2]
top_down_cn2_list <- top_down_cn2_list[!top_down_cn2_list %in% bad_words_list & nchar(top_down_cn2_list) > 2]


top_up_cn1_list_df <- as.data.frame(table(top_up_cn1_list))
top_up_cn2_list_df <- as.data.frame(table(top_up_cn2_list))
top_down_cn1_list_df <- as.data.frame(table(top_down_cn1_list))
top_down_cn2_list_df <- as.data.frame(table(top_down_cn2_list))




pdf("../results/11_edgeR_wordcloud_cn1_up.pdf")
wordcloud(
  words = top_up_cn1_list_df$top_up_cn1_list, freq = top_up_cn1_list_df$Freq, min.freq = 1,
  max.words = 200, random.order = FALSE, rot.per = 0.35,
  colors = brewer.pal(8, "Set1")
)
dev.off()

pdf("../results/11_edgeR_wordcloud_cn2_up.pdf")
wordcloud(
  words = top_up_cn2_list_df$top_up_cn2_list, freq = top_up_cn2_list_df$Freq, min.freq = 1,
  max.words = 200, random.order = FALSE, rot.per = 0.35,
  colors = brewer.pal(8, "Set1")
)
dev.off()

pdf("../results/11_edgeR_wordcloud_cn1_down.pdf")
wordcloud(
  words = top_down_cn1_list_df$top_down_cn1_list, freq = top_down_cn1_list_df$Freq, min.freq = 1,
  max.words = 200, random.order = FALSE, rot.per = 0.35,
  colors = brewer.pal(8, "Set1")
)
dev.off()

pdf("../results/11_edgeR_wordcloud_cn2_down.pdf")
wordcloud(
  words = top_down_cn2_list_df$top_down_cn2_list, freq = top_down_cn2_list_df$Freq, min.freq = 1,
  max.words = 200, random.order = FALSE, rot.per = 0.35,
  colors = brewer.pal(8, "Set1")
)
dev.off()

# table

edgeR_table <- data.frame(
  gene_id = sort(row.names(res_cond1))
)

create_tpm_column_and_merge <- function(files, condition, tbl = edgeR_table) {
  for (i in seq_along(files)) {
    df <- read.delim(files[[i]])[, c("target_id", "tpm")]
    df$target_id <- ifelse(
      is.na(str_match(df$target_id, "(.*)_.*")[, 2]),
      df$target_id,
      str_match(df$target_id, "(.*)_.*")[, 2]
    )
    names(df)[[1]] <- "gene_id"
    names(df)[[2]] <- paste0("TPM-Repl-", i, "-", condition)
    tbl <- merge(tbl, df, by = "gene_id")
  }
  return(tbl)
}

edgeR_table <- create_tpm_column_and_merge(control_files, "Control")
edgeR_table <- create_tpm_column_and_merge(condition1_files, "Condition1")
edgeR_table <- create_tpm_column_and_merge(condition2_files, "Condition2")

res_cond1$gene_id <- row.names(res_cond1)
edgeR_table <- merge(edgeR_table, res_cond1[, c("gene_id", "logFC")], by = "gene_id")
names(edgeR_table)[[11]] <- "logFC-Control-vs-Condition1"

res_cond2$gene_id <- row.names(res_cond2)
edgeR_table <- merge(edgeR_table, res_cond2[, c("gene_id", "logFC")], by = "gene_id")
names(edgeR_table)[[12]] <- "logFC-Control-vs-Condition2"

edgeR_table <- merge(edgeR_table, res_cond1[, c("gene_id", "FDR")], by = "gene_id")
names(edgeR_table)[[13]] <- "Adjusted p-value-Control-vs-Condition1"

edgeR_table <- merge(edgeR_table, res_cond2[, c("gene_id", "FDR")], by = "gene_id")
names(edgeR_table)[[14]] <- "Adjusted p-value-Control-vs-Condition2"

# head(edgeR_table)

# grouped barchart

top_genes_table <- results_control_vs_cond1$table
top_positive <- top_genes_table[top_genes_table$logFC > 0, ]
top_negative <- top_genes_table[top_genes_table$logFC < 0, ]
top_positive_sorted <- top_positive[order(top_positive$logFC,
  decreasing = TRUE
), ]
top_negative_sorted <- top_negative[order(top_negative$logFC), ]
top_positive_3 <- head(top_positive_sorted, 3)
top_negative_3 <- head(top_negative_sorted, 3)
print("Top 3 Positive Genes:")
print(top_positive_3)
print("Top 3 Negative Genes:")
print(top_negative_3)

# condition1 positive

top3_genes <- rownames(top_positive_3)
ab_comb <- cbind(
  txi_control$abundance[top3_genes, ],
  txi_condition1$abundance[top3_genes, ]
)
colnames(ab_comb) <- c(
  "control_rep1", "control_rep2", "control_rep3",
  "condition1_rep1", "condition1_rep2", "condition1_rep3"
)
abundance_long <- as.data.frame(ab_comb)
abundance_long$GeneID <- rownames(abundance_long)
abundance_long <- abundance_long %>%
  pivot_longer(
    cols = -GeneID,
    names_to = "Condition_Replicate",
    values_to = "Expression"
  )
abundance_long$Condition <- gsub(
  "_.*", "",
  abundance_long$Condition_Replicate
)
abundance_long$Replicate <- gsub(
  ".*_", "",
  abundance_long$Condition_Replicate
)

p <- ggplot(abundance_long, aes(
  x = Condition, y = Expression,
  fill = Replicate
)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~GeneID) +
  labs(
    x = "Condition",
    y = "Normalized Expression"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("../results/11_edgeR_top_positive_genes_expression_by_condition1.png",
  plot = p,
  width = 10, height = 8, dpi = 300
)
# conditon1 negative

top3_genes <- rownames(top_negative_3)
ab_comb <- cbind(
  txi_control$abundance[top3_genes, ],
  txi_condition1$abundance[top3_genes, ]
)
colnames(ab_comb) <- c(
  "control_rep1", "control_rep2", "control_rep3",
  "condition1_rep1", "condition1_rep2", "condition1_rep3"
)
abundance_long <- as.data.frame(ab_comb)
abundance_long$GeneID <- rownames(abundance_long)
abundance_long <- abundance_long %>%
  pivot_longer(
    cols = -GeneID,
    names_to = "Condition_Replicate",
    values_to = "Expression"
  )
abundance_long$Condition <- gsub(
  "_.*", "",
  abundance_long$Condition_Replicate
)
abundance_long$Replicate <- gsub(
  ".*_", "",
  abundance_long$Condition_Replicate
)

p <- ggplot(abundance_long, aes(
  x = Condition, y = Expression,
  fill = Replicate
)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~GeneID) +
  labs(
    x = "Condition",
    y = "Normalized Expression"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("../results/11_edgeR_top_negative_genes_expression_by_condition1.png",
  plot = p,
  width = 10, height = 8, dpi = 300
)

# condition2 top_positive

top_genes_table <- results_control_vs_cond2$table
top_positive <- top_genes_table[top_genes_table$logFC > 0, ]
top_negative <- top_genes_table[top_genes_table$logFC < 0, ]
top_positive_sorted <- top_positive[order(top_positive$logFC,
  decreasing = TRUE
), ]
top_negative_sorted <- top_negative[order(top_negative$logFC), ]
top_positive_3 <- head(top_positive_sorted, 3)
top_negative_3 <- head(top_negative_sorted, 3)
print("Top 3 Positive Genes:")
print(top_positive_3)
print("Top 3 Negative Genes:")
print(top_negative_3)

top3_genes <- rownames(top_positive_3)
ab_comb <- cbind(
  txi_control$abundance[top3_genes, ],
  txi_condition1$abundance[top3_genes, ]
)
colnames(ab_comb) <- c(
  "control_rep1", "control_rep2", "control_rep3",
  "condition1_rep1", "condition1_rep2", "condition1_rep3"
)
abundance_long <- as.data.frame(ab_comb)
abundance_long$GeneID <- rownames(abundance_long)
abundance_long <- abundance_long %>%
  pivot_longer(
    cols = -GeneID,
    names_to = "Condition_Replicate",
    values_to = "Expression"
  )
abundance_long$Condition <- gsub(
  "_.*", "",
  abundance_long$Condition_Replicate
)
abundance_long$Replicate <- gsub(
  ".*_", "",
  abundance_long$Condition_Replicate
)

p <- ggplot(abundance_long, aes(
  x = Condition, y = Expression,
  fill = Replicate
)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~GeneID) +
  labs(
    x = "Condition",
    y = "Normalized Expression"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("../results/11_edgeR_top_positive_genes_expression_by_condition2.png",
  plot = p,
  width = 10, height = 8, dpi = 300
)
# conditon2 negative

top3_genes <- rownames(top_negative_3)
ab_comb <- cbind(
  txi_control$abundance[top3_genes, ],
  txi_condition1$abundance[top3_genes, ]
)
colnames(ab_comb) <- c(
  "control_rep1", "control_rep2", "control_rep3",
  "condition1_rep1", "condition1_rep2", "condition1_rep3"
)
abundance_long <- as.data.frame(ab_comb)
abundance_long$GeneID <- rownames(abundance_long)
abundance_long <- abundance_long %>%
  pivot_longer(
    cols = -GeneID,
    names_to = "Condition_Replicate",
    values_to = "Expression"
  )
abundance_long$Condition <- gsub(
  "_.*", "",
  abundance_long$Condition_Replicate
)
abundance_long$Replicate <- gsub(
  ".*_", "",
  abundance_long$Condition_Replicate
)

p <- ggplot(abundance_long, aes(
  x = Condition, y = Expression,
  fill = Replicate
)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~GeneID) +
  labs(
    x = "Condition",
    y = "Normalized Expression"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("../results/11_edgeR_top_negative_genes_expression_by_condition2.png",
  plot = p,
  width = 10, height = 8, dpi = 300
)

# limma

dge <- DGEList(counts = counts)

keep <- filterByExpr(dge, design)
dge <- dge[keep, , keep.lib.sizes = FALSE]

dge <- calcNormFactors(dge)

png("../results/11_voom_plot1.png")
v <- voom(dge, design, plot = TRUE)
dev.off()

fit <- lmFit(v, design)
fit <- eBayes(fit)
control_vs_condition1_voom <- topTable(fit, coef = 2, number = Inf)
write.csv(control_vs_condition1_voom, "../results/11_limma_control_vs_condition1_voom.csv")

fit <- lmFit(v, design)
fit <- eBayes(fit)
control_vs_condition2_voom <- topTable(fit, coef = 3, number = Inf)
write.csv(control_vs_condition2_voom, "../results/11_limma_control_vs_condition2_voom.csv")

# volcano plots lima


control_vs_condition1_voom$sig <- ifelse(
  control_vs_condition1_voom$adj.P.Val < 0.05,
  "Significant",
  "Not Significant"
)

control_vs_condition2_voom$sig <- ifelse(
  control_vs_condition2_voom$adj.P.Val < 0.05,
  "Significant",
  "Not Significant"
)

limma_volcano_con1 <- ggplot(data = control_vs_condition1_voom, aes(x = logFC, y = -log10(P.Value), col = sig)) +
  geom_point() +
  scale_color_manual(values = c("Not Significant" = "grey", "Significant" = "red")) +
  theme_minimal() +
  labs(
    x = "logFC",
    y = "-log(PValue)",
    color = "Significance"
  )
ggsave("11_limma_volcano_con1.pdf", plot = limma_volcano_con1, device = "pdf", path = "../results")

limma_volcano_con2 <- ggplot(data = control_vs_condition2_voom, aes(x = logFC, y = -log10(P.Value), col = sig)) +
  geom_point() +
  scale_color_manual(values = c("Not Significant" = "grey", "Significant" = "red")) +
  theme_minimal() +
  labs(
    x = "logFC",
    y = "-log(PValue)",
    color = "Significance"
  )
ggsave("11_limma_volcano_con2.pdf", plot = limma_volcano_con2, device = "pdf", path = "../results")

# venn diagrams limma

lfc_ge_2_cn1 <- row.names(control_vs_condition1_voom[abs(control_vs_condition1_voom$logFC) >= 2, ])
fdr_le_0.05_cn1 <- row.names(control_vs_condition1_voom[control_vs_condition1_voom$adj.P.Val <= 0.05, ])
fdr_le_0.01_cn1 <- row.names(control_vs_condition1_voom[control_vs_condition1_voom$adj.P.Val <= 0.01, ])

deg_list_cn1 <- list(
  LFC_gr_2 = lfc_ge_2_cn1,
  FDR_le_0.05 = fdr_le_0.05_cn1,
  FDR_le_0.01 = fdr_le_0.01_cn1
)

limma_venn_diagram_cn1 <- ggVennDiagram(deg_list_cn1) +
  scale_fill_gradient(low = "white", high = "blue") +
  theme_minimal()
ggsave("../results/11_limma_venn_diagram_cn1.pdf", plot = limma_venn_diagram_cn1, device = "pdf")

lfc_ge_2_cn2 <- row.names(control_vs_condition2_voom[abs(control_vs_condition2_voom$logFC) >= 2, ])
fdr_le_0.05_cn2 <- row.names(control_vs_condition2_voom[control_vs_condition2_voom$adj.P.Val <= 0.05, ])
fdr_le_0.01_cn2 <- row.names(control_vs_condition2_voom[control_vs_condition2_voom$adj.P.Val <= 0.01, ])

deg_list_cn2 <- list(
  LFC_gr_2 = lfc_ge_2_cn2,
  FDR_le_0.05 = fdr_le_0.05_cn2,
  FDR_le_0.01 = fdr_le_0.01_cn2
)

limma_venn_diagram_cn2 <- ggVennDiagram(deg_list_cn2) +
  scale_fill_gradient(low = "white", high = "blue") +
  theme_minimal()
ggsave("../results/11_limma_venn_diagram_cn2.pdf", plot = limma_venn_diagram_cn2, device = "pdf")

# increased and decreased genes

limma_res_cond1_inc <- row.names(control_vs_condition1_voom[control_vs_condition1_voom$logFC > 0, ])
limma_res_cond1_dec <- row.names(control_vs_condition1_voom[control_vs_condition1_voom$logFC < 0, ])

limma_res_cond2_inc <- row.names(control_vs_condition2_voom[control_vs_condition2_voom$logFC > 0, ])
limma_res_cond2_dec <- row.names(control_vs_condition2_voom[control_vs_condition2_voom$logFC < 0, ])

# top 20 differentially expressed geens
top_up_cn1 <- head(control_vs_condition1_voom[order(control_vs_condition1_voom$logFC, decreasing = TRUE), ], 20)
top_up_cn2 <- head(control_vs_condition2_voom[order(control_vs_condition2_voom$logFC, decreasing = TRUE), ], 20)
top_do_cn1 <- head(control_vs_condition1_voom[order(control_vs_condition1_voom$logFC, decreasing = FALSE), ], 20)
top_do_cn2 <- head(control_vs_condition2_voom[order(control_vs_condition2_voom$logFC, decreasing = FALSE), ], 20)

top_up_cn1_limma <- row.names(top_up_cn1)
top_up_cn2_limma <- row.names(top_up_cn2)
top_down_cn1_limma <- row.names(top_do_cn1)
top_down_cn2_limma <- row.names(top_do_cn2)


top_up_cn1_desc <- subset(go_annotations, gene_id %in% rownames(top_up_cn1))
top_up_cn2_desc <- subset(go_annotations, gene_id %in% rownames(top_up_cn2))
top_do_cn1_desc <- subset(go_annotations, gene_id %in% rownames(top_do_cn1))
top_do_cn2_desc <- subset(go_annotations, gene_id %in% rownames(top_do_cn2))

top_up_cn1_list <- str_split(top_up_cn1_desc$protein_description, " ") %>%
  unlist() %>%
  tolower()
top_up_cn2_list <- str_split(top_up_cn2_desc$protein_description, " ") %>%
  unlist() %>%
  tolower()
top_do_cn1_list <- str_split(top_do_cn1_desc$protein_description, " ") %>%
  unlist() %>%
  tolower()
top_do_cn2_list <- str_split(top_do_cn2_desc$protein_description, " ") %>%
  unlist() %>%
  tolower()

top_up_cn1_list <- top_up_cn1_list[!top_up_cn1_list %in% bad_words_list & nchar(top_up_cn1_list) > 2]
top_up_cn2_list <- top_up_cn2_list[!top_up_cn2_list %in% bad_words_list & nchar(top_up_cn2_list) > 2]
top_do_cn1_list <- top_do_cn1_list[!top_do_cn1_list %in% bad_words_list & nchar(top_do_cn1_list) > 2]
top_do_cn2_list <- top_do_cn2_list[!top_do_cn2_list %in% bad_words_list & nchar(top_do_cn2_list) > 2]

top_up_cn1_list_df <- as.data.frame(table(top_up_cn1_list))
top_up_cn2_list_df <- as.data.frame(table(top_up_cn2_list))
top_do_cn1_list_df <- as.data.frame(table(top_do_cn1_list))
top_do_cn2_list_df <- as.data.frame(table(top_do_cn2_list))

pdf("../results/11_limma_wordcloud_cn1_up.pdf")
wordcloud(
  words = top_up_cn1_list_df$top_up_cn1_list, freq = top_up_cn1_list_df$Freq, min.freq = 1,
  max.words = 200, random.order = FALSE, rot.per = 0.35,
  colors = brewer.pal(8, "Set1")
)
dev.off()

pdf("../results/11_limma_wordcloud_cn2_up.pdf")
wordcloud(
  words = top_up_cn2_list_df$top_up_cn2_list, freq = top_up_cn2_list_df$Freq, min.freq = 1,
  max.words = 200, random.order = FALSE, rot.per = 0.35,
  colors = brewer.pal(8, "Set1")
)
dev.off()

pdf("../results/11_limma_wordcloud_cn1_down.pdf")
wordcloud(
  words = top_down_cn1_list_df$top_down_cn1_list, freq = top_down_cn1_list_df$Freq, min.freq = 1,
  max.words = 200, random.order = FALSE, rot.per = 0.35,
  colors = brewer.pal(8, "Set1")
)
dev.off()

pdf("../results/11_wordcloud_cn2_down.pdf")
wordcloud(
  words = top_down_cn2_list_df$top_down_cn2_list, freq = top_down_cn2_list_df$Freq, min.freq = 1,
  max.words = 200, random.order = FALSE, rot.per = 0.35,
  colors = brewer.pal(8, "Set1")
)
dev.off()

# table limma

limma_table <- edgeR_table[, c(
  "gene_id",
  "TPM-Repl-1-Control",
  "TPM-Repl-2-Control",
  "TPM-Repl-3-Control",
  "TPM-Repl-1-Condition1",
  "TPM-Repl-2-Condition1",
  "TPM-Repl-3-Condition1",
  "TPM-Repl-1-Condition2",
  "TPM-Repl-2-Condition2",
  "TPM-Repl-3-Condition2"
)]

control_vs_condition1_voom$gene_id <- row.names(control_vs_condition1_voom)
limma_table <- merge(limma_table, control_vs_condition1_voom[, c("gene_id", "logFC")], by = "gene_id")
names(limma_table)[[11]] <- "logFC-Control-vs-Condition1"

control_vs_condition2_voom$gene_id <- row.names(control_vs_condition2_voom)
limma_table <- merge(limma_table, control_vs_condition2_voom[, c("gene_id", "logFC")], by = "gene_id")
names(limma_table)[[12]] <- "logFC-Control-vs-Condition2"

limma_table <- merge(limma_table, control_vs_condition1_voom[, c("gene_id", "adj.P.Val")], by = "gene_id")
names(edgeR_table)[[13]] <- "Adjusted p-value-Control-vs-Condition1"

limma_table <- merge(limma_table, control_vs_condition2_voom[, c("gene_id", "adj.P.Val")], by = "gene_id")
names(edgeR_table)[[14]] <- "Adjusted p-value-Control-vs-Condition2"

# grouped barcharts limma

top_genes_table <- control_vs_condition1_voom
top_positive <- top_genes_table[top_genes_table$logFC > 0, ]
top_negative <- top_genes_table[top_genes_table$logFC < 0, ]
top_positive_sorted <- top_positive[order(top_positive$logFC,
  decreasing = TRUE
), ]
top_negative_sorted <- top_negative[order(top_negative$logFC), ]
top_positive_3 <- head(top_positive_sorted, 3)
top_negative_3 <- head(top_negative_sorted, 3)
print("Top 3 Positive Genes:")
print(top_positive_3)
print("Top 3 Negative Genes:")
print(top_negative_3)

# condition1 positive

top3_genes <- rownames(top_positive_3)
ab_comb <- cbind(
  txi_control$abundance[top3_genes, ],
  txi_condition1$abundance[top3_genes, ]
)
colnames(ab_comb) <- c(
  "control_rep1", "control_rep2", "control_rep3",
  "condition1_rep1", "condition1_rep2", "condition1_rep3"
)
abundance_long <- as.data.frame(ab_comb)
abundance_long$GeneID <- rownames(abundance_long)
abundance_long <- abundance_long %>%
  pivot_longer(
    cols = -GeneID,
    names_to = "Condition_Replicate",
    values_to = "Expression"
  )
abundance_long$Condition <- gsub(
  "_.*", "",
  abundance_long$Condition_Replicate
)
abundance_long$Replicate <- gsub(
  ".*_", "",
  abundance_long$Condition_Replicate
)

p <- ggplot(abundance_long, aes(
  x = Condition, y = Expression,
  fill = Replicate
)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~GeneID) +
  labs(
    title = "Top Genes Expression by Condition",
    x = "Condition",
    y = "Normalized Expression"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("../results/11_top_positive_genes_expression_by_condition1_limma.png",
  plot = p,
  width = 10, height = 8, dpi = 300
)
# conditon1 negative

top3_genes <- rownames(top_negative_3)
ab_comb <- cbind(
  txi_control$abundance[top3_genes, ],
  txi_condition1$abundance[top3_genes, ]
)
colnames(ab_comb) <- c(
  "control_rep1", "control_rep2", "control_rep3",
  "condition1_rep1", "condition1_rep2", "condition1_rep3"
)
abundance_long <- as.data.frame(ab_comb)
abundance_long$GeneID <- rownames(abundance_long)
abundance_long <- abundance_long %>%
  pivot_longer(
    cols = -GeneID,
    names_to = "Condition_Replicate",
    values_to = "Expression"
  )
abundance_long$Condition <- gsub(
  "_.*", "",
  abundance_long$Condition_Replicate
)
abundance_long$Replicate <- gsub(
  ".*_", "",
  abundance_long$Condition_Replicate
)

p <- ggplot(abundance_long, aes(
  x = Condition, y = Expression,
  fill = Replicate
)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~GeneID) +
  labs(
    title = "Top Genes Expression by Condition",
    x = "Condition",
    y = "Normalized Expression"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("../results/11_top_negative_genes_expression_by_condition1_limma.png",
  plot = p,
  width = 10, height = 8, dpi = 300
)

# condition2 top_positive

top_genes_table <- control_vs_condition2_voom
top_positive <- top_genes_table[top_genes_table$logFC > 0, ]
top_negative <- top_genes_table[top_genes_table$logFC < 0, ]
top_positive_sorted <- top_positive[order(top_positive$logFC,
  decreasing = TRUE
), ]
top_negative_sorted <- top_negative[order(top_negative$logFC), ]
top_positive_3 <- head(top_positive_sorted, 3)
top_negative_3 <- head(top_negative_sorted, 3)
print("Top 3 Positive Genes:")
print(top_positive_3)
print("Top 3 Negative Genes:")
print(top_negative_3)

top3_genes <- rownames(top_positive_3)
ab_comb <- cbind(
  txi_control$abundance[top3_genes, ],
  txi_condition1$abundance[top3_genes, ]
)
colnames(ab_comb) <- c(
  "control_rep1", "control_rep2", "control_rep3",
  "condition2_rep1", "condition2_rep2", "condition2_rep3"
)
abundance_long <- as.data.frame(ab_comb)
abundance_long$GeneID <- rownames(abundance_long)
abundance_long <- abundance_long %>%
  pivot_longer(
    cols = -GeneID,
    names_to = "Condition_Replicate",
    values_to = "Expression"
  )
abundance_long$Condition <- gsub(
  "_.*", "",
  abundance_long$Condition_Replicate
)
abundance_long$Replicate <- gsub(
  ".*_", "",
  abundance_long$Condition_Replicate
)

p <- ggplot(abundance_long, aes(
  x = Condition, y = Expression,
  fill = Replicate
)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~GeneID) +
  labs(
    title = "Top Genes Expression by Condition",
    x = "Condition",
    y = "Normalized Expression"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("../results/11_top_positive_genes_expression_by_condition2_limma.png",
  plot = p,
  width = 10, height = 8, dpi = 300
)
# conditon2 negative

top3_genes <- rownames(top_negative_3)
ab_comb <- cbind(
  txi_control$abundance[top3_genes, ],
  txi_condition1$abundance[top3_genes, ]
)
colnames(ab_comb) <- c(
  "control_rep1", "control_rep2", "control_rep3",
  "condition1_rep2", "condition2_rep2", "condition2_rep3"
)
abundance_long <- as.data.frame(ab_comb)
abundance_long$GeneID <- rownames(abundance_long)
abundance_long <- abundance_long %>%
  pivot_longer(
    cols = -GeneID,
    names_to = "Condition_Replicate",
    values_to = "Expression"
  )
abundance_long$Condition <- gsub(
  "_.*", "",
  abundance_long$Condition_Replicate
)
abundance_long$Replicate <- gsub(
  ".*_", "",
  abundance_long$Condition_Replicate
)

p <- ggplot(abundance_long, aes(
  x = Condition, y = Expression,
  fill = Replicate
)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~GeneID) +
  labs(
    title = "Top Genes Expression by Condition",
    x = "Condition",
    y = "Normalized Expression"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("../results/11_top_negative_genes_expression_by_condition2_limma.png",
  plot = p,
  width = 10, height = 8, dpi = 300
)

# jaccard sim
jacc_up_cn1 <- length(intersect(top_up_cn1_edgeR, top_up_cn1_limma)) / length(union(top_up_cn1_edgeR, top_up_cn1_limma))
jacc_up_cn2 <- length(intersect(top_up_cn2_edgeR, top_up_cn2_limma)) / length(union(top_up_cn2_edgeR, top_up_cn2_limma))
jacc_do_cn1 <- length(intersect(top_down_cn1_edgeR, top_down_cn1_limma)) / length(union(top_down_cn1_edgeR, top_down_cn1_limma))
jacc_do_cn2 <- length(intersect(top_down_cn2_edgeR, top_down_cn2_limma)) / length(union(top_down_cn2_edgeR, top_down_cn2_limma))

print("Jaccard similarity coefficient for up regulated genes in condition 1:")
print(jacc_up_cn1)
print("Jaccard similarity coefficient for up regulated genes in condition 2:")
print(jacc_up_cn2)
print("Jaccard similarity coefficient for down regulated genes in condition 1:")
print(jacc_do_cn1)
print("Jaccard similarity coefficient for down regulated genes in condition 2:")
print(jacc_do_cn2)
