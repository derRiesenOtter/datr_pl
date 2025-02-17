library(tidyverse)
library(ggplot2)
library(tximport)
library(WGCNA)
library(stringr)
library(wordcloud)

# csv import

tx2gene <- read_table("../results/tx2gene")

control_files <- list.files(path = "../material/control", pattern = "tsv$", recursive = TRUE, full.names = TRUE)
condition1_files <- list.files(path = "../material/condition1", pattern = "tsv$", recursive = TRUE, full.names = TRUE)
condition2_files <- list.files(path = "../material/condition2", pattern = "tsv$", recursive = TRUE, full.names = TRUE)

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

abundance <- cbind(
  txi_control$abundance,
  txi_condition1$abundance,
  txi_condition2$abundance
)

# data preparation

abundance_filtered <- abundance[rowSums(abundance > 0) > ncol(abundance) * 0.5, ]

datExpr <- t(abundance_filtered)
dim(datExpr)

# WGCNA

goodSamplesGenes <- goodSamplesGenes(datExpr, verbose = 3)
datExpr <- datExpr[goodSamplesGenes$goodSamples, goodSamplesGenes$goodGenes]

powers <- c(1:20)
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

softPower <- 2

net <- blockwiseModules(datExpr,
  power = softPower,
  TOMType = "unsigned", minModuleSize = 30,
  reassignThreshold = 0, mergeCutHeight = 0.25,
  numericLabels = TRUE, pamRespectsDendro = FALSE,
  saveTOMs = TRUE,
  saveTOMFileBase = "lunghemorrhageTOM",
  verbose = 3
)

# calculate correlation ...

# 9 samples, 3 replicates per condition
traits <- c(rep("control", 3), rep("hs90", 3), rep("hs120", 3))
traits_numeric <- as.numeric(factor(traits,
  levels = c("control", "hs90", "hs120")
))
# module_trait_correlations is the correlation between
# each module eigengene (ME1, ME2, ME3, etc.) and the trait values
# (traits_numeric). The result is a vector of correlation
# coefficients for each module
module_trait_correlations <- cor(net$MEs, traits_numeric, use = "p")
module_trait_pvalues <- corPvalueStudent(module_trait_correlations,
  nSamples = 9
)

print("Only ME2 seems to be significant")
print("Correlations:")
print(module_trait_correlations)
print("p-values:")
print(module_trait_pvalues)


# gene function analysis

ME2_genes <- names(net$colors[net$colors == 2])

go_annotations <- read_delim("../material/go_annotations", skip = 4, col_names = FALSE)

go_annotations$gene_id <- ifelse(
  is.na(str_match(go_annotations[[11]], ".*\\|(.*)")[, 2]),
  go_annotations[[11]],
  str_match(go_annotations[[11]], ".*\\|(.*)")[, 2]
)

names(go_annotations)[[3]] <- "protein_id"
names(go_annotations)[[10]] <- "protein_description"

words_me2 <- subset(go_annotations, gene_id %in% ME2_genes)

words_me2_list <- list()

for (desc in words_me2$protein_description) {
  words_me2_list <- c(words_me2_list, str_split(desc, " ")[[1]])
}

words_me2_list <- unlist(words_me2_list)
words_me2_list <- tolower((words_me2_list))

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

words_me2_list <- words_me2_list[!words_me2_list %in% bad_words_list]
words_me2_list <- words_me2_list[nchar(words_me2_list) >= 3]

words_me2_df <- as.data.frame(table(words_me2_list))
set.seed(1234)

pdf("../results/10_6_wordcloud_me2.pdf")
wordcloud(
  words = words_me2_df$words_me2_list, freq = words_me2_df$Freq, min.freq = 1,
  max.words = 200, random.order = FALSE, rot.per = 0.35,
  colors = brewer.pal(8, "Dark2")
)
dev.off()
