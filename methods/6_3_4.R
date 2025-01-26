library(tidyverse)
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(corrplot)

# 6_3

# kallisto

abundance_files <- list.files(
  path = "../material",
  pattern = "abundance.tsv",
  recursive = TRUE,
  full.names = TRUE
)

abundance_list <- list()

for (file in abundance_files) {
  seq_name <- str_match(file, "\\d{7}")
  abundance_list[[seq_name]] <- read_delim(file)
}

number_of_expressed_genes <- list()

for (i in seq_along(abundance_list)) {
  number_of_expressed_genes[[names(abundance_list)[i]]] <- sum(abundance_list[[i]]$est_counts > 0)
}

print("Number of the expressed genes: ")
print(unlist(number_of_expressed_genes))

fraction_of_expressed_genes <- as.numeric(number_of_expressed_genes) / 11599
names(fraction_of_expressed_genes) <- names(abundance_list)

print("Fraction of expressed genes: ")
print(fraction_of_expressed_genes)

abundance_list[1]

print("Summary of the TPM values for each replicate")
for (i in seq_along(abundance_list)) {
  print(names(abundance_list)[i])
  print(summary(abundance_list[[i]]$tpm))
}

for (i in seq_along(abundance_list)) {
  abundance_list[[i]]$replicate <- names(abundance_list)[i]
}

names(abundance_list)

abundance_list[["9929263"]]$group <- "control"
abundance_list[["9929264"]]$group <- "control"
abundance_list[["9929273"]]$group <- "control"
abundance_list[["9929265"]]$group <- "condition1"
abundance_list[["9929271"]]$group <- "condition1"
abundance_list[["9929274"]]$group <- "condition1"
abundance_list[["9929272"]]$group <- "condition2"
abundance_list[["9929279"]]$group <- "condition2"
abundance_list[["9929281"]]$group <- "condition2"

abundance_list[["9929263"]]$replicate <- "9929263.ctl"
abundance_list[["9929264"]]$replicate <- "9929264.ctl"
abundance_list[["9929273"]]$replicate <- "9929273.ctl"
abundance_list[["9929265"]]$replicate <- "9929265.cn1"
abundance_list[["9929271"]]$replicate <- "9929271.cn1"
abundance_list[["9929274"]]$replicate <- "9929274.cn1"
abundance_list[["9929272"]]$replicate <- "9929272.cn2"
abundance_list[["9929279"]]$replicate <- "9929279.cn2"
abundance_list[["9929281"]]$replicate <- "9929281.cn2"

abundance_long <- bind_rows(abundance_list)

kallisto_boxplot <- ggplot(abundance_long, aes(x = group, y = log(tpm + 1), fill = interaction(replicate, group))) +
  geom_boxplot() +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1") +
  labs(y = "log(TPM + 1)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("../results/6_3_kallisto_boxplot.pdf", device = "pdf", plot = kallisto_boxplot)

abundance_wide <- abundance_long[, c("target_id", "tpm", "replicate")] %>%
  pivot_wider(names_from = replicate, values_from = tpm)

numeric_data <- abundance_wide %>%
  select(-target_id)

cor_matrix <- cor(numeric_data, use = "complete.obs")

pdf("../results/6_3_kallisto_corr_matrix.pdf")
corrplot(cor_matrix,
  method = "color",
  type = "full",
  order = "hclust",
  hclust.method = "ward.D2",
  addrect = 5,
  tl.col = "black"
)
dev.off()

# bowtie

names <- c(
  "9929263.ctl",
  "9929264.ctl",
  "9929273.ctl",
  "9929265.cn1",
  "9929271.cn1",
  "9929274.cn1",
  "9929272.cn2",
  "9929279.cn2",
  "9929281.cn2"
)

uniquely_aligned_reads <- c(
  20501750,
  21853000,
  18269738,
  23078461,
  21422437,
  18774347,
  21072416,
  20559264,
  22626535
)

multi_mapped_reads <- c(
  11177080,
  7309758,
  6808188,
  6726601,
  6147476,
  5941791,
  6134765,
  5740451,
  7423736
)

not_aligned_reads <- c(
  1534102,
  2375543,
  1583826,
  2422843,
  1581190,
  2287348,
  2564081,
  1298564,
  1429119
)

bowtie_alignments <- tibble(
  names,
  uniquely_aligned_reads,
  multi_mapped_reads,
  not_aligned_reads
)

bowtie_alignments_long <- bowtie_alignments %>%
  pivot_longer(-names, names_to = "category", values_to = "alignments")

bowtie_alignment_bar <- ggplot(bowtie_alignments_long, aes(x = names, y = alignments, fill = category)) +
  geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal() +
  labs(x = "Replicate", y = "Number of Sequences") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("../results/6_3_bowtie_alignment_bar.pdf", device = "pdf", plot = bowtie_alignment_bar)

# 6_4

bowtie_files <- list.files(path = "../material", pattern = "gene_abundance", full.names = TRUE)

bowtie_list <- list()

for (i in seq_along(bowtie_files)) {
  replicate <- str_match(bowtie_files[[i]], "\\d{7}")
  bowtie_list[[replicate]] <- read_delim(bowtie_files[[i]])
}

bowtie_tpm_list <- list()

for (i in seq_along(bowtie_list)) {
  df <- tibble(
    gene_id = str_match(bowtie_list[[i]]$"Gene ID", "^[^_ ]+"),
    tpm = bowtie_list[[i]]$TPM
  )
  bowtie_tpm_list[[names(bowtie_list)[[i]]]] <- df
}

bowtie_tpm_agg_list <- map2(
  bowtie_tpm_list,
  names(bowtie_tpm_list),
  ~ {
    agg_data <- aggregate(tpm ~ gene_id, data = .x, FUN = sum) %>%
      as_tibble()
    agg_data <- agg_data %>%
      mutate(replicate = .y)
    return(agg_data)
  }
)

bowtie_tpm_all <- bind_rows(bowtie_tpm_agg_list)

bowtie_tpm_all <- bowtie_tpm_all %>%
  mutate(group = case_when(
    replicate %in% c("9929263", "9929264", "9929273") ~ "control",
    replicate %in% c("9929265", "9929271", "9929274") ~ "condition1",
    replicate %in% c("9929272", "9929279", "9929281") ~ "condition2",
  ))

bowtie_tpm_all <- bowtie_tpm_all %>%
  mutate(algorithm = "Bowtie2")

abundance_long$target_id <- str_match(abundance_long$target_id, "^[^_ ]+")

abundance_long$replicate <- str_match(abundance_long$replicate, "^[^_ ]+")

kallisto_tpm_all <- abundance_long %>%
  select(target_id, tpm, replicate, group)

colnames(kallisto_tpm_all) <- names(bowtie_tpm_all)

kallisto_tpm_all <- aggregate(tpm ~ gene_id + replicate + group, data = kallisto_tpm_all, FUN = sum) %>%
  as_tibble()

kallisto_tpm_all <- kallisto_tpm_all %>%
  mutate(algorithm = "Kallisto")

tpm_all <- bind_rows(list(bowtie_tpm_all, kallisto_tpm_all))

# post_quantification

all_tpm_boxplot <- ggplot(data = tpm_all, aes(x = interaction(group, algorithm), y = log(tpm + 1), fill = algorithm)) +
  geom_boxplot() +
  labs(x = "Gene expression per condition and algorithm", y = expression("log(TPM + 1)")) +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("../results/6_post_all_tpm_boxplot.pdf", plot = all_tpm_boxplot, device = "pdf")

tpm_all_agg <- aggregate(tpm ~ gene_id + group + algorithm, tpm_all, FUN = mean) %>%
  as_tibble()

tpm_all_agg_wide <- tpm_all_agg %>%
  pivot_wider(names_from = group, values_from = tpm)

tpm_all_lfc <- tpm_all_agg_wide %>%
  mutate(
    lfc1 = log2(condition1 / control),
    lfc2 = log2(condition2 / control)
  )

tpm_all_lfc_long <- tpm_all_lfc %>%
  select(gene_id, algorithm, lfc1, lfc2) %>%
  pivot_longer(-c(gene_id, algorithm))

density_lfc <- ggplot(tpm_all_lfc_long, aes(x = value, fill = interaction(name, algorithm))) +
  geom_density(alpha = 0.3) +
  labs(x = "TPM", fill = "Legend") +
  xlim(-8, 8) +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal()

ggsave("../results/6_post_all_lfc_density.pdf", plot = density_lfc, device = "pdf")

# comparison kallisto and bowtie2

tpm_comp <- tpm_all_agg %>%
  pivot_wider(names_from = c(group, algorithm), values_from = tpm)

tpm_comp <- tpm_comp %>%
  mutate(
    diff_ctl = control_Kallisto - control_Bowtie2,
    diff_cn1 = condition1_Kallisto - condition1_Bowtie2,
    diff_cn2 = condition2_Kallisto - condition2_Bowtie2,
    frac_ctl = control_Kallisto / control_Bowtie2,
    frac_cn1 = condition1_Kallisto / condition1_Bowtie2,
    frac_cn2 = condition2_Kallisto / condition2_Bowtie2
  )

tpm_comp_dif_long <- tpm_comp %>%
  select(diff_ctl, diff_cn1, diff_cn2, gene_id) %>%
  pivot_longer(-gene_id, names_to = "comp", values_to = "diff")

density_diff <- ggplot(tpm_comp_dif_long, aes(x = diff, fill = comp)) +
  geom_density(alpha = 0.3) +
  labs(x = "Difference in TPM", fill = "Legend") +
  xlim(-200, 100) +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal()

ggsave("../results/6_post_all_diff_density.pdf", plot = density_diff, device = "pdf")

tpm_comp_frac_long <- tpm_comp %>%
  select(frac_ctl, frac_cn1, frac_cn2, gene_id) %>%
  pivot_longer(-gene_id, names_to = "comp", values_to = "frac")

density_frac <- ggplot(tpm_comp_frac_long, aes(x = frac, fill = comp)) +
  geom_density(alpha = 0.3) +
  labs(x = "Fraction of TPM", fill = "Legend") +
  xlim(-0.5, 3) +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal()

ggsave("../results/6_post_all_frac_density.pdf", plot = density_frac, device = "pdf")

lfc_k_1 <- tpm_all_lfc_long %>%
  filter(value >= 2) %>%
  filter(algorithm == "Kallisto" & name == "lfc1") %>%
  select(gene_id)


lfc_k_2 <- tpm_all_lfc_long %>%
  filter(value >= 2) %>%
  filter(algorithm == "Kallisto" & name == "lfc2") %>%
  select(gene_id)

lfc_b_1 <- tpm_all_lfc_long %>%
  filter(value >= 2) %>%
  filter(algorithm == "Bowtie2" & name == "lfc1") %>%
  select(gene_id)

lfc_b_2 <- tpm_all_lfc_long %>%
  filter(value >= 2) %>%
  filter(algorithm == "Bowtie2" & name == "lfc2") %>%
  select(gene_id)

jaccard_sim_1 <- nrow(intersect(lfc_b_1, lfc_k_1)) / nrow(union(lfc_b_1, lfc_k_1))
jaccard_sim_2 <- nrow(intersect(lfc_b_2, lfc_k_2)) / nrow(union(lfc_b_2, lfc_k_2))

print("The jaccard similarity coefficient for condition 1 is:")
print(jaccard_sim_1)
print("The jaccard similarity coefficient for condition 2 is:")
print(jaccard_sim_2)
