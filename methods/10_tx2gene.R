library(tidyverse)
library(stringr)

abundance <- read_delim("../material/control/9929263/abundance.tsv")

abundance$gene_id <- ifelse(
  is.na(str_match(abundance$target_id, "(.*)_.*")[, 2]),
  abundance$target_id,
  str_match(abundance$target_id, "(.*)_.*")[, 2]
)

tx2gene <- tibble(
  gene_id = abundance$gene_id
)

tx2gene <- merge(tx2gene, abundance[, c("gene_id", "target_id")], by = "gene_id", all = TRUE)

tx2gene <- tx2gene[, c("target_id", "gene_id")]

names(tx2gene)[[1]] <- "transcript_id"

tx2gene <- unique(tx2gene)

write.table(tx2gene, "../results/tx2gene", quote = FALSE, row.names = FALSE)
