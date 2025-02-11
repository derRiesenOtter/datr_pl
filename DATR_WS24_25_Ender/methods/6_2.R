library(ggplot2)
library(tidyr)
library(RColorBrewer)

replicate <- c(
  "SRR9929263",
  "SRR9929264",
  "SRR9929273",
  "SRR9929265",
  "SRR9929271",
  "SRR9929274",
  "SRR9929272",
  "SRR9929279",
  "SRR9929281"
)

total_seq <- c(
  33341110,
  31683088,
  26814316,
  32357568,
  29250018,
  27342241,
  30046891,
  27741299,
  31668729
)

total_seq_filtered <- c(
  33212932,
  31538301,
  26661752,
  32227905,
  29151103,
  27003486,
  29771262,
  27598279,
  31479390
)

filter_df <- data.frame(
  replicate = replicate,
  total_seq = total_seq,
  total_seq_filtered = total_seq_filtered,
  diff_seq = total_seq - total_seq_filtered
)

filter_df_long <- filter_df[, !(names(filter_df) == "total_seq")] %>%
  pivot_longer(!replicate, names_to = "category", values_to = "sequences")

bar_chart <- ggplot(filter_df_long, aes(x = replicate, y = sequences, fill = category)) +
  geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal() +
  labs(x = "Replicate", y = "Number of Sequences") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("../results/6_2_filter_bar_chart.pdf",
  plot = bar_chart,
  device = "pdf"
)
