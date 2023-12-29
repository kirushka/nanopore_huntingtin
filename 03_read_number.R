## Load libraries --------------------------------------------------------------
library(tidyverse)


## Plot number of reads per barcode group --------------------------------------
### Read in files with stats metrics -------------------------------------------
files_stats <- list.files(
  "fastq_pass/QC/nanoplot_raw_barcoded", 
  pattern = "NanoStats.txt",
  full.names = TRUE)

nanoplot_stats <- map(
  files_stats, read_tsv, 
  skip = 1, col_names = c("metric", "value"), id = "sample") %>% 
  list_rbind() %>% 
  mutate(
    sample = str_split_i(sample, pattern = "/", i = -1),
    sample = str_sub(sample, 1, 3))


### Barplot --------------------------------------------------------------------
nanoplot_stats %>% 
  filter(sample != "ALL") %>% 
  filter(metric == "number_of_reads") %>% 
  mutate(
    group = case_when(
      sample %in% c("B04", "B05", "B06") ~ "Q138",
      sample %in% c("B07", "B08") ~ "Q15",
      TRUE ~ "Unclassified"),
    group = factor(group, levels = c("Q15", "Q138", "Unclassified"))) %>% 
  group_by(group) %>% 
  summarise(n_reads = sum(as.numeric(value))) %>% 
  ggplot(aes(x = group, y = n_reads, fill = group)) +
  geom_col() +
  geom_text(aes(y = n_reads + n_reads*0.2, label = n_reads)) +
  scale_fill_brewer(palette = "Set2") +
  scale_y_log10() +
  labs(x = "", y = "Number of reads", fill = "") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 14),
    axis.text = element_text(color = "black", size = 12),
    aspect.ratio = 1
  )

ggsave("plots/02_reads_number.png", width = 5, height = 4, units = "cm", dpi = 300, scale = 2)
ggsave("plots/02_reads_number.pdf", width = 5, height = 4, units = "cm", scale = 2)

