## Load libraries --------------------------------------------------------------
library(tidyverse)


##  Plot read length distribution ----------------------------------------------
### Read in files with length ~ quality ----------------------------------------
files_lq <- list.files(
  "fastq_pass/QC/nanoplot_raw_barcoded",
  pattern = "NanoPlot-data.tsv.gz",
  full.names = TRUE)

nanoplot_lq <- map(files_lq, read_tsv, id = "sample") %>% 
  list_rbind() %>% 
  mutate(
    sample = str_split_i(sample, pattern = "/", i = -1),
    sample = str_sub(sample, 1, 3))

### Plot histogram of read length ----------------------------------------------
#### All samples
nanoplot_lq %>% 
  filter(sample == "ALL") %>% 
  ggplot(aes(x = lengths)) +
  geom_histogram(fill = "cornflowerblue", color = "black") +
  scale_x_continuous(labels = scales::number_format(scale = 0.001, suffix = " kb")) +
  labs(x = "Read length", y = "Number of reads") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(color = "black", size = 12),
    axis.text.x = element_text(hjust = 0.7),
    aspect.ratio = 2/3
  )

ggsave("plots/01_read_length.png", width = 6, height = 4, units = "cm", dpi = 300, scale = 2)
ggsave("plots/01_read_length.pdf", width = 6, height = 4, units = "cm", scale = 2)

