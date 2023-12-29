## Load libraries --------------------------------------------------------------
library(tidyverse)


##  Plot repeat length distribution --------------------------------------------
### Read in files with additional repeat lenhgth -------------------------------
files_rl <- list.files(
  "tandem_genotypes/plasmid_flt", 
  pattern = ".txt", 
  full.names = TRUE)


repeat_length <- map(
  files_rl, 
  read_tsv,
  comment = "#", 
  col_names = c("seq_name", "start", "end", "repeat", "smth1", "smth2", "fwd", "rev"),
  col_types = "cddccccc",
  id = "file") %>% 
  list_rbind() %>% 
  mutate(
    group = str_extract(file, "Q\\d+|UNC")) %>% 
  select(group, fwd, rev)

repeats_lng <- repeat_length %>% 
  pivot_longer(c(fwd, rev), names_to = "strand", values_to = "repeat_num") %>% 
  separate_rows(repeat_num, sep = ",") %>% 
  # Add 23 Qs present in reference plasmid
  mutate(repeat_num = as.numeric(repeat_num) + 23) 


### Plot only Q138 with valid barcodes -----------------------------------------
repeats_lng_q138 <- repeats_lng %>% 
  filter(group == "Q138")

repeats_lng_q138_mode <- repeats_lng_q138 %>% 
  count(repeat_num, sort = TRUE) %>% 
  slice_head(n = 1) %>% 
  pull(repeat_num)

repeats_lng_q138_med <- repeats_lng_q138 %>% 
  summarise(repeat_num = median(repeat_num)) %>% 
  pull()

repeats_lng_q138 %>% 
  ggplot(aes(x = repeat_num)) +
  geom_histogram(binwidth = 5, fill = "#fc8d62",  color = "black") +
  geom_vline(xintercept = repeats_lng_q138_mode, linetype = "dashed") +
  annotate("text", x = 122, y = 100, label = repeats_lng_q138_mode, size = 5) +
  labs(x = "CAG copy number", y = "Number of reads", fill = "") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(color = "black", size = 12),
    aspect.ratio = 2/3
  )

ggsave("plots/03_repeat_extension.png", width = 6, height = 4, units = "cm", dpi = 300, scale = 2)
ggsave("plots/03_repeat_extension.pdf", width = 6, height = 4, units = "cm", scale = 2)


### Plot only Q138 + long Unclassified -----------------------------------------
repeats_lng_q138_unc <- repeats_lng %>% 
  filter(group == "Q138" | (group == "UNC" & repeat_num > 100))

repeats_lng_q138_unc_mode <- repeats_lng_q138_unc %>% 
  count(repeat_num, sort = TRUE) %>% 
  slice_head(n = 2) %>% 
  pull(repeat_num)

repeats_lng_q138_unc_med <- repeats_lng_q138_unc %>% 
  summarise(repeat_num = median(repeat_num)) %>% 
  pull()

repeats_lng_q138_unc %>% 
  ggplot(aes(x = repeat_num)) +
  geom_histogram(binwidth = 5, fill = "#fc8d62",  color = "black") +
  # geom_vline(xintercept = repeats_lng_q138_mode, linetype = "dashed") +
  annotate("text", x = 115, y = 115, label = paste(repeats_lng_q138_unc_mode, collapse = ";"), size = 5) +
  scale_x_continuous(breaks = seq(0, 160, 20)) +
  labs(x = "CAG copy number", y = "Number of reads", fill = "") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(color = "black", size = 12),
    aspect.ratio = 2/3
  )

ggsave("plots/03_repeat_extension_unc.png", width = 6, height = 4, units = "cm", dpi = 300, scale = 2)
ggsave("plots/03_repeat_extension_unc.pdf", width = 6, height = 4, units = "cm", scale = 2)
