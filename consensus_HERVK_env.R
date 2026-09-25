library(dplyr)
library(ggplot2)
library(stringr)
library(readr)
library(tidyr)
library(ggeasy)
library(edgeR)
library(ggrepel)

align <- read.delim("/data2/lackey_lab/DownloadedSequenceData/randazza/nerve_root/rrna_gone_herv_con/HERVK_env_isoform_counts.txt", sep = "\t")
align$sample <- c("ATRTct", "ATRTkd","N6", "N7", "N8", "P10", "P11", "C1", "C2", "C3", "C4", "P4b", "P4a", "P5", "P6", "P7", "P8", "P9")
counts <- align[,c("sample", "Total_HERVK", "Total_library_reads", "HERVK_to_library_ratio")]
addition <- data.frame(Sample = c("N6", "N7", "N8", "P10", "P11", "C1", "C2", "2TC", "C3", "C4",
                                  "P4b", "P4a", "P5", "P6", "P7", "P8", "P9", "ATRTct", "ATRTkd"),
                       Category = c("patient", "patient", "patient", "patient", "patient", "control",
                                    "control", "trauma", "control", "control", "patient", "patient", "patient",
                                    "patient", "patient", "patient", "patient", "ATRT", "ATRT"),
                       Sex = c("F", "F", "F", "M", "F", "M", "M", "M", "F", "M",
                               "F", "F", "F", "M", "F", "F", "F", "M", "M"))
counts <- merge(counts, addition, by.x = "sample", by.y = "Sample")
pos <- position_jitterdodge(jitter.width = 0.15, dodge.width = 0.9)

#HERVK_to_library_ratio = ratio of reads aligned to HERVK env to total reads mapped to human genome
counts %>%
  ggplot(aes(x = Category, y = HERVK_to_library_ratio, fill = Category)) +
  geom_boxplot(position = position_dodge(width = 0.9), outlier.shape = NA) +
  geom_jitter(
    aes(color = Sex, group = Category),
    position = pos,
    size = 3,
    alpha = 0.7
  ) +
  geom_text_repel(
    aes(label = sample, color = Sex, group = Category),
    position = pos,
    size = 4,
    max.overlaps = 100,
    box.padding = 0.3,
    point.padding = 0.2,
    segment.color = "grey50"
  ) +
  labs(y = "HERVK_reads/total_reads") +
  scale_fill_manual(values = c("control" = "gray", "patient" = "purple")) +
  theme_minimal()

#Read depth at each position of the consensus envelope sequence (samtools depth)
dist <- read.delim("/data2/lackey_lab/DownloadedSequenceData/randazza/nerve_root/rrna_gone_herv_con/coverage_matrix.txt", sep = "\t", header = F)
dist$V1 <- NULL
colnames(dist) <- c("position", "ATRTct", "ATRTNMDPos", "ATRTkd", "N6", "N8", "N7", "P10", "P11", "C1",
                    "C2", "C3", "C4", "P4b", "P4a", "P5", "P6", "P8", "P9")
dist <- dist[, c("ATRTct", "ATRTkd","N6", "N7", "N8", "P10", "P11", "C1", "C2", "C3", "C4",
                 "P4b", "P4a", "P5", "P6", "P8", "P9", "position")]
dist_long <- dist %>%
  pivot_longer(cols = c("ATRTct", "ATRTkd","N6", "N7", "N8", "P10", "P11", "C1", "C2", "C3", "C4",
                        "P4b", "P4a", "P5", "P6", "P8", "P9"), names_to = "sample", values_to = "coverage")
##RPM coverage = coverage per position/(mapped reads/10^6)
combo <- dist_long %>%
  left_join(counts, by = "sample")
combo <- combo %>%
  mutate(RPM_totlib = coverage / (Total_library_reads / 1e6)) 
combo <- combo %>%
  mutate(position = as.numeric(position))

combo %>%
  arrange(sample, position) %>%
  ggplot(aes(position, RPM_totlib, color = sample)) +
  geom_line() +
  theme_minimal()

combo_avg <- combo %>%
  filter(!sample %in% c("NIH", "P7")) %>%
  group_by(Category, position) %>%
  summarise(
    med = median(RPM_totlib, na.rm = TRUE),
    mean = mean(RPM_totlib, na.rm = TRUE),
    sd = sd(RPM_totlib, na.rm = TRUE),
    .groups = "drop"
  )

combo_avg %>%
  ggplot() +
  annotate("rect",
           xmin = 483, xmax = 540,
           ymin = -Inf, ymax = Inf,
           fill = "yellow", alpha = 0.3) +
  annotate("rect",
           xmin = 843, xmax = 930,
           ymin = -Inf, ymax = Inf,
           fill = "yellow", alpha = 0.3) +
  annotate("rect",
           xmin = 1263, xmax = 1350,
           ymin = -Inf, ymax = Inf,
           fill = "yellow", alpha = 0.3) +
  geom_ribbon(
    aes(x = position, ymin = mean - sd, ymax = mean + sd, fill = Category),
    alpha = 0.25
  ) +
  geom_line(
    aes(x = position, y = mean, color = Category),
    linewidth = 1
  ) +
  theme_minimal() +
  scale_fill_manual(values = c(
    patient = "purple",
    control = "orange"
  )) +
  scale_color_manual(values = c(
    patient = "purple",
    control = "orange"
  ))

##RPM coverage at epitope locations
combo_1 <- combo %>%
  filter(position >= 483 & position <= 540)
combo_2 <- combo %>%
  filter(position >= 843 & position <= 930)
combo_3 <- combo %>%
  filter(position >= 1263 & position <= 1350)
mega <- rbind(combo_1, combo_2, combo_3)

mega %>%
  filter(!sample %in% c("NIH", "P7")) %>%
  ggplot(aes(x = Category, y = RPM_totlib, fill = Category)) +
  geom_boxplot() +
  theme_minimal() +
  scale_fill_manual(values = c("control" = "gray", "patient" = "purple"))

mega_avg <- mega %>%
  group_by(sample) %>%
  summarise(
    med = median(RPM_totlib, na.rm = TRUE),
    mean = mean(RPM_totlib, na.rm = TRUE),
    .groups = "drop"
  )
mega_avg <- merge(mega_avg, addition, by.x = "sample", by.y = "Sample")
pos <- position_jitterdodge(jitter.width = 0.15, dodge.width = 0.9)

mega_avg %>%
  filter(!sample %in% c("NIH", "P7")) %>%
  ggplot(aes(x = Category, y = med, fill = Category)) +
  geom_boxplot(position = position_dodge(width = 0.9), outlier.shape = NA) +
  geom_jitter(
    aes(color = Sex, group = Category),
    position = pos,
    size = 3,
    alpha = 0.7
  ) +
  geom_text_repel(
    aes(label = sample, color = Sex, group = Category),
    position = pos,
    size = 4,
    max.overlaps = 100,
    box.padding = 0.3,
    point.padding = 0.2,
    segment.color = "grey50"
  ) +
  labs(y = "Median RPM") +
  scale_fill_manual(values = c("control" = "gray", "patient" = "purple")) +
  theme_minimal()
