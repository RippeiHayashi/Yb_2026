### nucleotide composition meta lineplot
library(Biostrings)
library(dplyr)
library(ggplot2)
library(tidyr)

setwd("${references}/dm6")

# Load sequences
fasta_file <- "RepBase_Drosophila_Gypsy_512ORFs_nucleotides.Gypsy_17.6_412_mdg1.curated.fasta"
seqs <- readDNAStringSet(fasta_file)

# Parameters
n_bins <- 50  # number of bins along each sequence

seq_to_bins <- function(seq, n_bins) {
  len <- length(seq)
  bases <- unlist(strsplit(as.character(seq), split=""))
  df <- data.frame(pos = 1:len, base = bases)
  df$bin <- cut(df$pos,
                breaks = seq(0, len, length.out = n_bins + 1),
                include.lowest = TRUE,
                labels = FALSE)
  df
}

# Build per-sequence, per-bin data
all_df <- do.call(rbind, lapply(names(seqs), function(nm) {
  df <- seq_to_bins(seqs[[nm]], n_bins)
  df$seq_name <- nm
  df
}))

# Compute frequency per sequence per bin
freq_per_seq <- all_df %>%
  filter(base %in% c("A","C","G","T")) %>%
  group_by(seq_name, bin, base) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(seq_name, bin) %>%
  mutate(freq = count / sum(count)) %>%
  ungroup()

# Compute mean and SD across sequences
stats_df <- freq_per_seq %>%
  group_by(bin, base) %>%
  summarise(
    mean_freq = mean(freq),
    sd_freq = sd(freq),
    .groups = "drop"
  )

stats_df$base <- factor(stats_df$base, levels = c("A","C","G","T"))

# Plot with standard deviation ribbon
#three_prime_utr <- ggplot(stats_df, aes(x = bin, y = mean_freq, color = base, fill = base)) +
#cds <- ggplot(stats_df, aes(x = bin, y = mean_freq, color = base, fill = base)) +
Gypsy_ORFs <- ggplot(stats_df, aes(x = bin, y = mean_freq, color = base, fill = base)) +
  geom_ribbon(aes(ymin = mean_freq - sd_freq,
                  ymax = mean_freq + sd_freq),
              alpha = 0.2, color = NA) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c("A"="red","C"="green","G"="#33CCFF","T"="purple")) +
  scale_fill_manual(values = c("A"="red","C"="green","G"="#33CCFF","T"="purple")) +
  scale_y_continuous(limits = c(0, 0.6)) +
  labs(x="scaled 3UTR region",
       y="Average nucleotide frequency ± SD",
       title="RepBase Drosophila Gypsy 17.6 MDG1 512 ORFs") +
  theme_minimal()

## printing the plot
pdf(file="Drosophila_melanogaster.BDGP6.54.115.non-redundant_three_prime_utr.gt300nt.nucleotide_composition.50bins.pdf",width=8,height=8)
three_prime_utr
dev.off()

pdf(file="Drosophila_melanogaster.BDGP6.54.cds.non_redundant.gt600nt.noSTOP.noSTART.nucleotide_composition.50bins.pdf",width=8,height=8)
cds
dev.off()

pdf(file="RepBase_Drosophila_Gypsy_512ORFs_nucleotides.Gypsy_17.6_412_mdg1.nucleotide_composition.50bins.pdf",width=8,height=8)
Gypsy_ORFs
dev.off()

sessionInfo()
# R version 4.3.1 (2023-06-16)
# Platform: aarch64-apple-darwin20 (64-bit)
# Running under: macOS Ventura 13.7.4
