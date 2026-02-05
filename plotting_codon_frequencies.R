library(ggplot2)
library(tidyverse)

### Dmel mRNAs vs RepBase Drosophila Gypsy elements --- START ---
setwd("${references}")
df <- read.table("Drosophila_melanogaster.BDGP6.54.cds.non_redundant.gt600nt.RepBase_Drosophila_Gypsy_512ORFs.codon_frequencies.txt", header=TRUE)

# Convert to long format
df_long <- df %>%
  pivot_longer(cols = c(mRNA, ERV),
               names_to = "type",
               values_to = "count")

# Convert to frequencies within each amino acid
df_long <- df_long %>%
  group_by(AA, type) %>%
  mutate(freq = count / sum(count))

# Plot — faceted bar chart
per_AA <- ggplot(df_long, aes(x = codon, y = freq, fill = type)) +
  geom_col(position = "dodge") +
  facet_wrap(~AA, scales = "free_x") +
  theme_bw() +
  labs(y = "Codon frequency within AA",
       x = "Codon",
       fill = "Source") +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    strip.text = element_text(face = "bold")
  )

### pdf
pdf(file="RepBase_Drosophila_512Gypsy_17.7_MDG1_ORFs_and_Dmel-mRNAs.codon_frequencies_per_AA.pdf",width=8,height=8)
per_AA
dev.off()



# compute amino-acid frequencies
aa_counts <- df %>%
  group_by(AA) %>%
  summarise(
    mRNA = sum(mRNA),
    ERV  = sum(ERV)
  )

aa_freq <- aa_counts %>%
  mutate(
    mRNA_freq = mRNA / sum(mRNA),
    ERV_freq  = ERV  / sum(ERV)
  )


# compute difference metrics
aa_freq <- aa_freq %>%
  mutate(diff = ERV_freq - mRNA_freq)


# sorted bar chart
usage_difference <- ggplot(aa_freq,
       aes(x = reorder(AA, diff),
           y = diff,
           fill = diff > 0)) +
  geom_col() +
  coord_flip() +
  theme_bw() +
  labs(
    x = "Amino acid",
    y = "ERV - mRNA frequency",
    title = "Amino acid usage difference (ERV vs mRNA)"
  ) +
  scale_fill_manual(values = c("#3b6fb6", "#c0392b"),
                    guide = "none")

### pdf
pdf(file="RepBase_Drosophila_512Gypsy_17.7_MDG1_ORFs_and_Dmel-mRNAs.AA_usage.pdf",width=8,height=8)
usage_difference
dev.off()
