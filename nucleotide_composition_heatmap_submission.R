library(ggplot2)
library(dplyr)
library(tidyr)

setwd("${references}/dm6")

### Gypsy whole sequence internal regions --- START ---
df=read.table("RepBase_Drosophila_542.GYPSY_SPECIES_POL_RT_INT_As_Cs_Gs_Ts_TYPE.txt", header=F)
colnames(df) <- c("name","SPECIES","POL","RT","INT","As","Cs","Gs","Ts","TYPE")

# Long format + frequencies
df_long <- df %>%
  pivot_longer(cols = c(As, Cs, Gs, Ts),
               names_to = "Nucleotide",
               values_to = "Count") %>%
  group_by(name) %>%
  mutate(Freq = Count / sum(Count)) %>%
  ungroup()

### setting the colour scheme gradient
brbg <- RColorBrewer::brewer.pal(11, "BrBG")

# Heatmap grouped by TYPE
p<-ggplot(df_long, aes(x = name, y = Nucleotide, fill = Freq)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colours = brbg[11:1], limits = c(0, 0.5), oob = scales::squish, values = scales::rescale(seq(0, 0.5, length.out = 11))) +
  facet_wrap(~ TYPE, scales = "free_x") +
  labs(
    x = "Element",
    y = "Base",
    fill = "Frequency",
    title = "Nucleotide Composition (whole internal sequence w/o LTRs) by Type (all Drosophila)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold")
  )

pdf(file="RepBase_Drosophila_542.GYPSY_whole_noLTRs.nucleotide-compositions_by_GYPSY-TYPE.ver2.pdf",width=15,height=8)
p
dev.off()
### Gypsy whole sequence internal regions --- END ---

sessionInfo()
# R version 4.3.1 (2023-06-16)
# Platform: aarch64-apple-darwin20 (64-bit)
# Running under: macOS Ventura 13.7.4
