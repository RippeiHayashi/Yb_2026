library(ggplot2)

### Dmel mRNAs vs RepBase Drosophila Gypsy elements --- START ---
setwd("${references}")

table=read.table("Drosophila_melanogaster.BDGP6.54.cds.non_redundant.gt600nt.RepBase_Drosophila_Gypsy_512ORFs.expected-and-measured_A_bias.txt", header=F)
colnames(table) <- c("CDS","expected_As","observed_As","type")

table_mRNA = subset(table, type == "mRNA")
table_ERV = subset(table, type == "ERV")

Gypsy<-ggplot(table)+
  stat_density_2d(data=subset(table, type == "mRNA"),geom = "polygon",
                  aes(x=`expected_As`,
                      y=`observed_As`),alpha=0.2,colour="gray", fill="gray")+
  geom_vline(data=subset(table, type == "mRNA"),
             aes(xintercept=mean(expected_As)),
             color="gray", linetype="dashed", linewidth=0.5)+
  geom_text(data=table[1,], aes(x=0.25, y=0.2, label=median(table_mRNA$expected_As)), color="gray")+
  geom_hline(data=subset(table, type == "mRNA"),
             aes(yintercept=mean(observed_As)),
             color="gray", linetype="dashed", linewidth=0.5)+
  geom_text(data=table[1,], aes(x=0.2, y=0.25, label=median(table_mRNA$observed_As)), color="gray")+
  stat_density_2d(data=subset(table, type == "ERV"),geom = "polygon",
                  aes(x=`expected_As`,
                      y=`observed_As`),alpha=0.2,colour="red", fill="red")+
  geom_vline(data=subset(table, type == "ERV"),
             aes(xintercept=mean(expected_As)),
             color="red", linetype="dashed", linewidth=0.5)+
  geom_text(data=table[1,], aes(x=0.35, y=0.35, label=median(table_ERV$expected_As)), color="red")+
  geom_hline(data=subset(table, type == "ERV"),
             aes(yintercept=mean(observed_As)),
             color="red", linetype="dashed", linewidth=0.5)+
  geom_text(data=table[1,], aes(x=0.25, y=0.4, label=median(table_ERV$observed_As)), color="red")+
  labs(title="A composition, gray (mRNAs) red (512 GYPSY 17.6 MDG1 ORFs)",
       x="expected_As", y="observed_As")+
  xlim(0.1,0.52)+
  ylim(0.1,0.52)+
  coord_fixed()+
  theme_bw()

### pdf
pdf(file="RepBase_Drosophila_512Gypsy_17.7_MDG1_ORFs_and_Dmel-mRNAs.As.contour.pdf",width=8,height=8)
Gypsy
dev.off()
