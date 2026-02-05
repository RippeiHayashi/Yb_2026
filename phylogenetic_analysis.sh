### the code was used to infer the phylogenetic relationship of Gypsy POL RT domains.

### software used
### mafft-7.505
https://mafft.cbrc.jp/alignment/software/source.html
### iqtree-1.6.12-Linux
https://iqtree.github.io/release/v1.6.12

references="<directory where the POL RT fasta file is kept>"

mafft --auto --thread 4 ${references}/RepBase_Drosophila_Gypsy.POL_RT.fasta  > ${references}/RepBase_Drosophila_Gypsy.POL_RT.mafft.fasta
iqtree -redo -s ${references}/RepBase_Drosophila_Gypsy.POL_RT.mafft.fasta -m rtREV+R4 -bb 1000 -nt 8
