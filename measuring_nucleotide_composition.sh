### the code was used to infer the phylogenetic relationship of Gypsy POL RT domains.

### software used
### fastx_toolkit
http://hannonlab.cshl.edu/fastx_toolkit/
### EMBOSS-6.6.0
https://emboss.sourceforge.net/download/


### define the variables
references="<directory where the Ensembl and RepBase sequence files and the gtf file are kept>"

### fetching Drosophila melangoaster protein coding mRNA 3'UTRs and CDSs from Ensembl annotations --- START ---
mkdir -p ${references}/dm6/
wget --directory-prefix="${references}/dm6/" ftp://ftp.ensembl.org/pub/release-115/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.54.115.gtf.gz
wget --directory-prefix="${references}/dm6/" ftp://ftp.ensembl.org/pub/release-115/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.54.dna_sm.toplevel.fa.gz
wget --directory-prefix="${references}/dm6/" ftp://ftp.ensembl.org/pub/release-115/fasta/drosophila_melanogaster/cds/Drosophila_melanogaster.BDGP6.54.cds.all.fa.gz

### decompress the genome fasta file
zcat ${references}/dm6/Drosophila_melanogaster.BDGP6.54.dna_sm.toplevel.fa.gz > ${references}/dm6/Drosophila_melanogaster.BDGP6.54.dna_sm.toplevel.fa

### take only the longest three_prime_utr per gene that are longer than 300nt
zcat ${references}/dm6/Drosophila_melanogaster.BDGP6.54.115.gtf.gz | tr -d '"' | grep "three_prime_utr" |\
awk -F'\t' '{split($9,a,";"); split(a[1],b," "); split(a[2],c," "); if($5-$4+1>300) print $1,$4-1,$5,b[2],c[2],$7,$5-$4+1}' | sort -k7,7nr | awk '!seen[$4]++' |\
awk '{print $1,$2,$3,$4":"$5,".",$6}' | tr ' ' '\t' |\
bedtools getfasta -s -name -fi ${references}/dm6/Drosophila_melanogaster.BDGP6.54.dna_sm.toplevel.fa -bed - > ${references}/dm6/Drosophila_melanogaster.BDGP6.54.115.non-redundant_three_prime_utr.gt300nt.fasta

### converting lowercase to uppercase for downstream analysis
fasta_formatter -i ${references}/dm6/Drosophila_melanogaster.BDGP6.54.115.non-redundant_three_prime_utr.gt300nt.fasta -t |\
awk '{print ">"$1"\n"toupper($NF)}' > ${references}/dm6/Drosophila_melanogaster.BDGP6.54.115.non-redundant_three_prime_utr.gt300nt.mod.fasta
### 5653 3'UTRs

### CDS of >600nt only
### excluding start and stop codons for the line plot
zcat ${references}/dm6/Drosophila_melanogaster.BDGP6.54.cds.all.fa.gz | fasta_formatter - -t |\
awk '{print $1,$4,toupper($NF),length($NF)}' | sort -k4,4nr | awk '!seen[$2]++' |\
awk '{split($2,a,":"); if($4>600) print ">"$1":"a[2]"\n"substr($3,4,length($3)-6)}' > ${references}/dm6/Drosophila_melanogaster.BDGP6.54.cds.non_redundant.gt600nt.noSTOP.noSTART.fasta
### 10835 CDSs
### from start to stop codons for measuring the codon frequencies
zcat ${references}/dm6/Drosophila_melanogaster.BDGP6.54.cds.all.fa.gz | fasta_formatter - -t |\
awk '{print $1,$4,toupper($NF),length($NF)}' | sort -k4,4nr | awk '!seen[$2]++' |\
awk '{split($2,a,":"); if($4>600) print ">"$1":"a[2]"\n"$3}' > ${references}/dm6/Drosophila_melanogaster.BDGP6.54.cds.non_redundant.gt600nt.fasta
### fetching Drosophila melangoaster protein coding mRNA 3'UTRs and CDSs from Ensembl annotations --- END ---


### measure expected and observed A compositions of protein coding sequences of Drosophila melanogaster --- START ---

### step 1. measure codon frequencies
fasta_formatter -i ${references}/dm6/Drosophila_melanogaster.BDGP6.54.cds.non_redundant.gt600nt.fasta -t |\
awk '{for (i=0; i <= length($2)/3-2; i++) {CODON[substr($2,i*3+1,3)]++}} END {for(var in CODON) {if(var != "TAA" && var != "TAG" && var != "TGA") print var,CODON[var]}}' |\
sort -k1,1 | join -1 1 - ${references}/codon_table.sorted.txt > ${references}/dm6/Drosophila_melanogaster.BDGP6.54.gt600nt.codon_frequencies.txt

### ${references}/dm6/Drosophila_melanogaster.BDGP6.54.gt600nt.codon_frequencies.txt looks as follows:
AAA 124427 K
AAC 177446 N
AAG 269730 K
AAT 154814 N
ACA 83576 T
ACC 148348 T
ACG 102370 T
ACT 72550 T
AGA 39834 R
AGC 144501 S
AGG 47590 R
AGT 86745 S
ATA 71580 I
ATC 155347 I
ATG 163814 M
ATT 120310 I
CAA 115817 Q
CAC 111115 H
CAG 253385 Q
CAT 76634 H
CCA 101011 P
CCC 125453 P
CCG 112215 P
CCT 52169 P
CGA 62269 R
CGC 120749 R
CGG 59370 R
CGT 60389 R
CTA 61827 L
CTC 97661 L
CTG 267866 L
CTT 66915 L
GAA 159445 E
GAC 167680 D
GAG 296543 E
GAT 198258 D
GCA 92288 A
GCC 226144 A
GCG 97984 A
GCT 101444 A
GGA 124524 G
GGC 176481 G
GGG 33502 G
GGT 91430 G
GTA 46946 V
GTC 94165 V
GTG 194019 V
GTT 79919 V
TAC 128228 Y
TAT 80934 Y
TCA 59347 S
TCC 138119 S
TCG 115878 S
TCT 52312 S
TGC 92489 C
TGG 70308 W
TGT 40389 C
TTA 34913 F
TTC 150805 F
TTG 118606 F
TTT 99511 F

### step 2. calculate nucleotide frequencies per amino acids using the codon frequency table
### step 2-1. print out codons per AA
cat ${references}/dm6/Drosophila_melanogaster.BDGP6.54.gt600nt.codon_frequencies.txt | while read CODON NUM AA; do
for i in $(seq 1 ${NUM}); do
printf ${CODON} >> ${references}/dm6/Drosophila_melanogaster.BDGP6.54_${AA}.txt
done
done

for NUC in A T G C; do
### step 2-2. calculate A_bias per amino acid
cat ${references}/dm6/Drosophila_melanogaster.BDGP6.54.gt600nt.codon_frequencies.txt | awk '!seen[$3]++' | while read CODON NUM AA; do
awk -v AA=${AA} -v NUC=${NUC} '{print AA,gsub(NUC,NUC,$1)/length($1)}' ${references}/dm6/Drosophila_melanogaster.BDGP6.54_${AA}.txt >> ${references}/dm6/Drosophila_melanogaster.BDGP6.54.codon_frequencies.${NUC}_bias.txt
done

cat ${references}/dm6/Drosophila_melanogaster.BDGP6.54.codon_frequencies.${NUC}_bias.txt | sort -k1,1 > ${references}/dm6/Drosophila_melanogaster.BDGP6.54.codon_frequencies.${NUC}_bias.sorted.txt
done

### step 3. print out amino acid sequences using the CDS
cat ${references}/dm6/Drosophila_melanogaster.BDGP6.54.cds.non_redundant.gt600nt.fasta | tr ':' '@' > ${references}/dm6/Drosophila_melanogaster.BDGP6.54.cds.non_redundant.gt600nt.mod.fasta
transeq ${references}/dm6/Drosophila_melanogaster.BDGP6.54.cds.non_redundant.gt600nt.mod.fasta ${references}/dm6/Drosophila_melanogaster.BDGP6.54.cds.non_redundant.gt600nt.pep
fasta_formatter -i ${references}/dm6/Drosophila_melanogaster.BDGP6.54.cds.non_redundant.gt600nt.pep -t | tr -d "*" |\
awk '{split($1,a,"_"); print a[1],$NF}' > ${references}/dm6/Drosophila_melanogaster.BDGP6.54.cds.non_redundant.gt600nt.pep2

### step 4. measure expected A_bias per peptides
NUC="A"
cat ${references}/dm6/Drosophila_melanogaster.BDGP6.54.cds.non_redundant.gt600nt.pep2 | while read PEP AA; do
cat ${references}/dm6/Drosophila_melanogaster.BDGP6.54.cds.non_redundant.gt600nt.pep2 | grep ${PEP} |\
awk '{for (i=1; i <= length($NF); i++) {CODON[substr($2,i,1)]++}} END {for(var in CODON) print var,CODON[var]}' |\
sort -k1,1 | join -1 1 - ${references}/dm6/Drosophila_melanogaster.BDGP6.54.codon_frequencies.${NUC}_bias.sorted.txt |\
awk -v PEP=${PEP} '{SIZE+=$2; AA+=$2*$3} END {print PEP,AA/SIZE}' >> ${references}/dm6/Drosophila_melanogaster.BDGP6.54.cds.non_redundant.gt600nt.expected_${NUC}_bias.txt
done
cat ${references}/dm6/Drosophila_melanogaster.BDGP6.54.cds.non_redundant.gt600nt.expected_${NUC}_bias.txt | sort -k1,1 > ${references}/dm6/Drosophila_melanogaster.BDGP6.54.cds.non_redundant.gt600nt.expected_${NUC}_bias.sorted.txt

### step 5. measure observed A_bias per peptides excluding the last three bases
fasta_formatter -i ${references}/dm6/Drosophila_melanogaster.BDGP6.54.cds.non_redundant.gt600nt.mod.fasta -t |\
awk '{print $1,substr($NF,1,length($NF)-3)}' |\
awk -v NUC=${NUC} '{print $1,gsub(NUC,NUC,$2)/length($2)}' |\
sort -k1,1 > ${references}/dm6/Drosophila_melanogaster.BDGP6.54.cds.non_redundant.gt600nt.observed_${NUC}_bias.txt

join -1 1 ${references}/dm6/Drosophila_melanogaster.BDGP6.54.cds.non_redundant.gt600nt.expected_${NUC}_bias.sorted.txt ${references}/dm6/Drosophila_melanogaster.BDGP6.54.cds.non_redundant.gt600nt.observed_${NUC}_bias.txt > ${references}/dm6/Drosophila_melanogaster.BDGP6.54.cds.non_redundant.gt600nt.expected-and-observed_${NUC}_bias.txt

### ${references}/dm6/Drosophila_melanogaster.BDGP6.54.cds.non_redundant.gt600nt.expected-and-observed_${NUC}_bias.txt looks as follows:
<CDS_name> <expected A frequency> <observed A frequency>
FBtr0006151@FBgn0000056 0.255724 0.269608
FBtr0070000@FBgn0031081 0.289947 0.263783
FBtr0070002@FBgn0031085 0.233951 0.209524
FBtr0070003@FBgn0062565 0.222633 0.197244
FBtr0070006@FBgn0031089 0.210012 0.177627
FBtr0070007@FBgn0031092 0.245445 0.21047
...
### measure expected and observed A compositions of protein coding sequences of Drosophila melanogaster --- END ---


### measure expected and observed A compositions of 512 GYPSY, 17_6 and MDG1 ORFs --- START ---
### step 1. measure expected A_bias per peptides
NUC="A"
fasta_formatter -i ${references}/RepBase_Drosophila/RepBase_Drosophila_Gypsy_512ORFs_AA-sequences.GYPSY_17.6_MDG1.curated.fasta -t | while read PEP AA; do
fasta_formatter -i ${references}/RepBase_Drosophila/RepBase_Drosophila_Gypsy_512ORFs_AA-sequences.GYPSY_17.6_MDG1.curated.fasta -t | grep ${PEP} |\
awk '{for (i=1; i <= length($NF); i++) {if(substr($2,i,1) != "X") CODON[substr($2,i,1)]++}} END {for(var in CODON) print var,CODON[var]}' |\
sort -k1,1 | join -1 1 - ${references}/dm6/Drosophila_melanogaster.BDGP6.54.codon_frequencies.${NUC}_bias.sorted.txt |\
awk -v PEP=${PEP} '{SIZE+=$2; AA+=$2*$3} END {print PEP,AA/SIZE}' >> ${references}/RepBase_Drosophila/RepBase_Drosophila_Gypsy_512ORFs_AA-sequences.GYPSY_17.6_MDG1.curated.expected_${NUC}_bias.txt
done
cat ${references}/RepBase_Drosophila/RepBase_Drosophila_Gypsy_512ORFs_AA-sequences.GYPSY_17.6_MDG1.curated.expected_${NUC}_bias.txt | sort -k1,1 > ${references}/RepBase_Drosophila/RepBase_Drosophila_Gypsy_512ORFs_AA-sequences.GYPSY_17.6_MDG1.curated.expected_${NUC}_bias.sorted.txt

### step 2. measure observed A_bias
fasta_formatter -i ${references}/RepBase_Drosophila/RepBase_Drosophila_Gypsy_512ORFs_nucleotide-sequences.GYPSY_17.6_MDG1.curated.fasta -t |\
awk '{print $1,toupper($2)}' |\
awk -v NUC=${NUC} '{print $1,gsub(NUC,NUC,$2)/length($2)}' | sort -k1,1 > ${references}/RepBase_Drosophila/RepBase_Drosophila_Gypsy_512ORFs_AA-sequences.GYPSY_17.6_MDG1.curated.observed_${NUC}_bias.txt

join -1 1 ${references}/RepBase_Drosophila/RepBase_Drosophila_Gypsy_512ORFs_AA-sequences.GYPSY_17.6_MDG1.curated.expected_${NUC}_bias.sorted.txt ${references}/RepBase_Drosophila/RepBase_Drosophila_Gypsy_512ORFs_AA-sequences.GYPSY_17.6_MDG1.curated.observed_${NUC}_bias.txt > ${references}/RepBase_Drosophila/RepBase_Drosophila_Gypsy_512ORFs_AA-sequences.GYPSY_17.6_MDG1.curated.expected-and-observed_${NUC}_bias.txt
### measure expected and observed A compositions of 512 GYPSY, 17_6 and MDG1 ORFs --- END ---

### combine host ORFs and ERV ORFs
NUC="A"
awk '{print $0,"mRNA"}' ${references}/dm6/Drosophila_melanogaster.BDGP6.54.cds.non_redundant.gt600nt.expected-and-observed_${NUC}_bias.txt >> ${references}/Drosophila_melanogaster.BDGP6.54.cds.non_redundant.gt600nt.RepBase_Drosophila_Gypsy_512ORFs.expected-and-measured_${NUC}_bias.txt
awk '{print $0,"ERV"}' ${references}/RepBase_Drosophila/RepBase_Drosophila_Gypsy_512ORFs_AA-sequences.GYPSY_17.6_MDG1.curated.expected-and-observed_${NUC}_bias.txt >> ${references}/Drosophila_melanogaster.BDGP6.54.cds.non_redundant.gt600nt.RepBase_Drosophila_Gypsy_512ORFs.expected-and-measured_${NUC}_bias.txt
### use expected_vs_observed_contour.R to make a contour plot

### measure codon frequencies of 512 GYPSY, 17_6 and MDG1 ORFs and combine them to those of mRNAs
fasta_formatter -i ${references}/RepBase_Drosophila/RepBase_Drosophila_Gypsy_512ORFs_nucleotide-sequences.GYPSY_17.6_MDG1.curated.fasta -t |\
awk '{for (i=0; i <= length($2)/3-2; i++) {CODON[substr($2,i*3+1,3)]++}} END {for(var in CODON) {if(var != "TAA" && var != "TAG" && var != "TGA") print var,CODON[var]}}' |\
sort -k1,1 | join -1 1 - ${references}/codon_table.sorted.txt > ${references}/RepBase_Drosophila/RepBase_Drosophila_Gypsy_512ORFs_nucleotide-sequences.GYPSY_17.6_MDG1.curated.codon_frequencies.txt

join -1 1 ${references}/dm6/Drosophila_melanogaster.BDGP6.54.gt600nt.codon_frequencies.txt ${references}/RepBase_Drosophila/RepBase_Drosophila_Gypsy_512ORFs_nucleotide-sequences.GYPSY_17.6_MDG1.curated.codon_frequencies.txt |\
awk 'BEGIN{print "codon AA mRNA ERV"} {print $1,$3,$2,$4}' > ${references}/Drosophila_melanogaster.BDGP6.54.cds.non_redundant.gt600nt.RepBase_Drosophila_Gypsy_512ORFs.codon_frequencies.txt

### ${references}/Drosophila_melanogaster.BDGP6.54.cds.non_redundant.gt600nt.RepBase_Drosophila_Gypsy_512ORFs.codon_frequencies.txt looks as follows:
codon AA mRNA ERV
AAA K 124427 20599
AAC N 177446 11318
AAG K 269730 7939
AAT N 154814 14630
ACA T 83576 9277
ACC T 148348 5401
ACG T 102370 3120
ACT T 72550 5434
...
### use plotting_codon_frequencies.R to make barcharts
