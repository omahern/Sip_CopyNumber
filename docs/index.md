---
title: "Trophic Cascades - Copy Number Comparison"
output: 
  html_document:
    keep_md: true
    toc: true
    toc_float: true
    toc_depth: 4
    code_folding: hide
    number_sections: true
    theme: cosmo

knit: (function(input_file, encoding) {
  out_dir <- 'docs';
  rmarkdown::render(input_file,
 encoding=encoding,
 output_file=file.path(dirname(input_file), out_dir, 'index.html'))})
---


# Workspace Setup


```r
library(phyloseq)
library(reshape2)
library(tidyverse)
library(vegan)
library(HTSSIP)
library(ape)
library(CoDaSeq)
library(philr)
library(ggtree)
library(cowplot)
library(ggplot2)
library(viridis)
library(treeio)
library(dplyr)
library(microbiome)
library(treeio)
```



# code for 16S sequencing data

```
#!/bin/bash
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --mem=20gb
#SBATCH --cpus-per-task=20
#SBATCH --time=04:35:00

module load anaconda/5.1
source activate qiime2-2018.8

echo merging and clustering
qiime dada2 denoise-paired \
    --i-demultiplexed-seqs paired-end-demux.qza \
    --p-n-threads 20 \
    --p-trunc-q 2 \
    --p-trim-left-f 19 \
    --p-trim-left-r 20 \
    --p-trunc-len-f 200 \
    --p-trunc-len-r 160 \
    --p-max-ee 5 \
    --p-n-reads-learn 1000000 \
    --p-chimera-method consensus \
    --o-table asv_table.qza \
    --o-representative-sequences rep-seqs.qza \
    --o-denoising-stats ASVs3/stats-dada2.qza

echo classifying taxa
qiime feature-classifier classify-sklearn \
    --i-classifier ../21Oct21/tagseq-qiime2-snakemake/silva_all.qza \
    --i-reads rep-seqs.qza \
    --o-classification asv_tax_sklearn.qza


qiime tools export --input-path asv_table.qza \
    --output-path asv_table

biom convert -i asv_table/feature-table.biom -o asv_table/asv-table.tsv --to-tsv

qiime tools export  --input-path asv_tax_sklearn.qza --output-path asv_tax_dir

qiime phylogeny align-to-tree-mafft-fasttree \
    --i-sequences rep-seqs.qza \
    --output-dir mafft-fasttree-output

```

# State Data

## Chemostat Stability



### SI Figure 1


```r
setwd("/Users/oliviaahern/Documents/GitHub/SIP_CopyNumber")


{
  
# file contains info on MC2, Chemostat, measurements taken every 15 minutes
data=read.csv(file="datafiles/m2_do_02_edit.csv",header=T)
# subset to get the right time frame
read=subset(data, T<275)
time=read$T
ph= read$pH.M2
do=read$DO.M2..uM.
heado=read$O2.M2....
heado2=read$CO2.M2....

  par(mar=c(2,7,1,5), mfrow=c(3,1),xpd=F)
  plot(time, do,type='l',col='blue',xlab=" ", ylab = "", xlim=c(-450,240), xaxt='n')
  rect(xleft=0, xright=240, ybottom=111, ytop=300, col= rgb(0.7,0.7,0.7,alpha=0.2),
    border=NA)
  mtext(side=2,"Dissolved Oxygen (uM)", col = 'blue', padj=-3)
    axis(side=1, at=c(-450,-144, 0, 24, 48, 72, 120, 240))
  par(new = TRUE)
  plot(time, ph,col='red',type = "l", axes = FALSE, bty = "n", xlab = "", ylab = "")
  mtext(side=4, "pH", col ='red', line=3)
  axis(side=4, at = c(6.5, 6.7, 7, 7.2, 7.4))
  mtext("a", adj=-.10,padj=-10,side=1,srt=-90,cex=1.5,font=1)
  legend(25,7.4,legend=c("12C Glucose", "13C Glucose"), pch=22, pt.bg=c("white",'gray70'), bt='n')
  
  
  
  plot(time, heado,type='l',col='navy',xlab=" ", ylab = "", xlim=c(-450,240), xaxt='n')
  rect(xleft=0, xright=240, ybottom=0, ytop=300, col= rgb(0.7,0.7,0.7,alpha=0.2),border=NA)
  mtext(side=2,"Oxygen (%)", col = 'navy', padj=-3)
  axis(side=1, at=c(-450,-144, 0, 24,  48, 72, 120, 240))
  par(new = TRUE)
  plot(time, heado2,col='forestgreen',type = "l", axes = FALSE, bty = "n", xlab = "", ylab = "")
  mtext(side=4, "CO2 (%)", col ='forestgreen', line=3)
  axis(side=4, at = c(0,0.1,0.2,0.3,0.4))
  mtext("b", adj=-.1,padj=-10,side=1,srt=-90,cex=1.5,font=1)
 
# file contains info on MC2, Chemostat, data was taken daily
read=read.csv(file="datafiles/states_gluc.csv",header=T)
glucose=read$Glucose_um
time_glucose=read$Time
  
   par(mar=c(5,30,1,5))
  plot(read$Time,read$Glucose_um,type='o',ylim=c(0,12),xlab= " ",
       ylab = " ",col='gray50',bg='gray50',pch=21,cex.axis=1,cex.lab=1,cex=1.2,lwd=1.4,
       yaxt='n', xaxt='n') #xaxt=c(-144,0,48,120,240))
  rect(xleft=0, xright=240, ybottom=0, ytop=12, col= rgb(0.7,0.7,0.7,alpha=0.2),
    border=NA)
  mtext(side=2, "Glucose (um)", padj=-3, col="gray50")
  axis(side=2,at=c(0,4,8,12), cex=1.3)
  axis(side=1, at=c(-144, 0, 24,  48, 72, 120, 240))
  par(new = TRUE)
  plot(read$Time,read$Atom..13C,pch=23, 
       bg='black',type = "o", axes = FALSE, bty = "n", xlab = "", ylab = "", ylim=c(0,30), cex=1.2)
  mtext(side=4, "Atomic 13C (%)", line=3)
  axis(side=4, at=c(0, 10,20,30))
  mtext(side=1, "Time (hours)", line=3)
 # legend('topleft',legend=c("Glucose (uM)","Atomic 13C"), pch =c(21,23), pt.bg= c('gray70','black'), bty='n')
  mtext("c", adj=-.15,padj=-9,side=1,srt=-90,cex=1.5,font=1)
}
```

![](/Users/oliviaahern/Documents/GitHub/Sip_CopyNumber/docs/index_files/figure-html/si-fig-1-1.png)<!-- -->

## qPCR Buoyant Density Gradients

### SI Figure 2


```r
data<-read.csv("/Users/oliviaahern/Documents/MBL_WHOI/Trophic_Cascades/Experiment3/sequencing_14Jan26/edit_qc/exp3_bd.csv",header=T,row.names=1)
data=data.frame(data)
par(mfrow=c(1,3),
    mar=c(4,4,1,1))
plot(
     subset(data, Culture== "Batch" & Substrate =="Single" & Isotope=="C12")$R_MAX,
     subset(data, Culture== "Batch" & Substrate =="Single" & Isotope=="C12")$BD, 
     type='o', pch=21, col='gray70', 
     bg=subset(data, Culture== "Batch" & Substrate =="Single" & Isotope=="C12")$col,
     ylab="Buoyant Density (g/mL)",
     xlab="Ratio of Maximum Quantity",
     ylim=c(1.71,1.84),
     xaxt='n')
axis(1, at=c(0,0.5,1))
lines(
     subset(data, Culture== "Batch" & Substrate =="Single" & Isotope=="C13")$R_MAX,
     subset(data, Culture== "Batch" & Substrate =="Single" & Isotope=="C13")$BD,
      type='o', pch=21, col='black', 
          bg=subset(data, Culture== "Batch" & Substrate =="Single" & Isotope=="C13")$col)



plot( 
     subset(data, Culture== "Batch" & Substrate =="Multi" & Isotope=="C12")$R_MAX,
     subset(data, Culture== "Batch" & Substrate =="Multi" & Isotope=="C12")$BD,
     type='o', pch=21, col='gray70', 
bg=subset(data, Culture== "Batch" & Substrate =="Multi" & Isotope=="C12")$col,
ylab="",
xlab="Ratio of Maximum Quantity",
     ylim=c(1.71,1.84),
     xaxt='n')
axis(1, at=c(0,0.5,1))
lines(
     subset(data, Culture== "Batch" & Substrate =="Multi" & Isotope=="C13")$R_MAX,
     subset(data, Culture== "Batch" & Substrate =="Multi" & Isotope=="C13")$BD,
      type='o', pch=21, col='black', bg=subset(data, Culture== "Batch" & Substrate =="Multi" & Isotope=="C13")$col)


plot(
     subset(data, Culture== "Chemostat" & Substrate =="Multi" & Isotope=="C12")$R_MAX,
     subset(data, Culture== "Chemostat" & Substrate =="Multi" & Isotope=="C12")$BD, 
     type='o', pch=21, col='gray70', 
     bg=subset(data, Culture== "Chemostat" & Substrate =="Multi" & Isotope=="C12")$col,
     ylab="",xlab="Ratio of Maximum Quantity",
     ylim=c(1.71,1.84),
     xaxt='n')
axis(1, at=c(0,0.5,1))
lines(
     subset(data, Culture== "Chemostat" & Substrate =="Multi" & Isotope=="C13")$R_MAX,
     subset(data, Culture== "Chemostat" & Substrate =="Multi" & Isotope=="C13")$BD,
      type='o', pch=21, col='black',
     bg=subset(data, Culture== "Chemostat" & Substrate =="Multi" & Isotope=="C13")$col)
```

![](/Users/oliviaahern/Documents/GitHub/Sip_CopyNumber/docs/index_files/figure-html/si-fig2a-1.png)<!-- -->


```r
data<-read.csv("datafiles/exp3_bd.csv",header=T,row.names=1)
qPCR_data=data.frame(data)

desired_order=c("Single", "Multi")
qPCR_data$Substrate <- factor(qPCR_data$Substrate, levels = desired_order)



ggplot(data=qPCR_data, aes(x=BD, y=R_MAX, col=Isotope))+
  geom_point() +
  geom_line() +
  theme_bw() +
    facet_wrap(Culture ~ Substrate)  +
  coord_flip() +
  scale_color_manual(values=c("gray","black")) + 
  xlab("Buoyant Density (g/mL)") +
  ylab("Ratio of Maximum Quantity")
```

![](/Users/oliviaahern/Documents/GitHub/Sip_CopyNumber/docs/index_files/figure-html/si-fig2b-1.png)<!-- -->


# Import 16S Data


```r
numbs=read.csv(file='datafiles/copy_number.csv',
               header=T)
number=numbs$x
x<-read.csv(file='datafiles/asv-table.csv',
            header=TRUE,row.names=1)
x2=sweep(x, 1, number, "/") # divide asvs by 16S gene copy number
x2=ceiling(x2) # round it

OTU = otu_table(x2, taxa_are_rows=T)
taxa<-read.csv(file='datafiles/taxonomy.csv',
               header=TRUE,row.names=1)
t<-as.matrix(taxa)
tax2<-tax_table(t)
map<-import_qiime_sample_data("datafiles/sampling_batch_exp3_combined.txt")
tree=read.newick("datafiles/tree.nwk")
phyo1 = phyloseq(OTU, tax2,map,tree)
 # 1238 ASVs
phyo1 = subset_taxa(phyo1, !Order=="Chloroplast")
 #  1125, removed 113 that were chloroplasts
phyo1 = subset_taxa(phyo1, !Family=="Mitochondria")
 # 1104, removed 21 mitochondria
phyo1=subset_taxa(phyo1, !Phylum =="Cyanobacteria")
 # 1089, removed 15 cyanobacteria
phyo1=subset_taxa(phyo1, !Phylum =="")
 # 1078, removed 19 no assignment to phyla level

phyo_frac=subset_samples(phyo1, Isotope !="NA")
phyo_frac <- prune_samples(sample_sums(phyo_frac) > 0, phyo_frac)
```


# Code for running qSIP

## 1-normalize reads to qPCR copy number concentrations


```r
phyo_frac=subset_samples(phyo1, Isotope !="NA")
ASVs_to_keep <- taxa_sums(phyo_frac) > 0
phyo_frac <- prune_taxa(ASVs_to_keep, phyo_frac)


map<-read.csv("datafiles/exp3_qpcr_combined.csv",
              header=T,row.names=1)
physeq_rep3_t = OTU_qPCR_trans(phyo_frac, map)
```

## 2-batch-single


```r
batch_glu=
  subset_samples(physeq_rep3_t, Culture=="Batch" & Treatment=="Single")

otu=otu_table(batch_glu)
deco=decostand(otu, 'pa')
frac2=merge_phyloseq(deco, sample_data(batch_glu), phy_tree(batch_glu), tax_table(batch_glu))
## Subset by tfrac1## Subset by treatment
ps1 <- prune_samples(frac2@sam_data$Isotope == "C12", frac2)
ps1 <- filter_taxa(ps1, function(x) sum(x) > 1, prune = TRUE)
ps2 <- prune_samples(frac2@sam_data$Isotope == "C13", frac2)
ps2 <- filter_taxa(ps2, function(x) sum(x) > 1, prune = TRUE)
## Get vectors of numbered OTUs/ASVs
treatment1 <- rownames(otu_table(ps1))
treatment2 <- rownames(otu_table(ps2))
## Get the intersection
shared <- intersect(treatment1, treatment2) # 97 ASVs shared between the C12 and C13 treatment
frac1=batch_glu
## Subset phyloseq object to shared taxa
ps.s <- subset(otu_table(frac1), rownames(otu_table(frac1)) %in% shared)
merged=merge_phyloseq(ps.s, tax_table(frac1), sample_data(frac1), phy_tree(frac1)) # 53 taxa
merged
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 98 taxa and 15 samples ]
## sample_data() Sample Data:       [ 15 samples by 13 sample variables ]
## tax_table()   Taxonomy Table:    [ 98 taxa by 8 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 98 tips and 97 internal nodes ]
```

```r
### run qSIP
atomX = qSIP_atom_excess(merged, control_expr='Isotope=="C12"', treatment_rep='Isotope')
df_atomX_boot = qSIP_bootstrap(atomX, n_boot=3,
                               isotope="13C")
df_atomX_boot2=na.omit(df_atomX_boot)
CI_threshold = 0.011
df_atomX_boot = df_atomX_boot %>%
  mutate(Incorporator = A_CI_low > CI_threshold,
         OTU = reorder(OTU, -A))
n_incorp = df_atomX_boot %>%
  filter(Incorporator == TRUE) %>%
  nrow
cat('Number of incorporators:', n_incorp, '\n')
```

```
## Number of incorporators: 84
```

```r
## 83 incorp
df_dBD = delta_BD(merged, control_expr='Isotope=="C12"',
                  n=12)
stopifnot(nrow(df_atomX_boot) == nrow(df_dBD))
df_j = dplyr::inner_join(df_atomX_boot, df_dBD, c('OTU'='OTU'))
stopifnot(nrow(df_atomX_boot) == nrow(df_j))
df_j = df_j %>%
  dplyr::mutate(OTU = reorder(OTU, -delta_BD))
batch_glu=df_j

batch_glui=
  subset(batch_glu, Incorporator==TRUE)
```


## 3-batch-multi


```r
batch_carb=
  subset_samples(physeq_rep3_t, Culture=="Batch" & Treatment=="Multi")

otu=otu_table(batch_carb)
deco=decostand(otu, 'pa')
frac2=merge_phyloseq(deco, sample_data(batch_carb), phy_tree(batch_carb), tax_table(batch_carb))
## Subset by tfrac1## Subset by treatment
ps1 <- prune_samples(frac2@sam_data$Isotope == "C12", frac2)
ps1 <- filter_taxa(ps1, function(x) sum(x) > 1, prune = TRUE)
ps2 <- prune_samples(frac2@sam_data$Isotope == "C13", frac2)
ps2 <- filter_taxa(ps2, function(x) sum(x) > 1, prune = TRUE)
## Get vectors of numbered OTUs/ASVs
treatment1 <- rownames(otu_table(ps1))
treatment2 <- rownames(otu_table(ps2))
## Get the intersection
shared <- intersect(treatment1, treatment2) # 64 ASVs shared between the C12 and C13 treatment
frac1=batch_carb
## Subset phyloseq object to shared taxa
ps.s <- subset(otu_table(frac1), rownames(otu_table(frac1)) %in% shared)
merged=merge_phyloseq(ps.s, tax_table(frac1), sample_data(frac1), phy_tree(frac1)) # 53 taxa
merged
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 63 taxa and 12 samples ]
## sample_data() Sample Data:       [ 12 samples by 13 sample variables ]
## tax_table()   Taxonomy Table:    [ 63 taxa by 8 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 63 tips and 62 internal nodes ]
```

```r
### run qSIP
atomX = qSIP_atom_excess(merged, control_expr='Isotope=="C12"', treatment_rep='Isotope')
df_atomX_boot = qSIP_bootstrap(atomX, n_boot=1,isotope="13C")
df_atomX_boot2=na.omit(df_atomX_boot)
CI_threshold = 0.011
df_atomX_boot = df_atomX_boot %>%
  mutate(Incorporator = A_CI_low > CI_threshold,
         OTU = reorder(OTU, -A))
n_incorp = df_atomX_boot %>%
  filter(Incorporator == TRUE) %>%
  nrow
cat('Number of incorporators:', n_incorp, '\n')
```

```
## Number of incorporators: 60
```

```r
## 61 incorp
df_dBD = delta_BD(merged, control_expr='Isotope=="C12"',
                  n=12)
stopifnot(nrow(df_atomX_boot) == nrow(df_dBD))
df_j = dplyr::inner_join(df_atomX_boot, df_dBD, c('OTU'='OTU'))
stopifnot(nrow(df_atomX_boot) == nrow(df_j))
df_j = df_j %>%
  dplyr::mutate(OTU = reorder(OTU, -delta_BD))
batch_carb=df_j

batch_carbi=
  subset(batch_carb, Incorporator==TRUE)
```

## 4-chemostat-multi


```r
chemo=
  subset_samples(physeq_rep3_t, Culture=="Chemostat" & Treatment=="Multi")

otu=otu_table(chemo)
deco=decostand(otu, 'pa')
frac2=merge_phyloseq(deco, sample_data(chemo), phy_tree(chemo), tax_table(chemo))
## Subset by tfrac1## Subset by treatment
ps1 <- prune_samples(frac2@sam_data$Isotope == "C12", frac2)
ps1 <- filter_taxa(ps1, function(x) sum(x) > 1, prune = TRUE)
ps2 <- prune_samples(frac2@sam_data$Isotope == "C13", frac2)
ps2 <- filter_taxa(ps2, function(x) sum(x) > 1, prune = TRUE)
## Get vectors of numbered OTUs/ASVs
treatment1 <- rownames(otu_table(ps1))
treatment2 <- rownames(otu_table(ps2))
## Get the intersection
shared <- intersect(treatment1, treatment2) # 78 ASVs shared between the C12 and C13 treatment
frac1=chemo
## Subset phyloseq object to shared taxa
ps.s <- subset(otu_table(frac1), rownames(otu_table(frac1)) %in% shared)
merged=merge_phyloseq(ps.s, tax_table(frac1), sample_data(frac1), phy_tree(frac1)) # 53 taxa
merged
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 77 taxa and 23 samples ]
## sample_data() Sample Data:       [ 23 samples by 13 sample variables ]
## tax_table()   Taxonomy Table:    [ 77 taxa by 8 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 77 tips and 76 internal nodes ]
```

```r
### run qSIP
atomX = qSIP_atom_excess(merged, control_expr='Isotope=="C12"', treatment_rep='Isotope')
df_atomX_boot = qSIP_bootstrap(atomX, n_boot=1,isotope="13C")
df_atomX_boot2=na.omit(df_atomX_boot)
CI_threshold = 0.011
df_atomX_boot = df_atomX_boot %>%
  mutate(Incorporator = A_CI_low > CI_threshold,
         OTU = reorder(OTU, -A))
n_incorp = df_atomX_boot %>%
  filter(Incorporator == TRUE) %>%
  nrow
cat('Number of incorporators:', n_incorp, '\n')
```

```
## Number of incorporators: 53
```

```r
## 52 incorp
df_dBD = delta_BD(merged, control_expr='Isotope=="C12"',
                  n=12)
stopifnot(nrow(df_atomX_boot) == nrow(df_dBD))
df_j = dplyr::inner_join(df_atomX_boot, df_dBD, c('OTU'='OTU'))
stopifnot(nrow(df_atomX_boot) == nrow(df_j))
df_j = df_j %>%
  dplyr::mutate(OTU = reorder(OTU, -delta_BD))
chemo_carb=df_j

chemo_carbi=
  subset(chemo_carb, Incorporator==TRUE)
```


## 5-combine-qSIP-one-file


```r
numbs=read.csv(file='datafiles/copy_number.csv',header=T)

batch_glu$Culture=c(rep("Batch",98))
batch_glu$Substrate=c(rep("Single",98))
batch_glu$copy_number=(subset(numbs, label %in% batch_glu$OTU))$x
batch_glu$id=c(rep("batch_glu",98))

batch_carb$Culture=c(rep("Batch",63))
batch_carb$Substrate=c(rep("Multi", 63))
batch_carb$copy_number=(subset(numbs, label %in% batch_carb$OTU))$x
batch_carb$id=c(rep("batch_carb",63))

chemo_carb$Culture=c(rep("Chemostat",77))
chemo_carb$Substrate=c(rep("Multi",77))
chemo_carb$copy_number=(subset(numbs, label %in% chemo_carb$OTU))$x
chemo_carb$id=c(rep("chemo_carb",77))


s1=bind_rows(batch_glu, batch_carb)
all=bind_rows(s1, chemo_carb)
taxa=tax_table(phyo1)
colnames(taxa)
```

```
## [1] "Domain"  "Phylum"  "Class"   "Order"   "Family"  "Genus"   "Species"
## [8] "Strain"
```

```r
colnames(taxa)=c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "OTU")
qsip_data=merge(all, taxa, by="OTU")

write.csv(qsip_data, 'datafiles/qsip_data_exp3.csv')
```


# qSIP basic stats

## 16S rRNA copy number among treatments
looking at all ASVs that were in the buoyant density gradients (must be present at in at least 2 fractions of both enriched and unenriched) and just the incorporators 


```r
par(mfrow=c(2,2))

boxplot(qsip_data$copy_number~qsip_data$id, xlab=" ", ylab="16S rRNA copy number", main = "16S copy number vs. treatments, all ASVs")
a=aov(qsip_data$copy_number~qsip_data$id)
summary(a)
```

```
##               Df Sum Sq Mean Sq F value Pr(>F)
## qsip_data$id   2    2.0   1.011   0.404  0.668
## Residuals    235  588.4   2.504
```

```r
tuk=TukeyHSD(a)
plot(tuk,las=2)

incorp_only = subset(qsip_data, Incorporator == "TRUE")

boxplot(incorp_only$copy_number~incorp_only$id, xlab=" ", ylab="16S rRNA copy number", main = "16S copy number vs. treatments, only incorp")
a=aov(incorp_only$copy_number~incorp_only$id)
summary(a)
```

```
##                 Df Sum Sq Mean Sq F value Pr(>F)
## incorp_only$id   2    1.4  0.6989   0.256  0.774
## Residuals      194  528.8  2.7259
```

```r
tuk=TukeyHSD(a)
plot(tuk,las=2)
```

![](/Users/oliviaahern/Documents/GitHub/Sip_CopyNumber/docs/index_files/figure-html/qSIP-1-1.png)<!-- -->


## EAF among treatments 
looking at all ASVs that were in the buoyant density gradients (must be present at in at least 2 fractions of both enriched and unenriched) and just the incorporators 


```r
par(mfrow=c(2,2))
boxplot(qsip_data$A~qsip_data$id, main="EAF among groups, all ASVs")
a=aov(qsip_data$A~qsip_data$id)
summary(a)
```

```
##               Df Sum Sq Mean Sq F value   Pr(>F)    
## qsip_data$id   2  2.825  1.4123    22.3 1.35e-09 ***
## Residuals    235 14.880  0.0633                     
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
tuk=TukeyHSD(a)
plot(tuk,las=2)

incorp_only = subset(qsip_data, Incorporator == "TRUE")

boxplot(incorp_only$A~incorp_only$id, main="EAF among groups, only incorp")
a=aov(incorp_only$A~incorp_only$id)
summary(a)
```

```
##                 Df Sum Sq Mean Sq F value   Pr(>F)    
## incorp_only$id   2  1.199  0.5993   21.92 2.62e-09 ***
## Residuals      194  5.304  0.0273                     
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
tuk=TukeyHSD(a)
plot(tuk,las=2)
```

![](/Users/oliviaahern/Documents/GitHub/Sip_CopyNumber/docs/index_files/figure-html/qSIP-2-1.png)<!-- -->

```r
ag <- aggregate(A ~ id, incorp_only, function(x) c(mean = mean(x), sd = sd(x)))
```


## EAF and copy number among bacterial classes


```r
par(mfrow=c(2,2),mar=c(8,5,5,1))

incorp_only = subset(qsip_data, Incorporator == "TRUE")

boxplot(incorp_only$A~incorp_only$Class,las=2)
a=aov(incorp_only$A~incorp_only$Class)
summary(a)
```

```
##                    Df Sum Sq Mean Sq F value Pr(>F)  
## incorp_only$Class  12  0.748 0.06231   1.992 0.0271 *
## Residuals         184  5.755 0.03128                 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
tuk=TukeyHSD(a)
plot(tuk,las=2)

# bacteroidia had highest copy number 
gamma = subset(incorp_only, Class == "Bacteroidia")
mean(gamma$A)
```

```
## [1] 0.3342792
```

```r
sd(gamma$A)
```

```
## [1] 0.16447
```

```r
range(gamma$A)
```

```
## [1] 0.04734661 0.79093808
```

```r
boxplot(incorp_only$copy_number~incorp_only$Class,las=2)
a=aov(incorp_only$copy_number~incorp_only$Class)
summary(a)
```

```
##                    Df Sum Sq Mean Sq F value  Pr(>F)   
## incorp_only$Class  12   83.2   6.936   2.855 0.00126 **
## Residuals         184  447.0   2.429                   
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
tuk=TukeyHSD(a)
plot(tuk,las=2)
```

![](/Users/oliviaahern/Documents/GitHub/Sip_CopyNumber/docs/index_files/figure-html/qSIP-3-1.png)<!-- -->

```r
range(incorp_only$copy_number)
```

```
## [1] 1 9
```

```r
# gamma had the most 
gamma = subset(incorp_only, Class == "Gammaproteobacteria")
mean(gamma$copy_number)
```

```
## [1] 3.382979
```

```r
sd(gamma$copy_number)
```

```
## [1] 1.73898
```



## Comparison of incorporators 12C buoyant density 


```r
par(mfrow=c(1,2))
boxplot(incorp_only$Wlight~incorp_only$id)
a=aov(incorp_only$Wlight~incorp_only$id)
summary(a)
```

```
##                 Df   Sum Sq   Mean Sq F value Pr(>F)    
## incorp_only$id   2 0.006092 0.0030461   295.3 <2e-16 ***
## Residuals      194 0.002001 0.0000103                   
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
tuk=TukeyHSD(a)
plot(tuk,las=2)
```

![](/Users/oliviaahern/Documents/GitHub/Sip_CopyNumber/docs/index_files/figure-html/qSIP-4-1.png)<!-- -->

```r
tuk
```

```
##   Tukey multiple comparisons of means
##     95% family-wise confidence level
## 
## Fit: aov(formula = incorp_only$Wlight ~ incorp_only$id)
## 
## $`incorp_only$id`
##                               diff           lwr          upr     p adj
## batch_glu-batch_carb   0.001235454 -4.677719e-05  0.002517686 0.0616602
## chemo_carb-batch_carb -0.011764357 -1.319432e-02 -0.010334389 0.0000000
## chemo_carb-batch_glu  -0.012999811 -1.433052e-02 -0.011669102 0.0000000
```

```r
ag <- aggregate(Wlight ~ id, incorp_only, function(x) c(mean = mean(x), sd = sd(x)))
ag
```

```
##           id Wlight.mean   Wlight.sd
## 1 batch_carb 1.773856480 0.002784779
## 2  batch_glu 1.775091934 0.002456411
## 3 chemo_carb 1.762092123 0.004478331
```

## GC vs. EAF - incorporators all treatments 

```r
incorp_only = subset(qsip_data, Incorporator == "TRUE")

plot(incorp_only$Wlight, incorp_only$A)

ll=lm(incorp_only$A~incorp_only$Wlight)
summary(ll)
```

```
## 
## Call:
## lm(formula = incorp_only$A ~ incorp_only$Wlight)
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -0.31779 -0.12951 -0.03567  0.09235  0.49869 
## 
## Coefficients:
##                    Estimate Std. Error t value Pr(>|t|)    
## (Intercept)         -12.445      3.478  -3.578 0.000437 ***
## incorp_only$Wlight    7.176      1.964   3.654 0.000331 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.1767 on 195 degrees of freedom
## Multiple R-squared:  0.06409,	Adjusted R-squared:  0.05929 
## F-statistic: 13.35 on 1 and 195 DF,  p-value: 0.0003315
```

```r
abline(ll)
```

![](/Users/oliviaahern/Documents/GitHub/Sip_CopyNumber/docs/index_files/figure-html/gc-1-1.png)<!-- -->

```r
plot(incorp_only$Wlight, incorp_only$copy_number)

ll=lm(incorp_only$copy_number~incorp_only$Wlight)
summary(ll)
```

```
## 
## Call:
## lm(formula = incorp_only$copy_number ~ incorp_only$Wlight)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -1.7822 -0.7146 -0.5811  0.4135  6.4083 
## 
## Coefficients:
##                    Estimate Std. Error t value Pr(>|t|)
## (Intercept)          -16.82      32.44  -0.518    0.605
## incorp_only$Wlight    11.00      18.31   0.600    0.549
## 
## Residual standard error: 1.647 on 195 degrees of freedom
## Multiple R-squared:  0.001846,	Adjusted R-squared:  -0.003273 
## F-statistic: 0.3605 on 1 and 195 DF,  p-value: 0.5489
```

```r
abline(ll)
```

![](/Users/oliviaahern/Documents/GitHub/Sip_CopyNumber/docs/index_files/figure-html/gc-1-2.png)<!-- -->


## GC vs EAF - incorporators grouped by treatments 

```r
incorp_only = subset(qsip_data, Incorporator == "TRUE")


plot(incorp_only$Wlight, incorp_only$A)
```

![](/Users/oliviaahern/Documents/GitHub/Sip_CopyNumber/docs/index_files/figure-html/gc-2-1.png)<!-- -->

```r
par(mfrow=c(1,3),mar=c(5,5,1,1))
glu_only=subset(incorp_only, id=="batch_glu")
plot(glu_only$A, glu_only$Wlight, pch=21, cex=1.22, bg='black', xlab="AFE", ylab="Gene Copy Number")
ll=lm(glu_only$Wlight~glu_only$A)
summary(ll)
```

```
## 
## Call:
## lm(formula = glu_only$Wlight ~ glu_only$A)
## 
## Residuals:
##        Min         1Q     Median         3Q        Max 
## -0.0087071 -0.0010528 -0.0001573  0.0004101  0.0072168 
## 
## Coefficients:
##               Estimate Std. Error  t value Pr(>|t|)    
## (Intercept)  1.7757854  0.0004909 3617.753   <2e-16 ***
## glu_only$A  -0.0020478  0.0012198   -1.679    0.097 .  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.00243 on 82 degrees of freedom
## Multiple R-squared:  0.03323,	Adjusted R-squared:  0.02144 
## F-statistic: 2.818 on 1 and 82 DF,  p-value: 0.097
```

```r
abline(ll)


batch_carbon=subset(incorp_only, id=="batch_carb")
plot(batch_carbon$A, batch_carbon$Wlight, pch=21, cex=1.22, bg='black', xlab="AFE", ylab="Gene Copy Number")
ll=lm(batch_carbon$Wlight~batch_carbon$A)
summary(ll)
```

```
## 
## Call:
## lm(formula = batch_carbon$Wlight ~ batch_carbon$A)
## 
## Residuals:
##        Min         1Q     Median         3Q        Max 
## -0.0065542 -0.0009677  0.0001807  0.0008690  0.0052887 
## 
## Coefficients:
##                  Estimate Std. Error t value Pr(>|t|)    
## (Intercept)     1.7784295  0.0007308 2433.58  < 2e-16 ***
## batch_carbon$A -0.0170353  0.0025274   -6.74 7.99e-09 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.002103 on 58 degrees of freedom
## Multiple R-squared:  0.4392,	Adjusted R-squared:  0.4296 
## F-statistic: 45.43 on 1 and 58 DF,  p-value: 7.994e-09
```

```r
abline(ll)


chemo_carbon=subset(incorp_only, id=="chemo_carb")
plot(chemo_carbon$A, chemo_carbon$Wlight, pch=21, cex=1.22, bg='black', xlab="AFE", ylab="Gene Copy Number")
ll=lm(chemo_carbon$Wlight~chemo_carbon$A)
summary(ll)
```

```
## 
## Call:
## lm(formula = chemo_carbon$Wlight ~ chemo_carbon$A)
## 
## Residuals:
##        Min         1Q     Median         3Q        Max 
## -0.0140654 -0.0006758  0.0011694  0.0024761  0.0065914 
## 
## Coefficients:
##                 Estimate Std. Error  t value Pr(>|t|)    
## (Intercept)     1.762946   0.001023 1723.186   <2e-16 ***
## chemo_carbon$A -0.005825   0.005577   -1.044    0.301    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.004474 on 51 degrees of freedom
## Multiple R-squared:  0.02094,	Adjusted R-squared:  0.001744 
## F-statistic: 1.091 on 1 and 51 DF,  p-value: 0.3012
```

```r
abline(ll)
```

![](/Users/oliviaahern/Documents/GitHub/Sip_CopyNumber/docs/index_files/figure-html/gc-2-2.png)<!-- -->




## figure with afe and control buoyant density 

### SI Figure 3 


```r
#dim(incorp_only)
#incorp_only$Class
summary(as.factor(incorp_only$Class))
```

```
##      Acidimicrobiia      Actinobacteria Alphaproteobacteria         Bacteroidia 
##                   4                  11                  48                  54 
##        Chloroflexia Deltaproteobacteria Gammaproteobacteria               OM190 
##                   1                   3                  47                   1 
##    Planctomycetacia        Rhodothermia     Thermoleophilia           vadinHA49 
##                  12                   2                   1                   1 
##    Verrucomicrobiae 
##                  12
```

```r
custom_colors_class=c("Actinobacteria"="#ff6b92",
  "Alphaproteobacteria"="#aa003a",
  "Bacteroidia"="#792e00",
  "Cyanobacteriia"="#ffae55",
  "Deltaproteobacteria"="#63af33",
  "Gammaproteobacteria" = "#3eeb93",
  "OM190" = "#0060ad",
  "Planctomycetes"= "#9281fd",
  "Rhodothermia" = "#3f0e50", 
  "vadinHA49" = "#e747a0",
  "Verrucomicrobiae" = "#ff6eb1")

custom_colors_class=c("Acidimicrobiia" ="#f21277",
                      "Actinobacteria" ="#a02016",
                      "Alphaproteobacteria" ="#ffb659",
                      "Bacteroidia" = "#48591c",
                      "Chloroflexia"='gray',
                      "Deltaproteobacteria" ="#01a34a",
                      "Gammaproteobacteria" ="#7bc6ff",
                      "OM190" = 'gray',
                      "Planctomycetacia" ="#00539e",
                      "Rhodothermia"='gray',
                      "Thermoleophilia" = 'gray',
                      "vadinHA49" = "gray", 
                      "Verrucomicrobiae" ="#6868f2")
incorp_only$Class_ag=incorp_only$Class
incorp_only$Class_ag=as.factor(incorp_only$Class_ag)
levels(incorp_only$Class_ag)=c("Acidimicrobiia","Actinobacteria","Alphaproteobacteria","Bacteroidia","Other","Deltaproteobacteria","Gammaproteobacteria","Other","Planctomycetacia", "Other", "Other", "Other","Verrucomicrobiae")

custom_colors_class=c("Acidimicrobiia" ="#f21277",
                      "Actinobacteria" ="#a02016",
                      "Alphaproteobacteria" ="#ffb659",
                      "Bacteroidia" = "#48591c",
                                            "Other" = 'gray',

                      "Deltaproteobacteria" ="#01a34a",
                      "Gammaproteobacteria" ="#7bc6ff",
                      "Planctomycetacia" ="#00539e",
                      "Verrucomicrobiae" ="#6868f2")

desired_order=c("Single", "Multi")
incorp_only$Substrate <- factor(incorp_only$Substrate, levels = desired_order)


ggplot(incorp_only, aes(x=A, y=Wlight, fill=Class_ag, color=Class_ag)) + 
   geom_point(size=as.numeric(incorp_only$copy_number)) +
  #geom_point(size=3) +
  facet_wrap(Culture ~ Substrate) +
  theme_bw() +
  theme(legend.position='bottom') +
  scale_fill_manual(values=custom_colors_class) + 
  scale_color_manual(values=custom_colors_class) +
  xlab("Excess Atomic Fraction") +
  ylab(expression(" "*C^12*" Control Buoyant Density (g/mL)"))
```

![](/Users/oliviaahern/Documents/GitHub/Sip_CopyNumber/docs/index_files/figure-html/si-fig-3-1.png)<!-- -->




# 16S copy number vs EAF

## Figure 1


```r
formula = incorp_only$copy_number ~ incorp_only$A * incorp_only$Culture * incorp_only$Substrate

ggplot(incorp_only, aes(x=A, y=copy_number)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(Culture ~ Substrate) +
  theme_bw() +
  xlab("Excess Atomic Fraction") +
  ylab("16S Gene Copy Number") + 
  theme(axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())
```

![](/Users/oliviaahern/Documents/GitHub/Sip_CopyNumber/docs/index_files/figure-html/copy-2-1.png)<!-- -->

```r
  #tat_poly_line(formula = formula) +
# stat_poly_eq(formula = formula)
```



## Figure 1 alt


```r
custom_colors_class=c("Acidimicrobiia" ="#f21277",
                      "Actinobacteria" ="#a02016",
                      "Alphaproteobacteria" ="#ffb659",
                      "Bacteroidia" = "#48591c",
                      "Chloroflexia"='gray',
                      "Deltaproteobacteria" ="#01a34a",
                      "Gammaproteobacteria" ="#7bc6ff",
                      "OM190" = 'gray',
                      "Planctomycetacia" ="#00539e",
                      "Rhodothermia"='gray',
                      "Thermoleophilia" = 'gray',
                      "vadinHA49" = "gray", 
                      "Verrucomicrobiae" ="#6868f2")
incorp_only$Class_ag=incorp_only$Class
incorp_only$Class_ag=as.factor(incorp_only$Class_ag)
levels(incorp_only$Class_ag)=c("Acidimicrobiia","Actinobacteria","Alphaproteobacteria","Bacteroidia","Other","Deltaproteobacteria","Gammaproteobacteria","Other","Planctomycetacia", "Other", "Other", "Other","Verrucomicrobiae")

custom_colors_class=c("Acidimicrobiia" ="#f21277",
                      "Actinobacteria" ="#a02016",
                      "Alphaproteobacteria" ="#ffb659",
                      "Bacteroidia" = "#48591c",
                                            "Other" = 'gray',

                      "Deltaproteobacteria" ="#01a34a",
                      "Gammaproteobacteria" ="#7bc6ff",
                      "Planctomycetacia" ="#00539e",
                      "Verrucomicrobiae" ="#6868f2")

desired_order=c("Single", "Multi")
incorp_only$Substrate <- factor(incorp_only$Substrate, levels = desired_order)

formula = incorp_only$copy_number ~ incorp_only$A * incorp_only$Culture * incorp_only$Substrate

ggplot(incorp_only, aes(x=A, y=copy_number, color=Class_ag)) +
  #geom_point(size=as.numeric(incorp_only$copy_number)*0.7, alpha=0.8) +
  geom_point(size=3) +
  geom_smooth(method = "lm", se = FALSE, col='black')  +
  facet_wrap(Culture ~ Substrate) +
  theme_bw() +
  xlab("Excess Atomic Fraction") +
  ylab("16S Gene Copy Number") + 
  theme(axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
    scale_color_manual(values=custom_colors_class) 
```

![](/Users/oliviaahern/Documents/GitHub/Sip_CopyNumber/docs/index_files/figure-html/copy-2ag-1.png)<!-- -->


## Figure 1 alt2


```r
custom_colors_class=c("Acidimicrobiia" ="#f21277",
                      "Actinobacteria" ="#a02016",
                      "Alphaproteobacteria" ="#ffb659",
                      "Bacteroidia" = "#48591c",
                      "Chloroflexia"='gray',
                      "Deltaproteobacteria" ="#01a34a",
                      "Gammaproteobacteria" ="#7bc6ff",
                      "OM190" = 'gray',
                      "Planctomycetacia" ="#00539e",
                      "Rhodothermia"='gray',
                      "Thermoleophilia" = 'gray',
                      "vadinHA49" = "gray", 
                      "Verrucomicrobiae" ="#6868f2")
incorp_only$Class_ag=incorp_only$Class
incorp_only$Class_ag=as.factor(incorp_only$Class_ag)
levels(incorp_only$Class_ag)=c("Acidimicrobiia","Actinobacteria","Alphaproteobacteria","Bacteroidia","Other","Deltaproteobacteria","Gammaproteobacteria","Other","Planctomycetacia", "Other", "Other", "Other","Verrucomicrobiae")

custom_colors_class=c("Acidimicrobiia" ="#f21277",
                      "Actinobacteria" ="#a02016",
                      "Alphaproteobacteria" ="#ffb659",
                      "Bacteroidia" = "#48591c",
                                            "Other" = 'gray',

                      "Deltaproteobacteria" ="#01a34a",
                      "Gammaproteobacteria" ="#7bc6ff",
                      "Planctomycetacia" ="#00539e",
                      "Verrucomicrobiae" ="#6868f2")

desired_order=c("Single", "Multi")
incorp_only$Substrate <- factor(incorp_only$Substrate, levels = desired_order)

formula = incorp_only$copy_number ~ incorp_only$A * incorp_only$Culture * incorp_only$Substrate

ggplot(incorp_only, aes(x=A, y=copy_number, color=Class_ag)) +
  #geom_point(size=as.numeric(incorp_only$copy_number)*0.7, alpha=0.8) +
  geom_point(size=3) +
  geom_smooth(method = "lm", se = TRUE, col='black',fullrange=TRUE)  +
  facet_wrap(Culture ~ Substrate) +
  theme_bw() +
  xlab("Excess Atomic Fraction") +
  ylab("16S Gene Copy Number") + 
  theme(axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    strip.text = element_text(size = 14, face = "bold"),
    strip.background = element_blank(),
    panel.grid.minor = element_blank()) +
    scale_color_manual(values=custom_colors_class) 
```

![](/Users/oliviaahern/Documents/GitHub/Sip_CopyNumber/docs/index_files/figure-html/copy-2ag2-1.png)<!-- -->



### stats generated from here

```r
par(mfrow=c(1,3))
glu_only=subset(incorp_only, id=="batch_glu")
set.seed(124)

plot(glu_only$A, glu_only$copy_number, pch=21, cex=1.22, bg='black', xlab="AFE", ylab="Gene Copy Number")
ll=lm(glu_only$copy_number~glu_only$A)
summary(ll)
```

```
## 
## Call:
## lm(formula = glu_only$copy_number ~ glu_only$A)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -2.3094 -1.1104 -0.3466  0.7604  4.9769 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   1.3970     0.3158   4.424 2.96e-05 ***
## glu_only$A    3.9952     0.7848   5.091 2.23e-06 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 1.563 on 82 degrees of freedom
## Multiple R-squared:  0.2401,	Adjusted R-squared:  0.2309 
## F-statistic: 25.92 on 1 and 82 DF,  p-value: 2.23e-06
```

```r
abline(ll)
text(0.2, 8, "R2 = 0.23 
F = 25.92 
p = 2.23e-06")


batch_carbon=subset(incorp_only, id=="batch_carb")
plot(batch_carbon$A, batch_carbon$copy_number, pch=21, cex=1.22, bg='black', xlab="AFE", ylab="Gene Copy Number")
ll=lm(batch_carbon$copy_number~batch_carbon$A)
summary(ll)
```

```
## 
## Call:
## lm(formula = batch_carbon$copy_number ~ batch_carbon$A)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -1.6989 -1.0372 -0.3140  0.4577  6.5021 
## 
## Coefficients:
##                Estimate Std. Error t value Pr(>|t|)   
## (Intercept)       1.640      0.529   3.100  0.00298 **
## batch_carbon$A    3.700      1.830   2.022  0.04778 * 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 1.523 on 58 degrees of freedom
## Multiple R-squared:  0.06586,	Adjusted R-squared:  0.04975 
## F-statistic: 4.089 on 1 and 58 DF,  p-value: 0.04778
```

```r
abline(ll)
text(0.2, 8, "R2 = 0.04975
F = 4.089
p = 0.04778")

chemo_carbon=subset(incorp_only, id=="chemo_carb")
plot(chemo_carbon$A, chemo_carbon$copy_number, pch=21, cex=1.22, bg='black', xlab="AFE", ylab="Gene Copy Number")
ll=lm(chemo_carbon$copy_number~chemo_carbon$A)
summary(ll)
```

```
## 
## Call:
## lm(formula = chemo_carbon$copy_number ~ chemo_carbon$A)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -2.0684 -0.6138 -0.4826  0.5223  6.4817 
## 
## Coefficients:
##                Estimate Std. Error t value Pr(>|t|)    
## (Intercept)      2.4076     0.3515   6.849 9.43e-09 ***
## chemo_carbon$A   0.9518     1.9163   0.497    0.622    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 1.537 on 51 degrees of freedom
## Multiple R-squared:  0.004814,	Adjusted R-squared:  -0.0147 
## F-statistic: 0.2467 on 1 and 51 DF,  p-value: 0.6215
```

```r
abline(ll)
```

![](/Users/oliviaahern/Documents/GitHub/Sip_CopyNumber/docs/index_files/figure-html/copy-3-1.png)<!-- -->

# Beta diversity metrics 
These are to see if there are differences in commuinty between the batch and chemostat because there was different innoculum and different sequencing centers 

## philr PCoA 

### pcoa of all ASVs

```r
GP <- transform_sample_counts(phyo_frac, function(x) x+1)
phy_tree(GP) <- makeNodeLabel(phy_tree(GP), method="number", prefix='n')
name.balance(phy_tree(GP), tax_table(GP), 'n1')
```

```
## [1] "Class_Actinobacteria/Domain_Bacteria"
```

```r
otu.table <- t(otu_table(GP))
treefr <- phy_tree(GP)
metadata <- sample_data(GP)
tax <- tax_table(GP)

gp.philr <- philr(otu.table, treefr, 
                  part.weights='enorm.x.gm.counts', 
                  ilr.weights='blw.sqrt')


gp.dist <- dist(gp.philr, method="euclidean")
gp.pcoa <- ordinate(GP, 'PCoA', distance=gp.dist)

{par(mar=c(5,5,1,7))
  plot(gp.pcoa$vectors[,1],gp.pcoa$vectors[,2],
       pch=21,
       bg=as.factor(sample_data(GP)$Culture),
      xlab = "PCoA1 85.51%", ylab= "PCoA2 4.86%",xpd=F)
  ordiellipse(gp.pcoa$vectors,
              group=as.factor(sample_data(GP)$Culture),
              label=T,xpd=F,
              kind ='sd',conf=0.8)
  }
```

![](/Users/oliviaahern/Documents/GitHub/Sip_CopyNumber/docs/index_files/figure-html/pcoa-1-1.png)<!-- -->

```r
adonis2(gp.dist ~ as.factor(sample_data(GP)$Culture) +as.factor(sample_data(GP)$Substrate))
```

```
## Permutation test for adonis under reduced model
## Terms added sequentially (first to last)
## Permutation: free
## Number of permutations: 999
## 
## adonis2(formula = gp.dist ~ as.factor(sample_data(GP)$Culture) + as.factor(sample_data(GP)$Substrate))
##                                      Df SumOfSqs      R2        F Pr(>F)    
## as.factor(sample_data(GP)$Culture)    1  1178.81 0.81320 216.9951  0.001 ***
## as.factor(sample_data(GP)$Substrate)  3    26.33 0.01816   1.6157  0.164    
## Residual                             45   244.46 0.16864                    
## Total                                49  1449.60 1.00000                    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

### pcoa of just incorporators 

```r
# phyo_frac2=subset_taxa(phyo_frac, Strain %in% incorp_only$OTU)
# for some reason the above command doesn't work in my markdown so I just saved it 
phyo_frac2=readRDS('datafiles/phyo_frac2')
GP <- transform_sample_counts(phyo_frac2, function(x) x+1)
phy_tree(GP) <- makeNodeLabel(phy_tree(GP), method="number", prefix='n')
name.balance(phy_tree(GP), tax_table(GP), 'n1')
```

```
## [1] "Class_Actinobacteria/Domain_Bacteria"
```

```r
otu.table <- t(otu_table(GP))
treefr <- phy_tree(GP)
metadata <- sample_data(GP)
tax <- tax_table(GP)

gp.philr <- philr(otu.table, treefr, 
                  part.weights='enorm.x.gm.counts', 
                  ilr.weights='blw.sqrt')

gp.dist <- dist(gp.philr, method="euclidean")
gp.pcoa <- ordinate(GP, 'PCoA', distance=gp.dist)

{par(mar=c(5,5,1,7))
  plot(gp.pcoa$vectors[,1],gp.pcoa$vectors[,2],
       pch=21,
       bg=as.factor(sample_data(GP)$Culture),
      xlab = "PCoA1 89.28 %", ylab= "PCoA2 3.31 %",xpd=F) # from gp.pcoa$values$Relative_eig
  ordiellipse(gp.pcoa$vectors,
              group=as.factor(sample_data(GP)$Culture),
              label=T,xpd=F,
              kind ='sd',conf=0.8)
  par(xpd=T)
  }
```

![](/Users/oliviaahern/Documents/GitHub/Sip_CopyNumber/docs/index_files/figure-html/pcoa-2-1.png)<!-- -->

```r
adonis2(gp.dist ~ as.factor(sample_data(GP)$Culture) +as.factor(sample_data(GP)$Treatment), by='margin')
```

```
## Permutation test for adonis under reduced model
## Marginal effects of terms
## Permutation: free
## Number of permutations: 999
## 
## adonis2(formula = gp.dist ~ as.factor(sample_data(GP)$Culture) + as.factor(sample_data(GP)$Treatment), by = "margin")
##                                      Df SumOfSqs      R2       F Pr(>F)    
## as.factor(sample_data(GP)$Culture)    1  1324.73 0.48942 152.408  0.001 ***
## as.factor(sample_data(GP)$Treatment)  1    25.98 0.00960   2.989  0.090 .  
## Residual                             47   408.52 0.15093                   
## Total                                49  2706.75 1.00000                   
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

# shared asvs among treatments 


## ASVs present in all three treatmnents

I am not keeping this in because we used different CsTFA gradients for the batch and chemostat


```r
# summary(qsip_data$OTU)
# two ASVs that were in all treatments positively 
# ffd660f96ba6668b2d8d86fc4e150862
# 84ea3564fd34c5377c4305b340f0c16f

# af27577de5edb10c09f7fc19b9ca45ca true for batch but false for chemostat
# false for batch single e416d0916760d2fc17b616e2ac3ad855
```

### ASV 1

```r
par(mfrow=c(1,3))
asv2=subset_taxa(physeq_rep3_t, Strain=="ffd660f96ba6668b2d8d86fc4e150862")
tax_table(asv2)
```

```
## Taxonomy Table:     [1 taxa by 8 taxonomic ranks]:
##                                  Domain     Phylum          Class        
## ffd660f96ba6668b2d8d86fc4e150862 "Bacteria" "Bacteroidetes" "Bacteroidia"
##                                  Order                Family       Genus
## ffd660f96ba6668b2d8d86fc4e150862 "Sphingobacteriales" "env.OPS 17" ""   
##                                  Species Strain                            
## ffd660f96ba6668b2d8d86fc4e150862 ""      "ffd660f96ba6668b2d8d86fc4e150862"
```

```r
asv2_batch_single = subset_samples(asv2, Culture == "Batch" & Treatment == "Single")

plot( otu_table( subset_samples(asv2_batch_single, Isotope =="C12"))/max(otu_table( subset_samples(asv2_batch_single, Isotope =="C12"))),
      sample_data(subset_samples(asv2_batch_single, Isotope =="C12"))$diff_C12_peak,
     type="o", ylim=c(-0.05,0.051), col='gray', pch=21, bg='gray',
     ylab="Difference from C12 qPCR Peak (g/mL)", xlab="Ratio of Maximum Quantity of ASV",
     main="Batch Single")
lines( otu_table( subset_samples(asv2_batch_single, Isotope =="C13"))/max(otu_table( subset_samples(asv2_batch_single, Isotope =="C13"))),
       sample_data(subset_samples(asv2_batch_single, Isotope =="C13"))$diff_C12_peak,
     type="o", col='black', pch=21, bg='black')
text(0.6, -0.04, "EAF = 0.38")


asv2_batch_five = subset_samples(asv2, Culture == "Batch" & Treatment == "Multi")

plot( otu_table( subset_samples(asv2_batch_five, Isotope =="C12"))/max(otu_table( subset_samples(asv2_batch_five, Isotope =="C12"))),
      sample_data(subset_samples(asv2_batch_five, Isotope =="C12"))$diff_C12_peak,
     type="o", ylim=c(-0.05,0.051), col='gray', pch=21, bg='gray',
     ylab="Difference from C12 qPCR Peak (g/mL)", xlab="Ratio of Maximum Quantity of ASV",
     main="Batch Multi")
lines( otu_table( subset_samples(asv2_batch_five, Isotope =="C13"))/max(otu_table( subset_samples(asv2_batch_five, Isotope =="C13"))),
       sample_data(subset_samples(asv2_batch_five, Isotope =="C13"))$diff_C12_peak,
     type="o", col='black', pch=21, bg='black')
text(0.6, -0.04, "EAF = 0.14")


asv2_chemo = subset_samples(asv2, Culture == "Chemostat")

plot( otu_table( subset_samples(asv2_chemo, Isotope =="C12"))/max(otu_table( subset_samples(asv2_chemo, Isotope =="C12"))),
      sample_data(subset_samples(asv2_chemo, Isotope =="C12"))$diff_C12_peak,
     type="o", ylim=c(-0.05,0.051), col='gray', pch=21, bg='gray',
     ylab="Difference from C12 qPCR Peak(g/mL)", xlab="Ratio of Maximum Quantity of ASV", 
     main="Chemostat Multi")
lines(otu_table( subset_samples(asv2_chemo, Isotope =="C13"))/max(otu_table( subset_samples(asv2_chemo, Isotope =="C13"))),
      sample_data(subset_samples(asv2_chemo, Isotope =="C13"))$diff_C12_peak, 
     type="o", col='black', pch=21, bg='black')
text(0.6, -0.04, "EAF = 0.32")
```

![](/Users/oliviaahern/Documents/GitHub/Sip_CopyNumber/docs/index_files/figure-html/asv1-1.png)<!-- -->


### ASV 2

```r
par(mfrow=c(1,3))
asv2=subset_taxa(physeq_rep3_t, Strain=="84ea3564fd34c5377c4305b340f0c16f")
tax_table(asv2)
```

```
## Taxonomy Table:     [1 taxa by 8 taxonomic ranks]:
##                                  Domain     Phylum          
## 84ea3564fd34c5377c4305b340f0c16f "Bacteria" "Proteobacteria"
##                                  Class                 Order        
## 84ea3564fd34c5377c4305b340f0c16f "Alphaproteobacteria" "SAR11 clade"
##                                  Family      Genus Species
## 84ea3564fd34c5377c4305b340f0c16f "Clade III" ""    ""     
##                                  Strain                            
## 84ea3564fd34c5377c4305b340f0c16f "84ea3564fd34c5377c4305b340f0c16f"
```

```r
asv2_batch_single = subset_samples(asv2, Culture == "Batch" & Treatment == "Single")

plot( otu_table( subset_samples(asv2_batch_single, Isotope =="C12"))/max(otu_table( subset_samples(asv2_batch_single, Isotope =="C12"))),
      sample_data(subset_samples(asv2_batch_single, Isotope =="C12"))$diff_C12_peak,
     type="o", ylim=c(-0.05,0.051), col='gray', pch=21, bg='gray',
     ylab="Difference from C12 qPCR Peak (g/mL)", xlab="Ratio of Maximum Quantity of ASV",
     main="Batch Single")
lines( otu_table( subset_samples(asv2_batch_single, Isotope =="C13"))/max(otu_table( subset_samples(asv2_batch_single, Isotope =="C13"))),
       sample_data(subset_samples(asv2_batch_single, Isotope =="C13"))$diff_C12_peak,
     type="o", col='black', pch=21, bg='black')
text(0.6, -0.04, "EAF = 0.13")


asv2_batch_five = subset_samples(asv2, Culture == "Batch" & Treatment == "Multi")

plot( otu_table( subset_samples(asv2_batch_five, Isotope =="C12"))/max(otu_table( subset_samples(asv2_batch_five, Isotope =="C12"))),
      sample_data(subset_samples(asv2_batch_five, Isotope =="C12"))$diff_C12_peak,
     type="o", ylim=c(-0.05,0.051), col='gray', pch=21, bg='gray',
     ylab="Difference from C12 qPCR Peak (g/mL)", xlab="Ratio of Maximum Quantity of ASV",
     main="Batch Multi")
lines( otu_table( subset_samples(asv2_batch_five, Isotope =="C13"))/max(otu_table( subset_samples(asv2_batch_five, Isotope =="C13"))),
       sample_data(subset_samples(asv2_batch_five, Isotope =="C13"))$diff_C12_peak,
     type="o", col='black', pch=21, bg='black')
text(0.6, -0.04, "EAF = 0.15")


asv2_chemo = subset_samples(asv2, Culture == "Chemostat")

plot( otu_table( subset_samples(asv2_chemo, Isotope =="C12"))/max(otu_table( subset_samples(asv2_chemo, Isotope =="C12"))),
      sample_data(subset_samples(asv2_chemo, Isotope =="C12"))$diff_C12_peak,
     type="o", ylim=c(-0.05,0.051), col='gray', pch=21, bg='gray',
     ylab="Difference from C12 qPCR Peak (g/mL)", xlab="Ratio of Maximum Quantity of ASV", 
     main="Chemostat Multi")
lines(otu_table( subset_samples(asv2_chemo, Isotope =="C13"))/max(otu_table( subset_samples(asv2_chemo, Isotope =="C13"))),
      sample_data(subset_samples(asv2_chemo, Isotope =="C13"))$diff_C12_peak, 
     type="o", col='black', pch=21, bg='black')
text(0.6, -0.04, "EAF = 0.69")
```

![](/Users/oliviaahern/Documents/GitHub/Sip_CopyNumber/docs/index_files/figure-html/asv2-1.png)<!-- -->

### ASV 3 

```r
par(mfrow=c(1,3))
asv2=subset_taxa(physeq_rep3_t, Strain=="af27577de5edb10c09f7fc19b9ca45ca")
tax_table(asv2)
```

```
## Taxonomy Table:     [1 taxa by 8 taxonomic ranks]:
##                                  Domain     Phylum          
## af27577de5edb10c09f7fc19b9ca45ca "Bacteria" "Proteobacteria"
##                                  Class                 Order            
## af27577de5edb10c09f7fc19b9ca45ca "Alphaproteobacteria" "Rhodobacterales"
##                                  Family             Genus Species
## af27577de5edb10c09f7fc19b9ca45ca "Rhodobacteraceae" ""    ""     
##                                  Strain                            
## af27577de5edb10c09f7fc19b9ca45ca "af27577de5edb10c09f7fc19b9ca45ca"
```

```r
asv2_batch_single = subset_samples(asv2, Culture == "Batch" & Treatment == "Single")

plot( otu_table( subset_samples(asv2_batch_single, Isotope =="C12"))/max(otu_table( subset_samples(asv2_batch_single, Isotope =="C12"))),
      sample_data(subset_samples(asv2_batch_single, Isotope =="C12"))$diff_C12_peak,
     type="o", ylim=c(-0.05,0.051), col='gray', pch=21, bg='gray',
     ylab="Difference from C12 qPCR Peak (g/mL)", xlab="Ratio of Maximum Quantity of ASV",
     main="Batch Single")
lines( otu_table( subset_samples(asv2_batch_single, Isotope =="C13"))/max(otu_table( subset_samples(asv2_batch_single, Isotope =="C13"))),
       sample_data(subset_samples(asv2_batch_single, Isotope =="C13"))$diff_C12_peak,
     type="o", col='black', pch=21, bg='black')
text(0.6, -0.04, "EAF = 0.56")


asv2_batch_five = subset_samples(asv2, Culture == "Batch" & Treatment == "Multi")

plot( otu_table( subset_samples(asv2_batch_five, Isotope =="C12"))/max(otu_table( subset_samples(asv2_batch_five, Isotope =="C12"))),
      sample_data(subset_samples(asv2_batch_five, Isotope =="C12"))$diff_C12_peak,
     type="o", ylim=c(-0.05,0.051), col='gray', pch=21, bg='gray',
     ylab="Difference from C12 qPCR Peak (g/mL)", xlab="Ratio of Maximum Quantity of ASV",
     main="Batch Multi")
lines( otu_table( subset_samples(asv2_batch_five, Isotope =="C13"))/max(otu_table( subset_samples(asv2_batch_five, Isotope =="C13"))),
       sample_data(subset_samples(asv2_batch_five, Isotope =="C13"))$diff_C12_peak,
     type="o", col='black', pch=21, bg='black')
text(0.6, -0.04, "EAF = 0.31")


asv2_chemo = subset_samples(asv2, Culture == "Chemostat")

plot( otu_table( subset_samples(asv2_chemo, Isotope =="C12"))/max(otu_table( subset_samples(asv2_chemo, Isotope =="C12"))),
      sample_data(subset_samples(asv2_chemo, Isotope =="C12"))$diff_C12_peak,
     type="o", ylim=c(-0.05,0.051), col='gray', pch=21, bg='gray',
     ylab="Difference from C12 qPCR Peak (g/mL)", xlab="Ratio of Maximum Quantity of ASV", 
     main="Chemostat Multi")
lines(otu_table( subset_samples(asv2_chemo, Isotope =="C13"))/max(otu_table( subset_samples(asv2_chemo, Isotope =="C13"))),
      sample_data(subset_samples(asv2_chemo, Isotope =="C13"))$diff_C12_peak, 
     type="o", col='black', pch=21, bg='black')
text(0.6, -0.04, "EAF = NA")
```

![](/Users/oliviaahern/Documents/GitHub/Sip_CopyNumber/docs/index_files/figure-html/asv3-1.png)<!-- -->



### ASV 4

```r
par(mfrow=c(1,3))
asv2=subset_taxa(physeq_rep3_t, Strain=="e416d0916760d2fc17b616e2ac3ad855")
tax_table(asv2)
```

```
## Taxonomy Table:     [1 taxa by 8 taxonomic ranks]:
##                                  Domain     Phylum          
## e416d0916760d2fc17b616e2ac3ad855 "Bacteria" "Proteobacteria"
##                                  Class                 Order                  
## e416d0916760d2fc17b616e2ac3ad855 "Gammaproteobacteria" "Betaproteobacteriales"
##                                  Family             Genus           Species
## e416d0916760d2fc17b616e2ac3ad855 "Methylophilaceae" "Methylotenera" ""     
##                                  Strain                            
## e416d0916760d2fc17b616e2ac3ad855 "e416d0916760d2fc17b616e2ac3ad855"
```

```r
asv2_batch_single = subset_samples(asv2, Culture == "Batch" & Treatment == "Single")

plot( otu_table( subset_samples(asv2_batch_single, Isotope =="C12"))/max(otu_table( subset_samples(asv2_batch_single, Isotope =="C12"))),
      sample_data(subset_samples(asv2_batch_single, Isotope =="C12"))$diff_C12_peak,
     type="o", ylim=c(-0.05,0.051), col='gray', pch=21, bg='gray',
     ylab="Difference from C12 qPCR Peak (g/mL)", xlab="Ratio of Maximum Quantity of ASV",
     main="Batch Single")
lines( otu_table( subset_samples(asv2_batch_single, Isotope =="C13"))/max(otu_table( subset_samples(asv2_batch_single, Isotope =="C13"))),
       sample_data(subset_samples(asv2_batch_single, Isotope =="C13"))$diff_C12_peak,
     type="o", col='black', pch=21, bg='black')
text(0.6, -0.04, "EAF = NA")


asv2_batch_five = subset_samples(asv2, Culture == "Batch" & Treatment == "Multi")

plot( otu_table( subset_samples(asv2_batch_five, Isotope =="C12"))/max(otu_table( subset_samples(asv2_batch_five, Isotope =="C12"))),
      sample_data(subset_samples(asv2_batch_five, Isotope =="C12"))$diff_C12_peak,
     type="o", ylim=c(-0.05,0.051), col='gray', pch=21, bg='gray',
     ylab="Difference from C12 qPCR Peak (g/mL)", xlab="Ratio of Maximum Quantity of ASV",
     main="Batch Multi")
lines( otu_table( subset_samples(asv2_batch_five, Isotope =="C13"))/max(otu_table( subset_samples(asv2_batch_five, Isotope =="C13"))),
       sample_data(subset_samples(asv2_batch_five, Isotope =="C13"))$diff_C12_peak,
     type="o", col='black', pch=21, bg='black')
text(0.6, -0.04, "EAF = 0.22")


asv2_chemo = subset_samples(asv2, Culture == "Chemostat")

plot( otu_table( subset_samples(asv2_chemo, Isotope =="C12"))/max(otu_table( subset_samples(asv2_chemo, Isotope =="C12"))),
      sample_data(subset_samples(asv2_chemo, Isotope =="C12"))$diff_C12_peak,
     type="o", ylim=c(-0.05,0.051), col='gray', pch=21, bg='gray',
     ylab="Difference from C12 qPCR Peak (g/mL)", xlab="Ratio of Maximum Quantity of ASV", 
     main="Chemostat Multi")
lines(otu_table( subset_samples(asv2_chemo, Isotope =="C13"))/max(otu_table( subset_samples(asv2_chemo, Isotope =="C13"))),
      sample_data(subset_samples(asv2_chemo, Isotope =="C13"))$diff_C12_peak, 
     type="o", col='black', pch=21, bg='black')
text(0.6, -0.04, "EAF = 0.11")
```

![](/Users/oliviaahern/Documents/GitHub/Sip_CopyNumber/docs/index_files/figure-html/asv4-1.png)<!-- -->


## ASVs present in batch both multi and single 


```r
batch=subset(qsip_data, Culture=="Batch")
clam=(summary(batch$OTU))
clam=data.frame(clam)
clam2=row.names(subset(clam, clam==2))
length(clam2)
```

```
## [1] 63
```

```r
# 63 ASVs in batch both single and multi 

batch=subset(incorp_only, Culture=="Batch")
clam=(summary(batch$OTU))
clam=data.frame(clam)
clam2=row.names(subset(clam, clam==2))
length(clam2)
```

```
## [1] 57
```

```r
# 57 ASVs that are batch both single and multi and incorporators 
# batch_incorp=subset(batch, OTU %in% clam2) for some reason this doesn't work in my markdown %in% so I have to save it all 
# saveRDS(batch_incorp, 'datafiles/batch_incorp')
batch_incorp=readRDS('datafiles/batch_incorp')
# batch_incorp=data.frame(batch_incorp)

# summary(as.factor(batch_incorp$Class))
# summary(as.factor(batch_incorp$Family))
```



```r
multi=subset(batch_incorp, Substrate=="Multi")
single=subset(batch_incorp, Substrate=="Single")

plot(multi$A, single$A, xlab= "EAF Batch Multi Substrate", ylab ="EAF Batch Single Substrate", pch=21, ylim=c(0,1), 
     xlim=c(0,1), bg='gray',cex=(multi$copy_number*0.7),alpha = 0.5)
ll=lm(single$A~multi$A)
summary(ll)
```

```
## 
## Call:
## lm(formula = single$A ~ multi$A)
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -0.29924 -0.14486  0.01827  0.10574  0.45947 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  0.08308    0.06094   1.363    0.178    
## multi$A      1.05975    0.21035   5.038 5.41e-06 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.1727 on 55 degrees of freedom
## Multiple R-squared:  0.3158,	Adjusted R-squared:  0.3033 
## F-statistic: 25.38 on 1 and 55 DF,  p-value: 5.407e-06
```

```r
abline(ll)
text(0.8, 0.4, "R2 = 0.30 \n F = 25.38 \n p = 5.4 e-6")
abline(a=0,b=1,lty=2)
```

![](/Users/oliviaahern/Documents/GitHub/Sip_CopyNumber/docs/index_files/figure-html/batch-b-1.png)<!-- -->

### Figure 1

### batch multi eaf vs single eaf 

```r
both=cbind(single$A, multi$A, multi$Class, multi$copy_number)
colnames(both)=c("Single","Multi", "Class", "copy_number")
row.names(both)=multi$OTU
both=data.frame(both)


custom_colors_class=c("Acidimicrobiia" ="#f21277",
                      "Actinobacteria" ="#a02016",
                      "Alphaproteobacteria" ="#ffb659",
                      "Bacteroidia" = "#48591c",
                      "Chloroflexia"='gray',
                      "Deltaproteobacteria" ="#01a34a",
                      "Gammaproteobacteria" ="#7bc6ff",
                      "OM190" = 'gray',
                      "Planctomycetacia" ="#00539e",
                      "Rhodothermia"='gray',
                      "Thermoleophilia" = 'gray',
                      "vadinHA49" = "gray", 
                      "Verrucomicrobiae" ="#6868f2")

both$Class_ag=both$Class
both$Class_ag=as.factor(both$Class_ag)
levels(both$Class_ag)=c("Acidimicrobiia","Actinobacteria","Alphaproteobacteria","Bacteroidia","Other","Deltaproteobacteria","Gammaproteobacteria","Other","Planctomycetacia", "Other", "Other", "Other","Verrucomicrobiae")

custom_colors_class=c("Acidimicrobiia" ="#f21277",
                      "Actinobacteria" ="#a02016",
                      "Alphaproteobacteria" ="#ffb659",
                      "Bacteroidia" = "#48591c",
                                            "Other" = 'gray',

                      "Deltaproteobacteria" ="#01a34a",
                      "Gammaproteobacteria" ="#7bc6ff",
                      "Planctomycetacia" ="#00539e",
                      "Verrucomicrobiae" ="#6868f2")



ggplot(both, aes(x=as.numeric(Multi), y=as.numeric(Single),  color=Class)) +
  geom_point(size=as.numeric(both$copy_number), alpha=0.8) +
  theme_bw() +
  ylim(c(0,0.8)) + xlim(c(0,0.8))  +
  theme(axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
    geom_smooth(method = "lm", se = TRUE, col='black',fullrange=TRUE) + 
    geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed", linewidth = 1) +
  xlab("EAF Batch Multi Substrate") +
  ylab("EAF Batch Single-Substrate") +
  annotate("text", label="R2 = 0.30 \n F = 25.38 \n p = 5.4 e-6", x=0.6, y=0.4) +
 #   scale_fill_manual(values=custom_colors_class) + 
  scale_color_manual(values=custom_colors_class) 
```

![](/Users/oliviaahern/Documents/GitHub/Sip_CopyNumber/docs/index_files/figure-html/batch-bb-1.png)<!-- -->

#### heatmap 

```r
both=cbind(single$A, multi$A)
colnames(both)=c("Single","Multi")
row.names(both)=multi$OTU

otus=multi$OTU
# subset_batch=subset_taxa(phyo_frac, Strain %in% otus)
summary(as.factor(incorp_only$Class))
```

```
##      Acidimicrobiia      Actinobacteria Alphaproteobacteria         Bacteroidia 
##                   4                  11                  48                  54 
##        Chloroflexia Deltaproteobacteria Gammaproteobacteria               OM190 
##                   1                   3                  47                   1 
##    Planctomycetacia        Rhodothermia     Thermoleophilia           vadinHA49 
##                  12                   2                   1                   1 
##    Verrucomicrobiae 
##                  12
```

```r
# subset_batch_tree=phy_tree(subset_batch)
subset_batch_tree=readRDS('datafiles/subset_batch_tree')
tree2<-phylogram::as.dendrogram(subset_batch_tree)
xxx=RColorBrewer::brewer.pal(9,"Reds")
gplots::heatmap.2(both, Rowv=tree2,
                  dendrogram='row', trace='none',
                  scale = 'none',
                  col=xxx, sepcolor = 'gray',
                  sepwidth = c(0.005,0.05),
                  rowsep = c(0:57),
                  cexRow = 0.5,
                  cexCol = 0.8,
                  colsep = c(0,1,2),
                  margins=c(5,10))
```

![](/Users/oliviaahern/Documents/GitHub/Sip_CopyNumber/docs/index_files/figure-html/batch-c-1.png)<!-- -->
