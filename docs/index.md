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


``` r
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

### SI Figure 2 


``` r
setwd("/Users/oliviaahern/Documents/GitHub/SIP_CopyNumber")


{
  
# file contains info on MC2, Chemostat, measurements taken every 15 minutes
data=read.csv(file="datafiles/m2_do_02_edit.csv",header=T)
# subset to get the right time frame
data=data.frame(data)
data$T=as.numeric(data$T)
read=subset(data, T <= 24)
read=subset(read, T >= -100)

time=read$T
ph= read$pH.M2
do=read$DO.M2..uM.
heado=read$O2.M2....
heado2=read$CO2.M2....

  par(mar=c(4,7,1,5), mfrow=c(3,1),xpd=F)
  plot(time, do,type='l',col='blue',xlab=" ", ylab = "", xlim=c(-100,25), xaxt='n', xaxs = 'i',cex.axis=1.4, yaxt='n', ylim=c(200,285))
  axis(2, at=c(200,240,285),cex.axis=1.4)
  rect(xleft=0, xright=240, ybottom=111, ytop=300, col= rgb(0.7,0.7,0.7,alpha=0.2),
    border=NA)
  mtext(side=2,"Dissolved Oxygen (uM)", col = 'blue', padj=-3)
    axis(side=1, at=c(-100, -75, -50, -25, 0, 24),cex.axis=1.4,labels = FALSE)
  par(new = TRUE)
  plot(time, ph,col='red',type = "l", axes = FALSE, bty = "n", xlab = "", ylab = "",xlim=c(-100,25),xaxs = 'i',cex.axis=1.4, ylim=c(6.5,7.5))
  mtext(side=4, "pH", col ='red', line=3)
  axis(side=4, at = c(6.5, 7,7.5),cex.axis=1.4)
  mtext("a", adj=-.12,padj=-10,side=1,srt=-90,cex=1.5,font=1)
#  legend(25,7.4,legend=c("12C Glucose", "13C Glucose"), pch=22, pt.bg=c("white",'gray70'), bt='n')
  
  
  
  plot(time, heado,type='l',col='navy',xlab=" ", ylab = "", xlim=c(-100,25), xaxt='n', xaxs = 'i',cex.axis=1.4, yaxt='n', ylim=c(20.5,21))
  axis(2, at=c(20.5,20.75,21),cex.axis=1.4)
  rect(xleft=0, xright=240, ybottom=0, ytop=300, col= rgb(0.7,0.7,0.7,alpha=0.2),border=NA)
  mtext(side=2,"Oxygen (%)", col = 'navy', padj=-3)
    axis(side=1, at=c(-100, -75, -50, -25, 0, 24),cex.axis=1.4,labels = FALSE)
  par(new = TRUE)
  plot(time, heado2,col='forestgreen',type = "l", axes = FALSE, bty = "n", xlab = "", xlim=c(-100,24),ylab = "", xaxs = 'i',cex.axis=1.4, yaxt='n', ylim=c(0,0.4))
  mtext(side=4, "CO2 (%)", col ='forestgreen', line=3)
  axis(side=4, at = c(0,0.2,0.4),cex.axis=1.4)
  mtext("b", adj=-.12,padj=-10,side=1,srt=-90,cex=1.5,font=1)
 
# file contains info on MC2, Chemostat, data was taken daily
read=read.csv(file="datafiles/states_gluc.csv",header=T)
read=read[1:7,]
read=data.frame(read)
glucose=read$Glucose_um
time_glucose=read$Time
  
  # par(mar=c(5,41.5,1,5))
  plot(time_glucose,read$Glucose_um,type='o',ylim=c(0,12),xlab= " ",
       ylab = " ",col='gray50',bg='gray50',pch=21,cex.axis=1,cex.lab=1,cex=1.2,lwd=1.7, xlim=c(-100,24), xaxs = 'i',
       yaxt='n', xaxt='n',cex.axis=1.4) #xaxt=c(-144,0,48,120,240))
  rect(xleft=0, xright=24, ybottom=0, ytop=12, col= rgb(0.7,0.7,0.7,alpha=0.2),
    border=NA)
  mtext(side=2, "Glucose (um)", padj=-3, col="gray50")
  axis(side=2,at=c(0,4,8,12),cex.axis=1.4)
  par(new = TRUE)
  plot(time_glucose,read$Atom..13C,pch=23, lwd=1.7,
       bg='black',type = "o", axes = FALSE, bty = "n", xlab = "", ylab = "", ylim=c(0,30), cex=1.2,xlim=c(-100,25), xaxs = 'i',cex.axis=1.4)
  mtext(side=4, "Atomic 13C (%)", line=3)
  axis(side=4, at=c(0, 10,20,30),cex.axis=1.4)
    axis(side=1, at=c(-100, -75, -50, -25, 0, 24),cex.axis=1.4)

  mtext(side=1, "Time (hours)", line=3)
  #legend('topleft',legend=c("Glucose (uM)","Atomic 13C"), pch =c(21,23), pt.bg= c('gray70','black'), bty='n')
 # mtext("c", adj=-.15,padj=-9,side=1,srt=-90,cex=1.5,font=1)
    mtext("c", adj=-.12,padj=-10,side=1,srt=-90,cex=1.5,font=1)

}
```

![](/Users/oliviaahern/Documents/GitHub/Sip_CopyNumber/docs/index_files/figure-html/si-fig-1-again-1.png)<!-- -->

### SI Figure 2 (old but newer)


``` r
setwd("/Users/oliviaahern/Documents/GitHub/SIP_CopyNumber")


{
  
# file contains info on MC2, Chemostat, measurements taken every 15 minutes
data=read.csv(file="datafiles/m2_do_02_edit.csv",header=T)
# subset to get the right time frame
data=data.frame(data)
data$T=as.numeric(data$T)
read=subset(data, T <= 24)
read=subset(read, T >= -100)

time=read$T
ph= read$pH.M2
do=read$DO.M2..uM.
heado=read$O2.M2....
heado2=read$CO2.M2....

  par(mar=c(2,7,1,5), mfrow=c(3,1),xpd=F)
  plot(time, do,type='l',col='blue',xlab=" ", ylab = "", xlim=c(-100,25), xaxt='n', xaxs = 'i',cex.axis=1.4, yaxt='n')
  axis(2, at=c(200,210,220),cex.axis=1.4)
  rect(xleft=0, xright=240, ybottom=111, ytop=300, col= rgb(0.7,0.7,0.7,alpha=0.2),
    border=NA)
  mtext(side=2,"Dissolved Oxygen (uM)", col = 'blue', padj=-3)
    axis(side=1, at=c(-450,-144,-100, 0, 24, 48, 72, 120, 240),cex.axis=1.4,labels = FALSE)
  par(new = TRUE)
  plot(time, ph,col='red',type = "l", axes = FALSE, bty = "n", xlab = "", ylab = "",xlim=c(-100,25),xaxs = 'i',cex.axis=1.4)
  mtext(side=4, "pH", col ='red', line=3)
  axis(side=4, at = c(6.4, 6.55,6.6),cex.axis=1.4)
  mtext("a", adj=-.12,padj=-10,side=1,srt=-90,cex=1.5,font=1)
#  legend(25,7.4,legend=c("12C Glucose", "13C Glucose"), pch=22, pt.bg=c("white",'gray70'), bt='n')
  
  
  
  plot(time, heado,type='l',col='navy',xlab=" ", ylab = "", xlim=c(-100,25), xaxt='n', xaxs = 'i',cex.axis=1.4, yaxt='n')
  axis(2, at=c(20.52, 20.57,20.62),cex.axis=1.4)
  rect(xleft=0, xright=240, ybottom=0, ytop=300, col= rgb(0.7,0.7,0.7,alpha=0.2),border=NA)
  mtext(side=2,"Oxygen (%)", col = 'navy', padj=-3)
    axis(side=1, at=c(-450,-144,-100, 0, 24, 48, 72, 120, 240),cex.axis=1.4,labels=FALSE)
  par(new = TRUE)
  plot(time, heado2,col='forestgreen',type = "l", axes = FALSE, bty = "n", xlab = "", xlim=c(-100,24),ylab = "", xaxs = 'i',cex.axis=1.4, yaxt='n')
  mtext(side=4, "CO2 (%)", col ='forestgreen', line=3)
  axis(side=4, at = c(0.27,0.3,0.34),cex.axis=1.4)
  mtext("b", adj=-.12,padj=-10,side=1,srt=-90,cex=1.5,font=1)
 
# file contains info on MC2, Chemostat, data was taken daily
read=read.csv(file="datafiles/states_gluc.csv",header=T)
read=read[1:7,]
read=data.frame(read)
glucose=read$Glucose_um
time_glucose=read$Time
  
  # par(mar=c(5,41.5,1,5))
  plot(time_glucose,read$Glucose_um,type='o',ylim=c(0,12),xlab= " ",
       ylab = " ",col='gray50',bg='gray50',pch=21,cex.axis=1,cex.lab=1,cex=1.2,lwd=1.7, xlim=c(-100,24), xaxs = 'i',
       yaxt='n', xaxt='n',cex.axis=1.4) #xaxt=c(-144,0,48,120,240))
  rect(xleft=0, xright=24, ybottom=0, ytop=12, col= rgb(0.7,0.7,0.7,alpha=0.2),
    border=NA)
  mtext(side=2, "Glucose (um)", padj=-3, col="gray50")
  axis(side=2,at=c(0,4,8,12),cex.axis=1.4)
  par(new = TRUE)
  plot(time_glucose,read$Atom..13C,pch=23, lwd=1.7,
       bg='black',type = "o", axes = FALSE, bty = "n", xlab = "", ylab = "", ylim=c(0,30), cex=1.2,xlim=c(-100,25), xaxs = 'i',cex.axis=1.4)
  mtext(side=4, "Atomic 13C (%)", line=3)
  axis(side=4, at=c(0, 10,20,30),cex.axis=1.4)
    axis(side=1, at=c(-144,-100, 0, 24),cex.axis=1.4)

  mtext(side=1, "Time (hours)", line=3)
  #legend('topleft',legend=c("Glucose (uM)","Atomic 13C"), pch =c(21,23), pt.bg= c('gray70','black'), bty='n')
 # mtext("c", adj=-.15,padj=-9,side=1,srt=-90,cex=1.5,font=1)
    mtext("c", adj=-.12,padj=-10,side=1,srt=-90,cex=1.5,font=1)

}
```

![](/Users/oliviaahern/Documents/GitHub/Sip_CopyNumber/docs/index_files/figure-html/si-fig-1-good-1.png)<!-- -->

### SI Figure 2 (old)


``` r
setwd("/Users/oliviaahern/Documents/GitHub/SIP_CopyNumber")


{
  
# file contains info on MC2, Chemostat, measurements taken every 15 minutes
data=read.csv(file="datafiles/m2_do_02_edit.csv",header=T)
# subset to get the right time frame
data=data.frame(data)
data$T=as.numeric(data$T)
read=subset(data, T <= 24)

time=read$T
ph= read$pH.M2
do=read$DO.M2..uM.
heado=read$O2.M2....
heado2=read$CO2.M2....

  par(mar=c(2,7,1,5), mfrow=c(3,1),xpd=F)
  plot(time, do,type='l',col='blue',xlab=" ", ylab = "", xlim=c(-450,25), xaxt='n', xaxs = 'i')
  rect(xleft=0, xright=240, ybottom=111, ytop=300, col= rgb(0.7,0.7,0.7,alpha=0.2),
    border=NA)
  mtext(side=2,"Dissolved Oxygen (uM)", col = 'blue', padj=-3)
    axis(side=1, at=c(-450,-144,-100, 0, 24, 48, 72, 120, 240))
  par(new = TRUE)
  plot(time, ph,col='red',type = "l", axes = FALSE, bty = "n", xlab = "", ylab = "",xlim=c(-450,25),xaxs = 'i')
  mtext(side=4, "pH", col ='red', line=3)
  axis(side=4, at = c(6.5, 6.7, 7, 7.2, 7.4))
  mtext("a", adj=-.10,padj=-10,side=1,srt=-90,cex=1.5,font=1)
#  legend(25,7.4,legend=c("12C Glucose", "13C Glucose"), pch=22, pt.bg=c("white",'gray70'), bt='n')
  
  
  
  plot(time, heado,type='l',col='navy',xlab=" ", ylab = "", xlim=c(-450,25), xaxt='n', xaxs = 'i')
  rect(xleft=0, xright=240, ybottom=0, ytop=300, col= rgb(0.7,0.7,0.7,alpha=0.2),border=NA)
  mtext(side=2,"Oxygen (%)", col = 'navy', padj=-3)
    axis(side=1, at=c(-450,-144,-100, 0, 24, 48, 72, 120, 240))
  par(new = TRUE)
  plot(time, heado2,col='forestgreen',type = "l", axes = FALSE, bty = "n", xlab = "", xlim=c(-450,24),ylab = "", xaxs = 'i')
  mtext(side=4, "CO2 (%)", col ='forestgreen', line=3)
  axis(side=4, at = c(0,0.1,0.2,0.3,0.4))
  mtext("b", adj=-.1,padj=-10,side=1,srt=-90,cex=1.5,font=1)
 
# file contains info on MC2, Chemostat, data was taken daily
read=read.csv(file="datafiles/states_gluc.csv",header=T)
read=read[1:7,]
read=data.frame(read)
glucose=read$Glucose_um
time_glucose=read$Time
  
   par(mar=c(5,41.5,1,5))
  plot(time_glucose,read$Glucose_um,type='l',ylim=c(0,12),xlab= " ",
       ylab = " ",col='gray50',bg='gray50',pch=21,cex.axis=1,cex.lab=1,cex=1.2,lwd=1.7, xlim=c(-144,24), xaxs = 'i',
       yaxt='n', xaxt='n') #xaxt=c(-144,0,48,120,240))
  rect(xleft=0, xright=24, ybottom=0, ytop=12, col= rgb(0.7,0.7,0.7,alpha=0.2),
    border=NA)
  mtext(side=2, "Glucose (um)", padj=-3, col="gray50")
  axis(side=2,at=c(0,4,8,12), cex=1.3)
  par(new = TRUE)
  plot(time_glucose,read$Atom..13C,pch=23, lwd=1.7,
       bg='black',type = "l", axes = FALSE, bty = "n", xlab = "", ylab = "", ylim=c(0,30), cex=1.2,xlim=c(-144,25), xaxs = 'i')
  mtext(side=4, "Atomic 13C (%)", line=3)
  axis(side=4, at=c(0, 10,20,30))
    axis(side=1, at=c(-144,-100, 0, 24))

  mtext(side=1, "Time (hours)", line=3)
  #legend('topleft',legend=c("Glucose (uM)","Atomic 13C"), pch =c(21,23), pt.bg= c('gray70','black'), bty='n')
  mtext("c", adj=-.15,padj=-9,side=1,srt=-90,cex=1.5,font=1)
}
```

![](/Users/oliviaahern/Documents/GitHub/Sip_CopyNumber/docs/index_files/figure-html/si-fig-1-1.png)<!-- -->

## qPCR Buoyant Density Gradients

### SI Figure 2
old figure
```

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
### SI Figure 1 used


``` r
data<-read.csv("datafiles/exp3_bd.csv",header=T,row.names=1)
qPCR_data=data.frame(data)

desired_order=c("Single", "Multi")
qPCR_data$Substrate <- factor(qPCR_data$Substrate, levels = desired_order)



ggplot(data=qPCR_data, aes(x=BD, y=R_MAX, col=Isotope))+
    geom_line() +
  geom_point(aes(shape=Sequenced),size=2.5) +
  theme_bw() +
    facet_wrap(Culture ~ Substrate)  +
  coord_flip() +
  scale_color_manual(values=c("gray","black")) + 
#   scale_fill_manual(values=c("white","black")) + 
  xlab("Buoyant Density (g/mL)") +
  ylab("Ratio of Maximum Quantity") +
    scale_x_reverse(breaks=c(1.72,1.74,1.76,1.78,1.80, 1.82)) +
  # xlim(c(1.72, 1.841)) +
   theme(axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    strip.text = element_text(size = 14, face = "bold"),
    strip.background = element_blank(),
    panel.grid.minor = element_blank())
```

![](/Users/oliviaahern/Documents/GitHub/Sip_CopyNumber/docs/index_files/figure-html/si-fig2b-1.png)<!-- -->


# Import 16S Data


``` r
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


physeq_rep3_t=readRDS('datafiles/physeq_rep3_t')
```


# Code for running qSIP

## 1-normalize reads to qPCR copy number concentrations

```
phyo_frac=subset_samples(phyo1, Isotope !="NA")
ASVs_to_keep <- taxa_sums(phyo_frac) > 0
phyo_frac <- prune_taxa(ASVs_to_keep, phyo_frac)


map<-read.csv("datafiles/exp3_qpcr_combined.csv",
              header=T,row.names=1)
physeq_rep3_t = OTU_qPCR_trans(phyo_frac, map)
saveRDS(physeq_rep3_t, 'datafiles/physeq_rep3_t')
```

## 2-batch-single

```
batch_glu=subset_samples(physeq_rep3_t, Culture=="Batch" & Treatment=="Single")

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

```
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

```
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

```

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
colnames(taxa)=c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "OTU")
qsip_data=merge(all, taxa, by="OTU")

write.csv(qsip_data, 'datafiles/qsip_data_exp3.csv')

```


``` r
#saveRDS(qsip_data, 'datafiles/qsip_data')
qsip_data=readRDS('datafiles/qsip_data')
```

# qSIP basic stats

## 16S rRNA copy number among treatments
looking at all ASVs that were in the buoyant density gradients (must be present at in at least 2 fractions of both enriched and unenriched) and just the incorporators 


``` r
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

``` r
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

``` r
tuk=TukeyHSD(a)
plot(tuk,las=2)
```

![](/Users/oliviaahern/Documents/GitHub/Sip_CopyNumber/docs/index_files/figure-html/qSIP-1-1.png)<!-- -->

``` r
mean_by_group <- aggregate(incorp_only$copy_number, by=list(incorp_only$id), FUN=mean)
sd_by_group <- aggregate(incorp_only$copy_number, by=list(incorp_only$id), FUN=sd)
```


## EAF among treatments 
looking at all ASVs that were in the buoyant density gradients (must be present at in at least 2 fractions of both enriched and unenriched) and just the incorporators 


``` r
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

``` r
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

``` r
tuk=TukeyHSD(a)
plot(tuk,las=2)
```

![](/Users/oliviaahern/Documents/GitHub/Sip_CopyNumber/docs/index_files/figure-html/qSIP-2-1.png)<!-- -->

``` r
ag <- aggregate(A ~ id, incorp_only, function(x) c(mean = mean(x), sd = sd(x)))

mean_by_group <- aggregate(incorp_only$A, by=list(incorp_only$id), FUN=mean)
sd_by_group <- aggregate(incorp_only$A, by=list(incorp_only$id), FUN=sd)
```


## EAF and copy number among bacterial classes


``` r
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

``` r
tuk=TukeyHSD(a)
plot(tuk,las=2)

# bacteroidia had highest copy number 
gamma = subset(incorp_only, Class == "Bacteroidia")
mean(gamma$A)
```

```
## [1] 0.3342792
```

``` r
sd(gamma$A)
```

```
## [1] 0.16447
```

``` r
range(gamma$A)
```

```
## [1] 0.04734661 0.79093808
```

``` r
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

``` r
tuk=TukeyHSD(a)
plot(tuk,las=2)
```

![](/Users/oliviaahern/Documents/GitHub/Sip_CopyNumber/docs/index_files/figure-html/qSIP-3-1.png)<!-- -->

``` r
range(incorp_only$copy_number)
```

```
## [1] 1 9
```

``` r
# gamma had the most 
gamma = subset(incorp_only, Class == "Gammaproteobacteria")
mean(gamma$copy_number)
```

```
## [1] 3.382979
```

``` r
sd(gamma$copy_number)
```

```
## [1] 1.73898
```



## Comparison of incorporators 12C buoyant density 


``` r
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

``` r
tuk=TukeyHSD(a)
plot(tuk,las=2)
```

![](/Users/oliviaahern/Documents/GitHub/Sip_CopyNumber/docs/index_files/figure-html/qSIP-4-1.png)<!-- -->

``` r
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

``` r
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

``` r
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

``` r
abline(ll)
```

![](/Users/oliviaahern/Documents/GitHub/Sip_CopyNumber/docs/index_files/figure-html/gc-1-1.png)<!-- -->

``` r
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

``` r
abline(ll)
```

![](/Users/oliviaahern/Documents/GitHub/Sip_CopyNumber/docs/index_files/figure-html/gc-1-2.png)<!-- -->


## GC vs EAF - incorporators grouped by treatments 

``` r
incorp_only = subset(qsip_data, Incorporator == "TRUE")


plot(incorp_only$Wlight, incorp_only$A)
```

![](/Users/oliviaahern/Documents/GitHub/Sip_CopyNumber/docs/index_files/figure-html/gc-2-1.png)<!-- -->

``` r
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

``` r
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

``` r
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

``` r
abline(ll)
```

![](/Users/oliviaahern/Documents/GitHub/Sip_CopyNumber/docs/index_files/figure-html/gc-2-2.png)<!-- -->


## IOC

``` r
qsip_data
```

```
##                                  OTU     Wlab   Wlight             Z       Gi
## 1   00feb692689d706370ba89503e9ae97e 1.762555 1.763477 -0.0009223364 1.406131
## 2   02dd75f976cd05e6517f3d4da0a1a415 1.761812 1.761111  0.0007008874 1.377795
## 3   06eb32e29df10778b60c49356ef3b5e0 1.769951 1.779058 -0.0091063377 1.592710
## 4   07eb166d091e1aebe641105881f4f2d3 1.769075 1.764848  0.0042273031 1.422541
## 5   0893f5bbbc8d21217aa590338706fe49 1.763235 1.787805 -0.0245706139 1.697462
## 6   0ae5531a9802e71139c959faebdde10c 1.769420 1.764580  0.0048396941 1.419337
## 7   0b06017e8169c7ce93b22574992470dc 1.796788 1.773863  0.0229243248 1.530505
## 8   0b06017e8169c7ce93b22574992470dc 1.787391 1.773817  0.0135740469 1.529955
## 9   0e4971b9fdfe3044222a19f5b7f167dc 1.785297 1.775737  0.0095609406 1.552937
## 10  0e4971b9fdfe3044222a19f5b7f167dc 1.791636 1.773568  0.0180673787 1.526973
## 11  101ab1c0cf819a3e48e09b3356bb6f3d 1.752202 1.748803  0.0033981136 1.230408
## 12  1083eb67a7140dca982f34e3968117a3 1.767051 1.763152  0.0038991359 1.402234
## 13  11b50e80047289f589868d85531954f4 1.792017 1.782529  0.0094871625 1.634283
## 14  11b50e80047289f589868d85531954f4 1.791030 1.773821  0.0172093069 1.529995
## 15  123fa4546d88beaa1017c112c3c600cc 1.786109 1.770283  0.0158256723 1.487628
## 16  137dc661fa672b506fd805198510edf9 1.766126 1.762271  0.0038558164 1.391680
## 17  1583d8af9fd9cb078f88c0e1a8b8a888 1.769742 1.763430  0.0063122896 1.405565
## 18  177d2176a3368f43cc132f9e34a9405c 1.788832 1.773502  0.0153301453 1.526179
## 19  177d2176a3368f43cc132f9e34a9405c 1.786390 1.775346  0.0110433937 1.548262
## 20  1bccea8882fb44e316fa3626a2062514 1.785726 1.772959  0.0127673647 1.519676
## 21  1bccea8882fb44e316fa3626a2062514 1.780886 1.775177  0.0057088579 1.546236
## 22  1c2b8892a0d7f828d72808cd01cfd199 1.772908 1.758688  0.0142204801 1.348774
## 23  1e104b993416995498346089f57fd3bf 1.768826 1.763328  0.0054978378 1.404346
## 24  20ada5c0fa2722151f9039d1896592e8 1.793310 1.770065  0.0232448598 1.485017
## 25  211b7bfb480d7e0ba8224a99279d1055 1.797864 1.775203  0.0226604463 1.546552
## 26  211b7bfb480d7e0ba8224a99279d1055 1.793327 1.773499  0.0198285267 1.526137
## 27  21bdb6c8cb7ed339792d39838ea86d08 1.764588 1.752696  0.0118913624 1.277027
## 28  2461a8cf82261ae32270bfb8a9e60b2e 1.764662 1.796457 -0.0317948080 1.801070
## 29  27af011f1ebcc0867bb168193e6a9d27 1.769230 1.761423  0.0078069934 1.381526
## 30  27c7baf02d6ba15f77c56b6acb9434c0 1.788195 1.789411 -0.0012156375 1.716692
## 31  282d970c7e74112691085e63921dbcfd 1.787645 1.774325  0.0133190450 1.536039
## 32  282d970c7e74112691085e63921dbcfd 1.784796 1.782262  0.0025339953 1.631082
## 33  285bfcf58e6f5f2669ba8a2efb4ed898 1.773778 1.757220  0.0165571962 1.331202
## 34  29731c1ae123c2faa2b4184d96385388 1.791869 1.776170  0.0156987066 1.558128
## 35  29b2870e522fa2d0b60ae25278a57c5b 1.767132 1.763183  0.0039486908 1.402605
## 36  29f3083d2d8f993bec5a67e7af219321 1.760520 1.764415 -0.0038941413 1.417355
## 37  2aee17e668dd944e7792d625ec479fcf 1.767498 1.794593 -0.0270951663 1.778751
## 38  2b4e6ccd11a78681bceaffed2053907f 1.792437 1.773974  0.0184633256 1.531826
## 39  2b4e6ccd11a78681bceaffed2053907f 1.787443 1.774198  0.0132445131 1.534515
## 40  2c070f4880cc78c48ec292dd9c89b597 1.803086 1.774964  0.0281220537 1.543686
## 41  2c070f4880cc78c48ec292dd9c89b597 1.794294 1.772998  0.0212958434 1.520143
## 42  2da31e61296f678bb466e1ae0593e33d 1.784579 1.774809  0.0097694333 1.541835
## 43  304d1a8fe851969c2eb82c85f9f88c11 1.782208 1.776511  0.0056974075 1.562210
## 44  34327ad682b105ca4607c44af24c5dd2 1.776604 1.764243  0.0123605547 1.415305
## 45  36dae39c87c9e05e2a808f08e5b671ce 1.764084 1.761419  0.0026646388 1.381486
## 46  37762c5e9b1d5a5233d6caeef3660d6a 1.772797 1.776735 -0.0039371413 1.564889
## 47  3bb0bf5132b0f988573de9bd6059b22c 1.784227 1.781022  0.0032050898 1.616236
## 48  3bb0bf5132b0f988573de9bd6059b22c 1.778663 1.774012  0.0046511980 1.532286
## 49  3c8a9630951999d2b809868618b05848 1.779649 1.774431  0.0052174246 1.537308
## 50  3c8a9630951999d2b809868618b05848 1.786620 1.774179  0.0124407307 1.534286
## 51  3d2961d34551d101a985e74adfc1e128 1.759380 1.764609 -0.0052287420 1.419682
## 52  3e5cc4006667756b128e551a166d5861 1.775425 1.763367  0.0120579521 1.404813
## 53  3eef12286c142030d497dc480cc49332 1.774672 1.786885 -0.0122121387 1.686436
## 54  3fa2974b1f70109275811beb1ae321c9 1.778209 1.764465  0.0137438551 1.417955
## 55  3fc227a79c6a108319409a05f87103f3 1.786796 1.782337  0.0044589463 1.631983
## 56  448ee48054e6edfd36e703a6baecc8a8 1.816894 1.776828  0.0400661583 1.566004
## 57  4744c03db226d903293f9342f292e94d 1.786122 1.773967  0.0121550003 1.531751
## 58  4744c03db226d903293f9342f292e94d 1.782099 1.774197  0.0079012380 1.534506
## 59  4800a9c596272fb238a316d893b0871e 1.787925 1.773833  0.0140926541 1.530136
## 60  4964f339e90a20437a19d7ec3bcd8a1a 1.784294 1.770735  0.0135595246 1.493037
## 61  4964f339e90a20437a19d7ec3bcd8a1a 1.799206 1.776338  0.0228679512 1.560135
## 62  4d65f65496c3d0238ebf74738245e526 1.790367 1.773802  0.0165648206 1.529772
## 63  4d65f65496c3d0238ebf74738245e526 1.797326 1.774210  0.0231159094 1.534657
## 64  50db42495de3aa6b40560c8c73e032d5 1.764215 1.765185 -0.0009691262 1.426575
## 65  520c7e89d0c7c72e7d85eed29aaedb70 1.768716 1.766358  0.0023583566 1.440623
## 66  55ce73379b753a6d229a4febc0b1362d 1.799584 1.774241  0.0253429550 1.535025
## 67  55ce73379b753a6d229a4febc0b1362d 1.790627 1.793242 -0.0026143670 1.762563
## 68  5914e4ffa20a5fa0b95bf68e5c9777b6 1.786057 1.776796  0.0092606082 1.565624
## 69  5aad5881139a653bc5964f1d88d7aaee 1.772210 1.764595  0.0076148095 1.419517
## 70  5d2c9965718d1254c00a642f570eb604 1.767213 1.773690 -0.0064776296 1.528434
## 71  5ed8f8f3a98e669d81b30c56825e0870 1.774118 1.761786  0.0123320477 1.385880
## 72  6070669b932f1826dc8b8e61495fd27c 1.781305 1.775000  0.0063050977 1.544117
## 73  6267512d7ea98a0da430e5a6f96698e6 1.787283 1.775338  0.0119445373 1.548168
## 74  6267512d7ea98a0da430e5a6f96698e6 1.793656 1.774359  0.0192968968 1.536437
## 75  63413adbe8227f80412e824dffac1f5a 1.769147 1.762848  0.0062998418 1.398589
## 76  64223c9ef560d611f5b9e7e5e78eb7ff 1.767008 1.798028 -0.0310198279 1.819883
## 77  64edc344889f4279a1292afbe81dc043 1.781159 1.773803  0.0073564557 1.529777
## 78  64edc344889f4279a1292afbe81dc043 1.782433 1.772366  0.0100663019 1.512577
## 79  65a75a16ddf0882daec9e72eea3aa180 1.765954 1.760196  0.0057579979 1.366835
## 80  6650e53404090b74619be81139104295 1.751270 1.789306 -0.0380359986 1.715431
## 81  68262884f1021c896f8e1bf7348d773c 1.807936 1.773209  0.0347274107 1.522665
## 82  68262884f1021c896f8e1bf7348d773c 1.796562 1.770815  0.0257473898 1.494000
## 83  6a679587c40215f64acb6b0231b57581 1.784981 1.775760  0.0092217393 1.553212
## 84  6b94959874b693c3304e68c57cefd287 1.786691 1.773845  0.0128466194 1.530283
## 85  6b94959874b693c3304e68c57cefd287 1.794994 1.775069  0.0199255655 1.544938
## 86  6e6f306a337ba255a7f4a7f8469ff691 1.783041 1.773533  0.0095082690 1.526551
## 87  6e6f306a337ba255a7f4a7f8469ff691 1.813338 1.774201  0.0391377699 1.534544
## 88  70f086828189cd2b21299d41e077a74a 1.784336 1.776697  0.0076388043 1.564444
## 89  710dab4db94586d934436dc307f766c9 1.793562 1.771826  0.0217356207 1.506113
## 90  710dab4db94586d934436dc307f766c9 1.804188 1.774705  0.0294821418 1.540589
## 91  740d9b2ea5713119fec5849936b9128b 1.770284 1.765790  0.0044947194 1.433821
## 92  748289d833e671b4848c77b9443c2034 1.770982 1.763152  0.0078302485 1.402235
## 93  755afc7ffddcfd5e093ea2222f4affd6 1.784569 1.773753  0.0108152122 1.529188
## 94  755afc7ffddcfd5e093ea2222f4affd6 1.785687 1.772228  0.0134588963 1.510920
## 95  7580ff42ec296410fd9b751933c26164 1.780934 1.780109  0.0008249781 1.605295
## 96  7580ff42ec296410fd9b751933c26164 1.779532 1.773951  0.0055817754 1.531551
## 97  76a2e443105acc53d3461f176f38148c 1.782574 1.773641  0.0089325539 1.527846
## 98  76a2e443105acc53d3461f176f38148c 1.810104 1.774931  0.0351726898 1.543292
## 99  788f25619c77abb69cf199dd19233eb4 1.785784 1.780742  0.0050420397 1.612882
## 100 788f25619c77abb69cf199dd19233eb4 1.796246 1.770775  0.0254708964 1.493517
## 101 78c20e329a7b47d11c13171c92a1aad2 1.754046 1.768827 -0.0147805280 1.470189
## 102 7beab4f249af89a8aab102a03113ba62 1.809030 1.773820  0.0352102634 1.529988
## 103 7beab4f249af89a8aab102a03113ba62 1.781474 1.769049  0.0124255263 1.472848
## 104 7c338eeece15edbfcf8ed52c568f6fc2 1.781629 1.778692  0.0029375307 1.588327
## 105 7c338eeece15edbfcf8ed52c568f6fc2 1.783516 1.792552 -0.0090363452 1.754310
## 106 7dd273cfe05a5c5e224f299f760201cf 1.786252 1.775399  0.0108530812 1.548894
## 107 7dd273cfe05a5c5e224f299f760201cf 1.768129 1.794519 -0.0263905878 1.777863
## 108 7e777c0f972309bbe9f7781b653b0662 1.751916 1.766730 -0.0148140596 1.445078
## 109 7ebbb1e6f79a2910951c191a4bf1df66 1.765851 1.763174  0.0026775726 1.402496
## 110 801efb14d6457f799fd2becee6d07b72 1.774518 1.763539  0.0109782887 1.406875
## 111 8215c539c392e8fc631e1d5cddeac7da 1.768166 1.764324  0.0038414873 1.416273
## 112 83c246c4ed456dea2a15d054f16009b3 1.768577 1.765813  0.0027635571 1.434103
## 113 84ea3564fd34c5377c4305b340f0c16f 1.781911 1.773885  0.0080261956 1.530763
## 114 84ea3564fd34c5377c4305b340f0c16f 1.802414 1.765192  0.0372221913 1.426668
## 115 84ea3564fd34c5377c4305b340f0c16f 1.780486 1.773475  0.0070106776 1.525854
## 116 8735f92a7f18db944bc862777fbd3222 1.807232 1.808400 -0.0011672482 1.944083
## 117 8735f92a7f18db944bc862777fbd3222 1.785641 1.765464  0.0201772891 1.429923
## 118 87a918f2af94ff451dd1f9ada0c46ecb 1.761538 1.786707 -0.0251686114 1.684306
## 119 8b08139d6eb5d7a849dff1bd579727d7 1.790274 1.773244  0.0170302595 1.523086
## 120 8b08139d6eb5d7a849dff1bd579727d7 1.796092 1.774939  0.0211525772 1.543387
## 121 8b1b6b50ee8e9d52ac980948823b441f 1.772619 1.764889  0.0077300245 1.423030
## 122 8b2494a11c2f75ef0302082e2b542682 1.791926 1.772784  0.0191420239 1.517574
## 123 8b2494a11c2f75ef0302082e2b542682 1.805022 1.774278  0.0307436591 1.535470
## 124 8d6ed6c4e133f0c9730f92c02bed365a 1.765435 1.797482 -0.0320468331 1.813347
## 125 8f813280e8cc1c972cdf97b409fa6840 1.765379 1.778989 -0.0136092607 1.591880
## 126 9100c673ddb0128f63023d0c6d423792 1.773645 1.772276  0.0013692731 1.511497
## 127 9176a4011d0b5a7a3fe654f25319091e 1.777094 1.774907  0.0021870419 1.543002
## 128 9176a4011d0b5a7a3fe654f25319091e 1.781262 1.772193  0.0090687538 1.510506
## 129 9255e939f667b42a03cb99246b86cb32 1.795627 1.778161  0.0174659010 1.581975
## 130 9255e939f667b42a03cb99246b86cb32 1.808657 1.781982  0.0266759145 1.627722
## 131 92c0c75c6858840cbc21744f047e2a02 1.784690 1.775093  0.0095973258 1.545226
## 132 92c0c75c6858840cbc21744f047e2a02 1.785245 1.774293  0.0109526558 1.535645
## 133 936f76885e2ff79dde9adf042e7d4617 1.765607 1.750175  0.0154320355 1.246828
## 134 97526cc3bdb30750d172ec27041d3211 1.762310 1.757536  0.0047740477 1.334978
## 135 98f53519513214e968c569829a1471bb 1.773201 1.762250  0.0109509533 1.391436
## 136 9a64f65cc4a51bc87a27dbff7ecefdd6 1.775364 1.774554  0.0008102825 1.538773
## 137 9a64f65cc4a51bc87a27dbff7ecefdd6 1.782228 1.773876  0.0083518417 1.530662
## 138 9ac9a274b64edbcb2b5e1e1c81743b51 1.771461 1.766904  0.0045569804 1.447169
## 139 9af3467db68cf6063627304cecd46a65 1.807723 1.775198  0.0325253858 1.546483
## 140 9af3467db68cf6063627304cecd46a65 1.794090 1.769408  0.0246815307 1.477152
## 141 9c54c522a386d5ed75686001dae803d9 1.767993 1.764795  0.0031979222 1.421913
## 142 9c886d46813a724862c055b0ed46e86d 1.789728 1.771525  0.0182020868 1.502508
## 143 9c886d46813a724862c055b0ed46e86d 1.802888 1.773860  0.0290278633 1.530463
## 144 9cfa52be027c8fa57d197b21dd7a958c 1.781277 1.780503  0.0007738351 1.610018
## 145 9cfa52be027c8fa57d197b21dd7a958c 1.786797 1.773778  0.0130189915 1.529479
## 146 9e26dc7c6642e4bf724f8d376d56a973 1.768423 1.764294  0.0041290341 1.415914
## 147 a4e7c4029a5970db9ae5de9cb1b8ced9 1.783503 1.778363  0.0051396474 1.584393
## 148 a51911fb468aee7787ca8edf2ad66641 1.789132 1.773488  0.0156433060 1.526015
## 149 a57e8133deb4fc1a6cbeaac02d482d8f 1.763822 1.769437 -0.0056144679 1.477496
## 150 a60f0b4bc20e5a162b45b3d2ed7a0d6c 1.768771 1.764049  0.0047216399 1.412976
## 151 a6db2d87f1e7a32216e5041a1920f12e 1.766711 1.762825  0.0038855569 1.398323
## 152 a7b57a3495513cb6d0fae02b7d400ac7 1.815897 1.776859  0.0390371021 1.566383
## 153 a7d9ac24cf597a17e98e400931d90089 1.762792 1.765255 -0.0024634895 1.427424
## 154 aca28e4f440b01d56bc951ba98270a87 1.761361 1.786355 -0.0249941496 1.680100
## 155 acc9603a819dc23a606ac6c1b9105b6d 1.769655 1.764542  0.0051122452 1.418883
## 156 ae88f77525a7096df7c0d9c6ede9f846 1.768593 1.765037  0.0035560073 1.424802
## 157 aed191bac0b1c3264208a084af848bed 1.808058 1.774352  0.0337054229 1.536359
## 158 aed191bac0b1c3264208a084af848bed 1.795637 1.773101  0.0225359817 1.521379
## 159 aee83626bf99f96165a5c90cc9511595 1.763031 1.787428 -0.0243970880 1.692947
## 160 af27577de5edb10c09f7fc19b9ca45ca 1.804732 1.774662  0.0300701524 1.540068
## 161 af27577de5edb10c09f7fc19b9ca45ca 1.764538 1.764905 -0.0003669592 1.423221
## 162 af27577de5edb10c09f7fc19b9ca45ca 1.790391 1.773551  0.0168403654 1.526766
## 163 afc0eeec83a181be740331928d883362 1.781474 1.773309  0.0081646985 1.523872
## 164 afc0eeec83a181be740331928d883362 1.780635 1.776195  0.0044402071 1.558430
## 165 b396df71ff4bad780a6862032851bce5 1.775030 1.768868  0.0061621349 1.470686
## 166 b6a0ad25677cccdc5a6ada23f7e276b6 1.781882 1.775858  0.0060242111 1.554387
## 167 b6a0ad25677cccdc5a6ada23f7e276b6 1.786160 1.779289  0.0068709657 1.595473
## 168 b7a2bc450146ff938fb3354ecc36444a 1.772969 1.793736 -0.0207666299 1.768483
## 169 bca82bc447bbd86d4a6a93e596e3aabd 1.768143 1.790625 -0.0224813000 1.731224
## 170 bca82bc447bbd86d4a6a93e596e3aabd 1.772377 1.776073 -0.0036952474 1.556961
## 171 bd217e3552c0912e166951d39a120601 1.792241 1.774955  0.0172863043 1.543573
## 172 bd217e3552c0912e166951d39a120601 1.794915 1.774485  0.0204296923 1.537949
## 173 be45460bb19ae26f588454f3dd7388ec 1.781901 1.777927  0.0039730767 1.579174
## 174 be9d6317e49313c0a8885815491ce6f0 1.790216 1.773776  0.0164405694 1.529457
## 175 be9d6317e49313c0a8885815491ce6f0 1.805817 1.774922  0.0308953909 1.543183
## 176 c1aff5b39b852ad01c74231c921af653 1.752676 1.762226 -0.0095496101 1.391144
## 177 c1d9eada909d2d5dbe078f25f1038508 1.754988 1.751960  0.0030281949 1.268211
## 178 c36646fabf5978777e0e7b08044e9568 1.784394 1.775056  0.0093380789 1.544791
## 179 c4e16e2ca85ff9eb360fd895f75227ca 1.786890 1.775016  0.0118739834 1.544313
## 180 c4e16e2ca85ff9eb360fd895f75227ca 1.787019 1.774504  0.0125150445 1.538175
## 181 c52287a7efd223c578b51937b849bf91 1.789177 1.773547  0.0156298887 1.526716
## 182 c52287a7efd223c578b51937b849bf91 1.809374 1.774009  0.0353647681 1.532252
## 183 c530fa2c576c67ce632ce34474b8ea71 1.788133 1.774386  0.0137463879 1.536766
## 184 c530fa2c576c67ce632ce34474b8ea71 1.790456 1.774752  0.0157044890 1.541142
## 185 c56ae4705fd5cf31ed0b25fa86e22ed6 1.791941 1.772741  0.0192008494 1.517059
## 186 c56ae4705fd5cf31ed0b25fa86e22ed6 1.804785 1.774326  0.0304586917 1.536051
## 187 c67744b93c1aac731ad718596675f25e 1.753119 1.761966 -0.0088462348 1.388026
## 188 c69352f1a000317323fd8841f1007691 1.762309 1.766217 -0.0039084684 1.438944
## 189 c9eeb9179e405b1fd7846a8aefb864d1 1.785315 1.782424  0.0028906284 1.633022
## 190 ca2b23c1b70c9034ea6412255f88eb7f 1.769616 1.762772  0.0068442347 1.397685
## 191 cf3381318906e1e08f58bca905636dc0 1.774512 1.764772  0.0097402241 1.421629
## 192 cfb98993000a7864ecd250723d57b8db 1.794956 1.774393  0.0205626114 1.536848
## 193 cfb98993000a7864ecd250723d57b8db 1.788018 1.773736  0.0142819234 1.528976
## 194 d1e7069c930bb4790713681f06f84cf0 1.792255 1.772845  0.0194098710 1.518309
## 195 d1e7069c930bb4790713681f06f84cf0 1.807038 1.774447  0.0325913268 1.537494
## 196 d295b9a63b6148b9f4663819da1c914b 1.768186 1.764238  0.0039480252 1.415240
## 197 d2af52fc6f48ab4d05c96dad9527d8a2 1.778422 1.766810  0.0116120360 1.446034
## 198 d4cc0a095dd016823ef8c5ad567c91fd 1.789423 1.774229  0.0151944099 1.534879
## 199 d4cc0a095dd016823ef8c5ad567c91fd 1.795821 1.774224  0.0215967760 1.534827
## 200 d52d14abfccb539e79068c89454cbbbd 1.769616 1.798219 -0.0286026129 1.822167
## 201 d56e5ac31cfe9e5dea7d4a63d13f1232 1.788757 1.773962  0.0147959729 1.531681
## 202 d56e5ac31cfe9e5dea7d4a63d13f1232 1.791763 1.773404  0.0183585449 1.525006
## 203 d5b23bb0b75b67df5ecd5ad44760b312 1.787401 1.773806  0.0135946392 1.529824
## 204 d5b23bb0b75b67df5ecd5ad44760b312 1.788806 1.775666  0.0131395799 1.552097
## 205 d5c62883f6e56836e55aefb9a15623a5 1.771251 1.765014  0.0062369187 1.424538
## 206 d5eb9156cd0592eb1cd773db88550ef3 1.767431 1.795134 -0.0277031292 1.785230
## 207 d94c20c00db235bb6fccc68f06f100c7 1.766792 1.761363  0.0054287345 1.380814
## 208 dc3835d12dc9edc210bddff9c66bf356 1.793618 1.772862  0.0207565910 1.518510
## 209 dc3835d12dc9edc210bddff9c66bf356 1.795754 1.775298  0.0204557747 1.547688
## 210 de05e9f7072b8da3d4f510900a0b870a 1.764921 1.796920 -0.0319982416 1.806608
## 211 e177a451db2600efdeb8bde47ef88671 1.813774 1.774451  0.0393231035 1.537542
## 212 e33c5576113465ae725d3993fc43738e 1.755417 1.764183 -0.0087661601 1.414580
## 213 e416d0916760d2fc17b616e2ac3ad855 1.769735 1.764098  0.0056370160 1.413563
## 214 e416d0916760d2fc17b616e2ac3ad855 1.785579 1.773750  0.0118288290 1.529147
## 215 e416d0916760d2fc17b616e2ac3ad855 1.777238 1.777711 -0.0004724842 1.576579
## 216 e4278f7eb9e1d1595866d1c4440c78c9 1.790571 1.771723  0.0188484985 1.504871
## 217 e4278f7eb9e1d1595866d1c4440c78c9 1.806858 1.773703  0.0331551025 1.528581
## 218 e50cc2a3da1dd7fd34ca4daba09a32ad 1.776983 1.763414  0.0135688417 1.405372
## 219 e61c5d15739059361f8ac80a09a88441 1.761113 1.764812 -0.0036987992 1.422114
## 220 e657907312415da47ee5d5a4796a675a 1.763183 1.762810  0.0003723722 1.398143
## 221 e6c8443214bca5d5809bbb0c5644af52 1.787413 1.780674  0.0067391508 1.612058
## 222 e6c8443214bca5d5809bbb0c5644af52 1.793248 1.775177  0.0180716444 1.546233
## 223 e739dae8872120bd9c58ff3f4afc28b3 1.764412 1.759981  0.0044315287 1.364260
## 224 e747694a818ae82e58f977992c196fd0 1.788071 1.772984  0.0150874321 1.519969
## 225 e747694a818ae82e58f977992c196fd0 1.786183 1.776111  0.0100723760 1.557419
## 226 eb2625c1474f98f554b002fcd5d1088e 1.795977 1.771577  0.0243995787 1.503126
## 227 eb2625c1474f98f554b002fcd5d1088e 1.817294 1.774930  0.0423644294 1.543279
## 228 efbe1f58b1e2984ddc53a64f047d94ff 1.810441 1.776811  0.0336291603 1.565809
## 229 efbe1f58b1e2984ddc53a64f047d94ff 1.780885 1.763664  0.0172204274 1.408368
## 230 f1e34308d952529287ada31976faec95 1.771693 1.761525  0.0101683712 1.382748
## 231 fb2dab56bf8047ade734c2fceccb49e4 1.780099 1.773141  0.0069582593 1.521850
## 232 fb800d4d7cd0a8e0da07b414fcf21b8c 1.803453 1.774556  0.0288965213 1.538803
## 233 fb800d4d7cd0a8e0da07b414fcf21b8c 1.797563 1.771619  0.0259437620 1.503627
## 234 fd3acdbb6eaf9e3669c9f303a6fa0b28 1.810100 1.774382  0.0357178489 1.536721
## 235 fecd0294f7ecdb7bc5c54d01e9f541f0 1.789626 1.775619  0.0140067064 1.551534
## 236 ffd660f96ba6668b2d8d86fc4e150862 1.794859 1.774402  0.0204571963 1.536953
## 237 ffd660f96ba6668b2d8d86fc4e150862 1.764209 1.747016  0.0171932756 1.209005
## 238 ffd660f96ba6668b2d8d86fc4e150862 1.781663 1.774070  0.0075931388 1.532982
##       Mlight Mheavymax     Mlab            A     A_CI_low    A_CI_high
## 1   308.3884  317.6617 308.2271 -0.017200093 -0.017200093 -0.017200093
## 2   308.3744  317.6618 308.4971  0.013067477  0.013067477  0.013067477
## 3   308.4810  317.6612 306.9020 -0.170088621 -0.170088621 -0.170088621
## 4   308.3966  317.6617 309.1353  0.078842871  0.078842871  0.078842871
## 5   308.5329  317.6609 304.2926 -0.459376798 -0.459376798 -0.459376798
## 6   308.3950  317.6617 309.2408  0.090262148  0.090262148  0.090262148
## 7   308.4501  317.6614 312.4364  0.427946387  0.427946387  0.427946387
## 8   308.4499  317.6614 310.8103  0.253396191  0.253396191  0.253396191
## 9   308.4613  317.6613 310.1221  0.178516562  0.178516562  0.178516562
## 10  308.4484  317.6614 311.5905  0.337267633  0.337267633  0.337267633
## 11  308.3013  317.6622 308.9003  0.063284948  0.063284948  0.063284948
## 12  308.3865  317.6617 309.0685  0.072710354  0.072710354  0.072710354
## 13  308.5016  317.6611 310.1435  0.177268760  0.177268760  0.177268760
## 14  308.4499  317.6614 311.4424  0.321258235  0.321258235  0.321258235
## 15  308.4289  317.6615 311.1861  0.295321786  0.295321786  0.295321786
## 16  308.3813  317.6618 309.0560  0.071896477  0.071896477  0.071896477
## 17  308.3882  317.6617 309.4921  0.117713536  0.117713536  0.117713536
## 18  308.4480  317.6614 311.1142  0.286169145  0.286169145  0.286169145
## 19  308.4589  317.6613 310.3777  0.206187637  0.206187637  0.206187637
## 20  308.4448  317.6614 310.6659  0.238316113  0.238316113  0.238316113
## 21  308.4579  317.6613 309.4499  0.106586350  0.106586350  0.106586350
## 22  308.3600  317.6619 310.8533  0.265069160  0.265069160  0.265069160
## 23  308.3876  317.6617 309.3491  0.102524384  0.102524384  0.102524384
## 24  308.4276  317.6615 312.4779  0.433761165  0.433761165  0.433761165
## 25  308.4581  317.6613 312.3956  0.423079536  0.423079536  0.423079536
## 26  308.4480  317.6614 311.8966  0.370140685  0.370140685  0.370140685
## 27  308.3244  317.6621 310.4163  0.221534335  0.221534335  0.221534335
## 28  308.5843  317.6607 303.1228 -0.595045801 -0.595045801 -0.595045801
## 29  308.3762  317.6618 309.7430  0.145559340  0.145559340  0.145559340
## 30  308.5425  317.6609 308.3329 -0.022731976 -0.022731976 -0.022731976
## 31  308.4529  317.6614 310.7683  0.248649039  0.248649039  0.248649039
## 32  308.5000  317.6611 308.9386  0.047346608  0.047346608  0.047346608
## 33  308.3513  317.6619 311.2567  0.308583713  0.308583713  0.308583713
## 34  308.4638  317.6613 311.1902  0.293130909  0.293130909  0.293130909
## 35  308.3867  317.6617 309.0773  0.073634663  0.073634663  0.073634663
## 36  308.3940  317.6617 307.7134 -0.072626060 -0.072626060 -0.072626060
## 37  308.5733  317.6607 303.9144 -0.506977801 -0.506977801 -0.506977801
## 38  308.4508  317.6614 311.6611  0.344673362  0.344673362  0.344673362
## 39  308.4521  317.6614 310.7547  0.247254348  0.247254348  0.247254348
## 40  308.4567  317.6614 313.3438  0.525036749  0.525036749  0.525036749
## 41  308.4450  317.6614 312.1498  0.397510620  0.397510620  0.397510620
## 42  308.4557  317.6614 310.1536  0.182391677  0.182391677  0.182391677
## 43  308.4659  317.6613 309.4551  0.106387519  0.106387519  0.106387519
## 44  308.3930  317.6617 310.5536  0.230521556  0.230521556  0.230521556
## 45  308.3762  317.6618 308.8427  0.049681475  0.049681475  0.049681475
## 46  308.4672  317.6613 307.7836 -0.073519862 -0.073519862 -0.073519862
## 47  308.4927  317.6612 309.0478  0.059877581  0.059877581  0.059877581
## 48  308.4510  317.6614 309.2597  0.086828910  0.086828910  0.086828910
## 49  308.4535  317.6614 309.3605  0.097403527  0.097403527  0.097403527
## 50  308.4520  317.6614 310.6149  0.232248523  0.232248523  0.232248523
## 51  308.3952  317.6617 307.4814 -0.097518314 -0.097518314 -0.097518314
## 52  308.3878  317.6617 310.4966  0.224859077  0.224859077  0.224859077
## 53  308.5275  317.6610 306.4189 -0.228296497 -0.228296497 -0.228296497
## 54  308.3943  317.6617 310.7965  0.256325291  0.256325291  0.256325291
## 55  308.5005  317.6611 309.2723  0.083314178  0.083314178  0.083314178
## 56  308.4677  317.6613 315.4235  0.748179279  0.748179279  0.748179279
## 57  308.4507  317.6614 310.5642  0.226909397  0.226909397  0.226909397
## 58  308.4521  317.6614 309.8258  0.147503745  0.147503745  0.147503745
## 59  308.4499  317.6614 310.9005  0.263077807  0.263077807  0.263077807
## 60  308.4315  317.6615 310.7934  0.253044951  0.253044951  0.253044951
## 61  308.4648  317.6613 312.4359  0.427004734  0.427004734  0.427004734
## 62  308.4498  317.6614 311.3303  0.309226560  0.309226560  0.309226560
## 63  308.4522  317.6614 312.4710  0.431538412  0.431538412  0.431538412
## 64  308.3986  317.6617 308.2293 -0.018075637 -0.018075637 -0.018075637
## 65  308.4055  317.6616 308.8173  0.043991893  0.043991893  0.043991893
## 66  308.4524  317.6614 312.8583  0.473115436  0.473115436  0.473115436
## 67  308.5652  317.6608 308.1154 -0.048909584 -0.048909584 -0.048909584
## 68  308.4675  317.6613 310.0753  0.172928280  0.172928280  0.172928280
## 69  308.3951  317.6617 309.7259  0.142019321  0.142019321  0.142019321
## 70  308.4491  317.6614 307.3226 -0.120920831 -0.120920831 -0.120920831
## 71  308.3784  317.6618 310.5370  0.229935719  0.229935719  0.229935719
## 72  308.4569  317.6614 309.5526  0.117716170  0.117716170  0.117716170
## 73  308.4589  317.6613 310.5342  0.223012417  0.223012417  0.223012417
## 74  308.4531  317.6614 311.8076  0.360248911  0.360248911  0.360248911
## 75  308.3847  317.6617 309.4868  0.117474838  0.117474838  0.117474838
## 76  308.5937  317.6606 303.2698 -0.580652449 -0.580652449 -0.580652449
## 77  308.4498  317.6614 309.7290  0.137327875  0.137327875  0.137327875
## 78  308.4412  317.6614 310.1931  0.187886499  0.187886499  0.187886499
## 79  308.3690  317.6618 309.3777  0.107343942  0.107343942  0.107343942
## 80  308.5419  317.6609 301.9831 -0.711250608 -0.711250608 -0.711250608
## 81  308.4462  317.6614 314.4870  0.648239878  0.648239878  0.648239878
## 82  308.4320  317.6615 312.9166  0.480496212  0.480496212  0.480496212
## 83  308.4614  317.6613 310.0633  0.172183601  0.172183601  0.172183601
## 84  308.4500  317.6614 310.6839  0.239817480  0.239817480  0.239817480
## 85  308.4573  317.6613 311.9198  0.372012968  0.372012968  0.372012968
## 86  308.4482  317.6614 310.1018  0.177492248  0.177492248  0.177492248
## 87  308.4521  317.6614 315.2564  0.730641143  0.730641143  0.730641143
## 88  308.4670  317.6613 309.7932  0.142641980  0.142641980  0.142641980
## 89  308.4380  317.6615 312.2218  0.405670680  0.405670680  0.405670680
## 90  308.4551  317.6614 313.5793  0.550414609  0.550414609  0.550414609
## 91  308.4022  317.6617 309.1872  0.083838125  0.083838125  0.083838125
## 92  308.3865  317.6617 309.7561  0.146017004  0.146017004  0.146017004
## 93  308.4495  317.6614 310.3302  0.201893760  0.201893760  0.201893760
## 94  308.4404  317.6614 310.7828  0.251205357  0.251205357  0.251205357
## 95  308.4872  317.6612 308.6302  0.015410732  0.015410732  0.015410732
## 96  308.4506  317.6614 309.4212  0.104200335  0.104200335  0.104200335
## 97  308.4488  317.6614 310.0022  0.166747161  0.166747161  0.166747161
## 98  308.4565  317.6614 314.5690  0.656669367  0.656669367  0.656669367
## 99  308.4910  317.6612 309.3645  0.094192655  0.094192655  0.094192655
## 100 308.4318  317.6615 312.8683  0.475334359  0.475334359  0.475334359
## 101 308.4202  317.6616 305.8430 -0.275777968 -0.275777968 -0.275777968
## 102 308.4499  317.6614 314.5726  0.657294712  0.657294712  0.657294712
## 103 308.4215  317.6615 310.5878  0.231843051  0.231843051  0.231843051
## 104 308.4788  317.6612 308.9883  0.054865181  0.054865181  0.054865181
## 105 308.5611  317.6608 307.0057 -0.169038224 -0.169038224 -0.169038224
## 106 308.4593  317.6613 310.3449  0.202635498  0.202635498  0.202635498
## 107 308.5728  317.6607 304.0349 -0.493790061 -0.493790061 -0.493790061
## 108 308.4078  317.6616 305.8218 -0.276345998 -0.276345998 -0.276345998
## 109 308.3866  317.6617 308.8550  0.049930975  0.049930975  0.049930975
## 110 308.3888  317.6617 310.3086  0.204728693  0.204728693  0.204728693
## 111 308.3935  317.6617 309.0649  0.071643431  0.071643431  0.071643431
## 112 308.4023  317.6617 308.8850  0.051547596  0.051547596  0.051547596
## 113 308.4503  317.6614 309.8459  0.149831637  0.149831637  0.149831637
## 114 308.3986  317.6617 314.9018  0.694249463  0.694249463  0.694249463
## 115 308.4478  317.6614 309.6671  0.130868556  0.130868556  0.130868556
## 116 308.6553  317.6603 308.4560 -0.021877919 -0.021877919 -0.021877919
## 117 308.4002  317.6617 311.9249  0.376346565  0.376346565  0.376346565
## 118 308.5264  317.6610 304.1803 -0.470498251 -0.470498251 -0.470498251
## 119 308.4465  317.6614 311.4088  0.317896822  0.317896822  0.317896822
## 120 308.4565  317.6614 312.1325  0.394916084  0.394916084  0.394916084
## 121 308.3968  317.6617 309.7476  0.144172236  0.144172236  0.144172236
## 122 308.4437  317.6614 311.7742  0.357299274  0.357299274  0.357299274
## 123 308.4526  317.6614 313.7973  0.573940807  0.573940807  0.573940807
## 124 308.5904  317.6606 303.0886 -0.599836897 -0.599836897 -0.599836897
## 125 308.4806  317.6612 306.1207 -0.254192541 -0.254192541 -0.254192541
## 126 308.4407  317.6614 308.6790  0.025557105  0.025557105  0.025557105
## 127 308.4563  317.6614 308.8364  0.040831674  0.040831674  0.040831674
## 128 308.4402  317.6614 310.0186  0.169264358  0.169264358  0.169264358
## 129 308.4757  317.6612 311.5056  0.326197549  0.326197549  0.326197549
## 130 308.4984  317.6611 313.1165  0.498412569  0.498412569  0.498412569
## 131 308.4574  317.6613 310.1252  0.179183807  0.179183807  0.179183807
## 132 308.4527  317.6614 310.3567  0.204470965  0.204470965  0.204470965
## 133 308.3094  317.6622 311.0279  0.287433204  0.287433204  0.287433204
## 134 308.3531  317.6619 309.1907  0.088978591  0.088978591  0.088978591
## 135 308.3812  317.6618 310.2975  0.204193703  0.204193703  0.204193703
## 136 308.4542  317.6614 308.5951  0.015127268  0.015127268  0.015127268
## 137 308.4502  317.6614 309.9025  0.155910606  0.155910606  0.155910606
## 138 308.4088  317.6616 309.2042  0.085008770  0.085008770  0.085008770
## 139 308.4581  317.6613 314.1097  0.607261510  0.607261510  0.607261510
## 140 308.4237  317.6615 312.7259  0.460539686  0.460539686  0.460539686
## 141 308.3963  317.6617 308.9551  0.059643720  0.059643720  0.059643720
## 142 308.4362  317.6615 311.6054  0.339710824  0.339710824  0.339710824
## 143 308.4501  317.6614 313.4977  0.541885740  0.541885740  0.541885740
## 144 308.4896  317.6612 308.6236  0.014455991  0.014455991  0.014455991
## 145 308.4496  317.6614 310.7135  0.243033585  0.243033585  0.243033585
## 146 308.3933  317.6617 309.1150  0.077005933  0.077005933  0.077005933
## 147 308.4769  317.6612 309.3684  0.095991421  0.095991421  0.095991421
## 148 308.4479  317.6614 311.1686  0.292014528  0.292014528  0.292014528
## 149 308.4238  317.6615 307.4452 -0.104762253 -0.104762253 -0.104762253
## 150 308.3918  317.6617 309.2173  0.088055862  0.088055862  0.088055862
## 151 308.3846  317.6617 309.0643  0.072454866  0.072454866  0.072454866
## 152 308.4679  317.6613 315.2449  0.728965549  0.728965549  0.728965549
## 153 308.3990  317.6617 307.9686 -0.045948041 -0.045948041 -0.045948041
## 154 308.5243  317.6610 304.2075 -0.467218299 -0.467218299 -0.467218299
## 155 308.3948  317.6617 309.2882  0.095344980  0.095344980  0.095344980
## 156 308.3977  317.6617 309.0190  0.066323841  0.066323841  0.066323841
## 157 308.4530  317.6614 314.3124  0.629237632  0.629237632  0.629237632
## 158 308.4456  317.6614 312.3659  0.420663667  0.420663667  0.420663667
## 159 308.5307  317.6609 304.3195 -0.456112911 -0.456112911 -0.456112911
## 160 308.4549  317.6614 313.6814  0.561389873  0.561389873  0.561389873
## 161 308.3969  317.6617 308.3328 -0.006844146 -0.006844146 -0.006844146
## 162 308.4483  317.6614 311.3771  0.314362151  0.314362151  0.314362151
## 163 308.4468  317.6614 309.8670  0.152408090  0.152408090  0.152408090
## 164 308.4640  317.6613 309.2351  0.082909086  0.082909086  0.082909086
## 165 308.4205  317.6616 309.4949  0.114974790  0.114974790  0.114974790
## 166 308.4620  317.6613 309.5084  0.112482151  0.112482151  0.112482151
## 167 308.4824  317.6612 309.6736  0.128339426  0.128339426  0.128339426
## 168 308.5682  317.6607 304.9958 -0.388524904 -0.388524904 -0.388524904
## 169 308.5497  317.6608 304.6758 -0.420451177 -0.420451177 -0.420451177
## 170 308.4633  317.6613 307.8215 -0.068998047 -0.068998047 -0.068998047
## 171 308.4566  317.6614 311.4607  0.322733757  0.322733757  0.322733757
## 172 308.4538  317.6614 312.0051  0.381401807  0.381401807  0.381401807
## 173 308.4743  317.6613 309.1636  0.074200321  0.074200321  0.074200321
## 174 308.4496  317.6614 311.3085  0.306906242  0.306906242  0.306906242
## 175 308.4564  317.6614 313.8256  0.576812223  0.576812223  0.576812223
## 176 308.3810  317.6618 306.7099 -0.178063560 -0.178063560 -0.178063560
## 177 308.3200  317.6621 308.8530  0.056411169  0.056411169  0.056411169
## 178 308.4572  317.6613 310.0799  0.174342956  0.174342956  0.174342956
## 179 308.4570  317.6614 310.5204  0.221687644  0.221687644  0.221687644
## 180 308.4539  317.6614 310.6294  0.233643754  0.233643754  0.233643754
## 181 308.4483  317.6614 311.1665  0.291765837  0.291765837  0.291765837
## 182 308.4510  317.6614 314.5999  0.660191932  0.660191932  0.660191932
## 183 308.4532  317.6614 310.8429  0.256628591  0.256628591  0.256628591
## 184 308.4554  317.6614 311.1849  0.293195179  0.293195179  0.293195179
## 185 308.4435  317.6614 311.7843  0.358395705  0.358395705  0.358395705
## 186 308.4529  317.6614 313.7479  0.568623740  0.568623740  0.568623740
## 187 308.3795  317.6618 306.8312 -0.164944221 -0.164944221 -0.164944221
## 188 308.4047  317.6616 307.7222 -0.072906089 -0.072906089 -0.072906089
## 189 308.5010  317.6611 309.0013  0.054011109  0.054011109  0.054011109
## 190 308.3843  317.6618 309.5816  0.127625355  0.127625355  0.127625355
## 191 308.3961  317.6617 310.0982  0.181662286  0.181662286  0.181662286
## 192 308.4533  317.6614 312.0278  0.383879593  0.383879593  0.383879593
## 193 308.4494  317.6614 310.9330  0.266608355  0.266608355  0.266608355
## 194 308.4441  317.6614 311.8211  0.362301119  0.362301119  0.362301119
## 195 308.4536  317.6614 314.1190  0.608444882  0.608444882  0.608444882
## 196 308.3930  317.6617 309.0831  0.073629741  0.073629741  0.073629741
## 197 308.4082  317.6616 310.4352  0.216616177  0.216616177  0.216616177
## 198 308.4523  317.6614 311.0939  0.283656770  0.283656770  0.283656770
## 199 308.4523  317.6614 312.2069  0.403179129  0.403179129  0.403179129
## 200 308.5948  317.6606 303.6863 -0.535417648 -0.535417648 -0.535417648
## 201 308.4507  317.6614 311.0234  0.276210873  0.276210873  0.276210873
## 202 308.4474  317.6614 311.6405  0.342697063  0.342697063  0.342697063
## 203 308.4498  317.6614 310.8138  0.253780315  0.253780315  0.253780315
## 204 308.4608  317.6613 310.7434  0.245333112  0.245333112  0.245333112
## 205 308.3976  317.6617 309.4873  0.116325833  0.116325833  0.116325833
## 206 308.5765  317.6607 303.8144 -0.518386889 -0.518386889 -0.518386889
## 207 308.3759  317.6618 309.3263  0.101216758  0.101216758  0.101216758
## 208 308.4442  317.6614 312.0554  0.387439421  0.387439421  0.387439421
## 209 308.4587  317.6613 312.0129  0.381921243  0.381921243  0.381921243
## 210 308.5871  317.6606 303.0920 -0.598886551 -0.598886551 -0.598886551
## 211 308.4536  317.6614 315.2892  0.734120191  0.734120191  0.734120191
## 212 308.3926  317.6617 306.8602 -0.163485950 -0.163485950 -0.163485950
## 213 308.3921  317.6617 309.3776  0.105127597  0.105127597  0.105127597
## 214 308.4495  317.6614 310.5065  0.220815446  0.220815446  0.220815446
## 215 308.4730  317.6613 308.3910 -0.008823808 -0.008823808 -0.008823808
## 216 308.4374  317.6615 311.7187  0.351782088  0.351782088  0.351782088
## 217 308.4492  317.6614 314.2149  0.618922032  0.618922032  0.618922032
## 218 308.3881  317.6617 310.7610  0.253035582  0.253035582  0.253035582
## 219 308.3964  317.6617 307.7500 -0.068985576 -0.068985576 -0.068985576
## 220 308.3845  317.6617 308.4496  0.006943700  0.006943700  0.006943700
## 221 308.4906  317.6612 309.6581  0.125896221  0.125896221  0.125896221
## 222 308.4579  317.6613 311.5981  0.337403840  0.337403840  0.337403840
## 223 308.3677  317.6618 309.1441  0.082613467  0.082613467  0.082613467
## 224 308.4449  317.6614 311.0697  0.281623292  0.281623292  0.281623292
## 225 308.4635  317.6613 310.2128  0.188073218  0.188073218  0.188073218
## 226 308.4366  317.6615 312.6846  0.455378814  0.455378814  0.455378814
## 227 308.4565  317.6614 315.8188  0.790938080  0.790938080  0.790938080
## 228 308.4676  317.6613 314.3059  0.627976290  0.627976290  0.627976290
## 229 308.3896  317.6617 311.4007  0.321139116  0.321139116  0.321139116
## 230 308.3768  317.6618 310.1569  0.189588442  0.189588442  0.189588442
## 231 308.4458  317.6614 309.6563  0.129885567  0.129885567  0.129885567
## 232 308.4542  317.6614 313.4771  0.539473006  0.539473006  0.539473006
## 233 308.4368  317.6615 312.9536  0.484200575  0.484200575  0.484200575
## 234 308.4532  317.6614 314.6623  0.666809185  0.666809185  0.666809185
## 235 308.4606  317.6613 310.8938  0.261522202  0.261522202  0.261522202
## 236 308.4533  317.6614 312.0095  0.381911968  0.381911968  0.381911968
## 237 308.2907  317.6623 311.3247  0.320151348  0.320151348  0.320151348
## 238 308.4514  317.6614 309.7716  0.141750137  0.141750137  0.141750137
##     Incorporator CM_control CM_treatment      delta_BD   Culture Substrate
## 1          FALSE   1.781355     1.777140 -0.0042149381 Chemostat     Multi
## 2           TRUE   1.772052     1.750714 -0.0213385838 Chemostat     Multi
## 3          FALSE   1.792544     1.772681 -0.0198627995     Batch    Single
## 4           TRUE   1.780523     1.764715 -0.0158078012 Chemostat     Multi
## 5          FALSE   1.793283     1.752116 -0.0411665445 Chemostat     Multi
## 6           TRUE   1.776392     1.758338 -0.0180538388 Chemostat     Multi
## 7           TRUE   1.785237     1.790742  0.0055048312     Batch    Single
## 8           TRUE   1.783239     1.788393  0.0051542404     Batch     Multi
## 9           TRUE   1.776284     1.782520  0.0062363222     Batch    Single
## 10          TRUE   1.782327     1.793629  0.0113019899     Batch     Multi
## 11          TRUE   1.770617     1.746097 -0.0245197991 Chemostat     Multi
## 12          TRUE   1.787969     1.758159 -0.0298104716 Chemostat     Multi
## 13          TRUE   1.796046     1.783490 -0.0125564627     Batch    Single
## 14          TRUE   1.784867     1.792952  0.0080850861     Batch     Multi
## 15          TRUE   1.784749     1.782692 -0.0020578227     Batch    Single
## 16          TRUE   1.777254     1.753735 -0.0235188618 Chemostat     Multi
## 17          TRUE   1.778736     1.765697 -0.0130383156 Chemostat     Multi
## 18          TRUE   1.781755     1.788203  0.0064481010     Batch     Multi
## 19          TRUE   1.787545     1.781036 -0.0065095078     Batch    Single
## 20          TRUE   1.779155     1.786008  0.0068529619     Batch     Multi
## 21          TRUE   1.787909     1.776830 -0.0110789062     Batch    Single
## 22          TRUE   1.772173     1.763971 -0.0082021345 Chemostat     Multi
## 23          TRUE   1.777855     1.764825 -0.0130300034 Chemostat     Multi
## 24          TRUE   1.779516     1.791658  0.0121421321     Batch    Single
## 25          TRUE   1.786422     1.790671  0.0042487242     Batch    Single
## 26          TRUE   1.781695     1.797147  0.0154519136     Batch     Multi
## 27          TRUE   1.765118     1.753679 -0.0114389832 Chemostat     Multi
## 28         FALSE   1.796444     1.754948 -0.0414961944 Chemostat     Multi
## 29          TRUE   1.788572     1.755946 -0.0326259590 Chemostat     Multi
## 30         FALSE   1.794149     1.777893 -0.0162565754     Batch    Single
## 31          TRUE   1.784415     1.783751 -0.0006640894     Batch     Multi
## 32          TRUE   1.796279     1.779815 -0.0164635622     Batch    Single
## 33          TRUE   1.764754     1.750224 -0.0145299224 Chemostat     Multi
## 34          TRUE   1.783647     1.786307  0.0026603963     Batch    Single
## 35          TRUE   1.783209     1.750622 -0.0325877389 Chemostat     Multi
## 36         FALSE   1.774822     1.744891 -0.0299313492 Chemostat     Multi
## 37         FALSE   1.795878     1.769330 -0.0265482840     Batch    Single
## 38          TRUE   1.787487     1.785461 -0.0020263175     Batch    Single
## 39          TRUE   1.784531     1.787939  0.0034081058     Batch     Multi
## 40          TRUE   1.784224     1.797631  0.0134065336     Batch    Single
## 41          TRUE   1.777254     1.797444  0.0201901923     Batch     Multi
## 42          TRUE   1.781120     1.780547 -0.0005738725     Batch    Single
## 43          TRUE   1.787935     1.778952 -0.0089827126     Batch    Single
## 44          TRUE   1.778892     1.776283 -0.0026097535 Chemostat     Multi
## 45          TRUE   1.782746     1.761869 -0.0208770848 Chemostat     Multi
## 46         FALSE   1.784578     1.773320 -0.0112581821     Batch    Single
## 47          TRUE   1.787149     1.782443 -0.0047055683     Batch     Multi
## 48          TRUE   1.791147     1.776811 -0.0143361647     Batch    Single
## 49          TRUE   1.790076     1.775795 -0.0142811720     Batch    Single
## 50          TRUE   1.781785     1.787599  0.0058133443     Batch     Multi
## 51         FALSE   1.780527     1.747581 -0.0329459455 Chemostat     Multi
## 52          TRUE   1.777518     1.775084 -0.0024340538 Chemostat     Multi
## 53         FALSE   1.796335     1.795021 -0.0013138861 Chemostat     Multi
## 54          TRUE   1.781107     1.787033  0.0059262452 Chemostat     Multi
## 55          TRUE   1.790988     1.781254 -0.0097339338     Batch    Single
## 56          TRUE   1.785176     1.817487  0.0323107530     Batch    Single
## 57          TRUE   1.783552     1.785040  0.0014876890     Batch     Multi
## 58          TRUE   1.786030     1.778209 -0.0078214551     Batch    Single
## 59          TRUE   1.779702     1.782024  0.0023216370     Batch    Single
## 60          TRUE   1.779251     1.787014  0.0077628874     Batch     Multi
## 61          TRUE   1.789525     1.796341  0.0068155345     Batch    Single
## 62          TRUE   1.782533     1.791047  0.0085145883     Batch     Multi
## 63          TRUE   1.784804     1.791051  0.0062465371     Batch    Single
## 64         FALSE   1.781100     1.764631 -0.0164691843     Batch    Single
## 65          TRUE   1.785123     1.758150 -0.0269729693 Chemostat     Multi
## 66          TRUE   1.786616     1.794986  0.0083701132     Batch    Single
## 67         FALSE   1.791733     1.789143 -0.0025898227     Batch     Multi
## 68          TRUE   1.775313     1.783569  0.0082555807     Batch    Single
## 69          TRUE   1.776808     1.776205 -0.0006026151 Chemostat     Multi
## 70         FALSE   1.783179     1.769275 -0.0139043539     Batch    Single
## 71          TRUE   1.782767     1.772948 -0.0098193820 Chemostat     Multi
## 72          TRUE   1.784234     1.777807 -0.0064261993     Batch    Single
## 73          TRUE   1.783960     1.787048  0.0030876283     Batch     Multi
## 74          TRUE   1.787938     1.785996 -0.0019414832     Batch    Single
## 75          TRUE   1.778551     1.763311 -0.0152398154 Chemostat     Multi
## 76         FALSE   1.800367     1.768929 -0.0314379168     Batch    Single
## 77          TRUE   1.783351     1.776549 -0.0068024100     Batch    Single
## 78          TRUE   1.786654     1.779230 -0.0074236478     Batch     Multi
## 79          TRUE   1.767591     1.751257 -0.0163346548 Chemostat     Multi
## 80         FALSE   1.805978     1.741164 -0.0648134652 Chemostat     Multi
## 81          TRUE   1.779760     1.807556  0.0277953763     Batch    Single
## 82          TRUE   1.772820     1.800612  0.0277912391     Batch     Multi
## 83          TRUE   1.775695     1.782787  0.0070917804     Batch    Single
## 84          TRUE   1.780421     1.786710  0.0062890355     Batch     Multi
## 85          TRUE   1.788102     1.787357 -0.0007445122     Batch    Single
## 86          TRUE   1.779557     1.779835  0.0002785927     Batch     Multi
## 87          TRUE   1.787736     1.815195  0.0274593297     Batch    Single
## 88          TRUE   1.793053     1.778459 -0.0145948023     Batch    Single
## 89          TRUE   1.776519     1.798303  0.0217833384     Batch     Multi
## 90          TRUE   1.784982     1.799027  0.0140456175     Batch    Single
## 91          TRUE   1.786190     1.765287 -0.0209035952 Chemostat     Multi
## 92          TRUE   1.773521     1.768969 -0.0045522552 Chemostat     Multi
## 93          TRUE   1.786146     1.779344 -0.0068019443     Batch    Single
## 94          TRUE   1.784509     1.784348 -0.0001607221     Batch     Multi
## 95          TRUE   1.786376     1.781511 -0.0048645775     Batch     Multi
## 96          TRUE   1.782102     1.774840 -0.0072619975     Batch    Single
## 97          TRUE   1.782667     1.779330 -0.0033369969     Batch     Multi
## 98          TRUE   1.787641     1.812560  0.0249184619     Batch    Single
## 99          TRUE   1.787593     1.785388 -0.0022044568     Batch     Multi
## 100         TRUE   1.792043     1.789061 -0.0029815343     Batch    Single
## 101        FALSE   1.790170     1.746007 -0.0441630197 Chemostat     Multi
## 102         TRUE   1.789907     1.809756  0.0198481780     Batch    Single
## 103         TRUE   1.776997     1.778003  0.0010063278     Batch     Multi
## 104         TRUE   1.786447     1.779016 -0.0074310557     Batch    Single
## 105        FALSE   1.791312     1.783431 -0.0078808992     Batch     Multi
## 106         TRUE   1.778288     1.785682  0.0073938189     Batch     Multi
## 107        FALSE   1.795924     1.770350 -0.0255738759     Batch    Single
## 108        FALSE   1.786586     1.737276 -0.0493095114 Chemostat     Multi
## 109         TRUE   1.776728     1.751513 -0.0252145997 Chemostat     Multi
## 110         TRUE   1.776171     1.776970  0.0007993260 Chemostat     Multi
## 111         TRUE   1.788059     1.747369 -0.0406900681 Chemostat     Multi
## 112         TRUE   1.778421     1.762814 -0.0156072721 Chemostat     Multi
## 113         TRUE   1.779851     1.779600 -0.0002508908     Batch     Multi
## 114         TRUE   1.810029     1.804158 -0.0058713894 Chemostat     Multi
## 115         TRUE   1.788062     1.776587 -0.0114755067     Batch    Single
## 116        FALSE   1.807215     1.804697 -0.0025175425     Batch    Single
## 117         TRUE   1.769814     1.790160  0.0203461516     Batch     Multi
## 118        FALSE   1.790280     1.749496 -0.0407836642 Chemostat     Multi
## 119         TRUE   1.780908     1.791040  0.0101320148     Batch     Multi
## 120         TRUE   1.792001     1.789001 -0.0030001559     Batch    Single
## 121         TRUE   1.781883     1.773368 -0.0085147492 Chemostat     Multi
## 122         TRUE   1.778583     1.794079  0.0154959567     Batch     Multi
## 123         TRUE   1.788657     1.801720  0.0130634634     Batch    Single
## 124        FALSE   1.798169     1.755262 -0.0429070061 Chemostat     Multi
## 125        FALSE   1.787073     1.753903 -0.0331692393 Chemostat     Multi
## 126         TRUE   1.787316     1.775120 -0.0121954458     Batch    Single
## 127         TRUE   1.789063     1.775358 -0.0137047347     Batch    Single
## 128         TRUE   1.783147     1.778577 -0.0045701315     Batch     Multi
## 129         TRUE   1.786226     1.794300  0.0080739176     Batch     Multi
## 130         TRUE   1.787956     1.809827  0.0218705111     Batch    Single
## 131         TRUE   1.783086     1.781520 -0.0015658807     Batch     Multi
## 132         TRUE   1.786139     1.780522 -0.0056166150     Batch    Single
## 133         TRUE   1.799947     1.768688 -0.0312588990 Chemostat     Multi
## 134         TRUE   1.772325     1.748178 -0.0241467502 Chemostat     Multi
## 135         TRUE   1.785359     1.763863 -0.0214965534 Chemostat     Multi
## 136         TRUE   1.789464     1.773663 -0.0158009333     Batch    Single
## 137         TRUE   1.784293     1.779176 -0.0051174141     Batch     Multi
## 138         TRUE   1.790032     1.774310 -0.0157217061     Batch    Single
## 139         TRUE   1.793676     1.808503  0.0148271645     Batch    Single
## 140         TRUE   1.770182     1.797372  0.0271904831     Batch     Multi
## 141         TRUE   1.775722     1.759109 -0.0166122123 Chemostat     Multi
## 142         TRUE   1.774412     1.793323  0.0189103786     Batch     Multi
## 143         TRUE   1.789167     1.798519  0.0093516261     Batch    Single
## 144         TRUE   1.786440     1.778349 -0.0080911277     Batch     Multi
## 145         TRUE   1.789527     1.781026 -0.0085009646     Batch    Single
## 146         TRUE   1.777545     1.766641 -0.0109046269 Chemostat     Multi
## 147         TRUE   1.787969     1.780129 -0.0078406368     Batch    Single
## 148         TRUE   1.788156     1.782614 -0.0055414289     Batch    Single
## 149        FALSE   1.785005     1.755960 -0.0290451539 Chemostat     Multi
## 150         TRUE   1.769756     1.752840 -0.0169159274 Chemostat     Multi
## 151         TRUE   1.772975     1.756166 -0.0168090830 Chemostat     Multi
## 152         TRUE   1.790305     1.817349  0.0270441674     Batch    Single
## 153        FALSE   1.792623     1.752457 -0.0401660816 Chemostat     Multi
## 154        FALSE   1.789748     1.752115 -0.0376333665 Chemostat     Multi
## 155         TRUE   1.774695     1.770590 -0.0041046819 Chemostat     Multi
## 156         TRUE   1.787541     1.759497 -0.0280443249 Chemostat     Multi
## 157         TRUE   1.786079     1.807931  0.0218526517     Batch    Single
## 158         TRUE   1.779373     1.797552  0.0181785589     Batch     Multi
## 159        FALSE   1.796400     1.753465 -0.0429348047 Chemostat     Multi
## 160         TRUE   1.789688     1.801516  0.0118275925     Batch    Single
## 161        FALSE   1.803863     1.806958  0.0030947928 Chemostat     Multi
## 162         TRUE   1.782351     1.791155  0.0088043988     Batch     Multi
## 163         TRUE   1.786552     1.778003 -0.0085488556     Batch     Multi
## 164         TRUE   1.785471     1.776438 -0.0090323560     Batch    Single
## 165         TRUE   1.780259     1.772941 -0.0073181891 Chemostat     Multi
## 166         TRUE   1.792870     1.777858 -0.0150126847     Batch    Single
## 167         TRUE   1.784652     1.784428 -0.0002238086     Batch     Multi
## 168        FALSE   1.792191     1.771972 -0.0202191453     Batch    Single
## 169        FALSE   1.790041     1.775642 -0.0143994123     Batch     Multi
## 170        FALSE   1.783673     1.775068 -0.0086056125     Batch    Single
## 171         TRUE   1.784771     1.788113  0.0033417089     Batch    Single
## 172         TRUE   1.782045     1.799318  0.0172729901     Batch     Multi
## 173         TRUE   1.795316     1.776945 -0.0183704237     Batch    Single
## 174         TRUE   1.782669     1.791028  0.0083592054     Batch     Multi
## 175         TRUE   1.789783     1.802750  0.0129675492     Batch    Single
## 176        FALSE   1.772780     1.757567 -0.0152128875 Chemostat     Multi
## 177         TRUE   1.767191     1.751841 -0.0153504583 Chemostat     Multi
## 178         TRUE   1.787656     1.779601 -0.0080546032     Batch    Single
## 179         TRUE   1.784726     1.781452 -0.0032731897     Batch    Single
## 180         TRUE   1.781333     1.785610  0.0042768573     Batch     Multi
## 181         TRUE   1.781595     1.788668  0.0070730540     Batch     Multi
## 182         TRUE   1.790412     1.806991  0.0165784675     Batch    Single
## 183         TRUE   1.781001     1.787596  0.0065953421     Batch     Multi
## 184         TRUE   1.788698     1.784432 -0.0042653760     Batch    Single
## 185         TRUE   1.778309     1.794015  0.0157060464     Batch     Multi
## 186         TRUE   1.788036     1.801183  0.0131472411     Batch    Single
## 187        FALSE   1.760443     1.744900 -0.0155429923 Chemostat     Multi
## 188        FALSE   1.781185     1.753604 -0.0275805876 Chemostat     Multi
## 189         TRUE   1.797811     1.780988 -0.0168232461     Batch    Single
## 190         TRUE   1.782139     1.758609 -0.0235297829 Chemostat     Multi
## 191         TRUE   1.779609     1.779349 -0.0002593902 Chemostat     Multi
## 192         TRUE   1.786596     1.788459  0.0018627783     Batch    Single
## 193         TRUE   1.783408     1.788186  0.0047786819     Batch     Multi
## 194         TRUE   1.779201     1.796200  0.0169996161     Batch     Multi
## 195         TRUE   1.782945     1.804799  0.0218543643     Batch    Single
## 196         TRUE   1.777019     1.761700 -0.0153189377 Chemostat     Multi
## 197         TRUE   1.780599     1.782131  0.0015319846 Chemostat     Multi
## 198         TRUE   1.782309     1.789961  0.0076517758     Batch     Multi
## 199         TRUE   1.786598     1.788282  0.0016843229     Batch    Single
## 200        FALSE   1.801601     1.770815 -0.0307859955     Batch    Single
## 201         TRUE   1.789642     1.784168 -0.0054743307     Batch    Single
## 202         TRUE   1.779234     1.789903  0.0106688796     Batch     Multi
## 203         TRUE   1.781677     1.786339  0.0046620027     Batch     Multi
## 204         TRUE   1.786903     1.782719 -0.0041834525     Batch    Single
## 205         TRUE   1.782487     1.768879 -0.0136083249 Chemostat     Multi
## 206        FALSE   1.795987     1.769461 -0.0265265836     Batch    Single
## 207         TRUE   1.790803     1.752686 -0.0381170813 Chemostat     Multi
## 208         TRUE   1.778713     1.794673  0.0159595187     Batch     Multi
## 209         TRUE   1.786906     1.789242  0.0023364775     Batch    Single
## 210        FALSE   1.797244     1.749862 -0.0473818198 Chemostat     Multi
## 211         TRUE   1.774605     1.812899  0.0382938903     Batch    Single
## 212        FALSE   1.776889     1.748690 -0.0281992160 Chemostat     Multi
## 213         TRUE   1.776861     1.766028 -0.0108329396 Chemostat     Multi
## 214         TRUE   1.783373     1.783867  0.0004942590     Batch     Multi
## 215        FALSE   1.796538     1.775053 -0.0214855913     Batch    Single
## 216         TRUE   1.774101     1.795335  0.0212333139     Batch     Multi
## 217         TRUE   1.784807     1.805272  0.0204641614     Batch    Single
## 218         TRUE   1.756980     1.775042  0.0180627947 Chemostat     Multi
## 219        FALSE   1.781096     1.746040 -0.0350567362 Chemostat     Multi
## 220        FALSE   1.773166     1.749510 -0.0236564221 Chemostat     Multi
## 221         TRUE   1.787424     1.786712 -0.0007118287     Batch     Multi
## 222         TRUE   1.787900     1.786920 -0.0009796940     Batch    Single
## 223         TRUE   1.776137     1.758460 -0.0176762131 Chemostat     Multi
## 224         TRUE   1.779544     1.787511  0.0079672151     Batch     Multi
## 225         TRUE   1.792732     1.780136 -0.0125958718     Batch    Single
## 226         TRUE   1.775278     1.797230  0.0219525749     Batch     Multi
## 227         TRUE   1.775624     1.818268  0.0426439835     Batch    Single
## 228         TRUE   1.780994     1.811500  0.0305060468     Batch    Single
## 229         TRUE   1.775660     1.794871  0.0192105722 Chemostat     Multi
## 230         TRUE   1.770319     1.769095 -0.0012234239 Chemostat     Multi
## 231         TRUE   1.783029     1.776169 -0.0068604670     Batch    Single
## 232         TRUE   1.784822     1.799609  0.0147867074     Batch    Single
## 233         TRUE   1.774394     1.800351  0.0259566264     Batch     Multi
## 234         TRUE   1.774473     1.809001  0.0345283822     Batch    Single
## 235         TRUE   1.779876     1.781457  0.0015810025     Batch    Single
## 236         TRUE   1.787062     1.787530  0.0004676971     Batch    Single
## 237         TRUE   1.774165     1.756932 -0.0172328921 Chemostat     Multi
## 238         TRUE   1.785133     1.778820 -0.0063126145     Batch     Multi
##     copy_number         id   Domain             Phylum               Class
## 1             2 chemo_carb Bacteria     Proteobacteria Deltaproteobacteria
## 2             3 chemo_carb Bacteria     Proteobacteria Alphaproteobacteria
## 3             3  batch_glu Bacteria     Proteobacteria Gammaproteobacteria
## 4             2 chemo_carb Bacteria     Proteobacteria Gammaproteobacteria
## 5             1 chemo_carb Bacteria     Planctomycetes               OM190
## 6             2 chemo_carb Bacteria     Proteobacteria Alphaproteobacteria
## 7             2  batch_glu Bacteria      Bacteroidetes         Bacteroidia
## 8             2 batch_carb Bacteria      Bacteroidetes         Bacteroidia
## 9             3  batch_glu Bacteria     Actinobacteria      Actinobacteria
## 10            3 batch_carb Bacteria     Actinobacteria      Actinobacteria
## 11            1 chemo_carb Bacteria    Verrucomicrobia    Verrucomicrobiae
## 12            2 chemo_carb Bacteria     Proteobacteria Alphaproteobacteria
## 13            3  batch_glu Bacteria      Bacteroidetes         Bacteroidia
## 14            3 batch_carb Bacteria      Bacteroidetes         Bacteroidia
## 15            1  batch_glu Bacteria      Bacteroidetes         Bacteroidia
## 16            3 chemo_carb Bacteria     Proteobacteria Gammaproteobacteria
## 17            2 chemo_carb Bacteria     Proteobacteria Alphaproteobacteria
## 18            1 batch_carb Bacteria     Planctomycetes    Planctomycetacia
## 19            1  batch_glu Bacteria     Planctomycetes    Planctomycetacia
## 20            1 batch_carb Bacteria     Planctomycetes    Planctomycetacia
## 21            1  batch_glu Bacteria     Planctomycetes    Planctomycetacia
## 22            2 chemo_carb Bacteria     Proteobacteria Alphaproteobacteria
## 23            2 chemo_carb Bacteria     Proteobacteria Gammaproteobacteria
## 24            3  batch_glu Bacteria    Verrucomicrobia    Verrucomicrobiae
## 25            2  batch_glu Bacteria      Bacteroidetes         Bacteroidia
## 26            2 batch_carb Bacteria      Bacteroidetes         Bacteroidia
## 27            2 chemo_carb Bacteria     Proteobacteria Gammaproteobacteria
## 28            1 chemo_carb Bacteria     Proteobacteria Alphaproteobacteria
## 29            2 chemo_carb Bacteria     Proteobacteria Gammaproteobacteria
## 30            2  batch_glu Bacteria        Chloroflexi        Chloroflexia
## 31            1 batch_carb Bacteria      Bacteroidetes         Bacteroidia
## 32            1  batch_glu Bacteria      Bacteroidetes         Bacteroidia
## 33            3 chemo_carb Bacteria    Verrucomicrobia    Verrucomicrobiae
## 34            6  batch_glu Bacteria     Actinobacteria      Actinobacteria
## 35            3 chemo_carb Bacteria     Proteobacteria Gammaproteobacteria
## 36            1 chemo_carb Bacteria    Verrucomicrobia    Verrucomicrobiae
## 37            1  batch_glu Bacteria     Proteobacteria Alphaproteobacteria
## 38            3  batch_glu Bacteria      Bacteroidetes         Bacteroidia
## 39            3 batch_carb Bacteria      Bacteroidetes         Bacteroidia
## 40            2  batch_glu Bacteria     Proteobacteria Alphaproteobacteria
## 41            2 batch_carb Bacteria     Proteobacteria Alphaproteobacteria
## 42            1  batch_glu Bacteria     Planctomycetes    Planctomycetacia
## 43            2  batch_glu Bacteria     Proteobacteria Alphaproteobacteria
## 44            3 chemo_carb Bacteria      Bacteroidetes         Bacteroidia
## 45            1 chemo_carb Bacteria     Planctomycetes           vadinHA49
## 46            1  batch_glu Bacteria     Actinobacteria      Acidimicrobiia
## 47            2 batch_carb Bacteria     Proteobacteria Gammaproteobacteria
## 48            2  batch_glu Bacteria     Proteobacteria Gammaproteobacteria
## 49            2  batch_glu Bacteria     Proteobacteria Alphaproteobacteria
## 50            2 batch_carb Bacteria     Proteobacteria Alphaproteobacteria
## 51            2 chemo_carb Bacteria     Proteobacteria Alphaproteobacteria
## 52            4 chemo_carb Bacteria     Proteobacteria Gammaproteobacteria
## 53            4 chemo_carb Bacteria      Bacteroidetes         Bacteroidia
## 54            3 chemo_carb Bacteria     Proteobacteria Alphaproteobacteria
## 55            1  batch_glu Bacteria      Bacteroidetes        Rhodothermia
## 56            8  batch_glu Bacteria     Proteobacteria Gammaproteobacteria
## 57            3 batch_carb Bacteria      Bacteroidetes         Bacteroidia
## 58            3  batch_glu Bacteria      Bacteroidetes         Bacteroidia
## 59            3  batch_glu Bacteria      Bacteroidetes         Bacteroidia
## 60            4 batch_carb Bacteria      Bacteroidetes         Bacteroidia
## 61            4  batch_glu Bacteria      Bacteroidetes         Bacteroidia
## 62            2 batch_carb Bacteria      Bacteroidetes         Bacteroidia
## 63            2  batch_glu Bacteria      Bacteroidetes         Bacteroidia
## 64            3  batch_glu Bacteria     Actinobacteria      Actinobacteria
## 65            1 chemo_carb Bacteria     Planctomycetes               OM190
## 66            2  batch_glu Bacteria      Bacteroidetes         Bacteroidia
## 67            2 batch_carb Bacteria      Bacteroidetes         Bacteroidia
## 68            2  batch_glu Bacteria        Chloroflexi        Chloroflexia
## 69            4 chemo_carb Bacteria      Bacteroidetes         Bacteroidia
## 70            1  batch_glu Bacteria     Proteobacteria Deltaproteobacteria
## 71            1 chemo_carb Bacteria    Verrucomicrobia    Verrucomicrobiae
## 72            1  batch_glu Bacteria     Planctomycetes    Planctomycetacia
## 73            1 batch_carb Bacteria     Actinobacteria      Actinobacteria
## 74            1  batch_glu Bacteria     Actinobacteria      Actinobacteria
## 75            2 chemo_carb Bacteria      Bacteroidetes         Bacteroidia
## 76            1  batch_glu Bacteria     Proteobacteria                    
## 77            1  batch_glu Bacteria     Proteobacteria Alphaproteobacteria
## 78            1 batch_carb Bacteria     Proteobacteria Alphaproteobacteria
## 79            2 chemo_carb Bacteria     Proteobacteria Gammaproteobacteria
## 80            2 chemo_carb Bacteria     Proteobacteria Alphaproteobacteria
## 81            5  batch_glu Bacteria      Bacteroidetes         Bacteroidia
## 82            5 batch_carb Bacteria      Bacteroidetes         Bacteroidia
## 83            2  batch_glu Bacteria     Planctomycetes    Planctomycetacia
## 84            1 batch_carb Bacteria     Actinobacteria      Actinobacteria
## 85            1  batch_glu Bacteria     Actinobacteria      Actinobacteria
## 86            4 batch_carb Bacteria     Proteobacteria Gammaproteobacteria
## 87            4  batch_glu Bacteria     Proteobacteria Gammaproteobacteria
## 88            4  batch_glu Bacteria     Proteobacteria Alphaproteobacteria
## 89            4 batch_carb Bacteria      Bacteroidetes         Bacteroidia
## 90            4  batch_glu Bacteria      Bacteroidetes         Bacteroidia
## 91            1 chemo_carb Bacteria     Proteobacteria Alphaproteobacteria
## 92            1 chemo_carb Bacteria     Proteobacteria Alphaproteobacteria
## 93            1  batch_glu Bacteria    Verrucomicrobia    Verrucomicrobiae
## 94            1 batch_carb Bacteria    Verrucomicrobia    Verrucomicrobiae
## 95            2 batch_carb Bacteria     Actinobacteria      Acidimicrobiia
## 96            2  batch_glu Bacteria     Actinobacteria      Acidimicrobiia
## 97            4 batch_carb Bacteria     Proteobacteria Gammaproteobacteria
## 98            4  batch_glu Bacteria     Proteobacteria Gammaproteobacteria
## 99            1 batch_carb Bacteria    Verrucomicrobia    Verrucomicrobiae
## 100           1  batch_glu Bacteria    Verrucomicrobia    Verrucomicrobiae
## 101           2 chemo_carb Bacteria     Proteobacteria Deltaproteobacteria
## 102           9  batch_glu Bacteria     Proteobacteria Alphaproteobacteria
## 103           9 batch_carb Bacteria     Proteobacteria Alphaproteobacteria
## 104           5  batch_glu Bacteria     Proteobacteria Gammaproteobacteria
## 105           5 batch_carb Bacteria     Proteobacteria Gammaproteobacteria
## 106           4 batch_carb Bacteria     Proteobacteria Gammaproteobacteria
## 107           4  batch_glu Bacteria     Proteobacteria Gammaproteobacteria
## 108           3 chemo_carb Bacteria      Bacteroidetes         Bacteroidia
## 109           4 chemo_carb Bacteria     Proteobacteria Gammaproteobacteria
## 110           5 chemo_carb Bacteria     Proteobacteria Alphaproteobacteria
## 111           3 chemo_carb Bacteria      Bacteroidetes         Bacteroidia
## 112           2 chemo_carb Bacteria     Proteobacteria Gammaproteobacteria
## 113           1 batch_carb Bacteria     Proteobacteria Alphaproteobacteria
## 114           1 chemo_carb Bacteria     Proteobacteria Alphaproteobacteria
## 115           1  batch_glu Bacteria     Proteobacteria Alphaproteobacteria
## 116           2  batch_glu Bacteria     Proteobacteria Gammaproteobacteria
## 117           2 batch_carb Bacteria     Proteobacteria Gammaproteobacteria
## 118           4 chemo_carb Bacteria     Proteobacteria Gammaproteobacteria
## 119           2 batch_carb Bacteria     Proteobacteria Gammaproteobacteria
## 120           2  batch_glu Bacteria     Proteobacteria Gammaproteobacteria
## 121           3 chemo_carb Bacteria      Bacteroidetes        Rhodothermia
## 122           2 batch_carb Bacteria     Proteobacteria Gammaproteobacteria
## 123           2  batch_glu Bacteria     Proteobacteria Gammaproteobacteria
## 124           4 chemo_carb Bacteria     Proteobacteria Alphaproteobacteria
## 125           2 chemo_carb Bacteria    Verrucomicrobia    Verrucomicrobiae
## 126           4  batch_glu Bacteria     Proteobacteria Gammaproteobacteria
## 127           2  batch_glu Bacteria     Proteobacteria Alphaproteobacteria
## 128           2 batch_carb Bacteria     Proteobacteria Alphaproteobacteria
## 129           4 batch_carb Bacteria     Proteobacteria Alphaproteobacteria
## 130           4  batch_glu Bacteria     Proteobacteria Alphaproteobacteria
## 131           2 batch_carb Bacteria     Planctomycetes    Planctomycetacia
## 132           2  batch_glu Bacteria     Planctomycetes    Planctomycetacia
## 133           2 chemo_carb Bacteria     Proteobacteria Alphaproteobacteria
## 134           2 chemo_carb Bacteria     Proteobacteria Alphaproteobacteria
## 135           4 chemo_carb Bacteria     Proteobacteria Alphaproteobacteria
## 136           3  batch_glu Bacteria     Proteobacteria Gammaproteobacteria
## 137           3 batch_carb Bacteria     Proteobacteria Gammaproteobacteria
## 138           3  batch_glu Bacteria      Bacteroidetes         Bacteroidia
## 139           5  batch_glu Bacteria     Proteobacteria Gammaproteobacteria
## 140           5 batch_carb Bacteria     Proteobacteria Gammaproteobacteria
## 141           2 chemo_carb Bacteria     Proteobacteria Gammaproteobacteria
## 142           2 batch_carb Bacteria      Bacteroidetes         Bacteroidia
## 143           2  batch_glu Bacteria      Bacteroidetes         Bacteroidia
## 144           1 batch_carb Bacteria     Actinobacteria      Actinobacteria
## 145           1  batch_glu Bacteria     Actinobacteria      Actinobacteria
## 146           2 chemo_carb Bacteria     Proteobacteria Alphaproteobacteria
## 147           1  batch_glu Bacteria     Planctomycetes    Planctomycetacia
## 148           1  batch_glu Bacteria     Actinobacteria      Actinobacteria
## 149           2 chemo_carb Bacteria     Proteobacteria Deltaproteobacteria
## 150           2 chemo_carb Bacteria     Proteobacteria Alphaproteobacteria
## 151           2 chemo_carb Bacteria     Proteobacteria Deltaproteobacteria
## 152           2  batch_glu Bacteria     Proteobacteria Gammaproteobacteria
## 153           5 chemo_carb Bacteria     Proteobacteria Gammaproteobacteria
## 154           2 chemo_carb Bacteria     Actinobacteria      Actinobacteria
## 155           3 chemo_carb Bacteria     Proteobacteria Gammaproteobacteria
## 156           2 chemo_carb Bacteria    Verrucomicrobia    Verrucomicrobiae
## 157           2  batch_glu Bacteria     Proteobacteria Gammaproteobacteria
## 158           2 batch_carb Bacteria     Proteobacteria Gammaproteobacteria
## 159           2 chemo_carb Bacteria     Proteobacteria Gammaproteobacteria
## 160           3  batch_glu Bacteria     Proteobacteria Alphaproteobacteria
## 161           3 chemo_carb Bacteria     Proteobacteria Alphaproteobacteria
## 162           3 batch_carb Bacteria     Proteobacteria Alphaproteobacteria
## 163           2 batch_carb Bacteria     Actinobacteria      Acidimicrobiia
## 164           2  batch_glu Bacteria     Actinobacteria      Acidimicrobiia
## 165           4 chemo_carb Bacteria     Proteobacteria Gammaproteobacteria
## 166           1  batch_glu Bacteria     Proteobacteria Alphaproteobacteria
## 167           1 batch_carb Bacteria     Proteobacteria Alphaproteobacteria
## 168           2  batch_glu Bacteria     Proteobacteria Gammaproteobacteria
## 169           3 batch_carb Bacteria      Bacteroidetes         Bacteroidia
## 170           3  batch_glu Bacteria      Bacteroidetes         Bacteroidia
## 171           2  batch_glu Bacteria     Proteobacteria Alphaproteobacteria
## 172           2 batch_carb Bacteria     Proteobacteria Alphaproteobacteria
## 173           4  batch_glu Bacteria     Proteobacteria Gammaproteobacteria
## 174           3 batch_carb Bacteria     Proteobacteria Alphaproteobacteria
## 175           3  batch_glu Bacteria     Proteobacteria Alphaproteobacteria
## 176           2 chemo_carb Bacteria     Proteobacteria Gammaproteobacteria
## 177           3 chemo_carb Bacteria     Proteobacteria Gammaproteobacteria
## 178           1  batch_glu Bacteria    Verrucomicrobia    Verrucomicrobiae
## 179           3  batch_glu Bacteria     Planctomycetes    Planctomycetacia
## 180           3 batch_carb Bacteria     Planctomycetes    Planctomycetacia
## 181           3 batch_carb Bacteria     Proteobacteria Alphaproteobacteria
## 182           3  batch_glu Bacteria     Proteobacteria Alphaproteobacteria
## 183           1 batch_carb Bacteria    Verrucomicrobia    Verrucomicrobiae
## 184           1  batch_glu Bacteria    Verrucomicrobia    Verrucomicrobiae
## 185           2 batch_carb Bacteria     Proteobacteria Gammaproteobacteria
## 186           2  batch_glu Bacteria     Proteobacteria Gammaproteobacteria
## 187           1 chemo_carb Bacteria     Proteobacteria Alphaproteobacteria
## 188           4 chemo_carb Bacteria     Proteobacteria Alphaproteobacteria
## 189           3  batch_glu Bacteria      Bacteroidetes         Bacteroidia
## 190           2 chemo_carb Bacteria     Proteobacteria Deltaproteobacteria
## 191           4 chemo_carb Bacteria     Proteobacteria Alphaproteobacteria
## 192           3  batch_glu Bacteria      Bacteroidetes         Bacteroidia
## 193           3 batch_carb Bacteria      Bacteroidetes         Bacteroidia
## 194           6 batch_carb Bacteria     Proteobacteria Gammaproteobacteria
## 195           6  batch_glu Bacteria     Proteobacteria Gammaproteobacteria
## 196           2 chemo_carb Bacteria     Proteobacteria Gammaproteobacteria
## 197           2 chemo_carb Bacteria     Actinobacteria      Actinobacteria
## 198           3 batch_carb Bacteria      Bacteroidetes         Bacteroidia
## 199           3  batch_glu Bacteria      Bacteroidetes         Bacteroidia
## 200           5  batch_glu Bacteria Epsilonbacteraeota     Campylobacteria
## 201           3  batch_glu Bacteria      Bacteroidetes         Bacteroidia
## 202           3 batch_carb Bacteria      Bacteroidetes         Bacteroidia
## 203           2 batch_carb Bacteria      Bacteroidetes         Bacteroidia
## 204           2  batch_glu Bacteria      Bacteroidetes         Bacteroidia
## 205           9 chemo_carb Bacteria     Proteobacteria Alphaproteobacteria
## 206           2  batch_glu Bacteria     Proteobacteria Gammaproteobacteria
## 207           2 chemo_carb Bacteria     Proteobacteria Alphaproteobacteria
## 208           2 batch_carb Bacteria      Bacteroidetes         Bacteroidia
## 209           2  batch_glu Bacteria      Bacteroidetes         Bacteroidia
## 210           1 chemo_carb Bacteria     Proteobacteria Alphaproteobacteria
## 211           5  batch_glu Bacteria      Bacteroidetes         Bacteroidia
## 212           1 chemo_carb Bacteria     Proteobacteria Gammaproteobacteria
## 213           2 chemo_carb Bacteria     Proteobacteria Gammaproteobacteria
## 214           2 batch_carb Bacteria     Proteobacteria Gammaproteobacteria
## 215           2  batch_glu Bacteria     Proteobacteria Gammaproteobacteria
## 216           7 batch_carb Bacteria      Bacteroidetes         Bacteroidia
## 217           7  batch_glu Bacteria      Bacteroidetes         Bacteroidia
## 218           1 chemo_carb Bacteria     Proteobacteria Alphaproteobacteria
## 219           3 chemo_carb Bacteria      Bacteroidetes         Bacteroidia
## 220           2 chemo_carb Bacteria     Proteobacteria Deltaproteobacteria
## 221           5 batch_carb Bacteria      Bacteroidetes         Bacteroidia
## 222           5  batch_glu Bacteria      Bacteroidetes         Bacteroidia
## 223           1 chemo_carb Bacteria     Proteobacteria Deltaproteobacteria
## 224           1 batch_carb Bacteria     Proteobacteria Alphaproteobacteria
## 225           1  batch_glu Bacteria     Proteobacteria Alphaproteobacteria
## 226           3 batch_carb Bacteria      Bacteroidetes         Bacteroidia
## 227           3  batch_glu Bacteria      Bacteroidetes         Bacteroidia
## 228           8  batch_glu Bacteria     Proteobacteria Gammaproteobacteria
## 229           8 chemo_carb Bacteria     Proteobacteria Gammaproteobacteria
## 230           2 chemo_carb Bacteria      Bacteroidetes         Bacteroidia
## 231           1  batch_glu Bacteria     Proteobacteria Alphaproteobacteria
## 232           3  batch_glu Bacteria      Bacteroidetes         Bacteroidia
## 233           3 batch_carb Bacteria      Bacteroidetes         Bacteroidia
## 234           6  batch_glu Bacteria     Proteobacteria Gammaproteobacteria
## 235           1  batch_glu Bacteria     Actinobacteria     Thermoleophilia
## 236           2  batch_glu Bacteria      Bacteroidetes         Bacteroidia
## 237           2 chemo_carb Bacteria      Bacteroidetes         Bacteroidia
## 238           2 batch_carb Bacteria      Bacteroidetes         Bacteroidia
##                         Order                     Family
## 1           Bdellovibrionales         Bdellovibrionaceae
## 2                 Rhizobiales          Xanthobacteraceae
## 3              Francisellales            Francisellaceae
## 4       Betaproteobacteriales             Rhodocyclaceae
## 5                                                       
## 6               Reyranellales             Reyranellaceae
## 7            Flavobacteriales             Cryomorphaceae
## 8            Flavobacteriales             Cryomorphaceae
## 9                       PeM15       uncultured bacterium
## 10                      PeM15       uncultured bacterium
## 11                 Opitutales                Opitutaceae
## 12           Rhodospirillales        Magnetospirillaceae
## 13           Flavobacteriales          Crocinitomicaceae
## 14           Flavobacteriales          Crocinitomicaceae
## 15            Chitinophagales           Chitinophagaceae
## 16      Betaproteobacteriales           Burkholderiaceae
## 17                Rhizobiales          Xanthobacteraceae
## 18               Pirellulales              Pirellulaceae
## 19               Pirellulales              Pirellulaceae
## 20               Pirellulales              Pirellulaceae
## 21               Pirellulales              Pirellulaceae
## 22                Rhizobiales                Devosiaceae
## 23      Betaproteobacteriales           Methylophilaceae
## 24         Verrucomicrobiales        Verrucomicrobiaceae
## 25            Chitinophagales             Saprospiraceae
## 26            Chitinophagales             Saprospiraceae
## 27         Diplorickettsiales        Diplorickettsiaceae
## 28                                                      
## 29      Betaproteobacteriales           Burkholderiaceae
## 30             Chloroflexales             Roseiflexaceae
## 31            Chitinophagales           Chitinophagaceae
## 32            Chitinophagales           Chitinophagaceae
## 33         Verrucomicrobiales        Verrucomicrobiaceae
## 34                      PeM15                           
## 35      Betaproteobacteriales                           
## 36                 Opitutales                Opitutaceae
## 37           Rhodospirillales           Magnetospiraceae
## 38           Flavobacteriales             Cryomorphaceae
## 39           Flavobacteriales             Cryomorphaceae
## 40            Caulobacterales           Caulobacteraceae
## 41            Caulobacterales           Caulobacteraceae
## 42               Pirellulales              Pirellulaceae
## 43            Caulobacterales           Caulobacteraceae
## 44               Cytophagales          Cyclobacteriaceae
## 45  uncultured soil bacterium  uncultured soil bacterium
## 46             Microtrichales                 uncultured
## 47      Betaproteobacteriales           Burkholderiaceae
## 48      Betaproteobacteriales           Burkholderiaceae
## 49            Rhodobacterales           Rhodobacteraceae
## 50            Rhodobacterales           Rhodobacteraceae
## 51                Rhizobiales Rhizobiales Incertae Sedis
## 52          Oceanospirillales       Saccharospirillaceae
## 53               Cytophagales              Spirosomaceae
## 54                Rhizobiales          Xanthobacteraceae
## 55             Rhodothermales            Rhodothermaceae
## 56      Betaproteobacteriales         Chitinibacteraceae
## 57           Flavobacteriales          Crocinitomicaceae
## 58           Flavobacteriales          Crocinitomicaceae
## 59         Sphingobacteriales       NS11-12 marine group
## 60               Cytophagales          Cyclobacteriaceae
## 61               Cytophagales          Cyclobacteriaceae
## 62           Flavobacteriales             Cryomorphaceae
## 63           Flavobacteriales             Cryomorphaceae
## 64                      PeM15       uncultured bacterium
## 65                                                      
## 66           Flavobacteriales          Flavobacteriaceae
## 67           Flavobacteriales          Flavobacteriaceae
## 68             Chloroflexales             Roseiflexaceae
## 69               Cytophagales              Spirosomaceae
## 70          Bdellovibrionales         Bdellovibrionaceae
## 71                 Opitutales                Opitutaceae
## 72               Pirellulales              Pirellulaceae
## 73                 Frankiales            Sporichthyaceae
## 74                 Frankiales            Sporichthyaceae
## 75         Sphingobacteriales                           
## 76                                                      
## 77                SAR11 clade             Ambiguous_taxa
## 78                SAR11 clade             Ambiguous_taxa
## 79           Salinisphaerales             Solimonadaceae
## 80                Rhizobiales          Xanthobacteraceae
## 81           Flavobacteriales          Flavobacteriaceae
## 82           Flavobacteriales          Flavobacteriaceae
## 83           Planctomycetales          Rubinisphaeraceae
## 84                 Frankiales            Sporichthyaceae
## 85                 Frankiales            Sporichthyaceae
## 86            Cellvibrionales           Cellvibrionaceae
## 87            Cellvibrionales           Cellvibrionaceae
## 88           Rhodospirillales           Magnetospiraceae
## 89           Flavobacteriales          Flavobacteriaceae
## 90           Flavobacteriales          Flavobacteriaceae
## 91             Sneathiellales            Sneathiellaceae
## 92              Rickettsiales             Rickettsiaceae
## 93                 Opitutales                Opitutaceae
## 94                 Opitutales                Opitutaceae
## 95             Microtrichales         Ilumatobacteraceae
## 96             Microtrichales         Ilumatobacteraceae
## 97            Cellvibrionales           Cellvibrionaceae
## 98            Cellvibrionales           Cellvibrionaceae
## 99                 Opitutales                Opitutaceae
## 100                Opitutales                Opitutaceae
## 101         Bdellovibrionales         Bdellovibrionaceae
## 102            Azospirillales            Azospirillaceae
## 103            Azospirillales            Azospirillaceae
## 104               OM182 clade             Ambiguous_taxa
## 105               OM182 clade             Ambiguous_taxa
## 106           Pseudomonadales           Pseudomonadaceae
## 107           Pseudomonadales           Pseudomonadaceae
## 108        Sphingobacteriales       NS11-12 marine group
## 109     Betaproteobacteriales             Rhodocyclaceae
## 110          Rhodospirillales           Terasakiellaceae
## 111           Chitinophagales           Chitinophagaceae
## 112     Betaproteobacteriales             Rhodocyclaceae
## 113               SAR11 clade                  Clade III
## 114               SAR11 clade                  Clade III
## 115               SAR11 clade                  Clade III
## 116     Betaproteobacteriales           Burkholderiaceae
## 117     Betaproteobacteriales           Burkholderiaceae
## 118     Betaproteobacteriales             Rhodocyclaceae
## 119     Betaproteobacteriales           Burkholderiaceae
## 120     Betaproteobacteriales           Burkholderiaceae
## 121               Balneolales               Balneolaceae
## 122     Betaproteobacteriales           Burkholderiaceae
## 123     Betaproteobacteriales           Burkholderiaceae
## 124               Rhizobiales               Rhizobiaceae
## 125        Verrucomicrobiales        Verrucomicrobiaceae
## 126     Betaproteobacteriales           Burkholderiaceae
## 127       Paracaedibacterales       Paracaedibacteraceae
## 128       Paracaedibacterales       Paracaedibacteraceae
## 129           Rhodobacterales           Rhodobacteraceae
## 130           Rhodobacterales           Rhodobacteraceae
## 131              Pirellulales              Pirellulaceae
## 132              Pirellulales              Pirellulaceae
## 133          Rhodospirillales        Magnetospirillaceae
## 134             Reyranellales             Reyranellaceae
## 135               Rhizobiales               Rhizobiaceae
## 136     Betaproteobacteriales           Methylophilaceae
## 137     Betaproteobacteriales           Methylophilaceae
## 138          Flavobacteriales          Flavobacteriaceae
## 139           Pseudomonadales           Pseudomonadaceae
## 140           Pseudomonadales           Pseudomonadaceae
## 141     Betaproteobacteriales           Burkholderiaceae
## 142          Flavobacteriales          Flavobacteriaceae
## 143          Flavobacteriales          Flavobacteriaceae
## 144                Frankiales            Sporichthyaceae
## 145                Frankiales            Sporichthyaceae
## 146                                                     
## 147              Pirellulales              Pirellulaceae
## 148                Frankiales            Sporichthyaceae
## 149         Bdellovibrionales         Bdellovibrionaceae
## 150           Caulobacterales           Caulobacteraceae
## 151              Myxococcales             Nannocystaceae
## 152     Betaproteobacteriales           Burkholderiaceae
## 153           Pseudomonadales              Moraxellaceae
## 154         Corynebacteriales           Mycobacteriaceae
## 155     Betaproteobacteriales             Rhodocyclaceae
## 156        Verrucomicrobiales        Verrucomicrobiaceae
## 157     Betaproteobacteriales           Burkholderiaceae
## 158     Betaproteobacteriales           Burkholderiaceae
## 159     Betaproteobacteriales           Burkholderiaceae
## 160           Rhodobacterales           Rhodobacteraceae
## 161           Rhodobacterales           Rhodobacteraceae
## 162           Rhodobacterales           Rhodobacteraceae
## 163            Microtrichales         Ilumatobacteraceae
## 164            Microtrichales         Ilumatobacteraceae
## 165             Legionellales             Legionellaceae
## 166         Thalassobaculales         Thalassobaculaceae
## 167         Thalassobaculales         Thalassobaculaceae
## 168     Betaproteobacteriales           Methylophilaceae
## 169          Flavobacteriales           NS9 marine group
## 170          Flavobacteriales           NS9 marine group
## 171               Rhizobiales                Devosiaceae
## 172               Rhizobiales                Devosiaceae
## 173     Betaproteobacteriales           Burkholderiaceae
## 174           Rhodobacterales           Rhodobacteraceae
## 175           Rhodobacterales           Rhodobacteraceae
## 176     Betaproteobacteriales           Methylophilaceae
## 177     Betaproteobacteriales           Burkholderiaceae
## 178            Pedosphaerales            Pedosphaeraceae
## 179             Isosphaerales             Isosphaeraceae
## 180             Isosphaerales             Isosphaeraceae
## 181           Rhodobacterales           Rhodobacteraceae
## 182           Rhodobacterales           Rhodobacteraceae
## 183            Pedosphaerales            Pedosphaeraceae
## 184            Pedosphaerales            Pedosphaeraceae
## 185     Betaproteobacteriales           Burkholderiaceae
## 186     Betaproteobacteriales           Burkholderiaceae
## 187           Caulobacterales            Hyphomonadaceae
## 188               Rhizobiales               Rhizobiaceae
## 189          Flavobacteriales          Crocinitomicaceae
## 190         Bdellovibrionales         Bacteriovoracaceae
## 191               Rhizobiales               Rhizobiaceae
## 192        Sphingobacteriales       NS11-12 marine group
## 193        Sphingobacteriales       NS11-12 marine group
## 194           Alteromonadales           Alteromonadaceae
## 195           Alteromonadales           Alteromonadaceae
## 196     Betaproteobacteriales           Burkholderiaceae
## 197             Micrococcales          Microbacteriaceae
## 198          Flavobacteriales             Cryomorphaceae
## 199          Flavobacteriales             Cryomorphaceae
## 200         Campylobacterales            Arcobacteraceae
## 201          Flavobacteriales                           
## 202          Flavobacteriales                           
## 203           Chitinophagales             Saprospiraceae
## 204           Chitinophagales             Saprospiraceae
## 205                Elsterales                Elsteraceae
## 206           Methylococcales           Methylococcaceae
## 207          Rhodospirillales        Magnetospirillaceae
## 208           Chitinophagales             Saprospiraceae
## 209           Chitinophagales             Saprospiraceae
## 210                                                     
## 211          Flavobacteriales          Flavobacteriaceae
## 212           Xanthomonadales           Xanthomonadaceae
## 213     Betaproteobacteriales           Methylophilaceae
## 214     Betaproteobacteriales           Methylophilaceae
## 215     Betaproteobacteriales           Methylophilaceae
## 216          Flavobacteriales          Flavobacteriaceae
## 217          Flavobacteriales          Flavobacteriaceae
## 218          Rhodospirillales                 uncultured
## 219              Cytophagales              Cytophagaceae
## 220         Bdellovibrionales         Bdellovibrionaceae
## 221          Flavobacteriales          Flavobacteriaceae
## 222          Flavobacteriales          Flavobacteriaceae
## 223         Bdellovibrionales         Bdellovibrionaceae
## 224         Thalassobaculales         Thalassobaculaceae
## 225         Thalassobaculales         Thalassobaculaceae
## 226              Cytophagales              Spirosomaceae
## 227              Cytophagales              Spirosomaceae
## 228     Betaproteobacteriales           Burkholderiaceae
## 229     Betaproteobacteriales           Burkholderiaceae
## 230              Cytophagales            Microscillaceae
## 231               SAR11 clade                  Clade III
## 232              Cytophagales          Cyclobacteriaceae
## 233              Cytophagales          Cyclobacteriaceae
## 234           Alteromonadales           Alteromonadaceae
## 235                Gaiellales                 uncultured
## 236        Sphingobacteriales                 env.OPS 17
## 237        Sphingobacteriales                 env.OPS 17
## 238        Sphingobacteriales                 env.OPS 17
##                                                  Genus
## 1                                         Bdellovibrio
## 2                                                     
## 3                   [Caedibacter] taeniospiralis group
## 4                                             Zoogloea
## 5                                                     
## 6                                           Reyranella
## 7                                           uncultured
## 8                                           uncultured
## 9                                 uncultured bacterium
## 10                                uncultured bacterium
## 11                                      Ambiguous_taxa
## 12                                    Magnetospirillum
## 13                                          Fluviicola
## 14                                          Fluviicola
## 15                                   Sediminibacterium
## 16                                                    
## 17                                        Ancylobacter
## 18                                          uncultured
## 19                                          uncultured
## 20                                          uncultured
## 21                                          uncultured
## 22                                             Devosia
## 23                                       Methylophilus
## 24                                     Prosthecobacter
## 25                               Candidatus Aquirestis
## 26                               Candidatus Aquirestis
## 27                                           Aquicella
## 28                                                    
## 29                                                    
## 30                                          uncultured
## 31                                   Sediminibacterium
## 32                                   Sediminibacterium
## 33                                     Prosthecobacter
## 34                                                    
## 35                                                    
## 36                                       Lacunisphaera
## 37                                          uncultured
## 38                                          uncultured
## 39                                          uncultured
## 40                                         Caulobacter
## 41                                         Caulobacter
## 42                                          uncultured
## 43                                                    
## 44                                                    
## 45                           uncultured soil bacterium
## 46                                                    
## 47                                                    
## 48                                                    
## 49                                                    
## 50                                                    
## 51                                       Phreatobacter
## 52                                                    
## 53                                         Dyadobacter
## 54                                  Pseudoxanthobacter
## 55                                          uncultured
## 56                                             Deefgea
## 57                                          Fluviicola
## 58                                          Fluviicola
## 59                                                    
## 60                                          uncultured
## 61                                          uncultured
## 62                                          uncultured
## 63                                          uncultured
## 64                                uncultured bacterium
## 65                                                    
## 66                                                    
## 67                                                    
## 68                                          uncultured
## 69                                         Dyadobacter
## 70                                          OM27 clade
## 71                                        Diplosphaera
## 72                                          uncultured
## 73                                                    
## 74                                                    
## 75                                                    
## 76                                                    
## 77                                                    
## 78                                                    
## 79                                             Nevskia
## 80                                        Xanthobacter
## 81                                      Flavobacterium
## 82                                      Flavobacterium
## 83                                          uncultured
## 84                                      Ambiguous_taxa
## 85                                      Ambiguous_taxa
## 86                                          Cellvibrio
## 87                                          Cellvibrio
## 88                                          uncultured
## 89                                      Flavobacterium
## 90                                      Flavobacterium
## 91                                         Ferrovibrio
## 92                                          uncultured
## 93                                     Cephaloticoccus
## 94                                     Cephaloticoccus
## 95                               CL500-29 marine group
## 96                               CL500-29 marine group
## 97                                          Cellvibrio
## 98                                          Cellvibrio
## 99                                        Diplosphaera
## 100                                       Diplosphaera
## 101                                       Bdellovibrio
## 102                                       Azospirillum
## 103                                       Azospirillum
## 104                                     Ambiguous_taxa
## 105                                     Ambiguous_taxa
## 106                                        Pseudomonas
## 107                                        Pseudomonas
## 108                                                   
## 109                                  Methyloversatilis
## 110                                         uncultured
## 111                                         uncultured
## 112                                           Zoogloea
## 113                                                   
## 114                                                   
## 115                                                   
## 116                                                   
## 117                                                   
## 118                                                   
## 119                                  RS62 marine group
## 120                                  RS62 marine group
## 121                               uncultured bacterium
## 122                                                   
## 123                                                   
## 124                                                   
## 125                                    Prosthecobacter
## 126                            MWH-UniP1 aquatic group
## 127                                         uncultured
## 128                                         uncultured
## 129                                                   
## 130                                                   
## 131                                          Pirellula
## 132                                          Pirellula
## 133                                   Magnetospirillum
## 134                                         Reyranella
## 135 Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium
## 136                                         OM43 clade
## 137                                         OM43 clade
## 138                                     Flavobacterium
## 139                                        Pseudomonas
## 140                                        Pseudomonas
## 141                                                   
## 142                                  NS3a marine group
## 143                                  NS3a marine group
## 144                                         hgcI clade
## 145                                         hgcI clade
## 146                                                   
## 147                                         uncultured
## 148                                                   
## 149                                       Bdellovibrio
## 150                                        Caulobacter
## 151                                        Nannocystis
## 152                                                   
## 153                                       Alkanindiges
## 154                                      Mycobacterium
## 155                                                   
## 156                                    Prosthecobacter
## 157                                                   
## 158                                                   
## 159                                                   
## 160                                                   
## 161                                                   
## 162                                                   
## 163                              CL500-29 marine group
## 164                              CL500-29 marine group
## 165                                         Legionella
## 166                                    Thalassobaculum
## 167                                    Thalassobaculum
## 168                                     Ambiguous_taxa
## 169                                                   
## 170                                                   
## 171                                                   
## 172                                                   
## 173                            MWH-UniP1 aquatic group
## 174                                                   
## 175                                                   
## 176                                       Methylovorus
## 177                                                   
## 178                                                   
## 179                                         uncultured
## 180                                         uncultured
## 181                                                   
## 182                                                   
## 183                                                   
## 184                                                   
## 185                                                   
## 186                                                   
## 187                                         Hyphomonas
## 188 Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium
## 189                                       Salinirepens
## 190                                      Bacteriovorax
## 191 Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium
## 192                                                   
## 193                                                   
## 194                                       Rheinheimera
## 195                                       Rheinheimera
## 196                                                   
## 197                                                   
## 198                                         uncultured
## 199                                         uncultured
## 200                                         Arcobacter
## 201                                                   
## 202                                                   
## 203                                         uncultured
## 204                                         uncultured
## 205                                            Elstera
## 206                                  Methyloparacoccus
## 207                                   Magnetospirillum
## 208                              Candidatus Aquirestis
## 209                              Candidatus Aquirestis
## 210                                                   
## 211                                     Flavobacterium
## 212                                                   
## 213                                      Methylotenera
## 214                                      Methylotenera
## 215                                      Methylotenera
## 216                                     Flavobacterium
## 217                                     Flavobacterium
## 218                                         metagenome
## 219                                          Cytophaga
## 220                                       Bdellovibrio
## 221                                     Flavobacterium
## 222                                     Flavobacterium
## 223                                       Bdellovibrio
## 224                                    Thalassobaculum
## 225                                    Thalassobaculum
## 226                                     Flectobacillus
## 227                                     Flectobacillus
## 228                                  Janthinobacterium
## 229                                  Janthinobacterium
## 230                                              OLB12
## 231                                                   
## 232                                       Algoriphagus
## 233                                       Algoriphagus
## 234                                       Rheinheimera
## 235                                                   
## 236                                                   
## 237                                                   
## 238                                                   
##                                                    Species
## 1                                               metagenome
## 2                                                         
## 3                                                         
## 4                                                         
## 5                                                         
## 6                                                         
## 7                                                         
## 8                                                         
## 9                                     uncultured bacterium
## 10                                    uncultured bacterium
## 11                                          Ambiguous_taxa
## 12                                                        
## 13                                                        
## 14                                                        
## 15                                                        
## 16                                                        
## 17                                                        
## 18                                                        
## 19                                                        
## 20                                                        
## 21                                                        
## 22                                                        
## 23                                                        
## 24                                                        
## 25                                                        
## 26                                                        
## 27                                                        
## 28                                                        
## 29                                                        
## 30                              uncultured Roseiflexus sp.
## 31                                                        
## 32                                                        
## 33                                                        
## 34                                                        
## 35                                                        
## 36                                                        
## 37                                                        
## 38       uncultured Bacteroidetes/Chlorobi group bacterium
## 39       uncultured Bacteroidetes/Chlorobi group bacterium
## 40                                                        
## 41                                                        
## 42                                                        
## 43                                                        
## 44                                                        
## 45                               uncultured soil bacterium
## 46                                                        
## 47                                                        
## 48                                                        
## 49                                                        
## 50                                                        
## 51                                          Ambiguous_taxa
## 52                                                        
## 53                                                        
## 54                                          Ambiguous_taxa
## 55                                uncultured Cytophaga sp.
## 56                                          Ambiguous_taxa
## 57                                                        
## 58                                                        
## 59                                                        
## 60                                                        
## 61                                                        
## 62                                                        
## 63                                                        
## 64                                    uncultured bacterium
## 65                                                        
## 66                                                        
## 67                                                        
## 68                              uncultured Roseiflexus sp.
## 69                                                        
## 70                                                        
## 71                            uncultured Desulfococcus sp.
## 72                                                        
## 73                                                        
## 74                                                        
## 75                                                        
## 76                                                        
## 77                                                        
## 78                                                        
## 79                                                        
## 80                                                        
## 81                                                        
## 82                                                        
## 83                                                        
## 84                                                        
## 85                                                        
## 86                                                        
## 87                                                        
## 88                                                        
## 89                                                        
## 90                                                        
## 91                                          Ambiguous_taxa
## 92                                                        
## 93                                    uncultured bacterium
## 94                                    uncultured bacterium
## 95                                                        
## 96                                                        
## 97                                                        
## 98                                                        
## 99                                    uncultured bacterium
## 100                                   uncultured bacterium
## 101                                                       
## 102                                                       
## 103                                                       
## 104                                         Ambiguous_taxa
## 105                                         Ambiguous_taxa
## 106                                                       
## 107                                                       
## 108                                                       
## 109                                                       
## 110                                     Aestuariispira sp.
## 111                                             metagenome
## 112                                                       
## 113                                                       
## 114                                                       
## 115                                                       
## 116                                                       
## 117                                                       
## 118                                                       
## 119                                   uncultured bacterium
## 120                                   uncultured bacterium
## 121                                   uncultured bacterium
## 122                                                       
## 123                                                       
## 124                                                       
## 125                                                       
## 126                                                       
## 127                                   uncultured bacterium
## 128                                   uncultured bacterium
## 129                                                       
## 130                                                       
## 131                                    uncultured organism
## 132                                    uncultured organism
## 133                                                       
## 134                                                       
## 135                                                       
## 136                                                       
## 137                                                       
## 138                                                       
## 139                                                       
## 140                                                       
## 141                                                       
## 142                                      marine metagenome
## 143                                      marine metagenome
## 144                                   uncultured bacterium
## 145                                   uncultured bacterium
## 146                                                       
## 147                                                       
## 148                                                       
## 149                                                       
## 150                                                       
## 151                                        Nannocystis sp.
## 152                                                       
## 153                                                       
## 154                                                       
## 155                                                       
## 156                                                       
## 157                                                       
## 158                                                       
## 159                                                       
## 160                                                       
## 161                                                       
## 162                                                       
## 163                                                       
## 164                                                       
## 165                                                       
## 166                           alpha proteobacterium BAL199
## 167                           alpha proteobacterium BAL199
## 168                                                       
## 169                                                       
## 170                                                       
## 171                                                       
## 172                                                       
## 173                                                       
## 174                                                       
## 175                                                       
## 176                                                       
## 177                                                       
## 178                                                       
## 179                                                       
## 180                                                       
## 181                                                       
## 182                                                       
## 183                                                       
## 184                                                       
## 185                                                       
## 186                                                       
## 187                                                       
## 188                                                       
## 189                                         Ambiguous_taxa
## 190                                         Ambiguous_taxa
## 191                                                       
## 192                                                       
## 193                                                       
## 194                                                       
## 195                                                       
## 196                                                       
## 197                                                       
## 198                                                       
## 199                                                       
## 200                                                       
## 201                                                       
## 202                                                       
## 203                uncultured Sphingobacteriales bacterium
## 204                uncultured Sphingobacteriales bacterium
## 205                                         Ambiguous_taxa
## 206                                                       
## 207                                                       
## 208                                                       
## 209                                                       
## 210                                                       
## 211                                                       
## 212                                                       
## 213                                                       
## 214                                                       
## 215                                                       
## 216                                                       
## 217                                                       
## 218                                             metagenome
## 219                     uncultured Bacteroidetes bacterium
## 220                                   uncultured bacterium
## 221                                                       
## 222                                                       
## 223 Bdellovibrionales bacterium RIFCSPHIGHO2_01_FULL_40_29
## 224                           alpha proteobacterium BAL199
## 225                           alpha proteobacterium BAL199
## 226                                                       
## 227                                                       
## 228                                                       
## 229                                                       
## 230                                                       
## 231                                                       
## 232                               uncultured Hongiella sp.
## 233                               uncultured Hongiella sp.
## 234                                                       
## 235                                                       
## 236                                                       
## 237                                                       
## 238
```

``` r
ioc_data=read.csv("datafiles/gtdb_r220_ioc.csv",header=T)
ioc_data=read.csv("datafiles/gtdb_r220_ioc_25june26.csv",header=T)

qsip_data_ioc=merge(qsip_data, ioc_data, by='OTU')


qsip_data_ioc_aa = subset(qsip_data_ioc, Incorporator == "TRUE")
plot(qsip_data_ioc_aa$Wlight, qsip_data_ioc_aa$IoC_genus_mean)
```

![](/Users/oliviaahern/Documents/GitHub/Sip_CopyNumber/docs/index_files/figure-html/ioc_1-1.png)<!-- -->

``` r
boxplot(qsip_data_ioc$IoC_genus_mean~qsip_data_ioc$id)

par(mfrow=c(1,2))
boxplot(qsip_data_ioc$IoC_genus_mean~qsip_data_ioc$id)
a=aov(qsip_data_ioc$IoC_genus_mean~qsip_data_ioc$id)
summary(a)
```

```
##                   Df Sum Sq Mean Sq F value  Pr(>F)   
## qsip_data_ioc$id   2   7.84   3.919   6.357 0.00228 **
## Residuals        139  85.68   0.616                   
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

``` r
tuk=TukeyHSD(a)
plot(tuk,las=2)
```

![](/Users/oliviaahern/Documents/GitHub/Sip_CopyNumber/docs/index_files/figure-html/ioc_1-2.png)<!-- -->

``` r
tuk
```

```
##   Tukey multiple comparisons of means
##     95% family-wise confidence level
## 
## Fit: aov(formula = qsip_data_ioc$IoC_genus_mean ~ qsip_data_ioc$id)
## 
## $`qsip_data_ioc$id`
##                             diff        lwr       upr     p adj
## batch_glu-batch_carb  0.01668301 -0.3550374 0.3884034 0.9937869
## chemo_carb-batch_carb 0.54036449  0.1239280 0.9568010 0.0071315
## chemo_carb-batch_glu  0.52368149  0.1404752 0.9068878 0.0042772
```

``` r
ag <- aggregate(IoC_genus_mean ~ id, qsip_data_ioc_aa, function(x) c(mean = mean(x), sd = sd(x)))
ag
```

```
##           id IoC_genus_mean.mean IoC_genus_mean.sd
## 1 batch_carb           3.2413737         0.7666890
## 2  batch_glu           3.2403142         0.8225627
## 3 chemo_carb           3.7639957         0.7504043
```

``` r
ss=subset(incorp_only, id=='chemo_carb')
ss1=ss$OTU

incorp_w_ioc=readRDS('datafiles/incorp_w_ioc')

ss3=subset(incorp_w_ioc, id=="chemo_carb")
ss4=ss3$OTU

setdiff(ss1,ss4)
```

```
##  [1] "137dc661fa672b506fd805198510edf9" "1c2b8892a0d7f828d72808cd01cfd199"
##  [3] "1e104b993416995498346089f57fd3bf" "27af011f1ebcc0867bb168193e6a9d27"
##  [5] "34327ad682b105ca4607c44af24c5dd2" "63413adbe8227f80412e824dffac1f5a"
##  [7] "748289d833e671b4848c77b9443c2034" "801efb14d6457f799fd2becee6d07b72"
##  [9] "8b1b6b50ee8e9d52ac980948823b441f" "9c54c522a386d5ed75686001dae803d9"
## [11] "9e26dc7c6642e4bf724f8d376d56a973" "acc9603a819dc23a606ac6c1b9105b6d"
## [13] "d2af52fc6f48ab4d05c96dad9527d8a2" "d94c20c00db235bb6fccc68f06f100c7"
## [15] "e416d0916760d2fc17b616e2ac3ad855"
```


## figure with afe and control buoyant density 

### SI Figure 3 


``` r
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

``` r
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
  theme(legend.position='right') +
  scale_fill_manual(values=custom_colors_class) + 
  scale_color_manual(values=custom_colors_class) +
      scale_y_reverse() +
  xlab("EAF") +
  ylab(expression(" "*C^12*" Control Buoyant Density (g/mL)")) + 
  theme(axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    strip.text = element_text(size = 12),
    strip.background = element_blank(),
    panel.grid.minor = element_blank())
```

![](/Users/oliviaahern/Documents/GitHub/Sip_CopyNumber/docs/index_files/figure-html/si-fig-3-1.png)<!-- -->

``` r
# new
ggplot(incorp_only, aes(x=A, y=Wlight,  fill=Class_ag, color= Class_ag,size=copy_number, linewidth=0)) +
#  geom_smooth(method = "lm", se = TRUE, col='black',fullrange=TRUE,linewidth=1) + 
 # geom_point(size=as.numeric(both$copy_number), alpha=0.8) +
   geom_point(alpha=0.8,show.legend=TRUE)+
    theme_bw(base_size = 12) + 

  facet_wrap(Culture ~ Substrate) +
  theme_bw() +
  theme(legend.position='right') +
  scale_fill_manual(values=custom_colors_class) + 
  scale_color_manual(values=custom_colors_class) +
      scale_y_reverse() +
  xlab("EAF") +
  ylab(expression(" "*C^12*" Control Buoyant Density (g/mL)")) + 
  theme(axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    strip.text = element_text(size = 12),
    strip.background = element_blank(),
    panel.grid.minor = element_blank()) + scale_color_manual(values=custom_colors_class) + guides(color = guide_legend(override.aes = list(size = 4)))
```

![](/Users/oliviaahern/Documents/GitHub/Sip_CopyNumber/docs/index_files/figure-html/si-fig-3-2.png)<!-- -->





``` r
batch_incorp=readRDS('datafiles/batch_incorp')


multi=subset(batch_incorp, Substrate=="Multi")
single=subset(batch_incorp, Substrate=="Single")
both=cbind(single$A, multi$A, multi$Class, multi$copy_number)
colnames(both)=c("Single","Multi", "Class", "copy_number")
row.names(both)=multi$OTU
both=data.frame(both)

size_var <- abs(as.numeric(both$copy_number))
# Define labels for the legend
legend_labels <- c("1", "5", "9")
# Define the sizes to use in the legend (e.g., representative values)
legend_sizes <- c(0.5, 2, 4) 

# 2. Create the main plot
plot(subset(incorp_only, id=="chemo_carb")$Wlight, subset(incorp_only, id=="chemo_carb")$A, 
     pch = 16, # Use a solid point type
     cex = size_var, # Map the 'size_var' to the point size in the plot
     col = "blue",
     main = "Scatter plot with point size by value",
     xlab = "X values",
     ylab = "Y values")

# 3. Add the legend using pt.cex
legend("topleft", # Position of the legend
       legend = legend_labels, # Text labels
       pt.cex = legend_sizes, # Set the specific sizes for the legend points
       pt.bg = "white",
       pch=21,
       title = "16S GCN")
```

![](/Users/oliviaahern/Documents/GitHub/Sip_CopyNumber/docs/index_files/figure-html/legend-1c-1.png)<!-- -->

# 16S copy number vs EAF

## Figure 1


``` r
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
  ylab("16S GCN") + 
  theme(axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    strip.text = element_text(size = 12),
    strip.background = element_blank(),
    panel.grid.minor = element_blank()) +
    scale_color_manual(values=custom_colors_class) 
```

![](/Users/oliviaahern/Documents/GitHub/Sip_CopyNumber/docs/index_files/figure-html/copy-2ag2-1.png)<!-- -->

``` r
plot_16s=ggplot(incorp_only, aes(x=A, y=copy_number, color=Class_ag)) +
  #geom_point(size=as.numeric(incorp_only$copy_number)*0.7, alpha=0.8) +
  geom_point(size=3) +
  geom_smooth(method = "lm", se = TRUE, col='black',fullrange=TRUE)  +
  facet_wrap(Culture ~ Substrate) +
  theme_bw(base_size=12) +
  xlab("Excess Atomic Fraction") +
  ylab("16S GCN") + 
  theme(axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    strip.text = element_text(size = 12),
    strip.background = element_blank(),
    panel.grid.minor = element_blank(), legend.position='none') +
    scale_color_manual(values=custom_colors_class) 
```



### stats generated from here

``` r
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

``` r
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

``` r
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

``` r
abline(ll)
```

![](/Users/oliviaahern/Documents/GitHub/Sip_CopyNumber/docs/index_files/figure-html/copy-3-1.png)<!-- -->



## Figure IoC


``` r
ioc_data=read.csv("datafiles/gtdb_r220_ioc.csv",header=T)
ioc_data=read.csv("datafiles/gtdb_r220_ioc_25june26.csv",header=T)

incorp1=read.csv('datafiles/incorporators.csv',header = T)

incorp_w_ioc=merge(incorp1, ioc_data, by='OTU')

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
incorp_w_ioc$Substrate <- factor(incorp_w_ioc$Substrate, levels = desired_order)

formula = incorp_w_ioc$IoC_genus_mean ~ incorp_w_ioc$A * incorp_w_ioc$Culture * incorp_w_ioc$Substrate

ggplot(incorp_w_ioc, aes(x=A, y=IoC_genus_mean, color=Class_ag)) +
  #geom_point(size=as.numeric(incorp_only$copy_number)*0.7, alpha=0.8) +
  geom_point(size=3) +
  geom_smooth(method = "lm", se = TRUE, col='black',fullrange=TRUE)  +
  facet_wrap(Culture ~ Substrate) +
  theme_bw(base_size=12) +
  xlab("Excess Atomic Fraction") +
  ylab("Index of Copiotrophy") + 
  theme(axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    strip.text = element_text(size = 12),
    strip.background = element_blank(),
    panel.grid.minor = element_blank()) +
    scale_color_manual(values=custom_colors_class) 
```

![](/Users/oliviaahern/Documents/GitHub/Sip_CopyNumber/docs/index_files/figure-html/copy-ico-1.png)<!-- -->

``` r
copio_plot=ggplot(incorp_w_ioc, aes(x=A, y=IoC_genus_mean, color=Class_ag)) +
  #geom_point(size=as.numeric(incorp_only$copy_number)*0.7, alpha=0.8) +
  geom_point(size=3) +
  geom_smooth(method = "lm", se = TRUE, col='black',fullrange=TRUE)  +
  facet_wrap(Culture ~ Substrate) +
  theme_bw(base_size=12) +
  xlab("Excess Atomic Fraction") +
  ylab("Index of Copiotrophy") + 
  theme(axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    strip.text = element_text(size = 12),
    strip.background = element_blank(),
    panel.grid.minor = element_blank(), legend.position = 'none') +
    scale_color_manual(values=custom_colors_class) 
```




``` r
library(ggh4x)
plot_16sa= plot_16s+ force_panelsizes(rows = unit(2.5, "in"), 
                     cols = unit(2.5, "in"))
copio_plota=copio_plot+ force_panelsizes(rows = unit(2.5, "in"), 
                     cols = unit(2.5, "in"))
quahog=cowplot::plot_grid(plot_16sa,copio_plota,
                   align='h',
                   axis='l',
                   byrow=TRUE,
                   nrow=2)
ggsave('figure2_withioc.svg', plot=quahog, dpi=300)


copio_plot3=ggplot(incorp_w_ioc, aes(x=A, y=IoC_genus_mean, color=Class_ag)) +
  #geom_point(size=as.numeric(incorp_only$copy_number)*0.7, alpha=0.8) +
  geom_point(size=3) +
  geom_smooth(method = "lm", se = TRUE, col='black',fullrange=TRUE)  +
  facet_wrap(Culture ~ Substrate) +
  theme_bw(base_size=12) +
  xlab("Excess Atomic Fraction") +
  ylab("Index of Copiotrophy") + 
  theme(axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    strip.text = element_text(size = 12),
    strip.background = element_blank(),
    panel.grid.minor = element_blank(), legend.position = 'bottom') +
    scale_color_manual(values=custom_colors_class) + force_panelsizes(rows = unit(2.5, "in"), 
                     cols = unit(2.5, "in"))
copio_plot3
```

![](/Users/oliviaahern/Documents/GitHub/Sip_CopyNumber/docs/index_files/figure-html/copy-ico2-1.png)<!-- -->

``` r
ggsave('copioplot_legend.svg', plot=copio_plot3, dpi=300)

cowplot::plot_grid(copio_plot3,
                   align='h',
                   axis='l',
                   byrow=TRUE,
                   nrow=2)
```

![](/Users/oliviaahern/Documents/GitHub/Sip_CopyNumber/docs/index_files/figure-html/copy-ico2-2.png)<!-- -->

### stats generated from here

``` r
par(mfrow=c(1,3))
glu_only=subset(incorp_w_ioc, id=="batch_glu")
set.seed(124)

plot(glu_only$A, glu_only$IoC_genus_mean, pch=21, cex=1.22, bg='black', xlab="AFE", ylab="IoC_genus_mean")
ll=lm(glu_only$IoC_genus_mean~glu_only$A)
summary(ll)
```

```
## 
## Call:
## lm(formula = glu_only$IoC_genus_mean ~ glu_only$A)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -1.0585 -0.4802 -0.1767  0.4235  2.1441 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   2.5694     0.1723   14.91  < 2e-16 ***
## glu_only$A    1.9421     0.4240    4.58  2.4e-05 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.7139 on 60 degrees of freedom
## Multiple R-squared:  0.2591,	Adjusted R-squared:  0.2467 
## F-statistic: 20.98 on 1 and 60 DF,  p-value: 2.396e-05
```

``` r
abline(ll)
text(0.2, 6, "R2 = 0.2656 
F = 23.07 
p = 1.083e-05")


batch_carbon=subset(incorp_w_ioc, id=="batch_carb")
plot(batch_carbon$A, batch_carbon$IoC_genus_mean, pch=21, cex=1.22, bg='black', xlab="AFE", ylab="IoC_genus_mean")
ll=lm(batch_carbon$IoC_genus_mean~batch_carbon$A)
summary(ll)
```

```
## 
## Call:
## lm(formula = batch_carbon$IoC_genus_mean ~ batch_carbon$A)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -0.9643 -0.4808 -0.2557  0.4172  1.9726 
## 
## Coefficients:
##                Estimate Std. Error t value Pr(>|t|)    
## (Intercept)      2.8135     0.2947   9.546 1.22e-11 ***
## batch_carbon$A   1.6201     1.0210   1.587    0.121    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.7522 on 38 degrees of freedom
## Multiple R-squared:  0.06215,	Adjusted R-squared:  0.03747 
## F-statistic: 2.518 on 1 and 38 DF,  p-value: 0.1208
```

``` r
abline(ll)
text(0.2, 5, "R2 = 0.077
F = 4.257
p = 0.045")

chemo_carbon=subset(incorp_w_ioc, id=="chemo_carb")
plot(chemo_carbon$A, chemo_carbon$IoC_genus_mean, pch=21, cex=1.22, bg='black', xlab="AFE", ylab="IoC_genus_mean")
ll=lm(chemo_carbon$IoC_genus_mean~chemo_carbon$A)
summary(ll)
```

```
## 
## Call:
## lm(formula = chemo_carbon$IoC_genus_mean ~ chemo_carbon$A)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -0.9761 -0.4559 -0.2254  0.4608  2.3729 
## 
## Coefficients:
##                Estimate Std. Error t value Pr(>|t|)    
## (Intercept)      3.9420     0.1894  20.818   <2e-16 ***
## chemo_carbon$A  -1.1890     0.9732  -1.222     0.23    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.7455 on 36 degrees of freedom
## Multiple R-squared:  0.03981,	Adjusted R-squared:  0.01314 
## F-statistic: 1.493 on 1 and 36 DF,  p-value: 0.2298
```

``` r
abline(ll)
```

![](/Users/oliviaahern/Documents/GitHub/Sip_CopyNumber/docs/index_files/figure-html/copy-5-1.png)<!-- -->



# Beta diversity metrics 
These are to see if there are differences in commuinty between the batch and chemostat because there was different innoculum and different sequencing centers 

## philr PCoA 

### pcoa of all ASVs

``` r
GP <- transform_sample_counts(phyo_frac, function(x) x+1)
phy_tree(GP) <- makeNodeLabel(phy_tree(GP), method="number", prefix='n')
name.balance(phy_tree(GP), tax_table(GP), 'n1')
```

```
## [1] "Strain_7789a99e360388313085539ee5765e33/Domain_Bacteria"
```

``` r
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

``` r
adonis2(gp.dist ~ as.factor(sample_data(GP)$Culture) +as.factor(sample_data(GP)$Substrate))
```

```
## Permutation test for adonis under reduced model
## Permutation: free
## Number of permutations: 999
## 
## adonis2(formula = gp.dist ~ as.factor(sample_data(GP)$Culture) + as.factor(sample_data(GP)$Substrate))
##          Df SumOfSqs      R2      F Pr(>F)    
## Model     4  1102.28 0.83837 57.055  0.001 ***
## Residual 44   212.52 0.16163                  
## Total    48  1314.80 1.00000                  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

### pcoa of just incorporators 

``` r
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

``` r
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

``` r
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
## as.factor(sample_data(GP)$Treatment)  1    25.98 0.00960   2.989  0.089 .  
## Residual                             47   408.52 0.15093                   
## Total                                49  2706.75 1.00000                   
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

# shared asvs among treatments 


## ASVs present in all three treatmnents

I am not keeping this in because we used different CsTFA gradients for the batch and chemostat


``` r
# summary(qsip_data$OTU)
# two ASVs that were in all treatments positively 
# ffd660f96ba6668b2d8d86fc4e150862
# 84ea3564fd34c5377c4305b340f0c16f

# af27577de5edb10c09f7fc19b9ca45ca true for batch but false for chemostat
# false for batch single e416d0916760d2fc17b616e2ac3ad855
```

### ASV 1

``` r
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

``` r
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

``` r
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

``` r
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

``` r
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

``` r
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

``` r
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

``` r
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


``` r
batch=subset(qsip_data, Culture=="Batch")
clam=(summary(batch$OTU))
clam=data.frame(clam)
clam2=row.names(subset(clam, clam==2))
length(clam2)
```

```
## [1] 63
```

``` r
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

``` r
# 57 ASVs that are batch both single and multi and incorporators 
# batch_incorp=subset(batch, OTU %in% clam2) for some reason this doesn't work in my markdown %in% so I have to save it all 
# saveRDS(batch_incorp, 'datafiles/batch_incorp')
batch_incorp=readRDS('datafiles/batch_incorp')
# batch_incorp=data.frame(batch_incorp)

# summary(as.factor(batch_incorp$Class))
# summary(as.factor(batch_incorp$Family))
```





``` r
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

``` r
abline(ll)
text(0.8, 0.4, "R2 = 0.30 \n F = 25.38 \n p = 5.4 e-6")
abline(a=0,b=1,lty=2)
```

![](/Users/oliviaahern/Documents/GitHub/Sip_CopyNumber/docs/index_files/figure-html/batch-b-1.png)<!-- -->
### Figure 3

### batch multi eaf vs single eaf 

``` r
multi=subset(batch_incorp, Substrate=="Multi")
single=subset(batch_incorp, Substrate=="Single")


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



ggplot(both, aes(x=as.numeric(Multi), y=as.numeric(Single),  color=Class, linewidth=as.numeric(copy_number))) +
  geom_point(size=as.numeric(both$copy_number), alpha=0.8) +
  theme_bw() +
  ylim(c(0,0.8)) + xlim(c(0,0.8))  +
  theme(axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
    geom_smooth(method = "lm", se = TRUE, col='black',fullrange=TRUE) + 
    geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed", linewidth = 1) +
  xlab("EAF multi-substrate batch") +
  ylab("EAF single-substrate batch") +
#  annotate("text", label="R2 = 0.30 \n F = 25.38 \n p = 5.4 e-6", x=0.6, y=0.4) +
 #   scale_fill_manual(values=custom_colors_class) + 
  scale_color_manual(values=custom_colors_class) 
```

![](/Users/oliviaahern/Documents/GitHub/Sip_CopyNumber/docs/index_files/figure-html/batch-bb-1.png)<!-- -->

``` r
# new
ggplot(both, aes(x=as.numeric(Multi), y=as.numeric(Single),  color=Class, size=as.numeric(copy_number), linewidth=0)) +
  geom_smooth(method = "lm", se = TRUE, col='black',fullrange=TRUE,linewidth=1) + 
 # geom_point(size=as.numeric(both$copy_number), alpha=0.8) +
   geom_point(alpha=0.8,show.legend=TRUE)+
  theme_bw(base_size = 12) + 
 # scale_size(range = c(1, 10)) +
  ylim(c(0,0.8)) + xlim(c(0,0.8))  +
  theme(axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 12)) +
 guides(color = guide_legend(override.aes = list(linetype = 0))) +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  xlab("EAF multi-substrate batch") +
  ylab("EAF single-substrate batch") +
#  annotate("text", label="R2 = 0.30 \n F = 25.38 \n p = 5.4 e-6", x=0.6, y=0.4) +
 #   scale_fill_manual(values=custom_colors_class) + 
  scale_color_manual(values=custom_colors_class) + guides(color = guide_legend(override.aes = list(size = 4)))
```

![](/Users/oliviaahern/Documents/GitHub/Sip_CopyNumber/docs/index_files/figure-html/batch-bb-2.png)<!-- -->

``` r
new2=ggplot(both, aes(x=as.numeric(Single), y=as.numeric(Multi),  color=Class, size=as.numeric(copy_number), linewidth=0)) +
  geom_smooth(method = "lm", se = TRUE, col='black',fullrange=TRUE,linewidth=1) + 
 # geom_point(size=as.numeric(both$copy_number), alpha=0.8) +
   geom_point(alpha=0.8,show.legend=TRUE)+
  theme_bw(base_size = 12) + 
 # scale_size(range = c(1, 10)) +
  ylim(c(0,0.8)) + xlim(c(0,0.8))  +
  theme(axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 12)) +
 guides(color = guide_legend(override.aes = list(linetype = 0))) +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  xlab("EAF single-substrate batch") +
  ylab("EAF multi-substrate batch") +
#  annotate("text", label="R2 = 0.30 \n F = 25.38 \n p = 5.4 e-6", x=0.6, y=0.4) +
 #   scale_fill_manual(values=custom_colors_class) + 
  scale_color_manual(values=custom_colors_class) + guides(color = guide_legend(override.aes = list(size = 4)))
new2
```

![](/Users/oliviaahern/Documents/GitHub/Sip_CopyNumber/docs/index_files/figure-html/batch-bb-3.png)<!-- -->

``` r
ggsave(filename = 'figure3.svg', new2, dpi=300,height=3.7, width=6, units='in')
```



### batch multi eaf vs single eaf 

``` r
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



ggplot(both, aes(x=as.numeric(Multi), y=as.numeric(Single),  color=Class, linewidth=as.numeric(copy_number))) +
  geom_point(size=as.numeric(both$copy_number), alpha=0.8) +
  theme_bw() +
  ylim(c(0,0.8)) + xlim(c(0,0.8))  +
  theme(axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
    geom_smooth(method = "lm", se = TRUE, col='black',fullrange=TRUE) + 
    geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed", linewidth = 1) +
  xlab("EAF multi-substrate batch") +
  ylab("EAF single-substrate batch") +
#  annotate("text", label="R2 = 0.30 \n F = 25.38 \n p = 5.4 e-6", x=0.6, y=0.4) +
 #   scale_fill_manual(values=custom_colors_class) + 
  scale_color_manual(values=custom_colors_class) 
```

![](/Users/oliviaahern/Documents/GitHub/Sip_CopyNumber/docs/index_files/figure-html/batch-bbc-1.png)<!-- -->

``` r
# new
ggplot(both, aes(x=as.numeric(Multi), y=as.numeric(Single),  color=Class, size=as.numeric(copy_number), linewidth=0)) +
  geom_smooth(method = "lm", se = TRUE, col='black',fullrange=TRUE,linewidth=1) + 
 # geom_point(size=as.numeric(both$copy_number), alpha=0.8) +
   geom_point(alpha=0.8,show.legend=TRUE)+
  theme_bw(base_size = 12) + 
 # scale_size(range = c(1, 10)) +
  ylim(c(0,0.8)) + xlim(c(0,0.8))  +
  theme(axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 12)) +
 guides(color = guide_legend(override.aes = list(linetype = 0))) +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  xlab("EAF multi-substrate batch") +
  ylab("EAF single-substrate batch") +
#  annotate("text", label="R2 = 0.30 \n F = 25.38 \n p = 5.4 e-6", x=0.6, y=0.4) +
 #   scale_fill_manual(values=custom_colors_class) + 
  scale_color_manual(values=custom_colors_class) + guides(color = guide_legend(override.aes = list(size = 4)))
```

![](/Users/oliviaahern/Documents/GitHub/Sip_CopyNumber/docs/index_files/figure-html/batch-bbc-2.png)<!-- -->

``` r
new2=ggplot(both, aes(x=as.numeric(Multi), y=as.numeric(Single),  color=Class, size=as.numeric(copy_number), linewidth=0)) +
  geom_smooth(method = "lm", se = TRUE, col='black',fullrange=TRUE,linewidth=1) + 
 # geom_point(size=as.numeric(both$copy_number), alpha=0.8) +
   geom_point(alpha=0.8,show.legend=TRUE)+
  theme_bw(base_size = 12) + 
 # scale_size(range = c(1, 10)) +
  ylim(c(0,0.8)) + xlim(c(0,0.8))  +
  theme(axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 12)) +
 guides(color = guide_legend(override.aes = list(linetype = 0))) +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  xlab("EAF multi-substrate batch") +
  ylab("EAF single-substrate batch") +
#  annotate("text", label="R2 = 0.30 \n F = 25.38 \n p = 5.4 e-6", x=0.6, y=0.4) +
 #   scale_fill_manual(values=custom_colors_class) + 
  scale_color_manual(values=custom_colors_class) + guides(color = guide_legend(override.aes = list(size = 4)))
new2
```

![](/Users/oliviaahern/Documents/GitHub/Sip_CopyNumber/docs/index_files/figure-html/batch-bbc-3.png)<!-- -->

``` r
ggsave(filename = 'figure3.svg', new2, dpi=300,height=3.7, width=6, units='in')
```


##### try figure 3 with ioc data 


``` r
multi=subset(incorp_w_ioc, Substrate=="Multi")
single=subset(incorp_w_ioc, Substrate=="Single")

both_i=merge(multi, single, by="OTU")
clam1=subset(both_i, Substrate.x=="Multi")
clam_single=subset(both_i, Substrate.x=="Single")

plot(both_i$A.y, both_i$A.x, xlab= "EAF Batch Single Substrate", ylab ="EAF Batch Multi Substrate", pch=21, ylim=c(0,1), 
     xlim=c(0,1), bg='gray',cex=(both_i$IoC_genus_mean.x))
ll=lm(both_i$A.x~both_i$A.y)
summary(ll)
```

```
## 
## Call:
## lm(formula = both_i$A.x ~ both_i$A.y)
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -0.22970 -0.05500  0.00184  0.06645  0.47647 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  0.18700    0.04124   4.535 4.94e-05 ***
## both_i$A.y   0.23516    0.09594   2.451   0.0186 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.1244 on 41 degrees of freedom
## Multiple R-squared:  0.1278,	Adjusted R-squared:  0.1065 
## F-statistic: 6.008 on 1 and 41 DF,  p-value: 0.01859
```

``` r
abline(ll)
#text(0.2, 0.8, "R2 = 0.1065 \n F = 6.008 \n p = 0.01859")
abline(a=0,b=1,lty=2)
```

![](/Users/oliviaahern/Documents/GitHub/Sip_CopyNumber/docs/index_files/figure-html/batch-b1-1.png)<!-- -->





``` r
both2=data.frame(both)

plot(both2$Single, both2$Multi, ylim=c(0,0.8),xlim=c(0,0.8), pch=21, bg=as.factor(both2$Class),cex=2)
ll=lm(as.numeric(both2$Multi)~as.numeric(both2$Single))
abline(ll)
abline(a=0,b=1,lty=2)
```

![](/Users/oliviaahern/Documents/GitHub/Sip_CopyNumber/docs/index_files/figure-html/basic-line-1.png)<!-- -->

``` r
summary(ll)
```

```
## 
## Call:
## lm(formula = as.numeric(both2$Multi) ~ as.numeric(both2$Single))
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -0.21695 -0.03318  0.01027  0.06574  0.16447 
## 
## Coefficients:
##                          Estimate Std. Error t value Pr(>|t|)    
## (Intercept)               0.15899    0.02490   6.385 3.81e-08 ***
## as.numeric(both2$Single)  0.29796    0.05914   5.038 5.41e-06 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.09159 on 55 degrees of freedom
## Multiple R-squared:  0.3158,	Adjusted R-squared:  0.3033 
## F-statistic: 25.38 on 1 and 55 DF,  p-value: 5.407e-06
```


``` r
size_var <- abs(as.numeric(both$copy_number))
# Define labels for the legend
legend_labels <- c("1", "5", "9")
# Define the sizes to use in the legend (e.g., representative values)
legend_sizes <- c(0.5, 2, 4) 

par(mar=c(4,4,1,1))
# 2. Create the main plot
plot(both$Multi, both$Single, 
     pch = 16, # Use a solid point type
     cex = size_var, # Map the 'size_var' to the point size in the plot
     col = "blue",
     xlab = "X values",
     ylab = "Y values")

# 3. Add the legend using pt.cex
legend("topleft", # Position of the legend
       legend = legend_labels, # Text labels
       pt.cex = legend_sizes, # Set the specific sizes for the legend points
       pt.bg = "white",
       pch=21,
       title = "16S Gene Copy Number")
```

![](/Users/oliviaahern/Documents/GitHub/Sip_CopyNumber/docs/index_files/figure-html/legend-1-1.png)<!-- -->


``` r
size_var <- abs(as.numeric(both$copy_number))
# Define labels for the legend
legend_labels <- c("1", "5", "9")
# Define the sizes to use in the legend (e.g., representative values)
legend_sizes <- c(0.5, 2, 4) 

par(mar=c(4,4,1,1))
# 2. Create the main plot
plot(both$Multi, both$Single, 
     pch = 16, # Use a solid point type
     cex = size_var, # Map the 'size_var' to the point size in the plot
     col = "blue",
     xlab = "X values",
     ylab = "Y values")

# 3. Add the legend using pt.cex
legend("topleft", # Position of the legend
       legend = legend_labels, # Text labels
       pt.cex = legend_sizes, # Set the specific sizes for the legend points
       pt.bg = "white",
       pch=21,
       title = "16S Gene Copy Number")
```

![](/Users/oliviaahern/Documents/GitHub/Sip_CopyNumber/docs/index_files/figure-html/legend-2-1.png)<!-- -->


``` r
# greater than 3
chowder=subset(both, copy_number >= 3)
dim(both)
```

```
## [1] 57  5
```

``` r
dim(chowder)
```

```
## [1] 25  5
```

``` r
fact=c(chowder$Single, chowder$Multi)
factor2=c(rep("Single",25), rep("Multi",25))
a=aov(fact~factor2)
summary(a)
```

```
##             Df Sum Sq Mean Sq F value  Pr(>F)   
## factor2      1 0.3212  0.3212   12.13 0.00107 **
## Residuals   48 1.2709  0.0265                   
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

``` r
boxplot(as.numeric(fact)~factor(factor2))
```

![](/Users/oliviaahern/Documents/GitHub/Sip_CopyNumber/docs/index_files/figure-html/unnamed-chunk-4-1.png)<!-- -->

``` r
# greater than 3
chowder=subset(both, copy_number >= 4)
dim(both)
```

```
## [1] 57  5
```

``` r
dim(chowder)
```

```
## [1] 11  5
```

``` r
fact=c(chowder$Single, chowder$Multi)
factor2=c(rep("Single",11), rep("Multi",11))
a=aov(fact~factor2)
summary(a)
```

```
##             Df Sum Sq Mean Sq F value   Pr(>F)    
## factor2      1 0.4087  0.4087   29.36 2.65e-05 ***
## Residuals   20 0.2785  0.0139                     
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

``` r
boxplot(as.numeric(fact)~factor(factor2))
```

![](/Users/oliviaahern/Documents/GitHub/Sip_CopyNumber/docs/index_files/figure-html/unnamed-chunk-4-2.png)<!-- -->

``` r
a=t.test(as.numeric(fact)~factor2)
summary(a)
```

```
##             Length Class  Mode     
## statistic   1      -none- numeric  
## parameter   1      -none- numeric  
## p.value     1      -none- numeric  
## conf.int    2      -none- numeric  
## estimate    2      -none- numeric  
## null.value  1      -none- numeric  
## stderr      1      -none- numeric  
## alternative 1      -none- character
## method      1      -none- character
## data.name   1      -none- character
```

``` r
boxplot(as.numeric(fact)~factor(factor2))

chowder=subset(both, copy_number < 4)
dim(both)
```

```
## [1] 57  5
```

``` r
dim(chowder)
```

```
## [1] 46  5
```

``` r
fact=c(chowder$Single, chowder$Multi)
factor2=c(rep("Single",46), rep("Multi",46))
a=aov(fact~factor2)
summary(a)
```

```
##             Df Sum Sq Mean Sq F value Pr(>F)  
## factor2      1 0.0764 0.07640    3.15 0.0793 .
## Residuals   90 2.1830 0.02426                 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

``` r
boxplot(as.numeric(fact)~factor(factor2))
```

![](/Users/oliviaahern/Documents/GitHub/Sip_CopyNumber/docs/index_files/figure-html/unnamed-chunk-4-3.png)<!-- -->

``` r
a=t.test(as.numeric(fact)~factor2)
summary(a)
```

```
##             Length Class  Mode     
## statistic   1      -none- numeric  
## parameter   1      -none- numeric  
## p.value     1      -none- numeric  
## conf.int    2      -none- numeric  
## estimate    2      -none- numeric  
## null.value  1      -none- numeric  
## stderr      1      -none- numeric  
## alternative 1      -none- character
## method      1      -none- character
## data.name   1      -none- character
```

``` r
boxplot(as.numeric(fact)~factor(factor2))
```
#### heatmap 

 
