source("http://bioconductor.org/biocLite.R")
biocLite("DESeq2")

BiocManager::install(c("DESeq", "affycoretools", "arrayQualityMetrics", "genefilter", "ggplot2","dplyr","pheatmap", "vegan", "ggrepel", "tidyverse"))
BiocManager::install(c("arrayQualityMetrics"))

BiocManager::install("DESeq", version = "3.12")
devtools::install_version("DESeq", "1.38.0")
#set your working directory
#
setwd("/usr4/bi594/skoppara/Assignment2/Assign2_gr2") #you will need to change to your own directory

###conduct array quality metrics to detect and remove outliers
library(DESeq) 
library(affycoretools)
library(arrayQualityMetrics)
library(genefilter)

library("DESeq2")
library("DESeq")
packageVersion("affycoretools")
packageVersion("arrayQualityMetrics")
packageVersion("genefilter")
packageVersion("DESeq2")
packageVersion("ggplot2")
packageVersion("dplyr")
packageVersion("pheatmap")
packageVersion("vegan")
packageVersion("ggrepel")
packageVersion("tidyverse")

#https://kbroman.org/blog/2017/08/08/eof-within-quoted-string/
gg=read.table(file="davies_Ssid_iso2gene.tab",sep="\t", quote="")

#read in counts 
countData <- read.csv("2.Assign2_GE_Coral_OA.csv")
head(countData)
length(countData[,1])
#16931
##names(countData)=c( "Control_1", "Hot_1", "Control_2", "Hot_2")
## RUN ONCE TO SET ISOGROUP TO ROWNAME 
countData <- data.frame(countData, row.names = 1)
row.names(countData)=sub("", "isogroup", rownames(countData))
head(countData)


setwd("/usr4/bi594/skoppara/Assignment2/Assign2_gr2/outlier")
v=setwd("/usr4/bi594/skoppara/Assignment2/Assign2_gr2/outlier")

# # # #look for outliers
treat=c( "Control_1", "Hot_1", "Control_2", "Hot_2")
g=data.frame(treat)
g
colData = g

#some vars in design formula are characters converting to factors
dds=DESeqDataSetFromMatrix(countData=countData, colData=g, design=~treat)


vsd.ge=assay(vst(dds))
rl=vst(dds)
something = AnnotatedDataFrame(as.data.frame(colData(rl)))
a1 = assay(rl)
e=ExpressionSet(assay(rl), AnnotatedDataFrame(as.data.frame(colData(rl))))
arrayQualityMetrics(e,outdir=v,intgroup=c("treat"),force=T)
dev.off()
# double-click index.html

#####you only ever need to run the above code once. Outliers are decided at the beginning. 
## I like to close R and restart with packages etc
##So, please save your script and restart R
#no outliers detected!

setwd("/usr4/bi594/skoppara/Assignment2/Assign2_gr2")
library("DESeq2")
library("ggplot2")

#read in counts 
countData <- read.csv("2.Assign2_GE_Coral_OA.csv")

##run this only once to move col1 data to become row names
#countData <- data.frame(countData, row.names = 1)

head(countData)
length(countData[,1])
#16931

#names(countData)=sub(".fastq.trim.sam.counts","",names(countData))

names(countData)=c( "Control_1", "Hot_1", "Control_2", "Hot_2" )
#row.names(countData)=sub("", "isogroup", rownames(countData))
head(countData)

summedCounts=as.data.frame(colSums(countData))

barplot(totalCounts, col=c("#ff6666", "#6666ff", "#ffcc00", "#339900"), ylab="raw counts")

totalCounts = as.data.frame(treat)
totalCounts$raw = summedCounts[,1]
totalCounts$treat <- factor(totalCounts$treat, levels = totalCounts$treat)
p <- ggplot(totalCounts, aes(x = treat, y = raw))
p + geom_bar(aes(fill = treat), stat="identity") + xlab("Temperature") + ylab("Raw Counts")


totalCounts=colSums(countData)

# # pH7.5a  pH7.5b  pH7.5c  pH7.6a  pH7.6b  pH7.6c    pH8a    pH8b    pH8c 
# 789550  918366 1027861  926497  816054  770342  612258  651727  480153 
min(totalCounts) #4110385
max(totalCounts)  # 21708269

#colName=c( "Control_1", "Hot_1", "Control_2", "Hot_2")

temperature=c( "Control", "Hot", "Control", "Hot")
genotype=c( "a", "a", "b", "b")
g=data.frame(temperature, genotype)
g
colData<- g

dds<-DESeqDataSetFromMatrix(countData=countData, colData=data.frame(temperature), design=~temperature) #can only test for the main effects of site, pco2, temp

#one step DESeq
dds<-DESeq(dds)
# estimating size factors
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates
# fitting model and testing

head(dds)
res<- results(dds)

#Look at dispersions plot
plotDispEsts(dds, main="Dispersion plot temperature")



####################pH8 vs pH7.6 pairwise comparisons
# ORDER MATTERS treatment, control
colData$heat<-factor(colData$temperature, levels=c("Hot", "Control"))
##second term is the "control"
resheat <- results(dds, contrast=c("temperature","Hot","Control"))
#how many FDR < 10%?
resheat
table(resheat$padj<0.01)
# FALSE  TRUE 
# 14826  1660 

summary(resheat)
# out of 16812 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 379, 2.3%
# LFC < 0 (down)     : 1839, 11%
# outliers [1]       : 0, 0%
# low counts [2]     : 326, 1.9%

nrow(resheat[resheat$padj<0.05 & !is.na(resheat$padj),])  # Num significantly differentially expressed genes excluding the no/low count genes   #228
# 1987

dev.off()
plotMA(resheat, main="Hot vs Control")
plotMA(resheat, main="Hot vs Control", ylim=c(-2,2))

results <- as.data.frame(resheat)
head(results)

nrow(resheat[resheat$padj<0.1 & resheat$log2FoldChange > 0 & !is.na(resheat$padj),])
#379 upregulated 
nrow(resheat[resheat$padj<0.1 & resheat$log2FoldChange < 0 & !is.na(resheat$padj),])
#1839 downregulated 

write.table(resheat, file="resheat.txt", quote=F, sep="\t")

cd <- read.table("resheat.txt")
head(cd)
#########################################################################################################
##make the GO table for MWU
head(cd)

library(dplyr)
cd
go_input_heat = cd %>%
  tibble::rownames_to_column(var = "iso") %>%
  mutate(mutated_p = -log(pvalue)) %>%
  mutate(mutated_p_updown = ifelse(log2FoldChange < 0, mutated_p*-1, mutated_p*1)) %>%
  na.omit() %>%
  select(iso, mutated_p_updown)
head(go_input_heat)
colnames(go_input_heat) <- c("gene", "pval")
head(go_input_heat)
write.csv(go_input_heat, file="heat_GO.csv", quote=F, row.names=FALSE)
####################################################################################################


###Pco2 pH8 vs pH7.5
# summary(res)
# resph75 <- results(dds, contrast=c("treat","pH7.5", "pH8"))
# table(resph75$padj<0.05)
# # 0.1=118
# # 0.05=72
# # 0.01=43
# summary(resph75)
# 
# plotMA(resph75, main="pH8 vs pH7.5")
# plotMA(resph75, main="pH8 vs pH7.5", ylim=c(-2,2))
# 
# results <- as.data.frame(resph75)
# head(results)
# 
# nrow(resph75[resph75$padj<0.1 & resph75$log2FoldChange > 0 & !is.na(resph75$padj),])
# nrow(resph75[resph75$padj<0.1 & resph75$log2FoldChange < 0 & !is.na(resph75$padj),])
# #UP in 7.5 38
# #DOWN in 7.5 80
# 
# write.table(resph75, file="7.5_2016.txt", quote=F, sep="\t")
# 
# cd2 <- read.table("7.5_2016.txt")
# head(cd2)
# 
# ##make the GO table for MWU for 7.5
# head(cd2)
# 
# library(dplyr)
# go_input_7.5 = cd2 %>%
#   tibble::rownames_to_column(var = "iso") %>%
#   mutate(mutated_p = -log(pvalue)) %>%
#   mutate(mutated_p_updown = ifelse(log2FoldChange < 0, mutated_p*-1, mutated_p*1)) %>%
#   na.omit() %>%
#   select(iso, mutated_p_updown)
# head(go_input_7.5)
# colnames(go_input_7.5) <- c("gene", "pval")
# head(go_input_7.5)
# write.csv(go_input_7.5, file="7.5_GO.csv", quote=F, row.names=FALSE)
# 
# write.table(resph75, file="7.5_2016.txt", quote=F, sep="\t")


#--------------get pvals
valheat=cbind(resheat$pvalue, resheat$padj)
head(valheat)
colnames(valheat)=c("pval.heat", "padj.heat")
length(valheat[,1])
#16931
table(complete.cases(valheat))
# FALSE  TRUE 
# 445 16486
# 
# val75=cbind(resph75$pvalue, resph75$padj)
# head(val75)
# colnames(val75)=c("pval.75", "padj.75")
# length(val75[,1])
# table(complete.cases(val75))

######-------------make rlogdata and pvals table
rlog=rlogTransformation(dds, blind=TRUE) 
rld=assay(rlog)
head(rld)
#colnames(rld)=paste(colData$treat)
head(rld)
length(rld[,1])

# helpful cuz then you can sort by which are significant when plotting heatmap 
rldpvals=cbind(rld,valheat)
head(rldpvals)
dim(rldpvals)
# [1] 16931     6
table(complete.cases(rldpvals))
# FALSE  TRUE 
# 445 16486 

write.csv(rldpvals, "resheat_RLDandPVALS.csv", quote=F)


# DO WE NEED TO RUN ????
colnames(rld)=paste(colData$treat)
head(rld)

library(RColorBrewer)
# Sample distance heatmap
sampleDists <- as.matrix(dist(t(rld)))
library(gplots)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          margin=c(10, 10), main="Sample Distance Matrix")



#################################################################################
# VENN Diagram to include both up and down regulated genes in common for PC02
# library(VennDiagram)
# install.packages("VennDiagram")
# 
# heat_up=row.names(resheat[resheat$padj<0.1 & !is.na(resheat$padj) & resheat$log2FoldChange>0,])
# length(heat_up) #379
# heat_down=row.names(resheat[resheat$padj<0.1 & !is.na(resheat$padj) & resheat$log2FoldChange<0,])
# length(heat_down) #1839
# # p75_up=row.names(resph75[resph75$padj<0.1 & !is.na(resph75$padj) & resph75$log2FoldChange>0,])
# # length(p75_up) #38
# # p75_down=row.names(resph75[resph75$padj<0.1 & !is.na(resph75$padj) & resph75$log2FoldChange<0,])
# # length(p75_down) #80
# 
# pheat=row.names(resheat[resheat$padj<0.1 & !is.na(resheat$padj),])
# #p75=row.names(resph75[resph75$padj<0.1 & !is.na(resph75$padj),])
# 
# #UP
# pdegs_up=union(heat_up, heat_down)
# length(pdegs_up)
# #2218 # DID NOT REMOVE ANY ??
# 
# #DOWN
# pdegs_down=union(p76_down,p75_down)
# length(pdegs05_down)
# #432
# 
# #ALL
# pdegs05=union(p76,p75)
# length(pdegs05)
# #524
# 
# ###do UP, DOWN, ALL
# #change list to be whatever you want to plot
# candidates=list("7.6"=p76_up, "7.5"=p75_up)
# #quartz()
# prettyvenn=venn.diagram(
#   x = candidates,
#   filename=NULL,
#   col = "transparent",
#   fill = c("coral2", "forestgreen"),
#   alpha = 0.5,
#   # label.col = c("darkred", "white", "darkgreen", "white", "white", "white", "blue4"),
#   cex = 2.5,
#   fontfamily = "sans",
#   fontface = "bold",
#   cat.default.pos = "text",
#   cat.col = c("darkred", "darkgreen"),
#   cat.cex = 2.5,
#   cat.fontfamily = "sans",
#   cat.dist = c(0.08, 0.08),
#   cat.pos = 1
# );
# grid.draw(prettyvenn)

####################################################################################################


###########################heat map of sample distances for pco2
rldpvals <- read.csv(file="resheat_RLDandPVALS.csv", row.names=1)
head(rldpvals)
rld=rldpvals[,1:4]
head(rld)

sampleDists <- dist(t(rld))
sampleDistMatrix <- as.matrix( sampleDists )
treat=c( "Hot", "Control")
#colnames(sampleDistMatrix)=paste(treat)
#rownames(sampleDistMatrix)=paste(treat)

library("pheatmap")
heat.colors = colorRampPalette(rev(c("blue","yellow","red")),bias=0.3)(100)
pheatmap(sampleDistMatrix,color = heat.colors,cex=0.9,border_color=NA,cluster_rows=T,cluster_cols=T)

library(vegan)
library(ggplot2)
library(ggrepel)
library(tidyverse)

rld_t=t(rld)
#rld_t[ , which(apply(rld_t, 2, var) != 0)]
pca <- prcomp(rld_t,center = TRUE)
head(pca)
li <- pca$sdev^2 / sum(pca$sdev^2)
pc1v <- round(li[1] * 100, 1)
pc2v <- round(li[2] * 100, 1)
pca_s <- as.data.frame(pca$x)
head(pca_s)
pca_s <- pca_s[,c(1,2)]
pca_s$Samples = row.names(pca_s)
pca_s$treat=colData$temperature
head(pca_s)

cbPalette <- c("darkgoldenrod2",  "darkolivegreen3") #, "dodgerblue3", "red")
ggplot(pca_s, aes(PC1, PC2, color = treat, pch = treat)) +
  geom_point(size=3) +
  #  geom_text_repel(aes(label=Samples)) +
  scale_colour_manual(values=cbPalette)+
  theme_bw() +
  # geom_density2d(alpha=.5)+
  geom_polygon(alpha=.2)+
  xlab(paste0("PC1: ",pc1v,"% variance")) +
  ylab(paste0("PC2: ",pc2v,"% variance")) 
head(pca)
library(vegan)
adonis(pca$x ~ treat, data = pca_s, method='eu', na.rm = TRUE)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# treat      1   24166.6 24166.6  21.422 0.91461 0.3333
# Residuals  2    2256.3  1128.1         0.08539       
# Total      3   26422.9                 1.00000   


###################################heatmaps for genes NS vs FR
rldpvals <- read.csv(file="resheat_RLDandPVALS.csv", row.names=1)
head(rldpvals)
rld_site= rldpvals[,1:4]
head(rld_site)
gg=read.table("davies_Ssid_iso2gene.tab",sep="\t", row.names=1, quote ='')
head(gg)

nrow(rldpvals[rldpvals$padj.heat<0.01& !is.na(rldpvals$padj.heat),])
#1660

topnum= 100 # number of DEGS
head(rldpvals)
top100=head(rldpvals[order(rldpvals$padj.heat), ],topnum)
head(top100)
length(top100[,1])
summary(top100)
###
library(pheatmap)
head(top100)
p.val=0.1 # FDR cutoff
conds=top100[top100$padj.heat<=p.val & !is.na(top100$padj.heat),]
length(conds[,1])

exp=conds[,1:4] # change numbers to be your vsd data columns
means=apply(exp,1,mean) # means of rows
explc=exp-means # subtracting them
head(explc)

ccol=colorRampPalette(rev(c("red","chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)
col0=colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)

pheatmap(explc,cluster_cols=T,scale="row",color=ccol, show_rownames = F)

###################################Heatmap for the genes in common
rldpvals <- read.csv(file="resheat_RLDandPVALS.csv", row.names=1)
head(rldpvals)
p.val=0.1 # FDR cutoff
conds=rldpvals[rldpvals$padj.heat<=p.val & !is.na(rldpvals$padj.heat),]
rld_data= conds[,c(1:4)]
head(rld_data)
nrow(rld_data)
gg=read.table("davies_Ssid_iso2gene.tab",sep="\t", row.names=1, quote='')
library(pheatmap)
means=apply(rld_data,1,mean) # means of rows
explc=rld_data-means # subtracting them

ccol=colorRampPalette(rev(c("red","chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)
col0=colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)

pheatmap(explc,cluster_cols=T,scale="row",color=col0, show_rownames = F)

# Make annotation table for pheatmap
ann = data.frame(cond = c('Hot', 'Hot', 'Control', 'Control'))
rownames(ann) <- names(explc)

# Set colors
Var1        <- c("darkgoldenrod2",  "darkolivegreen3", "dodgerblue3")
names(Var1) <- c("Hot", "Control")
anno_colors <- list(cond = Var1)

pheatmap(as.matrix(explc),annotation_col=ann,annotation_colors=anno_colors,cex=1.2,color=col0,border_color=NA,clustering_distance_rows="correlation",clustering_distance_cols="correlation", show_rownames=T)
