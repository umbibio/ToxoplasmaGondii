library(limma)
library(openxlsx)
library(gplots)
library(dplyr)
library(tidyverse)
library(DESeq2)
library(sva)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(PoiClaClu)

############## Count matrix from feature count ################

fc.files.dir  <- "../Input/RNAseqCounts/"
fc.files  <- file.path(fc.files.dir, list.files(fc.files.dir)[grep("counts.txt$", list.files(fc.files.dir))])
f <- read.table(fc.files[1], header = T, sep = '\t', stringsAsFactors = F)
f <- f[,c(1,7)]
colnames(f) <- c("GeneName", "Count1")
for(i in 2:length(fc.files)){
  tmp <- read.table(fc.files[i], header = T, sep = '\t', stringsAsFactors = F)
  tmp <- tmp[,c(1,7)]
  colnames(tmp) <- c("GeneName", paste("Count", i, sep = ''))
  f = merge(f, tmp, by.x = 1, by.y = 1)
  
}

#write.xlsx(f, "../Output/count.table.xlsx")

############## Sample Info - experiment design ###############

## run5 are DNA data
RunInfo <- data.frame(Sample = gsub('_S.*', '', gsub('.counts.txt', '', basename(fc.files))), 
                      Count = paste('Count', 1:length(fc.files), sep = ''),
                      stringsAsFactors = F)
## This was manually constructed from RNAseq Info excell files
RNAseqInfo <- read.xlsx('../Input/RNAseqCounts/sampleInfo.xlsx')
RNAseqInfo <- left_join(RNAseqInfo, RunInfo, by = c('Sample.ID' = 'Sample'))

RNAseqInfo.B2 <- RNAseqInfo %>% dplyr::filter((Clone == 'B2' | Clone == 'RH' | Clone == 'B') & cond != 'fresh')
RNAseqInfo.B2$treatment <- paste(RNAseqInfo.B2$cond, RNAseqInfo.B2$passage, sep = '.')
RNAseqInfo.B2 <- RNAseqInfo.B2 %>% mutate(Names = paste(passage, cond, 'rep', Bio_rep, sep = '.'))

#write.xlsx(RNAseqInfo.B2, "../Output/sampleInfo.xlsx")
# count matrix befor filtering low expressions 
y <- f %>% dplyr::select(c('GeneName', RNAseqInfo.B2$Count)) ## Just B, RH and B2 samples
xx <- colnames(y[,2:ncol(y)])
yy <- RNAseqInfo.B2$Count

treatment <- RNAseqInfo.B2$treatment[match(xx,yy)]
treatment <- factor(treatment, levels = unique(treatment))

#write.xlsx(RNAseqInfo.B2, "../Output/SampleInfo.xlsx")
############## creating DEseq object ###############

## 1.count matrix
countDataMatrix <- as.matrix(y[ , -1])
rownames(countDataMatrix) <- y[ , 1]
colnames(countDataMatrix) <- RNAseqInfo.B2$Names

## 2.sample Info
coldata <- RNAseqInfo.B2

## 3. design -> treatment = passage + cond

## ?? have 2 columns batch , treatment to used in matrix modeling 

# DESeq object 
dds <- DESeqDataSetFromMatrix(countData = countDataMatrix,
                              colData = coldata,
                              design = ~treatment)
dds <- estimateSizeFactors(dds)
dds

# filtering the genes that have ab average of greater than 10 across all samples
dds <- dds[rowSums(counts(dds)) > 10, ]
nrow(dds)
########### transformed dds object #######################
# stabilizing the variance across the mean value using 
# transformation of raw counts table to plot PCA and MDS
# for statistical analysis we use the raw counts(not normalized one)

rld <- rlog(dds, blind=FALSE)
#head(assay(rld), 3)

# let's see the difference of log2 and rlog transformation 
par(mfrow = c( 1, 2 ))
plot(log2(counts(dds, normalized=TRUE)[,2:3] + 1),
     pch=16, cex=0.3, main = "log2 transformation")
plot(assay(rld)[,2:3],pch=16, cex=0.3, main = "rlog transformation")

# plotting first and second sample against each other 
# to show how the rlog transdormation staibilize the variance 
# the genes with low expression are getting close to mean values 
# the genes with high expression rlog gives similar result to the
#ordinary log2 transformation of normalized counts

############### similarity of samples ################

# sample distance -> euclidean method

sampleDists <- dist(t(assay(rld)))

sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")))(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors,
         main = "euclidean dist from DESeq bject before removing bad samples")


par(mfrow = c(1,1))
plotMDS(assay(rld), cex = 0.8)

clust <- hclust(dist(t(assay(rld))))
plot(clust)

clusterCut <- cutree(clust, k = 4)
CC <- data.frame(Samples = names(clusterCut), Batch = clusterCut, stringsAsFactors = F)
RNAseqInfo.B2 <- left_join(RNAseqInfo.B2, CC, by = c('Names' = 'Samples'))


# poisson distance
poisd <- PoissonDistance(t(counts(dds)))
samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- coldata$Names
colnames(samplePoisDistMatrix) <- coldata$Names
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows=poisd$dd,
         clustering_distance_cols=poisd$dd,
         col=colors, 
         main = "poisson distance from DESeq object before removing bad samples")

# PCA 

#plotPCA(rld, intgroup = c("cond", "passage"))
pcadata <- plotPCA(rld, intgroup = c("cond", "passage"), returnData=TRUE)
percentVar <- round(100 * attr(pcadata, "percentVar"))


ggplot(pcadata, aes(PC1, PC2, color= treatment, shape=cond)) + geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  labs(title="PCA before removing batch effect")

# MDS plot
mdsData <- data.frame(cmdscale(sampleDistMatrix))
mds <- cbind(mdsData, as.data.frame(colData(rld)))
ggplot(mds, aes(X1,X2,color=treatment,shape=cond)) + geom_point(size=3) +
  labs(title="MDS before removing batch effect")



## what are bad samples 
# 1.	P35.extra.rep.1
# 2.	P210.extra.rep.1
# 3.	RH.intra.rep.2
# 4.	RH.intra.rep.3
# 5.	P148.extra.rep.2
# 6.	Rh.extra.rep.2
# 7.	P85.extra.rep.3
# 8.	P7.intra.rep.2


########### removing hidden batch effect  using sva package ###########
library("sva")
dat <- counts(dds, normalized = TRUE)
idx <- rowMeans(dat) > 1
dat <- dat[idx,]
mod <- model.matrix(~ treatment, colData(dds))
mod0 <- model.matrix(~1, colData(dds))
svseq <- svaseq(dat, mod = mod, mod0 = mod0)
#n.sv = NULL svaseq package returns 3 surrogate variable.names()
svseq
#???? what should we specify as number of surrogate variables? 
# when it is not specified the sva package estimate 3 sv in this experiment 

par(mfrow=c(3,1),mar=c(3,5,3,1))
stripchart(svseq$sv[,1] ~ dds$Names,vertical=TRUE,main="SV1")
abline(h=0)
stripchart(svseq$sv[,2] ~ dds$Names,vertical=TRUE,main="SV2")
abline(h=0)
stripchart(svseq$sv[,3] ~ dds$Names,vertical=TRUE,main="SV3")
abline(h=0)

# DESeq object after removing batch
ddssva <- dds 
ddssva$SV1 <- svseq$sv[,1]
ddssva$SV2 <- svseq$sv[,2]
ddssva$SV3 <- svseq$sv[,3]
## add 3 surrogate variable to the dds design 
design(ddssva) <- ~ SV1 + SV2 + SV3 + treatment

# transformation of counts after sva
rldsva <- rlog(ddssva, blind=FALSE)
#head(assay(rldsva), 3)


# With the following plots still I get the same figures as before batch correction
# The batch effect has not been added to the design 

# dendogram 
par(mfrow = c(1,1))
plotMDS(assay(rldsva), cex = 0.8)
clust <- hclust(dist(t(assay(rldsva))))
plot(clust)

par(mfrow = c( 1, 2 ))
plot(log2(counts(ddssva, normalized=TRUE)[,2:3] + 1),
     pch=16, cex=0.3)
plot(assay(rldsva)[,2:3],pch=16, cex=0.3)


# sample dist after sva

sampleDists.sva <- dist( t( assay(rldsva) ) )
sampleDistMatrix.sva <- as.matrix( sampleDists.sva )
# rownames(sampleDistMatrix.sva) <- rldsva$Names
# colnames(sampleDistMatrix.sva) <- rldsva$Names
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix.sva,
         clustering_distance_rows=sampleDists.sva,
         clustering_distance_cols=sampleDists.sva,
         col=colors)

########################################################################################

### removing bad samples 
# # bad samples : 
bad.samples <- c("P210.extra.rep.1", "P35.extra.rep.1", "RH.intra.rep.2", "P148.extra.rep.2",
                 "RH.extra.rep.2","P85.extra.rep.3", "P7.intra.rep.2")

good.samples.indx <- !colnames(countDataMatrix) %in% bad.samples
sum(good.samples.indx == TRUE)

countDataMatrix <- countDataMatrix[, good.samples.indx]

coldata <- coldata[!(coldata$Names %in% bad.samples) ,]

dds <- DESeqDataSetFromMatrix(countData = countDataMatrix,
                              colData = coldata,
                              design = ~treatment)
dds <- estimateSizeFactors(dds)
dds


# filtering the genes that have ab average of greater than 10 across all samples

dds <- dds[rowSums(counts(dds)) > 10, ]
nrow(dds)

########### transformed dds object #######################
# stabilizing the variance across the mean value using rlog()

rld <- rlog(dds, blind=FALSE)

# let's see the difference of log2 and rlog transformation 
par(mfrow = c( 1, 2 ))
plot(log2(counts(dds, normalized=TRUE)[,2:3] + 1),
     pch=16, cex=0.3)
plot(assay(rld)[,2:3],pch=16, cex=0.3)

############### similarity of samples ################

# sample distance -> euclidean method

sampleDists <- dist(t(assay(rld)))

sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")))(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors,
         main = "euclidean dist from DESeq bject before removing bad samples")

par(mfrow = c(1,1))
plotMDS(assay(rld), cex = 0.8)
clust <- hclust(dist(t(assay(rld))))
plot(clust)

# poisson distance

poisd <- PoissonDistance(t(counts(dds)))
samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- coldata$Names
colnames(samplePoisDistMatrix) <- coldata$Names
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows=poisd$dd,
         clustering_distance_cols=poisd$dd,
         col=colors, 
         main = "poisson distance from DESeq object before removing bad samples")

# PCA 

plotPCA(rld, intgroup = c("cond", "passage"))
pcadata <- plotPCA(rld, intgroup = c("cond", "passage"), returnData=TRUE)
percentVar <- round(100 * attr(pcadata, "percentVar"))


# ggplot(pcadata, aes(PC1, PC2, color= treatment, shape=cond)) + geom_point(size=3) +
#   xlab(paste0("PC1: ",percentVar[1],"% variance")) +
#   ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
#   labs(title="PCA before removing batch effect")

# MDS plot
mdsData <- data.frame(cmdscale(sampleDistMatrix))
mds <- cbind(mdsData, as.data.frame(colData(rld)))
ggplot(mds, aes(X1,X2,color=treatment,shape=cond)) + geom_point(size=3) +
  labs(title="MDS before removing batch effect")




