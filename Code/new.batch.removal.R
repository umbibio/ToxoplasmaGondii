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
library(forcats)
library(ggfortify)
library(jcolors)
library(gridExtra)

source("./util_funcs.R")

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

############################################################################################
################################ Removing bad samples ######################################
################################ Creating DESeq Object #####################################
################################ Data transformation  ######################################
################################ Sva batch correction ######################################

############################################################################################

#bad.samples <- c("P210.extra.rep.1", "P35.extra.rep.1", "RH.intra.rep.2", "P148.extra.rep.2",
#                 "RH.extra.rep.2","P85.extra.rep.3", "P7.intra.rep.2")

bad.samples <- c("P210.extra.rep.1","P35.extra.rep.1", "P85.extra.rep.3", "P148.extra.rep.2", 
"RH.intra.rep.2")

## 1.count matrix
countDataMatrix <- as.matrix(y[ , -1])
rownames(countDataMatrix) <- y[ , 1]
colnames(countDataMatrix) <- RNAseqInfo.B2$Names
good.samples.indx <- !colnames(countDataMatrix) %in% bad.samples
sum(good.samples.indx == TRUE)
countDataMatrix <- countDataMatrix[, good.samples.indx]


## 2.sample Info
coldata <- RNAseqInfo.B2
coldata <- coldata[!(coldata$Names %in% bad.samples) ,]


######################################################
################### Gene Set selection ###############
######################################################


### keeping genes from specific gene sets ###

# 1  Brady In vivo, not in vitro
# 2. Brady In vitro, not in vivo
# 3. Tachyzoite, over bradyzoite in vivo
# 4. Bradyzoite in vivo, over tachyzoite
# 5. Tachyzoite, over bradyzoite in vitro
# 6. Compound1 Brady
# 7. Apicoplast
# 8. tgo01040_Biosynthesis_of_unsaturated_fatty_acids
# 9. All G1 (Monerri, Cell 2015)
# 10. All S/M (Monerri, Cell 2015)
# 11. Trafficking & Transport (Varberg, mBio 2018; Monerri, Cell 2015)
# 12. Invasion (Varberg, mBio 2018; ToxoDB)
# 13. Predicted Signal Peptide (stringent, ToxoDB)
# 14. TissueCysts/Tachyzoite
# 15. Can we do all 319 candidates, the ones that correlate with phenotype? These include hypotheticals
# 16. Can we do all genes that significantly and consistently trend up or trend down over time (regardless 
#     their correlation with phenotype)? 

FourWay <- read.xlsx("~/OneDrive - University of Massachusetts Boston/Toxo_yasaman_2019/may_10th_2019/FourWay_detailed_no_filter.xlsx")
#toxo_tab <- read.xlsx("~/OneDrive - University of Massachusetts Boston/Toxo_yasaman_2019/may_10th_2019/toxo_table_batch_corrected_logCPM_expression_edgeR_DEGs_ap2_targs_phenotype_correlations_smeV2.xlsx")
GeneSet.190 <- read.xlsx("~/OneDrive - University of Massachusetts Boston/Toxo_yasaman_2019/may_10th_2019/GSEA sytenic list 3.28.19.xlsx")
lifeStage.genes <- read.xlsx("~/OneDrive - University of Massachusetts Boston/Toxo_yasaman_2019/lifestage/lifestage.geses.3.4.quantiles.xlsx")

geneSets <- c("Brady.In.vivo.not.in.vitro", "Brady.In.vitro.not.in.vivo", "Brady.in.vivo.over.Tachy",
              "Tachy.over.Brady.in.vivo", "Tachy.over.Brady.in.vitro","Tachy.Ubiquitome.", "Tachy.RH","All.G1.", 
              "Compound1.Brady", "Apicoplast", "Biosynthesis.of.unsaturated.fatty.acids", "All.S/M.", 
              "Trafficking.&.Transport.", "Invasion.", "Predicted.Signal.Peptide.", "TissueCysts/Tachy.UPRegulated.", 
              "TissueCysts/Tachy.DownRegulated.")

GeneSet.190 <- prep_geneset.190(GeneSet.190)
lifeStage <- prep_geneset.190(lifeStage.genes)
GeneSet.190 <- rbind(GeneSet.190, lifeStage)
GeneSet.190 <- GeneSet.190 %>% unnest(genes) %>% dplyr::select(-c(total))
colnames(GeneSet.190) <- c('GeneSet', 'GeneName')

#tachy.over.brady.in.vivo <- GeneSet.190 %>% dplyr::filter(GeneSet == "Tachy.over.Brady.in.vivo") #ind 1
Hypo <- 
  read.xlsx("~/Dropbox/Toxoplasma Gondii/BatchCorrected/hypothetical_KO_candidates_04.30.19/hypotheticals_candidates.xlsx")

all.cand.GSE <- 
  read.xlsx("~/Dropbox/Toxoplasma Gondii/BatchCorrected/Files_KO_4_25_2019/all_candidates_geneset_enrichment.xlsx")


tab <- FourWay %>% dplyr::filter(quantile.cond1 == 4 | quantile.cond2 == 4) %>% dplyr::filter(qval < 0.05)

P11.extra <- Cont[grep("\\P11.extra$", Cont)]
tab <- subset(tab, tab$Contrast %in% P11.extra)

#g <- geneSets[7]
g <- "190_hypothetical_genes"

#GS <- tab %>% dplyr::filter(GeneSet %in% g)
genesss <- intersect(unique(tab$GeneName), Hypo$GeneName)
dim(GS)
#GS.dat <- subset(countDataMatrix, rownames(countDataMatrix) %in% unique(GS$GeneName))
GS.dat <- subset(countDataMatrix, rownames(countDataMatrix) %in% genesss)
dim(GS.dat)
Counts <- GS.dat
dim(Counts)

#countDataMatrix <- countDataMatrix[intersect(rownames(countDataMatrix), GS$GeneName),]

dds <- DESeqDataSetFromMatrix(countData = Counts,
                              colData = coldata,
                              design = ~treatment)
dds <- estimateSizeFactors(dds)
dds
dds <- dds[rowSums(counts(dds)) > 10, ] # filtering counts >10
nrow(dds)

# data transformation *before* batch correction  
rld <- rlog(dds, blind=FALSE, fitType='local')
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")))(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors,
         main = "Euclidean dist - 190 hypothetical genes -before batch correction")


# sample variation before and after transformation

par(mfrow = c( 1, 2 ))
plot(log2(counts(dds, normalized=TRUE)[,1:4] + 1),
     pch=16, cex=0.3)
plot(assay(rld)[,1:4],pch=16, cex=0.3)

# batch correction 

dat <- counts(dds, normalized = TRUE)
idx <- rowMeans(dat) > 1
dat <- dat[idx,]
mod <- model.matrix(~ treatment, colData(dds))
mod0 <- model.matrix(~1, colData(dds))
svseq <- svaseq(dat = dat, mod = mod, mod0 = mod0, n.sv = 1)
svseq


par(mfrow=c(3,1),mar=c(3,5,3,1))
stripchart(svseq$sv[,1] ~ dds$treatment,vertical=TRUE,main="SV1")
abline(h=0)

# DESeq object after identifying sva

ddssva <- dds
ddssva$SV1 <- svseq$sv[,1]
design(ddssva) <- ~ SV1  +  treatment

# new deseq transformation
rldsva <- rlog(ddssva, blind = FALSE)

# removing batch effects
# using limma package 

cov <- cbind(svseq$sv[,1:svseq$n.sv])
mat <- assay(rldsva)
mat <- limma::removeBatchEffect(mat, covariates = cov)
assay(rldsva) <- mat

############################ Plots #################################

# PCA plot
# sample PCA function for transformed data 
# input: DESeq transformed object assay(rldsva)

plotPCA(rldsva, intgroup = c("cond", "passage"))
pcadata <- plotPCA(rldsva, intgroup = c("cond", "passage"), returnData = T)
percentVar <- round(100 * attr(pcadata, "percentVar"))

plot1 <- ggplot(pcadata, aes(PC1, PC2, colour = passage, shape = cond)) + 
  geom_point(size = 3)+
  xlab(paste("PC1:", percentVar[1], "% Variance"))+
  ylab(paste("PC2:", percentVar[2], "% Variance"))+
  theme(plot.title = element_text(hjust = 0.5), title = element_text(size = 9))+
  ggtitle(g)
  
plot1

# MDS plot 
# multidimensional scaling (MDS) 
# when we don’t have a matrix of data, but only a matrix of distances.
# multidimentional scaling takes a set of dissimilarities and returns a set of points
#such that the distances between the points are approximately equal to the dissimilarities

sampleDists.sva <- dist(t(assay(rldsva)))
sampleDistMatrix.sva <- as.matrix(sampleDists.sva)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")))(255)
pheatmap(sampleDistMatrix.sva,
         clustering_distance_rows=sampleDists.sva,
         clustering_distance_cols=sampleDists.sva,
         col=colors,
         main = "Euclidean dist - good samples -after batch correction")

mdsData <- data.frame(cmdscale(sampleDistMatrix.sva))
mds <- cbind(mdsData, as.data.frame(colData(rldsva)))

plot2 <- ggplot(mds, aes(X1, X2, colour = passage, shape = cond))+ geom_point(size = 3)+
  theme(plot.title = element_text(hjust = 0.5), title = element_text(size = 9))+
  ggtitle(g)
plot2
 
# pdf(paste("../Output/June_26/", paste(g, ".pdf", sep = ""), sep = ""))
# grid.arrange(plot1)
# dev.off()


####################################################################################
###################### PCA on means of replicate ###################################
####################################################################################

# mean of replicates saved in countTable
countTable <- data.frame(assay(rldsva))
nm <- sub("\\d$", "", colnames(countTable))
countTable.mean <- sapply(split.default(countTable, nm), rowMeans)
colnames(countTable.mean) <- gsub("rep.", "mean", colnames(countTable.mean))

# extra & intra 
nn <- colnames(countTable.mean) <- gsub("RH", "P300", colnames(countTable.mean))
nn <- gsub(".mean", "", nn)
nn.passage <- sort(unique(as.numeric(gsub('P', '', gsub('\\.intra', '', gsub('\\.extra', '', nn))))), index.return = T)$ix
nn.extra <- grep('extra', nn)
nn.intra <- grep('intra', nn)
ind1 <- nn.extra[sort(as.numeric(gsub('\\.extra', '', gsub('P', '', nn[nn.extra]))), index.return = T)$ix]
ind2 <- nn.intra[sort(as.numeric(gsub('\\.intra', '', gsub('P', '', nn[nn.intra]))), index.return = T)$ix]
ind <- c(ind1, ind2)
colnames(countTable.mean) <- gsub("P300", "RH", colnames(countTable.mean))

df <- data.frame(t(countTable.mean))
#colnames(df) <- rownames(toxo.mean)
rownames(df) <- NULL
df_pca <- prcomp(df,center = T, scale. = T)
df_out <- as.data.frame(df_pca$x)
df_out <- df_out[, 1:2]
df_out$sample <- colnames(countTable.mean)
df_out$passage <- gsub('\\.intra.*', '', gsub('\\.extra.*', '', colnames(countTable.mean)))
df_out$condition <- gsub('RH\\.', '', gsub('P.*[[:digit:]]\\.', '', gsub('\\.rep.*', '', colnames(countTable.mean))))
df_out$sample <- factor(df_out$sample, levels = unique(df_out$sample)[ind])
df_out$passage <- factor(df_out$passage, levels = unique(df_out$passage)[nn.passage])
df_out$passage_condition <- paste(df_out$passage, df_out$condition, sep = '_')


plot3 <- ggplot(df_out, aes(PC1, PC2, colour = passage, shape = condition)) +
  geom_point(size = 3)+
  xlab(paste("PC1:", percentVar[1], "% Variance"))+
  ylab(paste("PC2:", percentVar[2], "% Variance"))+
  theme(plot.title = element_text(hjust = 0.5), title = element_text(size = 9))+
  ggtitle(g)
  
plot3
#ggtitle(expression(atop("Mean expression of replicates", atop(italic("190 hypo genes"), ""))))

###########################################################################
####################### MDS on means of Replicates ########################
###########################################################################

sampleDist.mean <- dist(t(countTable.mean))
sampledist.matrix.means <- as.matrix(sampleDist.mean)
mds.mean.data <- data.frame(cmdscale(sampledist.matrix.means))
mds.mean.data <- cbind(mds.mean.data,  df_out[,3:6])

plot4 <- ggplot(mds.mean.data, aes(X1, X2, colour = passage, shape = condition))+
  geom_point(size = 3)+
  theme(plot.title = element_text(hjust = 0.5),  title = element_text(size = 9))+
  ggtitle(g)

plot4
############################################################################
########################### Distance Original Coordinates ##################
############################################################################

sampleDist.orig <- dist(t(countTable.mean))
sampleDistMatrix.orig <- as.matrix(sampleDist.orig)
sampleDist.df.orig <- data.frame(sampleDistMatrix.orig)
colnames(sampleDist.df.orig) <- gsub(".mean", "", colnames(sampleDist.df.orig))
rownames(sampleDist.df.orig) <- gsub(".mean", "", rownames(sampleDist.df.orig))
# rownames(sampleDist.df.orig) <- gsub("P300" , "RH", rownames(sampleDist.df.orig))
# colnames(sampleDist.df.orig) <- gsub("P300" , "RH", colnames(sampleDist.df.orig))

extra.df <- sampleDist.df.orig[ind1,]
extra.df <- tibble::rownames_to_column(extra.df, "samples")
extra.df  <- extra.df[, c("samples", "RH.extra")]

intra.df <- sampleDist.df.orig[ind2,]
intra.df <- tibble::rownames_to_column(intra.df, "samples")
intra.df  <- intra.df[, c("samples", "RH.intra")]

plot5 <- ggplot(extra.df, aes(x =fct_inorder(samples), y = RH.extra.mean, group = 1)) + 
  geom_line(aes(x = fct_inorder(samples), y = RH.extra),colour = "blue")+
  theme(plot.title = element_text(hjust = 0.5), axis.text=element_text(size=7),
        title =element_text(size=10, face='bold'))+
  xlab("samples")+
  ylab("Distance to RH")+
  ggtitle(paste(g, "Original Coordinates", sep = " - "))

plot5


plot6 <- ggplot(intra.df, aes(x = fct_inorder(samples), y = RH.intra.mean, group = 1)) + 
  geom_line(aes(x = fct_inorder(samples), y = RH.intra), colour = "red") +
  theme(plot.title = element_text(hjust = 0.5), title =element_text(size=10, face='bold')) +
  xlab("samples")+
  ylab("Distance to RH")+
  ggtitle(paste(g, "Original Coordinates",  sep = " - "))

plot6



bar.extra <- ggplot(extra.df, aes(x = fct_inorder(samples), y = RH.extra.mean))+ 
  geom_bar(stat = "identity", colour = "darkblue", fill = "blue", width = 0.7)+
  theme(plot.title = element_text(hjust = 0.5))+
  xlab("Samples")+
  ggtitle("Euclidean distance To RH")
bar.extra

bar.intra <- ggplot(intra.df, aes(x = fct_inorder(samples), y = RH.intra.mean))+
  geom_bar(stat = "identity", width = 0.6,colour = "darkred", fill = "darkred")+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Euclidean Distance To RH")+
  xlab("samples")
bar.intra


############################################################################
########################### Distance  PCA coordinates ######################
############################################################################

PCA.df <- df_out[,1:2]
#colnames(PCA.df) <- df_out$sample
rownames(PCA.df) <- df_out$sample
PCA.dist <- as.matrix(dist(PCA.df))
colnames(PCA.dist) <- gsub(".mean", "", colnames(PCA.dist))
rownames(PCA.dist) <- gsub(".mean", "", rownames(PCA.dist))


extra.PCA.df <- data.frame(PCA.dist[ind1,])
extra.PCA.df <- tibble::rownames_to_column(extra.PCA.df, "samples")
extra.PCA.df  <- extra.PCA.df[, c("samples", "RH.extra")]

intra.PCA.df <- data.frame(PCA.dist[ind2,])
intra.PCA.df <- tibble::rownames_to_column(intra.PCA.df, "samples")
intra.PCA.df  <- intra.PCA.df[, c("samples", "RH.intra")]


plot7 <- ggplot(extra.PCA.df, aes(x =fct_inorder(samples), y = RH.extra.mean, group = 1)) +
  geom_line(aes(x = fct_inorder(samples), y = RH.extra),colour = "darkblue")+
  theme(plot.title = element_text(hjust = 0.5), title =element_text(size=10, face='bold'))+
  xlab("samples")+
  ylab("Distance to RH")+
  ggtitle(paste(g, "PCA Coordinates",  sep = " - "))

plot7


plot8 <- ggplot(intra.PCA.df, aes(x =fct_inorder(samples), y = RH.intra.mean, group = 1)) +
  geom_line(aes(x = fct_inorder(samples), y = RH.intra),colour = "red")+
  theme(plot.title = element_text(hjust = 0.5), title =element_text(size=10, face='bold'))+
  xlab("samples")+
  ylab("Distance to RH")+
  ggtitle(paste(g, "PCA Coordinates",  sep = " - "))

plot8



MyPlots = list(plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8)
pdf(paste("../Output/June_26/", paste(g, ".pdf", sep = ""), sep = ""))
MyPlots
dev.off()
#
bar.PCA.extra <- ggplot(extra.PCA.df, aes(x = fct_inorder(samples), y = RH.extra.mean))+ 
  geom_bar(stat = "identity", colour = "darkblue", fill = "darkblue")+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle(g)
  
bar.PCA.extra



bar.PCA.intra <- ggplot(intra.PCA.df, aes(x = fct_inorder(samples), y = RH.intra.mean))+ 
  geom_bar(stat = "identity", colour = "darkred", fill = "darkred")+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle(g)

bar.PCA.intra


################## MDS distance #######################

MDS.matrix <- sampledist.matrix.means
colnames(MDS.matrix) <- gsub(".mean", "", colnames(MDS.matrix))
rownames(MDS.matrix) <- gsub(".mean", "", rownames(MDS.matrix))


extra.mds.df <- data.frame(MDS.matrix[ind1,])
extra.mds.df <- tibble::rownames_to_column(extra.mds.df, "samples")
extra.mds.df  <- extra.mds.df[, c("samples", "RH.extra")]

intra.mds.df <- data.frame(MDS.matrix[ind2,])
intra.mds.df <- tibble::rownames_to_column(intra.mds.df, "samples")
intra.mds.df  <- intra.mds.df[, c("samples", "RH.intra")]

plot9 <- ggplot(extra.mds.df, aes(x =fct_inorder(samples), y = RH.extra, group = 1)) +
  geom_line(aes(x = fct_inorder(samples), y = RH.extra),colour = "darkblue")+
  theme(plot.title = element_text(hjust = 0.5), title =element_text(size=10, face='bold'))+
  xlab("samples")+
  ylab("Distance to RH")+
  ggtitle(paste(g, "mds Coordinates",  sep = " - "))

plot9


plot10 <- ggplot(intra.mds.df, aes(x =fct_inorder(samples), y = RH.intra, group = 1)) +
  geom_line(aes(x = fct_inorder(samples), y = RH.intra),colour = "red")+
  theme(plot.title = element_text(hjust = 0.5), title =element_text(size=10, face='bold'))+
  xlab("samples")+
  ylab("Distance to RH")+
  ggtitle(paste(g, "MDS Coordinates",  sep = " - "))

plot10

MyPlots = list(plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8, plot9, plot10)
pdf(paste("../Output/June_26/", paste(g, "_DEG_all_passages_vs_P11.extra.pdf", sep = ""), sep = ""))
MyPlots
dev.off()

# autoplot(df_pca, data = df_out,  colour = 'passage', shape = 'condition', size =4,
#          stroke = 4) +  scale_color_jcolors(palette = "rainbow") +
#   scale_shape_manual(values=c(19, 15)) +
#   theme_bw()+
#   labs(title="PCA - good samples + mean expression + SV1 ")



####################################################################################
################ PCA and MDS on replicates mean expression value ###################
####################################################################################

# mean of replicates saved in countTable
countTable <- data.frame(assay(rldsva))
nm <- sub("\\d$", "", colnames(countTable))
countTable.mean <- sapply(split.default(countTable, nm), rowMeans)
colnames(countTable.mean) <- gsub("rep.", "mean", colnames(countTable.mean))

# extra & intra 
nn <- colnames(countTable.mean) <- gsub("RH", "P300", colnames(countTable.mean))
nn <- gsub(".mean", "", nn)
nn.passage <- sort(unique(as.numeric(gsub('P', '', gsub('\\.intra', '', gsub('\\.extra', '', nn))))), index.return = T)$ix
nn.extra <- grep('extra', nn)
nn.intra <- grep('intra', nn)
ind1 <- nn.extra[sort(as.numeric(gsub('\\.extra', '', gsub('P', '', nn[nn.extra]))), index.return = T)$ix]
ind2 <- nn.intra[sort(as.numeric(gsub('\\.intra', '', gsub('P', '', nn[nn.intra]))), index.return = T)$ix]
ind <- c(ind1, ind2)

# distance 
# sampleDists <- dist(t(assay(rldsva)))
# sampleDistMatrix <- as.matrix(sampleDists)

sampleDist.orig <- dist(t(countTable.mean))
sampleDistMatrix.orig <- as.matrix(sampleDist.orig)
sampleDist.df.orig <- data.frame(sampleDistMatrix.orig)

extra.df <- sampleDist.df.orig[ind1,]
grep("p300", colnames(extra.df))
extra.df <- tibble::rownames_to_column(extra.df, "samples")
extra.df  <- extra.df[, c("samples", "RH.extra.mean")]

intra.df <- sampleDist.df[ind2,]
intra.df <- tibble::rownames_to_column(intra.df, "samples")
intra.df  <- extra.df[, c("samples", "RH.intra.mean")]


P.extra <- ggplot(extra.df, aes(x =fct_inorder(samples), y = RH.extra.mean, group = 1)) + 
  geom_line(aes(x = fct_inorder(samples), y = RH.extra.mean),colour = "red")+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Extracellular samples' progression towards RH" )

P.extra


P.intra <- ggplot(intra.df, aes(x = fct_inorder(samples), y = RH.intra.mean, group = 1)) + 
  geom_line(aes(x = fct_inorder(samples), y = RH.intra.mean), colour = "blue") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Intracellular samples progression towards RH")
P.intra


bar.extra <- ggplot(extra.df, aes(x = fct_inorder(samples), y = RH.extra.mean))+ 
  geom_bar(stat = "identity")
bar.extra

bar.intra <- ggplot(intra.df, aes(x = samples, y = RH.extra.mean))+
  geom_bar(stat = "identity", colour = "red", fill = "red")
bar.intra









