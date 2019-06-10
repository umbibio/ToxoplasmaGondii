library(sva)
library(bladderbatch)
data(bladderdata)
library(pamr)
library(limma)

pheno = pData(bladderEset)
edata = exprs(bladderEset)

mod = model.matrix(~as.factor(cancer), data=pheno)
mod0 = model.matrix(~1,data=pheno)
n.sv = num.sv(edata,mod,method="leek")
n.sv
svobj = sva(edata,mod,mod0,n.sv=n.sv)


batch = pheno$batch
modcombat = model.matrix(~1, data=pheno)
combat_edata = ComBat(dat=edata, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots = TRUE)


library("airway")
dir <- system.file("extdata", package="airway", mustWork=TRUE)
list.files(dir)
csvfile <- file.path(dir,"sample_table.csv")
(sampleTable <- read.csv(csvfile,row.names=1))
sampleTable <- read.csv(csvfile,row.names=1)
sampleTable
filenames <- file.path(dir, paste0(sampleTable$Run, "_subset.bam"))
file.exists(filenames)
library("Rsamtools")
bamfiles <- BamFileList(filenames, yieldSize=2000000)
seqinfo(bamfiles[1])
library("GenomicFeatures")
gtffile <- file.path(dir,"Homo_sapiens.GRCh37.75_subset.gtf")
(txdb <- makeTxDbFromGFF(gtffile, format="gtf", circ_seqs=character()))
(ebg <- exonsBy(txdb, by="gene"))
library("GenomicAlignments")
library("BiocParallel")
register(SerialParam())
se <- summarizeOverlaps(features=ebg, reads=bamfiles,
                        mode="Union",
                        singleEnd=FALSE,
                        ignore.strand=TRUE,
                        fragments=TRUE )
se
dim(se)
assayNames(se)
head(assay(se), 3)
colSums(assay(se))
rowRanges(se)
str(metadata(rowRanges(se)))
colData(se)
(colData(se) <- DataFrame(sampleTable))


se$cell
se$dex
se$dex <- relevel(se$dex, "untrt")
se$dex

## data is ready 
data("airway")
se <- airway
se$cell
se$dex
se$dex <- relevel(se$dex, "untrt") # factor to be reference level
se$dex
assay(se) # raw counts 
round( colSums(assay(se)) / 1e6, 1 ) # millions of fragments that uniquely aligned to the genes
colData(se) # SummarizedExperiment object 

library("DESeq2")
dds <- DESeqDataSet(se, design = ~cell + dex)

countdata <- assay(se)
head(countdata, 3)
coldata <- colData(se)

ddsMat <- DESeqDataSetFromMatrix(countData = countdata,
                                 colData = coldata,
                                 design = ~ cell + dex)

ddsMat

ddsMat <- ddsMat[ rowSums(counts(ddsMat)) > 1, ]


# transformation and visual PCA 
rld <- rlog(ddsMat, blind=FALSE)
head(assay(rld), 3)
# show the effect of transformation 
par( mfrow = c( 1, 2 ) )
ddsMat <- estimateSizeFactors(ddsMat)
plot(log2(counts(ddsMat, normalized=TRUE)[,1:2] + 1),
     pch=16, cex=0.3)
plot(assay(rld)[,1:2],
     pch=16, cex=0.3)
# sample dist 
sampleDists <- dist( t( assay(rld) ) )
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$dex, rld$cell, sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

# PCA plot 
plotPCA(rld, intgroup = c("dex", "cell"))

# MDS plot 
mdsData <- data.frame(cmdscale(sampleDistMatrix))
mds <- cbind(mdsData, as.data.frame(colData(rld)))
ggplot(mds, aes(X1,X2,color=dex,shape=cell)) + geom_point(size=3)

# Differential expression 

# building results

# plotting the results 

# removing batch effects



# Exploratory analysis and visualization
# Pre-filtering the dataset
# The rlog transformation
# Sample distances
# PCA 
# MDS
# Differential expression analysis
# Building the results table
# several other analysis 

## Removing hidden batch effects
library("sva")
dds <- estimateSizeFactors(dds)
dat <- counts(dds, normalized=TRUE)
idx <- rowMeans(dat) > 1
dat <- dat[idx,]
mod <- model.matrix(~dex, colData(dds))
mod0 <- model.matrix(~ 1, colData(dds))
svseq <- svaseq(dat = dat,mod = mod, mod0 = mod0, n.sv= 2)
svseq$sv
ddssva <- dds
ddssva$SV1 <- svseq$sv[,1]
ddssva$SV2 <- svseq$sv[,2]
design(ddssva) <- ~ SV1 + SV2 + dex
ddssva
ddssva <- DESeq(ddssva)
assay(ddssva)
