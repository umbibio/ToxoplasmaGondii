library(tidyverse)
library(openxlsx)
library(ggplot2)
library(sme)

source("./util_funcs.R")

## Intersting gene that correlate with Phenotypes for KO experiments

## Lets start by ploting the phenotypes
plaque_sep_rep <- read.xlsx('../Plaque_sep_rep.xlsx')
reinvasion_sep_rep <- read.xlsx('../Reinvasion_sep_rep.xlsx')
replication_sep_rep <- read.xlsx('../Replication_sep_repV2.xlsx')
survival_sep_rep <- read.xlsx('../survival_sep_rep.xlsx')
pheno.fits <- readRDS(file = "~/work/ToxoplasmaGondii/ToxoR/output/sme_fits_phenoV2.RData")

pohenotype.tc <- get_pohenotype_time_course(plaque_sep_rep, reinvasion_sep_rep,
                                            replication_sep_rep, survival_sep_rep)

my.pheno <- unique(pohenotype.tc$variable) ## 9 phenotypes

out.pic <- "~/work/ToxoplasmaGondii/KO_candidates/plots/phenotypes.pdf"
pdf(out.pic, width=8, height=8)
par(mfrow = c(3,3))
for(phe in my.pheno){
  fit.ph <- pheno.fits[[which(unique(pohenotype.tc$variable) == phe)]]
  plot.sme(fit.ph, phe)
}
dev.off()
## Judging by this graph, phenotypes seem to be very noisy. 
## Reinvasion and V8 replication may be more accurate?
## Consider plaque, reinv, V8, and survival for correlation analysis
cor.to.consider <- c("plaque", "reinv", "V8", "survival")

## Make a FourWay table (including all genes)
toxo_tab <- read.xlsx('../toxo_table_batch_corrected_logCPM_expression_edgeR_DEGs_ap2_targs_phenotype_correlations_smeV2.xlsx')
## Gene Sets
GeneSet.190 <- read.xlsx('../GSEA sytenic list 3.28.19.xlsx')
GeneSet.190 <- prep_geneset.190(GeneSet.190)
GeneSet.190 <- GeneSet.190 %>% unnest(genes) %>% dplyr::select(-c('total'))
colnames(GeneSet.190) <- c('GeneSet', 'GeneName')

## Checking APIAp2 targets
AP2.targets <- prep_ap2_targets(toxo_tab)
AP2.targets <- AP2.targets %>% unnest(genes) %>% dplyr::select(-c('total'))
colnames(AP2.targets) <- c('AP2', 'GeneName')

## Fold Change, Qvalue, Quantiles, Mean exprs
toxo.fc.qval.quant <- get_fc_qval_quant_tabs(toxo_tab)
FourWay <- make_four_way(toxo_tab, toxo.fc.qval.quant, GeneSet.190, AP2.targets)
FourWay <- FourWay %>% 
  dplyr::select(GeneSet:Product.Description,c(paste("core.s.",cor.to.consider, sep = '')))
colnames(FourWay) <- gsub('core.s.', '', colnames(FourWay))

write.xlsx(FourWay, "~/work/ToxoplasmaGondii/FourWay_detailed_no_filter.xlsx")

FourWay <- FourWay %>% gather(key = phenotype, value = correlation, plaque:survival)

## Phenotypes (plaque, reinvasion, and V8) are all trending up
## For each phenotype identify genes that
## 1) Highly and positively correlate ( > 0.7)
## 2) Are highly expressed at some passage
## 3) Are highly up-regulated compared to earlier passage P11 (early adaptation)
FourWay.pos <- FourWay %>% dplyr::filter(Category == 'extra.vs.extra' & Passage1 == 'P11') %>% 
  spread(key = phenotype, value = correlation, fill = NA) %>% 
  dplyr::filter(fc > 1 & qval < 0.05 & quantile.cond2 >= 2) %>%
  dplyr::filter(plaque > 0.7 | reinv > 0.7 | survival > 0.7 | V8 > 0.7)

keep.columns <- c("GeneSet", "GeneName", "AP2", "Product.Description")
FourWay.pos.select <- FourWay.pos %>% dplyr::select(keep.columns, plaque:V8) %>% distinct()
  
FourWay.pos.sum <- FourWay.pos.select %>% group_by(GeneName) %>% 
  summarise(Product.Description = unique(Product.Description), 
            GeneSet = list(unique(as.character(GeneSet))), 
            AP2 = list(unique(as.character(AP2))),
            plaque = unique(plaque), survival = unique(survival), reinv = unique(reinv), V8 = unique(V8))
FourWay.pos.sum <- FourWay.pos.sum %>% 
  dplyr::filter(!is.na(Product.Description) & Product.Description != "hypothetical protein")

write.xlsx(FourWay.pos.sum, "~/work/ToxoplasmaGondii/KO_candidates/files/positive_cor_summary.xlsx")
## plot the expressions
tc.expr.extra <- get_extra_time_course(toxo_tab)
expr.fits <- readRDS(file = "~/work/ToxoplasmaGondii/ToxoR/output/sme_fits_exprs.RData")

FourWay.pos.sum.g <- FourWay.pos.sum %>% gather(key = phenotype, value = correlation, plaque:V8) %>%
  dplyr::select(GeneName, Product.Description, phenotype, correlation) %>% mutate(correlation = round(correlation, 2))

for(ph in unique(FourWay.pos.sum.g$phenotype)){
  XX <- FourWay.pos.sum.g %>% dplyr::filter(phenotype == ph) %>% dplyr::filter(correlation > 0.7)
  num.plots <- ceiling(nrow(XX)/8)
  for(j in 1:num.plots){
    out.pic <- paste("~/work/ToxoplasmaGondii/KO_candidates/plots/trend_positive_cor_", ph, '_', j, '.pdf', sep = '')
    pdf(out.pic, width=8, height=8)
    par(mfrow = c(3,3))
    fit.ph <- pheno.fits[[which(unique(pohenotype.tc$variable) == ph)]]
    plot.sme(fit.ph, ph)
    for(i in ((j-1) * 8 + 1):min((j * 8), nrow(XX)) ){
      my.gene <- XX[i,c(1,2,4)]
      fit.gene <- expr.fits[[which(unique(tc.expr.extra$variable) == as.character(my.gene[1]))]]
      plot.sme(fit.gene, as.character(my.gene[1]),as.character(my.gene[2]))
      text(x = 150, min(fit.gene$data$y) + (-min(fit.gene$data$y) + max(fit.gene$data$y))/2, 
           paste('cor:', my.gene[3]), cex = 1)
      
    }
    dev.off()
  }
  
}


## Negative Correlation
## Phenotypes (plaque, reinvasion, and V8) are all trending up
## For each phenotype identify genes that
## 1) Highly and negatively correlate ( < -0.7)
## 2) Is highly expressed at P11 (early)
## 3) Are highly down-regulated compared to earlier passage P11 (early adaptation)
FourWay.neg <- FourWay %>% dplyr::filter(Category == 'extra.vs.extra' & Passage1 == 'P11') %>% 
  spread(key = phenotype, value = correlation, fill = NA) %>% 
  dplyr::filter(fc < -1 & qval < 0.05 & quantile.cond1 >= 2) %>%
  dplyr::filter(plaque < -0.7 | reinv < -0.7 | survival < -0.7 | V8 < -0.7)

keep.columns <- c("GeneSet", "GeneName", "AP2", "Product.Description")
FourWay.neg.select <- FourWay.neg %>% dplyr::select(keep.columns, plaque:V8) %>% distinct()

FourWay.neg.sum <- FourWay.neg.select %>% group_by(GeneName) %>% 
  summarise(Product.Description = unique(Product.Description), 
            GeneSet = list(unique(as.character(GeneSet))), 
            AP2 = list(unique(as.character(AP2))),
            plaque = unique(plaque), survival = unique(survival), reinv = unique(reinv), V8 = unique(V8))
FourWay.neg.sum <- FourWay.neg.sum %>% 
  dplyr::filter(!is.na(Product.Description) & Product.Description != "hypothetical protein")

write.xlsx(FourWay.neg.sum, "~/work/ToxoplasmaGondii/KO_candidates/files/negative_cor_summary.xlsx")

FourWay.neg.sum.g <- FourWay.neg.sum %>% gather(key = phenotype, value = correlation, plaque:V8) %>%
  dplyr::select(GeneName, Product.Description, phenotype, correlation) %>% 
  mutate(correlation = round(correlation, 2))

for(ph in unique(FourWay.neg.sum.g$phenotype)){
  XX <- FourWay.neg.sum.g %>% dplyr::filter(phenotype == ph) %>% dplyr::filter(correlation < -0.7)
  num.plots <- ceiling(nrow(XX)/8)
  for(j in 1:num.plots){
    out.pic <- paste("~/work/ToxoplasmaGondii/KO_candidates/plots/trend_negative_cor_", ph, '_', j, '.pdf', sep = '')
    pdf(out.pic, width=8, height=8)
    par(mfrow = c(3,3))
    fit.ph <- pheno.fits[[which(unique(pohenotype.tc$variable) == ph)]]
    plot.sme(fit.ph, ph)
    for(i in ((j-1) * 8 + 1):min((j * 8), nrow(XX)) ){
      my.gene <- XX[i,c(1,2,4)]
      fit.gene <- expr.fits[[which(unique(tc.expr.extra$variable) == as.character(my.gene[1]))]]
      plot.sme(fit.gene, as.character(my.gene[1]),as.character(my.gene[2]))
      text(x = 150, min(fit.gene$data$y) + (-min(fit.gene$data$y) + max(fit.gene$data$y))/2, 
           paste('cor:', my.gene[3]), cex = 1)
      
    }
    dev.off()
  }
  
}


## What are these genes enriched in
Candidates <- rbind(FourWay.pos.sum, FourWay.neg.sum)
Candidates.filt <- Candidates %>% 
  dplyr::filter(!is.na(Product.Description) & Product.Description != "hypothetical protein")
GeneSet.190 <- read.xlsx('../GSEA sytenic list 3.28.19.xlsx')
GeneSet.190 <- prep_geneset.190(GeneSet.190)

GeneSet.list <- Candidates.filt %>% dplyr::select(GeneName) %>% 
  summarise(genes = list(as.character(GeneName)), total = n())


candidates.enrich <- fisherEnrichment(GeneSet.190, GeneSet.list)
candidates.enrich <- candidates.enrich %>% dplyr::filter(pvalue < 0.01)

write.xlsx(candidates.enrich, "~/work/ToxoplasmaGondii/KO_candidates/files/all_candidates_geneset_enrichment.xlsx")

## look at hypotheticals
FourWay.pos.sum <- FourWay.pos.select %>% group_by(GeneName) %>% 
  summarise(Product.Description = unique(Product.Description), 
            GeneSet = list(unique(as.character(GeneSet))), 
            AP2 = list(unique(as.character(AP2))),
            plaque = unique(plaque), survival = unique(survival), reinv = unique(reinv), V8 = unique(V8))
FourWay.neg.sum <- FourWay.neg.select %>% group_by(GeneName) %>% 
  summarise(Product.Description = unique(Product.Description), 
            GeneSet = list(unique(as.character(GeneSet))), 
            AP2 = list(unique(as.character(AP2))),
            plaque = unique(plaque), survival = unique(survival), reinv = unique(reinv), V8 = unique(V8))
Candidates <- rbind(FourWay.pos.sum, FourWay.neg.sum)
Candidates.hyp <- Candidates %>% 
  dplyr::filter(is.na(Product.Description) | Product.Description == "hypothetical protein")

GeneSet.190 <- read.xlsx('../GSEA sytenic list 3.28.19.xlsx')
GeneSet.190 <- prep_geneset.190(GeneSet.190)

GeneSet.list <- Candidates.hyp %>% dplyr::select(GeneName) %>% 
  summarise(genes = list(as.character(GeneName)), total = n())


candidates.enrich <- fisherEnrichment(GeneSet.190, GeneSet.list)
candidates.enrich <- candidates.enrich %>% dplyr::filter(pvalue < 0.01)

write.xlsx(candidates.enrich, "~/work/ToxoplasmaGondii/KO_candidates/files/hypotheticals_candidates_geneset_enrichment.xlsx")
write.xlsx(Candidates.hyp, "~/work/ToxoplasmaGondii/KO_candidates/files/hypotheticals_candidates.xlsx")

### Quantify Tachy- and Brady-ness

## GeneSets corresponding to Tachy and Brady
Tachy.over.Brady <- c("Tachy.over.Brady.in.vivo", "Tachy.over.Brady.in.vitro") 
Brady.over.Tachy <- c("Brady.in.vitro.over.Tachy", "Brady.in.vivo.over.Tachy")


toxo_tab <- read.xlsx('../toxo_table_batch_corrected_logCPM_expression_edgeR_DEGs_ap2_targs_phenotype_correlations_smeV2.xlsx')
## Gene Sets
GeneSet.190 <- read.xlsx('../GSEA sytenic list 3.28.19.xlsx')
GeneSet.190 <- prep_geneset.190(GeneSet.190)
GeneSet.190.tachy <- GeneSet.190 %>% dplyr::filter(GeneSet %in% Tachy.over.Brady) %>% dplyr::select(-total) %>%
  spread(key = GeneSet, value = genes) %>% mutate(Tachy = list(Reduce(union, unlist(.)))) %>% 
  dplyr::select(Tachy) %>% gather(key = GeneSet, value = genes, Tachy) %>% mutate(total = length(unlist(genes)))
GeneSet.190.brady <- GeneSet.190 %>% dplyr::filter(GeneSet %in% Brady.over.Tachy) %>% dplyr::select(-total) %>%
  spread(key = GeneSet, value = genes) %>% mutate(Brady = list(Reduce(union, unlist(.)))) %>% 
  dplyr::select(Brady) %>% gather(key = GeneSet, value = genes, Brady) %>% mutate(total = length(unlist(genes)))

#### Get FC and Qvalue in tidy format
toxo.fc.qval <- get_fc_qval_tabs(toxo_tab)
toxo.fc.qval.sig <- toxo.fc.qval %>% dplyr::filter(abs(fc) > log2(2) & qval < 0.05)
toxo.fc.qval.up.reg <- toxo.fc.qval.sig %>% dplyr::filter(fc > 0)
toxo.fc.qval.down.reg <- toxo.fc.qval.sig %>% dplyr::filter(fc < 0)



GeneSet.list.up <- toxo.fc.qval.up.reg %>% group_by(Contrast) %>% 
  summarise(genes = list(as.character(GeneName)), total = n())

GeneSet.list.down <- toxo.fc.qval.down.reg %>% group_by(Contrast) %>% 
  summarise(genes = list(as.character(GeneName)), total = n())

enrich.up.tachy <- getEnrichment(GeneSet.190.tachy, GeneSet.list.up)
enrich.up.tachy$tachy.score <- -log(enrich.up.tachy$pvalue)
enrich.down.tachy <- getEnrichment(GeneSet.190.tachy, GeneSet.list.down)
enrich.down.tachy$tachy.score
enrich.up.tachy$direction <- 'upregulated'
enrich.down.tachy$direction <- 'downregulated'

enrich.up.brady <- getEnrichment(GeneSet.190.brady, GeneSet.list.up)
enrich.down.brady <- getEnrichment(GeneSet.190.brady, GeneSet.list.down)
enrich.up.brady$direction <- 'upregulated'
enrich.down.brady$direction <- 'downregulated'

enrich.tachy <- rbind(enrich.up.tachy, enrich.down.tachy)
enrich.brady <- rbind(enrich.up.brady, enrich.down.brady)
enrich <- rbind(enrich.tachy, enrich.brady)

XX1 <- enrich %>% 
  filter(Category == 'extra.vs.extra' & Passage1 == 'P11' & 
           direction == 'upregulated' & GeneSet == 'Tachy') 
XX2 <- enrich %>% 
  filter(Category == 'extra.vs.extra' & Passage1 == 'P11' & 
           direction == 'downregulated' & GeneSet == 'Brady')
tachy.score = -log(XX1$pvalue) - log(XX2$pvalue)

XX1 <- enrich %>% 
  filter(Category == 'extra.vs.extra' & Passage1 == 'P11' & 
           direction == 'upregulated' & GeneSet == 'Brady') 
XX2 <- enrich %>% 
  filter(Category == 'extra.vs.extra' & Passage1 == 'P11' & 
           direction == 'downregulated' & GeneSet == 'Tachy')
brady.score = -log(XX1$pvalue) - log(XX2$pvalue)

data.frame(Category = 'extra.vs.extra', passage = XX1$Passage2, tachy.score = tachy.score, brady.score = brady.score)

XX1 <- enrich %>% filter(Category == 'extra.vs.intra' & direction == 'upregulated' & GeneSet == 'Brady')
XX2 <- enrich %>% filter(Category == 'extra.vs.intra' & direction == 'downregulated' & GeneSet == 'Tachy')
-log(XX1$pvalue) -log(XX2$pvalue)


## PCA on gene sets
geneSets <- c("Brady.in.vitro.over.Tachy", "Brady.in.vivo.over.Tachy", "Tachy.over.Brady.in.vitro",       
              "Tachy.over.Brady.in.vivo", "Tachy.Ubiquitome.", "TissueCysts/Tachy.DownRegulated.", 
              "TissueCysts/Tachy.UPRegulated.")


FourWay.tachy <- FourWay %>% dplyr::filter(GeneSet %in% geneSets[4]) %>% 
  dplyr::select(GeneName,Condition1, Condition2,contains('expr'))
FourWay.tachy.cond1 <- FourWay.tachy %>% dplyr::select(GeneName,Condition1,log2.expr.cond1)
colnames(FourWay.tachy.cond1) <- c('GeneName', 'Condition', 'log_expr')
FourWay.tachy.cond2 <- FourWay.tachy %>% dplyr::select(GeneName,Condition2,log2.expr.cond2)
colnames(FourWay.tachy.cond2) <- c('GeneName', 'Condition', 'log_expr')
FourWay.tachy <- rbind(FourWay.tachy.cond1, FourWay.tachy.cond2) 
FourWay.tachy$Condition<- gsub('RH', 'P300', FourWay.tachy$Condition)
FourWay.tachy$Passage <- gsub('\\..*', '', FourWay.tachy$Condition)
FourWay.tachy$Condition <- gsub('P.*[[:digit:]]\\.', '', FourWay.tachy$Condition)
FourWay.tachy$Passage <- gsub('P300', 'RH', FourWay.tachy$Passage)
FourWay.tachy <- FourWay.tachy %>% distinct() %>% 
  spread(key = Condition, value = log_expr) 

FourWay.tachy <- FourWay.tachy %>% 
library(ggfortify)
df <- data.frame(t(FourWay.tachy[,2:ncol(FourWay.tachy)]))
colnames(df) <- FourWay.tachy$GeneName
rownames(df) = NULL
df$Condition <- colnames(FourWay.tachy)[2:ncol(FourWay.tachy)]
autoplot(prcomp(as.matrix(df[,1:(ncol(df)-1)])), data = df, colour = 'Condition')


spread(key = phenotype, value = correlation, fill = NA) %>% 
  dplyr::filter(fc > 1 & qval < 0.05 & quantile.cond2 >= 2) %>%
  dplyr::filter(plaque > 0.7 | reinv > 0.7 | survival > 0.7 | V8 > 0.7)

keep.columns <- c("GeneSet", "GeneName", "AP2", "Product.Description")
FourWay.pos.select <- FourWay.pos %>% dplyr::select(keep.columns, plaque:V8) %>% distinct()


#### Why the cell cylce arrest at G1 in extracellular parasites
## Make a FiveWay table (including all genes)
toxo_tab <- read.xlsx('../toxo_table_batch_corrected_logCPM_expression_edgeR_DEGs_ap2_targs_phenotype_correlations_smeV2.xlsx')
## Gene Sets
GeneSet.190 <- read.xlsx('../GSEA sytenic list 3.28.19.xlsx')
GeneSet.190 <- prep_geneset.190(GeneSet.190)

GeneSet.190 <- GeneSet.190 %>% unnest(genes) %>% dplyr::select(-c('total'))
colnames(GeneSet.190) <- c('GeneSet', 'GeneName')

## Checking APIAp2 targets
AP2.targets <- prep_ap2_targets(toxo_tab)
AP2.targets <- AP2.targets %>% unnest(genes) %>% dplyr::select(-c('total'))
colnames(AP2.targets) <- c('AP2', 'GeneName')

## Fold Change, Qvalue, Quantiles, Mean exprs
toxo.fc.qval.quant <- get_fc_qval_quant_tabs(toxo_tab)
FourWay <- make_four_way(toxo_tab, toxo.fc.qval.quant, GeneSet.190, AP2.targets)

cor.to.consider <- c("plaque", "reinv", "V8", "survival")
FourWay <- FourWay %>% 
  dplyr::select(GeneSet:Product.Description,c(paste("core.s.",cor.to.consider, sep = '')))
colnames(FourWay) <- gsub('core.s.', '', colnames(FourWay))


## Calculate significant overlpas between DEGs and genesets
toxo.fc.qval.quant.sig <- toxo.fc.qval.quant %>% dplyr::filter(abs(fc) > log2(1.5) & qval < 0.05)
toxo.fc.qval.quant.up.reg <- toxo.fc.qval.quant.sig %>% dplyr::filter(fc > 0)
toxo.fc.qval.quant.down.reg <- toxo.fc.qval.quant.sig %>% dplyr::filter(fc < 0)



GeneSet.list.up <- toxo.fc.qval.quant.up.reg %>% group_by(Contrast) %>% 
  summarise(genes = list(as.character(GeneName)), total = n())

GeneSet.list.down <- toxo.fc.qval.quant.down.reg %>% group_by(Contrast) %>% 
  summarise(genes = list(as.character(GeneName)), total = n())

GeneSet.190 <- read.xlsx('../GSEA sytenic list 3.28.19.xlsx')
GeneSet.190 <- prep_geneset.190(GeneSet.190)

enrich.up <- getEnrichment(GeneSet.190, GeneSet.list.up)
enrich.up <- enrich.up %>% dplyr::select(GeneSet, Contrast, pvalue)
colnames(enrich.up) <- c('GeneSet', 'Contrast', 'pval.overlap.up')
enrich.down <- getEnrichment(GeneSet.190, GeneSet.list.down)
enrich.down <- enrich.down %>% dplyr::select(GeneSet, Contrast, pvalue)
colnames(enrich.down) <- c('GeneSet', 'Contrast', 'pval.overlap.down')
enrich <- full_join(enrich.up, enrich.down, by = c('GeneSet', 'Contrast'))
FiveWay <- left_join(FourWay, enrich, by = c('GeneSet', 'Contrast'))

write.xlsx(FiveWay, "~/work/ToxoplasmaGondii/FiveWay_detailed_no_filter.xlsx")

## G phase and SM phase
All.Gs <- unique(FiveWay$GeneSet[grep('G1', FiveWay$GeneSet)])
early.G1 <- unique(FiveWay$GeneSet[grep('G1.R[0-1]',  FiveWay$GeneSet)])
middle.G1 <- unique(FiveWay$GeneSet[grep('G1.R[2-5]',  FiveWay$GeneSet)])
late.G1 <- unique(FiveWay$GeneSet[grep('G1.R[67]',  FiveWay$GeneSet)])
All.SMs <- unique(FiveWay$GeneSet[grep('SM', FiveWay$GeneSet)])

## Criteria:
## 1) Differentially up regulated in all passages between extra and intra
## 2) Is not in SM genes
## 3) Is in G genes
## 4) is highly expressed in at least one passage in extracellular
up.in.all <- toxo.fc.qval.quant %>% dplyr::filter(Category == 'extra.vs.intra') %>% group_by(GeneName) %>%
  summarise(sums = sum(fc > log(1.5) & qval < 0.05 & quantile.cond2 >= 2)) %>% dplyr::filter(sums == 4)

#trends.up <- toxo.fc.qval.quant %>% dplyr::filter(Category == 'extra.vs.extra' & Passage1 == 'P11') %>% 
#  dplyr::filter(fc > log(1.5) & qval < 0.05) %>% dplyr::select(GeneName) %>% distinct()

FiveWay.G.pos <- FiveWay %>% 
  dplyr::filter((GeneName %in% up.in.all$GeneName) )  %>% 
  dplyr::filter(Category == 'extra.vs.intra' & quantile.cond2 >= 2) %>% 
  dplyr::filter(GeneSet %in% All.Gs) %>%
  dplyr::filter(!(GeneSet %in% All.SMs))
keep.columns <- c("GeneSet", "GeneName", "AP2", "Product.Description")
FiveWay.G.pos.select <- FiveWay.G.pos %>% dplyr::select(keep.columns, plaque:survival) %>% distinct() 

FiveWay.G.pos.sum <- FiveWay.G.pos.select %>% group_by(GeneName) %>% 
  summarise(Product.Description = unique(Product.Description), 
            GeneSet = list(unique(as.character(GeneSet))), 
            AP2 = list(unique(as.character(AP2))),plaque = unique(plaque),
            survival = unique(survival), reinv = unique(reinv), V8 = unique(V8))
FiveWay.G.pos.sum <- FiveWay.G.pos.sum %>% rowwise() %>%
  mutate(is.in.EarlyG = ifelse(any(unlist(GeneSet) %in% early.G1), 'yes', 'no'),
         is.in.MiddleG = ifelse(any(unlist(GeneSet) %in% middle.G1), 'yes', 'no'),
         is.in.LateG = ifelse(any(unlist(GeneSet) %in% late.G1), 'yes', 'no')) %>%
  dplyr::filter(!is.na(Product.Description) & Product.Description != "hypothetical protein")


write.xlsx(FiveWay.G.pos.sum, "~/work/ToxoplasmaGondii/G1_candidates/files/G1_up_regulated_in_extra_summary.xlsx")

## Criteria:
## 1) Differentially down regulated in all passages between extra and intra
## 2) Is loweley expressed accross all extra cellular passage.
## 3) Is not in SM genes
## 4) In in G genes
down.in.all <- toxo.fc.qval.quant %>% dplyr::filter(Category == 'extra.vs.intra') %>% group_by(GeneName) %>%
  summarise(sums = sum(fc < -log(1.5) & qval < 0.05)) %>% dplyr::filter(sums == 4)

#trends.down <- toxo.fc.qval.quant %>% dplyr::filter(Category == 'extra.vs.extra' & Passage1 == 'P11') %>% 
#  dplyr::filter(fc < -log(1.2) & qval < 0.1) %>% dplyr::select(GeneName) %>% distinct()

FiveWay.G.neg <- FiveWay %>% dplyr::filter(Category == 'extra.vs.intra' & GeneName %in% down.in.all$GeneName) %>% 
  dplyr::filter(GeneSet %in% All.Gs & quantile.cond1 >= 2) %>%
  dplyr::filter(!(GeneSet %in% All.SMs))
keep.columns <- c("GeneSet", "GeneName", "AP2", "Product.Description")
FiveWay.G.neg.select <- FiveWay.G.neg %>% dplyr::select(keep.columns, plaque:survival) %>% distinct() 

FiveWay.G.neg.sum <- FiveWay.G.neg.select %>% group_by(GeneName) %>% 
  summarise(Product.Description = unique(Product.Description), 
            GeneSet = list(unique(as.character(GeneSet))), 
            AP2 = list(unique(as.character(AP2))),plaque = unique(plaque),
            survival = unique(survival), reinv = unique(reinv), V8 = unique(V8))
FiveWay.G.neg.sum <- FiveWay.G.neg.sum %>% rowwise() %>%
  mutate(is.in.EarlyG = ifelse(any(unlist(GeneSet) %in% early.G1), 'yes', 'no'),
         is.in.MiddleG = ifelse(any(unlist(GeneSet) %in% middle.G1), 'yes', 'no'),
         is.in.LateG = ifelse(any(unlist(GeneSet) %in% late.G1), 'yes', 'no')) %>%
  dplyr::filter(!is.na(Product.Description) & Product.Description != "hypothetical protein")


write.xlsx(FiveWay.G.neg.sum, "~/work/ToxoplasmaGondii/G1_candidates/files/G1_down_regulated_in_extra_summary.xlsx")


### Trends
## plot the expressions
tc.expr.extra <- get_extra_time_course(toxo_tab)
expr.fits <- readRDS(file = "~/work/ToxoplasmaGondii/ToxoR/output/sme_fits_exprs.RData")

XX <- FiveWay.G.pos.sum
num.plots <- ceiling(nrow(XX)/9)
for(j in 1:num.plots){
  out.pic <- paste("~/work/ToxoplasmaGondii/G1_candidates/plots/trend_G1_up_reg_", j, '.pdf', sep = '')
  pdf(out.pic, width=8, height=8)
  par(mfrow = c(3,3))
  for(i in ((j-1) * 9 + 1):min((j * 9), nrow(XX)) ){
    my.gene <- XX[i,c(1,2,4)]
    fit.gene <- expr.fits[[which(unique(tc.expr.extra$variable) == as.character(my.gene[1]))]]
    plot.sme(fit.gene, as.character(my.gene[1]),as.character(my.gene[2]))
  }
  dev.off()
}

XX <- FiveWay.G.neg.sum
num.plots <- ceiling(nrow(XX)/9)
for(j in 1:num.plots){
  out.pic <- paste("~/work/ToxoplasmaGondii/G1_candidates/plots/trend_G1_down_reg_", j, '.pdf', sep = '')
  pdf(out.pic, width=8, height=8)
  par(mfrow = c(3,3))
  for(i in ((j-1) * 9 + 1):min((j * 9), nrow(XX)) ){
    my.gene <- XX[i,c(1,2,4)]
    fit.gene <- expr.fits[[which(unique(tc.expr.extra$variable) == as.character(my.gene[1]))]]
    plot.sme(fit.gene, as.character(my.gene[1]),as.character(my.gene[2]))
  }
  dev.off()
}




###########
## Older stuff
###########
enrich.up.G <- enrich.up %>% 
  dplyr::filter(GeneSet %in% All.Gs & Category == 'extra.vs.intra' & GeneSet != 'All.G1.')
enrich.down.G <- enrich.down %>% 
  dplyr::filter(GeneSet %in% All.Gs & Category == 'extra.vs.intra' & GeneSet != 'All.G1.')

enrich.up.G$direction = 'Upregulated'
enrich.down.G$direction = 'Downregulated'

enrich.up.down.G <- rbind(enrich.up.G, enrich.down.G)


p1 <- ggplot(enrich.up.down.G, aes(x = GeneSet, y = Contrast)) + 
  geom_point(aes(color = pvalue, size = -log(pvalue))) +
  theme_bw(base_size = 12) +
  scale_colour_gradient(limits=c(0, 0.01), low="red") +
  ylab(NULL) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + facet_grid(Category ~ ., scales='free') + 
  theme(strip.background = element_rect(colour="black", fill="white", 
                                        size=1, linetype="solid")) +
  facet_grid(direction~.)

plot(p1)

isDownOverlapsLateNotMiddle
XX <- ThreeWay %>% 
  dplyr::filter(GeneSet %in% All.Gs & Category == 'extra.vs.intra' & GeneSet != 'All.G1.' & 
                  !is.na(Product.Description) & Product.Description != 'hypothetical protein')

isDownOverlapsEarly <- XX %>% dplyr::filter(GeneSet %in% early.G1 & fc < 0) %>% 
  dplyr::select(GeneName) %>% distinct()
isDownNoOverlapsMiddle <- XX %>% dplyr::filter(!(GeneSet %in% middle.G1) & fc < 0) %>% 
  dplyr::select(GeneName) %>% distinct()
isDownOverlapsLate <- XX %>% dplyr::filter((GeneSet %in% late.G1) & fc < 0) %>% 
  dplyr::select(GeneName) %>% distinct()

inOverlap <- intersect(isDownNoOverlapsMiddle$GeneName,isDownOverlapsLate$GeneName)
YY <- XX %>% dplyr::filter(GeneName %in% inOverlap)


isUpNoOverlapsEarly <- XX %>% dplyr::filter(!(GeneSet %in% early.G1) & fc > 0) %>% 
  dplyr::select(GeneName) %>% distinct()
isUpOverlapsMiddle <- XX %>% dplyr::filter((GeneSet %in% middle.G1) & fc > 0) %>% 
  dplyr::select(GeneName) %>% distinct()
isUpNoOverlapsLate <- XX %>% dplyr::filter(!(GeneSet %in% late.G1) & fc > 0) %>% 
  dplyr::select(GeneName) %>% distinct()

inOverlapUp <- intersect(isUpNoOverlapsEarly$GeneName,
                         intersect(isUpOverlapsMiddle$GeneName,isUpNoOverlapsLate$GeneName))

ZZ <- XX %>% dplyr::filter(GeneName %in% inOverlapUp)

transmute(GeneName = GeneName, 
          isDownInEarly = ifelse((GeneSet %in% early.G1) & fc < 0, T, F))
isDownInEarly <- unique(isDownInEarly$GeneName[isDownInEarly$isDownInEarly])

isDownInMiddle <- XX %>% transmute(GeneName = GeneName, 
                                   isDownInMiddle = ifelse((GeneSet %in% middle.G1) & fc < 0, T, F))
isDownInMiddle <- unique(isDownInMiddle$GeneName[isDownInMiddle$isDownInMiddle])
isDownInLate <- XX %>% transmute(GeneName = GeneName, 
                                 isDownInLate = ifelse((GeneSet %in% late.G1) & fc < 0, T, F))
isDownInLate <- unique(isDownInLate$GeneName[isDownInLate$isDownInLate])

YY <- XX %>% dplyr::filter((GeneName %in% isDownInLate) &  !(GeneName %in% isDownInMiddle))


#dplyr::filter((GeneSet %in% middle.G1) & fc >= 0)

## Make a FourWay table (including all genes)
toxo_tab <- read.xlsx('../toxo_table_batch_corrected_logCPM_expression_edgeR_DEGs_ap2_targs_phenotype_correlations_smeV2.xlsx')
## Gene Sets
GeneSet.190 <- read.xlsx('../GSEA sytenic list 3.28.19.xlsx')
GeneSet.190 <- prep_geneset.190(GeneSet.190)
GeneSet.190 <- GeneSet.190 %>% unnest(genes) %>% dplyr::select(-c('total'))
colnames(GeneSet.190) <- c('GeneSet', 'GeneName')

## Checking APIAp2 targets
AP2.targets <- prep_ap2_targets(toxo_tab)
AP2.targets <- AP2.targets %>% unnest(genes) %>% dplyr::select(-c('total'))
colnames(AP2.targets) <- c('AP2', 'GeneName')

## Fold Change, Qvalue, Quantiles, Mean exprs
toxo.fc.qval.quant <- get_fc_qval_quant_tabs(toxo_tab)
cor.to.consider <- c("plaque", "reinv", "V8", "survival")
FourWay <- make_four_way(toxo_tab, toxo.fc.qval.quant, GeneSet.190, AP2.targets)
FourWay <- FourWay %>% 
  dplyr::select(GeneSet:Product.Description,c(paste("core.s.",cor.to.consider, sep = '')))
colnames(FourWay) <- gsub('core.s.', '', colnames(FourWay))

FourWay <- FourWay %>% gather(key = phenotype, value = correlation, plaque:survival)

