library(fda)
library(openxlsx)
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(tidyverse)
library(sme)


fitSingleSme <-function(my.tc, v){
  tryCatch(
    expr = {
      fit <- sme(my.tc[my.tc$variable==v,c("y","tme","ind")], criteria = 'AIC')
      return(fit)
    },
    error = function(v){
      message(paste('error:', v))
    },
    warning = function(w){
      message(paste('error:', v))
    },
    finally = {
      fit <- sme(my.tc[my.tc$variable==v,c("y","tme","ind")])
      return(fit)
    }
  )    
}

## spline the fitted values to get the means
splineSmeFits <- function(fits, variables){
  mus <- lapply(fits, function(x) spline(x = as.numeric(colnames(x$coefficients)),
                                         y = x$coefficients[1,], 
                                         n = as.numeric(colnames(x$coefficients))[ncol(x$coefficients)] - 
                                           as.numeric(colnames(x$coefficients))[1] + 1, 
                                         method = "natural")) 
  mus.y <- unlist(lapply(mus, `[[`, 2))
  mus.x <- unlist(lapply(mus, `[[`, 1))
  lens  <- unlist(lapply(lapply(mus, `[[`, 1), length))
  
  mus <- data.frame(variable = rep(variables, times = lens),
                    tme = mus.x,
                    y = mus.y)
  
  return(mus)
  
}

plot.sme <-function(fit, v){
  mu <- spline(x = as.numeric(colnames(fit$coefficients)), 
               y = fit$coefficients[1,], n = 500, 
               method = "natural")
  fs <- lapply(2:nrow(fit$coefficients), function(i) {
    spline(x = as.numeric(colnames(fit$coefficients)), 
           y = fit$coefficients[1,] + fit$coefficients[i, ], 
           method = "natural", 
           n = 500)})
  
  
  ylim <- range(fit$data$y, mu$y, sapply(fs, "[[", "y"))
  xlim <- range(as.numeric(colnames(fit$coefficients)))
  mu.variance <- diag(vcov(fit))
  upper.band <- spline(x = as.numeric(colnames(fit$coefficients)), 
                       y = fit$coefficients[1, ] + 1.96 * sqrt(mu.variance), 
                       method = "natural", n = 500)
  lower.band <- spline(x = as.numeric(colnames(fit$coefficients)), 
                       y = fit$coefficients[1, ] - 1.96 * sqrt(mu.variance), 
                       method = "natural", n = 500)
  ylim <- range(ylim, upper.band$y, lower.band$y)
  
  plot(x = fit$data$tme, y = fit$data$y, ylim = ylim, xlim = xlim, 
       xlab = 'Passage', ylab = 'y', col = 'black', cex = 0.8, main = v)
  
  for (i in 1:length(fs)) {
    lines(fs[[i]], lty = "dashed", col = 'black', lwd = 0.8)
  }
  
  lines(mu, lwd = 2, col = 'red')
  
  col.meanCurve.rgb <- col2rgb('red')
  polygon(x = c(upper.band$x, rev(lower.band$x)), 
          y = c(upper.band$y,rev(lower.band$y)), 
          col = rgb(col.meanCurve.rgb[1], 
                    col.meanCurve.rgb[2], col.meanCurve.rgb[3], alpha = 125, 
                    maxColorValue = 255), border = NA)
}

## Read in the new Phenotypes with separate replicate info
plaque_sep_rep <- read.xlsx('../Plaque_sep_rep.xlsx')
plaque_sep_rep <- plaque_sep_rep %>% gather(key = Replicate, value = y, -Passage)
plaque_sep_rep$Passage <- as.numeric(gsub('B2 P', '', plaque_sep_rep$Passage))
plaque_sep_rep$Replicate <- as.numeric(gsub('Rep', '', plaque_sep_rep$Replicate))
plaque_sep_rep$variable <- 'plaque'
plaque_sep_rep <- plaque_sep_rep[,c('y', 'Passage', 'Replicate', 'variable')]
colnames(plaque_sep_rep) <- c('y', 'tme', 'ind', 'variable')

reinvasion_sep_rep <- read.xlsx('../Reinvasion_sep_rep.xlsx')
reinvasion_sep_rep <- reinvasion_sep_rep %>% gather(key = Replicate, value = y, -Passage) 
reinvasion_sep_rep$Passage <- as.numeric(gsub('B2 P', '', reinvasion_sep_rep$Passage))
reinvasion_sep_rep$Replicate <- as.numeric(gsub('Rep', '', reinvasion_sep_rep$Replicate))
reinvasion_sep_rep$variable <- 'reinv'
reinvasion_sep_rep <- reinvasion_sep_rep[,c('y', 'Passage', 'Replicate', 'variable')]
colnames(reinvasion_sep_rep) <- c('y', 'tme', 'ind', 'variable')

#replication_sep_rep <- read.xlsx('../Replication_sep_rep.xlsx')
replication_sep_rep <- read.xlsx('../Replication_sep_repV2.xlsx')
replication_sep_rep$Passage <- as.numeric(gsub('B2 P', '', replication_sep_rep$Passage))
replication_sep_rep <- replication_sep_rep %>% gather(key = variable, value = y, -Passage,-Replicate)
replication_sep_rep <- replication_sep_rep %>% arrange(variable, Replicate, Passage)
replication_sep_rep <- replication_sep_rep[,c('y', 'Passage', 'Replicate', 'variable')]
colnames(replication_sep_rep) <- c('y', 'tme', 'ind', 'variable')

survival_sep_rep <- read.xlsx('../survival_sep_rep.xlsx')
survival_sep_rep <- survival_sep_rep %>% gather(key = Replicate, value = y, -Passage)
survival_sep_rep$Passage <- as.numeric(gsub('B2 P', '', survival_sep_rep$Passage))
survival_sep_rep$Replicate <- as.numeric(gsub('Rep', '', survival_sep_rep$Replicate))
survival_sep_rep$variable <- 'survival'
survival_sep_rep <- survival_sep_rep[,c('y', 'Passage', 'Replicate', 'variable')]
colnames(survival_sep_rep) <- c('y', 'tme', 'ind', 'variable')

phenotypes_sep_rep <- rbind(plaque_sep_rep, reinvasion_sep_rep, replication_sep_rep, survival_sep_rep)
tc <- phenotypes_sep_rep %>% na.omit()

## Fit a Smoothing Spline Mixed Effect Model.
pheno.fits <- lapply(unique(tc$variable),function(v) 
  fitSingleSme(tc, v))
#saveRDS(object = pheno.fits ,file = "~/work/ToxoplasmaGondii/ToxoR/output/sme_fits_pheno.RData")
saveRDS(object = pheno.fits ,file = "~/work/ToxoplasmaGondii/ToxoR/output/sme_fits_phenoV2.RData")
#pheno.fits <- readRDS(file = "~/work/ToxoplasmaGondii/ToxoR/output/sme_fits_phenoV2.RData")

## Read in expression values for each replicate
expr_sep_rep <- read.xlsx('../toxo_table_batch_corrected_logCPM_expression_edgeR_DEGs_ap2_targs.xlsx')

## Exclude P7 and RH and look at extra only
expr_sep_rep.extra <- expr_sep_rep %>% 
  dplyr::select(matches('extra.*rep')) %>% dplyr::select(-contains('RH')) %>%
  dplyr::select(-contains('P7'))
colnames(expr_sep_rep.extra) <- gsub('.extra', '', colnames(expr_sep_rep.extra))
expr_sep_rep.extra$GeneName <- expr_sep_rep$GeneName

expr_sep_rep.extra <- expr_sep_rep.extra %>% gather(key=PassageRep, value = Expr, -GeneName)
expr_sep_rep.extra <- expr_sep_rep.extra %>% 
  mutate(Passage = gsub("\\.rep.*", '', PassageRep), Rep = gsub(".*rep\\.", '', PassageRep)) %>%
  dplyr::select(-c('PassageRep'))

Passage <- as.numeric(gsub('P', '', expr_sep_rep.extra$Passage))
expr_sep_rep.extra$Passage <- Passage
expr_sep_rep.extra$Rep <- as.numeric(expr_sep_rep.extra$Rep)

expr_sep_rep.extra <- expr_sep_rep.extra %>% arrange(GeneName, Rep, Passage)
expr_sep_rep.extra <- expr_sep_rep.extra[,c('Expr', 'Passage', 'Rep', 'GeneName')]
colnames(expr_sep_rep.extra) <- c('y', 'tme', 'ind', 'variable')

## Fit a smoothing spline mixed effect model to each gene
tc.expr <- expr_sep_rep.extra %>% na.omit()
expr.fits <- lapply(unique(tc.expr$variable),function(v) 
  fitSingleSme(tc.expr, v))
saveRDS(object = expr.fits,file = "~/work/ToxoplasmaGondii/ToxoR/output/sme_fits_exprs.RData")
#expr.fits <- readRDS(file = "~/work/ToxoplasmaGondii/ToxoR/output/sme_fits_exprs.RData")

# expr.fits <- list()
# counter = 1
# for(v in unique(tc.expr$variable)){
#   expr.fits[[counter]] <- sme(tc.expr[tc.expr$variable==v,c("y","tme","ind")], 
#                               criteria = 'AIC')
#   counter = counter + 1
#   print(counter)
# }
# 

pheno.mus <- splineSmeFits(pheno.fits, unique(tc$variable))
colnames(pheno.mus) <- c('phenotype', 'tme', 'y') 
expr.mus <- splineSmeFits(expr.fits, unique(tc.expr$variable))
colnames(expr.mus) <- c('GeneName', 'tme', 'y') 

## Calculate Spearman Rank Correlation
mus <- inner_join(expr.mus, pheno.mus, by = 'tme')

cor.mat <- mus %>% group_by(GeneName, phenotype) %>% 
  summarise(cor.p = cor(y.x, y.y, method = "pearson"),
            cor.s = cor(y.x, y.y, method = "spearman"))


cor.mat.p <- cor.mat %>% dplyr::select(-c('cor.s')) %>% spread(key = phenotype, value = cor.p)
cor.mat.s <- cor.mat %>% dplyr::select(-c('cor.p')) %>% spread(key = phenotype, value = cor.s)
colnames(cor.mat.p) <- c('GeneName', paste('core.p.',colnames(cor.mat.p)[2:ncol(cor.mat.p)], sep = ''))
colnames(cor.mat.s) <- c('GeneName', paste('core.s.',colnames(cor.mat.s)[2:ncol(cor.mat.p)], sep = ''))

cor.mat <- left_join(cor.mat.p, cor.mat.s, by = 'GeneName')

sep_rep <- left_join(expr_sep_rep, cor.mat, by = 'GeneName')
write.xlsx(sep_rep, '../toxo_table_batch_corrected_logCPM_expression_edgeR_DEGs_ap2_targs_phenotype_correlations_smeV2.xlsx')


### Make a 4 Way matirx with correlations, AP2 targets, gene sets and DEGs.
ThreeWay <- read.xlsx('../ThreeWayDEGoverlaps.xlsx')
FourWay <- left_join(ThreeWay, cor.mat, by = 'GeneName')
#write.xlsx(FourWay, '../FourWayDEGoverlaps_sme.xlsx')
write.xlsx(FourWay, '../FourWayDEGoverlaps_smeV2.xlsx')


AP2.extra <- sep_rep %>% dplyr::filter(!is.na(Product.Description)) %>% 
  dplyr::filter(Product.Description != 'hypothetical protein') %>% 
  dplyr::filter(str_detect(Product.Description, 'AP2')) %>%
  dplyr::select(c('GeneName', 'Product.Description') ) %>% distinct()

my.AP2s <- c('AP2XI-5', 'AP2X-5','AP2IX-9', 'AP2IX-4', 'AP2IV-4', 'AP2IV-3', 'AP2XI-4', 'AP2Ib-1', 'AP2XII-6')
my.AP2s <- data.frame(AP2 = my.AP2s, 
           GeneName = AP2.extra$GeneName[unlist(lapply(my.AP2s, function(x) grep(x, AP2.extra$Product.Description)))])

pdf("~/work/ToxoplasmaGondii/AP2s_sme_fits.pdf", width=8, height=8)
par(mfrow=c(3,3))
for(i in 1:nrow(my.AP2s)){
  fit.gene <- expr.fits[[which(unique(tc.expr$variable) == my.AP2s[i,2])]]
  plot.sme(fit.gene, my.AP2s[i,1])
  
}
dev.off()

pheno.fits <- lapply(unique(tc$variable),function(v) 
  fitSingleSme(tc, v))

pdf("~/work/ToxoplasmaGondii/phenotypes_sme_fitsV2.pdf", width=8, height=8)
par(mfrow=c(3,3))
phenotypes <- unique(tc$variable)
Vs <- grep('V', phenotypes)
others <- c(1:length(phenotypes))[-Vs]
others <- c(1,2,9,3) ## send appv after survival
phenotypes <- phenotypes[c(others, Vs[sort(as.numeric(gsub('V', '', phenotypes[Vs])), index.return = T)$ix])]
for(i in 1:length(phenotypes)){
  j <- which(unique(tc$variable) == phenotypes[i])
  fit.ph <- pheno.fits[[j]]
  plot.sme(fit.ph, phenotypes[i])
  
}
dev.off()
