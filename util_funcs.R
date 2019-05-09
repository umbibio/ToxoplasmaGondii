library(dplyr)
library(ggplot2)
library(tidyverse)


get_fc_qval_tabs <- function(toxo_tab){
  #### Get FC and Qvalue in tidy format
  toxo.fc   <- toxo_tab %>% dplyr::select(c(1, contains('fc')))
  colnames(toxo.fc) <- gsub("\\._log_fc", "", colnames(toxo.fc))
  toxo.qval <- toxo_tab %>% dplyr::select(c(1, contains('qval')))
  colnames(toxo.qval) <- gsub("\\._qval", "", colnames(toxo.qval))
  
  toxo.fc <- toxo.fc  %>% gather(key = Contrast, value = fc, -GeneName)
  toxo.qval <- toxo.qval  %>% gather(key = Contrast, value = qval, -GeneName)
  
  toxo.fc$fc <- as.numeric(toxo.fc$fc)
  toxo.qval$qval <- as.numeric(toxo.qval$qval)
  
  toxo.fc.qval <- inner_join(toxo.fc, toxo.qval, by = c('GeneName', 'Contrast'))
  
  Passages <- strsplit(as.character(toxo.fc.qval$Contrast), split='\\.vs\\.')
  Passage2 <- unlist(lapply(strsplit(unlist(lapply(Passages, `[[`, 1)), split='\\.'), `[[`,1))
  Passage1 <- unlist(lapply(strsplit(unlist(lapply(Passages, `[[`, 2)), split='\\.'), `[[`,1))
  
  Category <- gsub("RH.intra", "RH",
                   gsub("RH.extra", "RH", 
                        gsub("^P.*[[:digit:]]\\.", "",
                             gsub("\\.P.*[[:digit:]]\\.", "\\.", 
                                  as.character(toxo.fc$Contrast)))))
  
  sorted.ind1 <- sort(as.numeric(gsub("P", "", gsub('RH', 'P300', unique(Passage1)))), index.return = T)$ix
  sorted.ind2 <- sort(as.numeric(gsub("P", "", gsub('RH', 'P300', unique(Passage2)))), index.return = T)$ix
  Passage1 <- factor((Passage1), levels = unique(Passage1)[sorted.ind1])
  Passage2 <- factor((Passage2), levels = unique(Passage2)[sorted.ind2])
  Category <- factor(Category, levels = sort(unique(Category)))
  
  toxo.fc.qval$Category <- Category
  toxo.fc.qval$Passage1 <- Passage1
  toxo.fc.qval$Passage2 <- Passage2
  toxo.fc.qval$Contrast <- factor(toxo.fc.qval$Contrast, levels = unique(toxo.fc.qval$Contrast))
  
  toxo.fc.qval <- toxo.fc.qval %>% na.omit()
  
  return(toxo.fc.qval)
}

get_cor.mat <- function(toxo_tab){
  #### Get FC and Qvalue in tidy format
  toxo.core   <- toxo_tab %>% dplyr::select(c(1, contains('core')))
  return(toxo.core)
}

get_fc_qval_quant_tabs <- function(toxo_tab){
  toxo.fc   <- toxo_tab %>% dplyr::select(c(1, contains('fc')))
  colnames(toxo.fc) <- gsub("\\._log_fc", "", colnames(toxo.fc))
  toxo.qval <- toxo_tab %>% dplyr::select(c(1, contains('qval')))
  colnames(toxo.qval) <- gsub("\\._qval", "", colnames(toxo.qval))
  
  toxo.fc <- toxo.fc  %>% gather(key = Contrast, value = fc, -GeneName)
  toxo.qval <- toxo.qval  %>% gather(key = Contrast, value = qval, -GeneName)
  
  toxo.fc$fc <- as.numeric(toxo.fc$fc)
  toxo.qval$qval <- as.numeric(toxo.qval$qval)
  
  toxo.fc.qval <- inner_join(toxo.fc, toxo.qval, by = c('GeneName', 'Contrast'))
  
  Passages <- strsplit(as.character(toxo.fc.qval$Contrast), split='\\.vs\\.')
  Passage2 <- unlist(lapply(strsplit(unlist(lapply(Passages, `[[`, 1)), split='\\.'), `[[`,1))
  Passage1 <- unlist(lapply(strsplit(unlist(lapply(Passages, `[[`, 2)), split='\\.'), `[[`,1))
  Condition2 <- unlist(lapply(Passages, `[[`, 1))
  Condition1 <- unlist(lapply(Passages, `[[`, 2))
  
  
  Category <- gsub("RH.intra", "RH",
                   gsub("RH.extra", "RH", 
                        gsub("^P.*[[:digit:]]\\.", "",
                             gsub("\\.P.*[[:digit:]]\\.", "\\.", 
                                  as.character(toxo.fc$Contrast)))))
  
  sorted.ind1 <- sort(as.numeric(gsub("P", "", gsub('RH', 'P300', unique(Passage1)))), index.return = T)$ix
  sorted.ind2 <- sort(as.numeric(gsub("P", "", gsub('RH', 'P300', unique(Passage2)))), index.return = T)$ix
  Passage1 <- factor((Passage1), levels = unique(Passage1)[sorted.ind1])
  Passage2 <- factor((Passage2), levels = unique(Passage2)[sorted.ind2])
  
  #Condition1 <- factor((Condition1), levels = unique(Condition1))
  #Condition2 <- factor((Condition2), levels = unique(Condition1))
  
  Category <- factor(Category, levels = sort(unique(Category)))
  
  
  toxo.fc.qval$Category <- Category
  toxo.fc.qval$Passage1 <- Passage1
  toxo.fc.qval$Passage2 <- Passage2
  toxo.fc.qval$Condition1 <- Condition1
  toxo.fc.qval$Condition2 <- Condition2
  toxo.fc.qval$Contrast <- factor(toxo.fc.qval$Contrast, levels = unique(toxo.fc.qval$Contrast))
  
  
  toxo.fc.qval <- toxo.fc.qval %>% na.omit()
  
  toxo.exprs <- toxo_tab %>% dplyr::select(c(1, contains('_mean')))
  colnames(toxo.exprs) <- gsub("_mean", "", colnames(toxo.exprs))
  toxo.exprs <- toxo.exprs %>% gather(key = Condition, value = expr, -GeneName)
  #toxo.exprs$Condition <- factor(toxo.exprs$Condition, levels = unique(Condition1))
  toxo.exprs <- toxo.exprs %>% group_by(Condition) %>% 
    mutate(quantile = ifelse(expr < quantile(expr)[2], 1, 
                          ifelse(expr < quantile(expr)[3], 2, 
                                 ifelse(expr < quantile(expr)[4], 3,4))))
  
  toxo.sds <- toxo_tab %>% dplyr::select(c(1, contains('_sd')))
  colnames(toxo.sds) <- gsub("_sd", "", colnames(toxo.sds))
  toxo.sds <- toxo.sds %>% gather(key = Condition, value = sd, -GeneName)
  toxo.exprs.sd <- left_join(toxo.exprs, toxo.sds, by = c("GeneName", "Condition"))
  
  toxo.fc.qval.expr <- left_join(toxo.fc.qval, toxo.exprs.sd, by = c('GeneName', 'Condition1' = 'Condition'))
  toxo.fc.qval.expr <- left_join(toxo.fc.qval.expr, toxo.exprs.sd, by = c('GeneName', 'Condition2' = 'Condition'))
  colnames(toxo.fc.qval.expr) <- gsub('expr', 'log2.expr', 
                                      gsub('\\.y', '\\.cond2', 
                                           gsub('\\.x', '\\.cond1', colnames(toxo.fc.qval.expr))))
  
  return(toxo.fc.qval.expr)
}


get_mean_sd_tabs <- function(toxo_tab, cond = 'extra'){
  ## Get expressions for extra cellular
  toxo.tab.cond <- toxo_tab %>% 
    dplyr::select(contains(paste(cond, '_', sep = '')))
  colnames(toxo.tab.cond) <- gsub(paste(cond, '_', sep = ''), '', colnames(toxo.tab.cond))
  
  toxo.tab.cond.mean <- toxo.tab.cond %>% dplyr::select(contains('mean'))
  colnames(toxo.tab.cond.mean) <- gsub('\\.mean', '', colnames(toxo.tab.cond.mean))
  toxo.tab.cond.mean$GeneName <- toxo_tab$GeneName
  toxo.tab.cond.mean <- toxo.tab.cond.mean %>% gather(key=Passage, value= mean, -GeneName)
  
  toxo.tab.cond.sd <- toxo.tab.cond %>% dplyr::select(contains('sd'))
  colnames(toxo.tab.cond.sd) <- gsub('\\.sd', '', colnames(toxo.tab.cond.sd))
  toxo.tab.cond.sd$GeneName <- toxo_tab$GeneName
  toxo.tab.cond.sd <- toxo.tab.cond.sd %>% gather(key=Passage, value= sd, -GeneName)
  
  toxo.tab.cond.mean.sd <- left_join(toxo.tab.cond.mean, toxo.tab.cond.sd, by = c('GeneName', 'Passage'))
  
  toxo.tab.cond.mean.sd$cond <- cond
  
  return(toxo.tab.cond.mean.sd)
  
}

prep_geneset.190 <- function(GeneSet.190){
  ## Make less lenghty names
  colnames(GeneSet.190) <- 
    gsub("development", "devel", 
         gsub('tgo.*[[:digit:]]\\.', '', 
              gsub('Tachyzoite', 'Tachy', 
                   gsub('Bradyzoite', 'Brady', 
                        gsub('tachyzoite', 'Tachy',  
                             gsub('bradyzoite', 'Brady',
                                  gsub('\\,', '', 
                                       gsub('\\(.*', '', 
                                            gsub('_', '.', 
                                                 gsub('-', '.', colnames(GeneSet.190)))))))))))
  
  GeneSet.190 <- GeneSet.190 %>% gather(key = GeneSet, value = GeneName) %>% na.omit() %>% 
    group_by(GeneSet) %>% summarise(genes = list(as.character(GeneName)), total = n())
  return(GeneSet.190)
}

prep_ap2_targets <- function(toxo_tab){
  AP2.targets <- data.frame(GeneName = toxo_tab$GeneName, toxo_tab[,grep('AP2.*_genes', colnames(toxo_tab))])
  #AP2.targets <- AP2.targets %>% dplyr::select(-c('AP2X5_ikD_upregulated_genes')) %>%
  #colnames(AP2.targets) <- gsub('XI5', 'XI-5',gsub('\\.', '-',gsub("_downregulated_genes", '', 
  #                                                                 gsub('_target_genes', '',colnames(AP2.targets)))))
  
  AP2.targets <- AP2.targets %>% 
    dplyr::select(-c('AP2X5_ikD_upregulated_genes')) %>%
    dplyr::select(-c('AP2X5_ikD_downregulated_genes')) %>%
    dplyr::select(-c('AP2IX.9_target_genes')) %>%
    dplyr::select(-c("AP2IX9_target_genes_up")) ## they play no role
  colnames(AP2.targets) <- gsub("AP2IX9_down", "AP2IX9", 
                                gsub('\\.', '', gsub('_target_genes', '',colnames(AP2.targets))))
  
  
  AP2.targets <- AP2.targets %>% gather(key = GeneSet, value = targets, -GeneName)
  AP2.targets <- AP2.targets %>% group_by(GeneSet) %>% 
    summarise(genes = list(as.character(GeneName[!is.na(targets)])), total = sum(!is.na(targets)))
  
  AP2.targets$GeneSet <- factor(AP2.targets$GeneSet, levels = unique(AP2.targets$GeneSet))
  
  return(AP2.targets)
  
}

make_three_way <- function(toxo_tab, toxo.fc.qval, GeneSet.190, AP2.targets){
  ### Perform a Three Way overlap between DEGs, AP2 targets, and GeneSets
  XX <- full_join(GeneSet.190, toxo.fc.qval, by = 'GeneName')
  
  ThreeWay <- full_join(XX, AP2.targets, by = 'GeneName')
  ThreeWay$Product.Description <- toxo_tab$Product.Description[match(ThreeWay$GeneName, 
                                                                          toxo_tab$GeneName)]
  ThreeWay <- ThreeWay[!is.na(ThreeWay$fc), ]
  return(ThreeWay)
}


make_four_way <- function(toxo_tab, toxo.fc.qval, GeneSet.190, AP2.targets){
  ThreeWay <- make_three_way(toxo_tab, toxo.fc.qval, GeneSet.190, AP2.targets)
  toxo.core <- get_cor.mat(toxo_tab)
  
  FourWay <- left_join(ThreeWay, toxo.core, by = 'GeneName')
  return(FourWay)
}

plot.sme <-function(fit, v, xlab = 'Passage'){
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
       xlab = xlab, ylab = 'y', col = 'black', cex = 0.8, main = v, xaxt = "n")
  axis(1, at=unique(fit$data$tme), labels = paste('P', unique(fit$data$tme), sep = ''))
  #text(unique(fit$data$tme), par("usr")[3] - 0.2, labels = paste('P', unique(fit$data$tme), sep = '')
  #     , srt = 90, pos = 1, xpd = TRUE,cex = 0.8)
  
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

get_extra_time_course <- function(toxo_tab){
  ## Exclude P7 and RH and look at extra only
  expr_sep_rep.extra <- toxo_tab %>% 
    dplyr::select(matches('extra.*rep')) %>% dplyr::select(-contains('RH')) %>%
    dplyr::select(-contains('P7'))
  colnames(expr_sep_rep.extra) <- gsub('.extra', '', colnames(expr_sep_rep.extra))
  expr_sep_rep.extra$GeneName <- toxo_tab$GeneName
  
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
  
  return(tc.expr)
}

get_pohenotype_time_course <- function(plaque_sep_rep, reinvasion_sep_rep, 
                                       replication_sep_rep, survival_sep_rep){
  plaque_sep_rep <- plaque_sep_rep %>% gather(key = Replicate, value = y, -Passage)
  plaque_sep_rep$Passage <- as.numeric(gsub('B2 P', '', plaque_sep_rep$Passage))
  plaque_sep_rep$Replicate <- as.numeric(gsub('Rep', '', plaque_sep_rep$Replicate))
  plaque_sep_rep$variable <- 'plaque'
  plaque_sep_rep <- plaque_sep_rep[,c('y', 'Passage', 'Replicate', 'variable')]
  colnames(plaque_sep_rep) <- c('y', 'tme', 'ind', 'variable')
  
  reinvasion_sep_rep <- reinvasion_sep_rep %>% gather(key = Replicate, value = y, -Passage) 
  reinvasion_sep_rep$Passage <- as.numeric(gsub('B2 P', '', reinvasion_sep_rep$Passage))
  reinvasion_sep_rep$Replicate <- as.numeric(gsub('Rep', '', reinvasion_sep_rep$Replicate))
  reinvasion_sep_rep$variable <- 'reinv'
  reinvasion_sep_rep <- reinvasion_sep_rep[,c('y', 'Passage', 'Replicate', 'variable')]
  colnames(reinvasion_sep_rep) <- c('y', 'tme', 'ind', 'variable')
  
  replication_sep_rep$Passage <- as.numeric(gsub('B2 P', '', replication_sep_rep$Passage))
  replication_sep_rep <- replication_sep_rep %>% gather(key = variable, value = y, -Passage,-Replicate)
  replication_sep_rep <- replication_sep_rep %>% arrange(variable, Replicate, Passage)
  replication_sep_rep <- replication_sep_rep[,c('y', 'Passage', 'Replicate', 'variable')]
  colnames(replication_sep_rep) <- c('y', 'tme', 'ind', 'variable')
  
  survival_sep_rep <- survival_sep_rep %>% gather(key = Replicate, value = y, -Passage)
  survival_sep_rep$Passage <- as.numeric(gsub('B2 P', '', survival_sep_rep$Passage))
  survival_sep_rep$Replicate <- as.numeric(gsub('Rep', '', survival_sep_rep$Replicate))
  survival_sep_rep$variable <- 'survival'
  survival_sep_rep <- survival_sep_rep[,c('y', 'Passage', 'Replicate', 'variable')]
  colnames(survival_sep_rep) <- c('y', 'tme', 'ind', 'variable')
  
  phenotypes_sep_rep <- rbind(plaque_sep_rep, reinvasion_sep_rep, replication_sep_rep, survival_sep_rep)
  tc <- phenotypes_sep_rep %>% na.omit()
 
  return(tc) 
}

fisherEnrichment <- function(GeneSet, GeneSet.list){
  ## This is a Hack to avoid nested for loops
  GeneSet$all <- 8637
  GeneSet.list$all <- 8637
  
  XX <- full_join(GeneSet, GeneSet.list, by = 'all')
  XX <- XX %>% rowwise() %>% mutate(overlap = length(intersect(c(genes.x), c(genes.y))),
                                    overlap.genes = list(intersect(c(genes.x), c(genes.y))))
  XX <- XX %>% rowwise() %>% 
    mutate(pvalue = fisher.test(matrix(c(overlap, total.x - overlap, total.y - overlap, 
                                         all - (total.x + total.y - overlap) ), byrow = T, ncol = 2, nrow = 2),
                                alternative = "greater")$p.value)
  
  return(XX)
}
  
getEnrichment <- function(GeneSet, GeneSet.list){
  ## This is a Hack to avoid nested for loops
  GeneSet$all <- 8637
  GeneSet.list$all <- 8637
  
  XX <- full_join(GeneSet, GeneSet.list, by = 'all')
  XX <- XX %>% rowwise() %>% mutate(overlap = length(intersect(c(genes.x), c(genes.y))),
                                    overlap.genes = list(intersect(c(genes.x), c(genes.y))))
  XX <- XX %>% rowwise() %>% 
    mutate(pvalue = fisher.test(matrix(c(overlap, total.x - overlap, total.y - overlap, 
                                         all - (total.x + total.y - overlap) ), byrow = T, ncol = 2, nrow = 2),
                                alternative = "greater")$p.value)
  
  contrasts <- unique(as.character(XX$Contrast))
  
  extra.vs.extra.ind <- grep('(P.*extra)(.*P.*extra)',contrasts)
  sorted.extra.vs.extra.ind <- extra.vs.extra.ind[
    sort(as.numeric(gsub('P', '', gsub('\\.extra.*', '',contrasts[extra.vs.extra.ind]))), index.return = T)$ix]
  
  intra.vs.intra.ind <- grep('(P.*intra)(.*P.*intra)',contrasts)
  sorted.intra.vs.intra.ind <- intra.vs.intra.ind[
    sort(as.numeric(gsub('P', '', gsub('\\.intra.*', '',contrasts[intra.vs.intra.ind]))), index.return = T)$ix]
  
  extra.vs.intra.ind <- grep('(P.*extra)(.*P.*intra)',contrasts)
  sorted.extra.vs.intra.ind <- extra.vs.intra.ind[
    sort(as.numeric(gsub('P', '', gsub('\\.extra.*', '',contrasts[extra.vs.intra.ind]))), index.return = T)$ix]
  
  
  RH.vs.RH.ind <- grep('(.*RH)(.*RH)',contrasts)
  
  RH.vs.intra.ind <- grep('(.*RH.intra)(.*intra)',contrasts)
  sorted.RH.vs.intra.ind <- RH.vs.intra.ind[
    sort(as.numeric(gsub('P', '', 
                         gsub('\\.intra.*', '',
                              gsub('RH.intra.vs.', '',contrasts[RH.vs.intra.ind])))), index.return = T)$ix]
  
  RH.vs.extra.ind <- grep('(.*RH.extra)(.*extra)',contrasts)
  sorted.RH.vs.extra.ind <- RH.vs.extra.ind[
    sort(as.numeric(gsub('P', '', gsub('\\.extra.*', '',
                                       gsub('RH.extra.vs.', '', contrasts[RH.vs.extra.ind])))), index.return = T)$ix]
  
  my.inds <- c(sorted.extra.vs.extra.ind, sorted.intra.vs.intra.ind, sorted.extra.vs.intra.ind, 
               sorted.RH.vs.extra.ind, sorted.RH.vs.intra.ind, RH.vs.RH.ind )
  XX$Contrast <- factor(contrasts, levels = contrasts[my.inds])
  
  
  Passages <- strsplit(as.character(XX$Contrast), split='\\.vs\\.')
  Passage2 <- unlist(lapply(strsplit(unlist(lapply(Passages, `[[`, 1)), split='\\.'), `[[`,1))
  Passage1 <- unlist(lapply(strsplit(unlist(lapply(Passages, `[[`, 2)), split='\\.'), `[[`,1))
  
  Category <- gsub("RH.intra", "RH",
                   gsub("RH.extra", "RH", 
                        gsub("^P.*[[:digit:]]\\.", "",
                             gsub("\\.P.*[[:digit:]]\\.", "\\.", 
                                  as.character(XX$Contrast)))))
  
  sorted.ind1 <- sort(as.numeric(gsub("P", "", gsub('RH', 'P300', unique(Passage1)))), index.return = T)$ix
  sorted.ind2 <- sort(as.numeric(gsub("P", "", gsub('RH', 'P300', unique(Passage2)))), index.return = T)$ix
  Passage1 <- factor((Passage1), levels = unique(Passage1)[sorted.ind1])
  Passage2 <- factor((Passage2), levels = unique(Passage2)[sorted.ind2])
  Category <- factor(Category, levels = sort(unique(Category)))
  
  XX$Category <- Category
  XX$Passage1 <- Passage1
  XX$Passage2 <- Passage2
  
  GeneSets <- XX$GeneSet
  
  XX$GeneSet <- factor(GeneSets, levels = unique(GeneSets))
  
  
  #XX <- XX %>% dplyr::filter(pvalue < 0.01) %>% dplyr::filter(Category != 'RH.vs.RH') %>% arrange(Contrast, GeneSet)
  #XX <- XX %>% dplyr::filter(pvalue < 0.01) %>% arrange(Contrast, GeneSet)
  XX <- XX %>%  arrange(Contrast, GeneSet)
  XX$pvalue <- as.numeric(XX$pvalue)
  colnames(XX) <- c("GeneSet", "genes.in.set", "total.genes.in.set", "all.genes",
                    "Contrast", "genes.in.contrast", "total.genes.in.contrast", "overlap",
                    "overlap.genes","pvalue","Category", "Passage1", "Passage2")
  
  return(XX)
}
