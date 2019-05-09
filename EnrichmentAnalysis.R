library(openxlsx)
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(forcats)

count.tab  <- read.xlsx('~/work/ToxoplasmaGondii/Enrichment_Counts.xlsx')
pvalue.tab <- read.xlsx('~/work/ToxoplasmaGondii/Enrichment_pValues.xlsx')

colnames(pvalue.tab) <- colnames(count.tab)[1:ncol(pvalue.tab)]
colnames(pvalue.tab)[1] <- 'Contrast' 

pvalue.tab <- pvalue.tab[-grep('B4', pvalue.tab$Contrast),]
pvalue.tab <- pvalue.tab[-grep('freshlyse', pvalue.tab$Contrast),]
pvalue.tab <- pvalue.tab[-grep("B2\\.P145\\.intra\\.over\\.B2\\.P11\\.intra" , pvalue.tab$Contrast),]
## this one is duplicated
pvalue.tab <- pvalue.tab[-grep("RH\\.intra\\.over\\.B2\\.P84\\.intra" , pvalue.tab$Contrast),] 

## Names are not consistant. Sometimes RH comes first and sometimes second. Same with intra over extra
pvalue.tab <- pvalue.tab[-grep("B2\\.P148\\.6hes\\.over\\.B2\\.P11\\.6hes" , pvalue.tab$Contrast),]
pvalue.tab$Contrast[pvalue.tab$Contrast == "B.P7.intra.over.B2.P11.intra"] <- "B2.P11.intra.over.B.P7.intra"
pvalue.tab$Contrast[pvalue.tab$Contrast == "B.P7.intra.over.B.P7.6hes"] <- "B.P7.6hes.over.B.P7.intra"
pvalue.tab$Contrast[pvalue.tab$Contrast == "RH.intra.over.B.P7.intra"] <- "B.P7.intra.over.RH.intra"


contrasts <- gsub('B\\.', '',gsub('over', 'vs',gsub('6hes', 'extra', gsub('B2\\.', '',pvalue.tab$Contrast))))

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

intra.vs.RH.ind <- grep('(.*intra)(.*RH)',contrasts)
sorted.intra.vs.RH.ind <- intra.vs.RH.ind[
  sort(as.numeric(gsub('P', '', gsub('\\.intra.*', '',contrasts[intra.vs.RH.ind]))), index.return = T)$ix]

extra.vs.RH.ind <- grep('(.*extra)(.*RH.extra)',contrasts)
sorted.extra.vs.RH.ind <- extra.vs.RH.ind[
  sort(as.numeric(gsub('P', '', gsub('\\.extra.*', '',contrasts[extra.vs.RH.ind]))), index.return = T)$ix]

my.inds <- c(sorted.extra.vs.extra.ind, sorted.intra.vs.intra.ind, sorted.extra.vs.intra.ind, 
             sorted.extra.vs.RH.ind, sorted.intra.vs.RH.ind, RH.vs.RH.ind )
contrasts[my.inds]
pvalue.tab$Contrast <- factor(contrasts, levels = contrasts[my.inds])
pvalue.tab <- pvalue.tab[my.inds, ]
XX <- pvalue.tab  %>% gather(key = GO, value = pvalue, -Contrast)
GO <- gsub("tgo.*[[:digit:]]_", "", gsub("development", "devel", gsub("tachyzoite", "tachy", gsub("bradyzoite", "brady", gsub("Bradyzoite", "Brady", 
                                 gsub("Tachyzoite", "Tachy", gsub("\\.\\.", '\\.', XX$GO)))))))

categories <- c("RH.on.RH", "Intra.on.RH", "Extra.on.RH", "Extra.on.Intra", "Intra.on.Intra", "Extra.on.Extra")
categories <- gsub("RH.intra", "RH",gsub("RH.extra", "RH", gsub("^P.*[[:digit:]]\\.", "", 
                                                                gsub("\\.P.*[[:digit:]]\\.", "\\.", as.character(XX$Contrast)))))
categories <- factor(categories, levels = unique(categories))
XX$Category <- categories

XX$GO <- factor(GO, levels = unique(GO))
XX <- XX %>% dplyr::filter(pvalue < 0.01) %>% dplyr::filter(Category != 'RH.vs.RH') %>% arrange(Contrast)

XX$pvalue <- as.numeric(XX$pvalue)

## Category of contrasts
p <- ggplot(XX, aes(x = GO, y = Contrast)) + 
  geom_point(aes(color = pvalue, size = -log(pvalue))) +
  theme_bw(base_size = 14) +
  scale_colour_gradient(limits=c(0, 0.01), low="red") +
  ylab(NULL) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + facet_grid(Category ~ ., scales='free') + 
  theme(strip.background = element_rect(colour="black", fill="white", 
                                        size=1, linetype="solid"))
  
plot(p)

ggsave(filename="~/work/ToxoplasmaGondii/ToxoR/output/enrichment.pdf", plot=p,
       width = 10, height = 10, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)


#### Repeat the same for AP2's
AP2.count.tab  <- read.xlsx('~/work/ToxoplasmaGondii/AP2_overlap_counts.xlsx')
AP2.pvalue.tab <- read.xlsx('~/work/ToxoplasmaGondii/AP2_overlap_pvalues.xlsx')

colnames(AP2.pvalue.tab)[1] <- 'Contrast' 

AP2.pvalue.tab <- AP2.pvalue.tab[-grep('B4', AP2.pvalue.tab$Contrast),]
AP2.pvalue.tab <- AP2.pvalue.tab[-grep('freshlyse', AP2.pvalue.tab$Contrast),]
AP2.pvalue.tab <- AP2.pvalue.tab[-grep("B2\\.P145\\.intra\\.over\\.B2\\.P11\\.intra" , AP2.pvalue.tab$Contrast),]
AP2.pvalue.tab <- AP2.pvalue.tab[-grep("B2\\.P148\\.6hes\\.over\\.B2\\.P11\\.6hes" , AP2.pvalue.tab$Contrast),]
## this one is duplicated
AP2.pvalue.tab <- AP2.pvalue.tab[-grep("RH\\.intra\\.over\\.B2\\.P84\\.intra" , AP2.pvalue.tab$Contrast),] 


## Names are not consistant. Sometimes RH comes first and sometimes second. Same with intra over extra
AP2.pvalue.tab$Contrast[AP2.pvalue.tab$Contrast == "B.P7.intra.over.B2.P11.intra"] <- "B2.P11.intra.over.B.P7.intra"
AP2.pvalue.tab$Contrast[AP2.pvalue.tab$Contrast == "B.P7.intra.over.B.P7.6hes"] <- "B.P7.6hes.over.B.P7.intra"
AP2.pvalue.tab$Contrast[AP2.pvalue.tab$Contrast == "RH.intra.over.B.P7.intra"] <- "B.P7.intra.over.RH.intra"


contrasts <- gsub('B\\.', '',gsub('over', 'vs',gsub('6hes', 'extra', gsub('B2\\.', '',AP2.pvalue.tab$Contrast))))

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

intra.vs.RH.ind <- grep('(.*intra)(.*RH)',contrasts)
sorted.intra.vs.RH.ind <- intra.vs.RH.ind[
  sort(as.numeric(gsub('P', '', gsub('\\.intra.*', '',contrasts[intra.vs.RH.ind]))), index.return = T)$ix]

extra.vs.RH.ind <- grep('(.*extra)(.*RH.extra)',contrasts)
sorted.extra.vs.RH.ind <- extra.vs.RH.ind[
  sort(as.numeric(gsub('P', '', gsub('\\.extra.*', '',contrasts[extra.vs.RH.ind]))), index.return = T)$ix]

my.inds <- c(sorted.extra.vs.extra.ind, sorted.intra.vs.intra.ind, sorted.extra.vs.intra.ind, 
             sorted.extra.vs.RH.ind, sorted.intra.vs.RH.ind, RH.vs.RH.ind )
contrasts[my.inds]
AP2.pvalue.tab$Contrast <- factor(contrasts, levels = contrasts[my.inds])
AP2.pvalue.tab <- AP2.pvalue.tab[my.inds, ]

## Only look at AP2's
## For some reason, there are duplicated columns!
AP2.pvalue.tab <- AP2.pvalue.tab[,!duplicated(colnames(AP2.pvalue.tab))]
AP2.pvalue.tab <- AP2.pvalue.tab %>% dplyr::select(c(1,starts_with('AP2')))
colnames(AP2.pvalue.tab) <- gsub('\\.', '-', colnames(AP2.pvalue.tab))
XX <- AP2.pvalue.tab  %>% gather(key = AP2, value = pvalue, -Contrast)

XX$AP2 <- factor(XX$AP2, levels = unique(XX$AP2))
XX <- XX %>% dplyr::filter(pvalue < 0.05) %>% arrange(Contrast)

XX$pvalue <- as.numeric(XX$pvalue)
p <- ggplot(XX, aes(x = AP2, y = Contrast)) + 
  geom_point(aes(color = pvalue, size = -log(pvalue))) +
  theme_bw(base_size = 14) +
  scale_colour_gradient(limits=c(0, 0.01), low="red") +
  ylab(NULL) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot(p)

ggsave(filename="~/work/ToxoplasmaGondii/ToxoR/output/AP2s.pdf", plot=p,
       width = 8, height = 8, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)
