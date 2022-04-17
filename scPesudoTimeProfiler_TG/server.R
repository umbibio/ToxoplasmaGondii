library(shiny)
library(DT)
library(Seurat)
library(ggplot2)
library(rsconnect)
library(sme)
library(tidyverse)
library(igraph)
library(ggraph)
library(graphlayouts)

source("util_funcs.R")
S.O.tg.gam <- readRDS('S.O.tg.gam.RData')
sc.tc.fits <- readRDS("tg_sme_fits_sc_tc_20min.RData")
sc.tc.df <- readRDS('tg_sc_tc_df.RData')
prod.desc <- readRDS('tg_prod_desc.RData')
tg.clust.net <- readRDS('tg_opt_net.RData')
tg.clust.net.d <- readRDS('tg_opt_net_d.RData') 
#tg.clust.net <- readRDS('tg_clust_net.RData')
#tg.clust.net.d <- readRDS('tg_clust_net_d.RData') 


tg.clust.net.d <- tg.clust.net.d %>% dplyr::select(GeneID, degree, phase)
tg.e.l <- readRDS('tg_opt_e_l.RData') %>% dplyr::select(src, trg, trgDesc)
#tg.e.l <- readRDS('tg_e_l.RData') %>% dplyr::select(src, trg, trgDesc)


genes <- gsub('-', '_', rownames(S.O.tg.gam))
prod.desc <- prod.desc[prod.desc$GeneID %in% genes, ]
prod.desc <- prod.desc[prod.desc$GeneID %in% sc.tc.df$variable, ]
prod.desc <- prod.desc[prod.desc$GeneID %in% tg.clust.net.d$GeneID, ]
prod.desc <- left_join(prod.desc, tg.clust.net.d, by = 'GeneID')

sc.tc.mus <- smoothSplineSmeFits(sc.tc.fits, unique(sc.tc.df$variable), extend = F)
colnames(sc.tc.mus) <- c('GeneID', 't', 'y')
sc.tc.mus <- sc.tc.mus%>%
  pivot_wider(names_from = 'GeneID', values_from = 'y') %>% 
  as.data.frame()

sc.tc.mus.scale <- sc.tc.mus
sc.tc.mus.scale[,2:ncol(sc.tc.mus.scale)] <- scale(sc.tc.mus.scale[,2:ncol(sc.tc.mus.scale)],
                                                   center = T,scale = T)
sc.tc.mus.scale <- sc.tc.mus.scale %>%  as.data.frame() %>% 
  pivot_longer(starts_with('TG'), names_to = 'GeneID', values_to = 'y')

got_palette <- c("#1A5878", "#C44237", "#AD8941", "#E99093", "#50594B")

server <- function(input, output, session) {
  
  df <- as.data.frame(prod.desc) %>% arrange(desc(degree))
  options(DT.options = list(pageLength = 5))
  
  # row selection
  output$x4 = DT::renderDataTable(df, server = FALSE, selection = 'single')
  
  ## expression  
  output$x1 = renderPlot({
    s = input$x4_rows_selected
    selected.Gene <- gsub('_', '-', df$GeneID[s])
    if(length(s)){
      if(!is.na(selected.Gene)){
        p1 <- FeaturePlot(object = S.O.tg.gam, features = selected.Gene, 
                         label = T, pt.size = 1, label.size = 5,
                         cols = c("lightgrey", "blue"), reduction = "pca") + 
          theme_minimal() + 
          labs(title="T. Gondii",x ="PC1", y = "PC2")+
          theme(
            plot.title = element_text(size=12, face="bold.italic"),
            axis.title.x = element_text(size=10, face="bold"),
            axis.title.y = element_text(size=10, face="bold")
          )
        
        plot(p1)
      }
    }
  })
  
  output$x2 = renderPlot({
    s = input$x4_rows_selected
    selected.Gene <- df$GeneID[s]
    if(length(s)){
      if(!is.na(selected.Gene)){
        ind <- which(unique(sc.tc.df$variable) == selected.Gene )
        plot.sme(sc.tc.fits[[ind]], paste('sc', selected.Gene))
      }
    }
  })
  

  # output$x3 = renderPlot({
  #   s = input$x4_rows_selected
  #   selected.Gene <- df$GeneID[s]
  #   if(length(s)){
  #     if(!is.na(selected.Gene)){
  #       ## Node networks
  #       nds <- gsub('_', '-', selected.Gene)
  #       tg.clust.net.nds <- set_edge_attr(tg.clust.net, 'subnet', E(tg.clust.net)[from(nds)], nds)
  #       tg.clust.net.nds <- set_edge_attr(tg.clust.net.nds, 'subnet.weight', E(tg.clust.net.nds)[from(nds)], 2)
  #       tg.clust.net.nds <- set_edge_attr(tg.clust.net.nds, 'subnet.weight', E(tg.clust.net.nds)[!from(nds)], 1)
  #       
  #       ggraph(tg.clust.net.nds,layout = "stress") +
  #         geom_edge_link(aes(colour = subnet, width = subnet.weight)) + 
  #         #geom_edge_link0(aes(colour = subnet))+
  #         #geom_edge_link0(edge_colour = "grey66")+
  #         geom_node_point(aes(fill = phase, size = size),shape=21)+
  #         geom_node_text(aes(filter = name == nds, label = name),family="serif")+
  #         scale_fill_manual(values = got_palette)+
  #         scale_edge_width(range = c(0.2,3))+
  #         scale_size(range = c(1,6))+
  #         theme_graph()+
  #         theme(legend.position = "none")
  #       
  #     }
  #   }
  # })
  
  observeEvent(input$x4_rows_selected,{
    # observe element change and render table
    s <- input$x4_rows_selected
    selected.Gene <- df$GeneID[s]
    df2 <- tg.e.l[tg.e.l$src == selected.Gene,] %>% distinct()
    output$x3 = DT::renderDataTable(df2, selection = 'none', options = list(pageLength = 5, 
                                                                            searching = FALSE, 
                                                                            lengthChange = FALSE, 
                                                                            rownames= FALSE))
    #output$x3 = DT::renderDataTable(df2, options = list(dom = 'ft'))
  })
  # # row selection
  # s = input$x4_rows_selected
  # selected.Gene <- gsub('_', '-', df$GeneID[s])
  # df2 <- tg.e.l[tg.e.l$src == selected.Gene,]
  
  
  
  
  # output$x3 = renderPlot({
  #   s = input$x4_rows_selected
  #   selected.Gene <- df$GeneID[s]
  #   if(length(s)){
  #     if(!is.na(selected.Gene)){
  #       tmp <- sc.tc.mus.scale %>% dplyr::filter(GeneID == selected.Gene)
  #       v.ind <- which(V(tg.clust.net) == selected.Gene)
  #       V(tg.clust.net)$color[v.ind] <- 'red'
  #       p <- ggraph(tg.clust.net,layout = "stress")+
  #         #geom_edge_link0(aes(edge_width = weight),edge_colour = "grey66")+
  #         geom_edge_link0(edge_colour = "grey66")+
  #         #geom_node_point(aes(fill = phase, size = size),shape=21)+
  #         geom_node_point(aes(size = size),shape=21)+
  #         geom_node_text(aes(filter = size>30, label = name),family="serif")+
  #         scale_fill_manual(values = got_palette)+
  #         scale_edge_width(range = c(0.2,3))+
  #         scale_size(range = c(1,6))+
  #         theme_graph()+
  #         theme(legend.position = "right")
  #       
  #       # p <- ggplot(tmp, aes(x = t, y = y)) + 
  #       #   geom_line(color = 'red') + 
  #       #   theme_minimal() + 
  #       #   labs(title=selected.Gene,x ="Time (h)", y = "z-score(log2(expr))")+
  #       #   theme(
  #       #     plot.title = element_text(size=12, face="bold.italic"),
  #       #     axis.title.x = element_text(size=10, face="bold"),
  #       #     axis.title.y = element_text(size=10, face="bold")
  #       #   )
  #       
  #       plot(p)
  #     }
  #   }
  # })
  
  
  session$allowReconnect(TRUE)
}