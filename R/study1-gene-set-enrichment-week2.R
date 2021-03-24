##-------------------------------------
## study1-gene-set-enrichment-week2.R
##
## Title:
## Purpose: 
## Author:
##
##
##
##-------------------------------------
## Notes:
#
#
#
#
#
#
#
#
## ------------------------------------

source("./R/libraries.R")
source("./R/themes.R")



new_line_fun <- function(x){
  
  l <- sapply(strsplit(x, " "), length)
  
  if(l > 2){
    
    paste(paste(strsplit(x, " ")[[1]][1:2], collapse = " "), 
          "<br>", 
          paste(strsplit(x, " ")[[1]][3:l], collapse = " "), 
          collaspe = "")
    
  } else {
    x
  }
  
}

new_line_fun_n <- function(x){
  
  l <- sapply(strsplit(x, " "), length)
  
  if(l > 2){
    
    paste(paste(strsplit(x, " ")[[1]][1:2], collapse = " "), 
          "\n", 
          paste(strsplit(x, " ")[[1]][3:l], collapse = " "), 
          collaspe = "")
    
  } else {
    x
  }
  
}


firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}



########### Gene set enrichement analysis ###############################

# Load data 

# Resting biopsy models
## Interaction model 
mm_interaction <- readRDS("./data/study-1b/data/derivedData/DE/mixedmodel2_results1.RDS")

## Time model
mm_time <- readRDS("./data/study-1b/data/derivedData/DE/mixedmodel2_results2.RDS")

# Acute models 
acute.interaction <- readRDS("./data/study-1b/data/derivedData/DE/mixedmodel2_acutemodel.RDS")

acute.time <- readRDS(file = "./data/study-1b/data/derivedData/DE/mixedmodel2_acutemodel_timeonly.RDS")

## Combine data frames 
mm_dat <- mm_interaction %>%
        dplyr::select(gene, method, model, coef, estimate, se, z.val, p.val) %>%
        mutate(interaction = TRUE) %>%
        rbind(mm_time %>%
                      dplyr::select(gene, method, model, coef, estimate, se, z.val, p.val) %>%
                      mutate(interaction = FALSE)) %>% 
        rbind(acute.interaction %>%
                      dplyr::select(gene, method, model, coef, estimate, se, z.val, p.val) %>%
                      mutate(interaction = TRUE)) %>%
        rbind(acute.time %>%
                      dplyr::select(gene, method, model, coef, estimate, se, z.val, p.val) %>%
                      mutate(interaction = FALSE)) %>%
        group_by(model, interaction, coef) %>%
        mutate(ciu = estimate + qnorm(0.975) * se, 
               cil = estimate - qnorm(0.975) * se, 
               msd = if_else(estimate > 0, cil, -ciu), 
               log2fc = estimate/log(2)) %>%  # msd calculation Wald 95% CI
        
        filter((coef %in% c("timew2pre", "timew12", "timew2post") & interaction == FALSE) | 
                       coef %in% c("timew2pre:setsmultiple", 
                                   "timew12:setsmultiple", 
                                   "timew2post:setsmultiple")) %>%
        mutate(mod_coef = paste0(model,".", coef)) 



### Load enrichment analysis results #####

# Read go terms from go-analysis.R

go_terms <- readRDS("./data/study-1b/data/derivedData/GOanalysis/go_terms.RDS")

## Extracts gene ontology ids that were identified with ORA
ORAID <- go_terms %>%
        
        filter(coef == "timew2pre:setsmultiple",
               model == "tissue_offset_lib_size_normalized",
               type == "full") %>% 
        mutate(ora = "IN_ORA") %>%
        pull(ID)



# fgsea
full_gsea <- readRDS(file = "./data/study-1b/data/derivedData/GOanalysis/fgsea_collapse.RDS")
collapse_gsea <- readRDS(file = "./data/study-1b/data/derivedData/GOanalysis/fgsea_full.RDS")

# cerno
cerno.test <- readRDS("./data/study-1b/data/derivedData/cerno.RDS")
cerno.test <- cerno.test$cerno.test


### Combine information from fgsea and cerno tests 

comb_gsea <- cerno.test %>%
        mutate(name = gsub("Go ", "", Title)) %>%
        dplyr::select(ID, name, model, coef, go.cat, cerno, N1, cerno.auc = AUC, 
                      cerno.p = P.Value, cerno.padj = adj.P.Val) %>%
        
        inner_join(full_gsea %>%
                           dplyr::select(-leadingEdge) %>%
                           mutate(name = tolower(gsub("GO_", "", pathway)),
                                  name = gsub("_", " ", name)) %>%
                           dplyr::select(name, model, coef, go.cat, fgsea.p = pval, fgsea.padj = padj, 
                                         ES, NES, size)) %>%
        
        print()


comb_data_w2pre <- comb_gsea %>%
        filter(coef == "timew2pre:setsmultiple", 
               model == "tissue_offset_lib_size_normalized") %>%
        filter(cerno.p < 0.000001) %>%   ### Strict filter for plotting
        rowwise() %>%
        mutate(name = firstup(name), 
               name = new_line_fun_n(name)) %>%
        ungroup() %>%
        mutate(ora = if_else(ID %in% ORAID, "Identified\nin ORA", "Not identified\nin ORA")) %>%
        mutate(go.cat = factor(go.cat, levels = c("bp", "cc", "mf"), 
                               labels = c("Biological processes", 
                                          "Cellular component", 
                                          "Molecular function"))) %>%
        print()





### Gene set enrichment week 2 pre bubble graph ########

sets_select <- comb_gsea %>%
        filter(coef == "timew2pre:setsmultiple", 
               model == "tissue_offset_lib_size_normalized",
               cerno.p < 0.0000000000000001) %>%
        pull(name)







top_go_week2 <- comb_gsea %>%
        filter(coef == "timew2pre:setsmultiple",
               model != "naive",
               name %in% sets_select) %>%
        rowwise() %>%
        mutate(name = firstup(name), 
               name = new_line_fun(name)) %>%
        ungroup() %>%
        mutate(go.cat = factor(go.cat, levels = c("bp", "cc", "mf"), 
                               labels = c("Biological processes", 
                                          "Cellular component", 
                                          "Molecular function")), 
               model = factor(model, levels = c("lib_size_normalized", 
                                                "tissue_offset_lib_size_normalized"), 
                              labels = c("Effective library size", 
                                         "Tissue offset"))) %>%

        
        ggplot(aes(-log10(cerno.padj), reorder(name, -log10(cerno.padj)),
                   fill = model)) + 
        
        
        scale_fill_manual(values = c(group.study.color[4], group.study.color[2]) )+
        
        labs(x = "Rank based enrichment test<br>(-Log<sub>10</sub> <i>P</i>-values)") +
        
        geom_point(size = 3, shape = 21)  +
        theme_bw() +
        theme(panel.grid = element_blank(), 
              panel.border = element_blank(), 
              axis.line.x  = element_line(size = line.size), 
              axis.line.y = element_blank(), 
              axis.ticks.x = element_line(size = line.size), 
              axis.ticks.y = element_blank(), 
              axis.text.x =  element_markdown(color = "black", size = 7), 
              axis.title.x = element_markdown(color = "black", size = 7),
              axis.text.y =  element_markdown(color = "black", size = 7), 
              axis.title.y = element_blank(),
              strip.text.x = element_markdown(size = 8, face = "bold"),
              legend.title = element_blank(), 
              legend.background = element_rect(fill = "white"),
              legend.margin = margin(t = 0, r = 1, b = 1, l = 1, unit = "pt"),
              legend.key = element_rect(fill = "white"),
              legend.text = element_text(size = 7),
              legend.direction = "vertical",
              legend.position = "bottom")



w2pre_tissue_bubble <-  comb_data_w2pre %>%
        
        ggplot(aes(-log10(fgsea.padj), NES, label = name, fill = ora)) +
        geom_point(alpha = 0.4, shape = 21, size = 3) + 
        geom_text_repel(data = comb_data_w2pre[-log10(comb_data_w2pre$cerno.padj) > 10,], 
                        size = 2.5, 
                        lineheight = 0.8) +
        facet_wrap(~ go.cat, ncol = 1) +
        
        scale_fill_manual(values = c(group.study.color[4],  group.study.color[2]), guide = guide_legend(ncol = 2)) +

        labs(size = "Set size", 
             fill = "", 
             x = "-Log<sub>10</sub>(Adjusted <em>P</em>-values)", 
             y = "Normalized enrichment score (NES)") +
        
        dissertation_theme() +
        
        theme(strip.background = element_rect(fill = "white", color = "white"),
              panel.spacing = unit(1, "lines"),
              strip.text = element_text(size = 7, hjust = 0),
              legend.text = element_text(size = 7),
              axis.title.x = element_markdown(size = 7),
              legend.title = element_blank(),
              legend.spacing.x = unit(0.2, 'cm'),
              legend.spacing.y = unit(0.3, 'cm'),
              legend.key.size = unit(0.2, "cm"),
              strip.text.x = element_markdown(size = 7),
             legend.direction  = "horizontal",
             #legend.box = "horizontal",
             #legend.margin = margin(t = -0.5,r = -0.2,  b = 0, l = -0.3, unit = "pt"),
              # panel.border = element_rect(colour = "black"),
              legend.position = "bottom") 


# A figure that shows global shift in log2 fold changes and 
# subsequent go categories 


#### Location of "leading edge" plots

gene_name_convert <- clusterProfiler::bitr(unique(mm_dat$gene), fromType = "ENSEMBL",
                          toType = c( "ENTREZID", "SYMBOL"),
                          OrgDb =org.Hs.eg.db::org.Hs.eg.db)


go_set <- full_gsea %>%
        filter(model == "tissue_offset_lib_size_normalized", 
               coef == "timew2pre:setsmultiple") %>% 
        filter(pathway == "GO_COLLAGEN_CONTAINING_EXTRACELLULAR_MATRIX") %>%
        pull(leadingEdge) %>%
        unlist()


### extract gene sets that are called as significant in DE
sig_to <- mm_dat %>%
        group_by(model, coef) %>%
        mutate(padj = p.adjust(p.val, method = "fdr")) %>%
        ungroup() %>%
        inner_join(gene_name_convert %>% mutate(gene = ENSEMBL)) %>%
        filter(SYMBOL %in% go_set, 
               model == "tissue_offset_lib_size_normalized", 
               coef == "timew2pre:setsmultiple", 
               log2fc > abs(0.5), padj < 0.05) %>%
        dplyr::select(SYMBOL, log2fc, padj) %>%
        
        pull(SYMBOL)

sig_lib <- mm_dat %>%
        group_by(model, coef) %>%
        mutate(padj = p.adjust(p.val, method = "fdr")) %>%
        ungroup() %>%
        inner_join(gene_name_convert %>% mutate(gene = ENSEMBL)) %>%
        filter(SYMBOL %in% go_set, 
               model == "lib_size_normalized", 
               coef == "timew2pre:setsmultiple", 
               log2fc > abs(0.5), padj < 0.05) %>%
        dplyr::select(SYMBOL, log2fc, padj) %>%
        
        pull(SYMBOL)

# No genes in naive model differentially expressed in the go category 


#### Stats for enrichment analysis

# Cerno 
cerno.p <- cerno.test %>%
        filter(Title == "Go collagen containing extracellular matrix", 
               coef == "timew2pre:setsmultiple") %>%
        pull(adj.P.Val)

names(cerno.p) <- c("Effective library size", "Tissue offset")



tissue_offset_text <- paste0("Tissue offset")
lib_size_text <- paste0("Effective library-size")
naive_text <- paste0("Naive")


# tissue_offset_text <- paste0("Tissue offset\nRank p-value = ", formatC(cerno.p[2], format = "e", digits = 2))
# lib_size_text <- paste0("Effective library-size\nRank p-value = ", formatC(cerno.p[1], format = "e", digits = 2))
# naive_text <- paste0("Naive")


# Nes
gsea.stats <- full_gsea %>%
        filter(coef == "timew2pre:setsmultiple") %>% 
        filter(pathway == "GO_COLLAGEN_CONTAINING_EXTRACELLULAR_MATRIX") %>%
        dplyr::select(- leadingEdge) %>%
        dplyr::select(model, pathway, padj, ES, NES)

density.dat <- mm_dat %>%
        group_by(model, coef) %>%
        mutate(padj = p.adjust(p.val)) %>%
        inner_join(gene_name_convert %>% mutate(gene = ENSEMBL)) %>%
        group_by(SYMBOL, model, coef) %>%
        summarise(log2fc = mean(log2fc)) %>%
        ungroup() %>%
        mutate(SYMBOL = fct_reorder(SYMBOL, log2fc), 
               highlight = if_else(SYMBOL %in% go_set, "high", "background")) %>%
        mutate(symbol = paste0("italic(", SYMBOL, ")")) %>%
        print()


### Tissue offset model 

tissue_offset_rug <- density.dat %>%  
        filter(model %in% c("tissue_offset_lib_size_normalized"),  
               coef == "timew12:setsmultiple") %>%   
        
        ggplot(aes(log2fc)) + 
        
        geom_vline(xintercept = 0, lty = 2, color = "gray20") +
        
        geom_density(alpha = 0.5, fill = group.study.color[4]) +
        
        # scale for log2 fc
        
        annotate("segment", 
                 x = c(2.7, 2.65, 2.65), 
                 xend = c(2.7, 2.75, 2.75), 
                 y = c(1, 1, 2), 
                 yend = c(2, 1, 2), size = 0.4) +
        
        annotate("text", x = c(2.83, 2.83), y = c(1, 2), label = c(1, 2), size = 2.2, color = "black") +
        annotate("richtext", x = 3.2, y = 1.5, 
                 hjust = 0.5, 
                 label = "Log<sub>2</sub><br>(Moderate-/Low-volume)", 
                 size = 2.5, 
                 angle = -90, 
                 color = "black", 
                 fill = NA, 
                 label.color = NA,
                 label.padding = unit(0, "pt")) +
        
        
        annotate("text", 
                 x = -3.4, 
                 y = 3.8, 
                 label = "Collagen containing\nextracellular matrix", 
                 lineheight = 0.8,
                 size = 3, 
                 hjust = 0) +

       annotate("text", 
                x = -3.4, 
                y = 2, 
                label = tissue_offset_text, 
                size = 2.5, 
                hjust = 0, 
                lineheight = 0.8) +
      
        
        # Plot all log2 fold changes
        geom_segment(aes(x = log2fc, xend = log2fc, y =  0, yend = 0 + log2fc), 
                     color = "gray80") +
        
        # Plot log2 fold changes contained in gene set 
        geom_segment(data =  density.dat %>%
                             filter(model %in% c("tissue_offset_lib_size_normalized"),  
                                    coef == "timew2pre:setsmultiple", 
                                    SYMBOL %in% go_set),
                     aes(x = log2fc, xend = log2fc, y =  0, yend = 0 + log2fc), color = "gray20") +
        
        # Plot gene symbols connected to "top leading edge" genes
        
        geom_text_repel(data = density.dat %>%
                                filter(model %in% c("tissue_offset_lib_size_normalized"),  
                                       coef == "timew2pre:setsmultiple", 
                                       SYMBOL %in% sig_to), 
                        aes(x = log2fc, 
                            y = 0, 
                            label = symbol, 
                            color = NULL), 
                        nudge_y   = -3.6,
                        nudge_x = 0.6,
                        direction    = "x",
                        segment.color = "black",
                        box.padding = 0.01,
                        segment.alpha = 0.3,
                        parse = TRUE,
                        angle = 90,
                        vjust        = 1,
                        segment.size = 0.2,
                        size = 2.2)  +
        
        
        scale_x_continuous(limits = c( -3.5, 3.5)) + 
        theme(legend.position = "none") +
        theme_void()


#### Library size rug

lib_size_rug <- density.dat %>%  
        filter(model %in% c("lib_size_normalized"),  
               coef == "timew12:setsmultiple") %>%   
        
        ggplot(aes(log2fc)) + 
        
        geom_vline(xintercept = 0, lty = 2, color = "gray20") +
        
        geom_density(alpha = 0.5, fill = group.study.color[2]) +
        
        
        
      annotate("text", 
               x = -3.4, 
               y = 2, 
               label = lib_size_text, 
               size = 2.5, 
               hjust = 0, 
               lineheight = 0.8) +
        
        
        # Plot all log2 fold changes
        geom_segment(aes(x = log2fc, xend = log2fc, y =  0, yend = 0 + log2fc), 
                     color = "gray80") +
        
        # Plot log2 fold changes contained in gene set 
        geom_segment(data =  density.dat %>%
                             filter(model %in% c("lib_size_normalized"),  
                                    coef == "timew2pre:setsmultiple", 
                                    SYMBOL %in% go_set), 
                     aes(x = log2fc, xend = log2fc, y =  0, yend = 0 + log2fc), color = "gray20") +
        
        # Plot gene symbols connected to "top leading edge" genes
        
        geom_text_repel(data = density.dat %>%
                                filter(model %in% c("lib_size_normalized"),  
                                       coef == "timew2pre:setsmultiple", 
                                       SYMBOL %in% sig_lib), 
                        aes(x = log2fc, 
                            y = 0, 
                            label = symbol, 
                            color = NULL), 
                        nudge_y   = -3.6,
                        nudge_x = 0.6,
                        color = "black",
                        direction    = "x",
                        segment.color = "black",
                        box.padding = 0.01,
                        segment.alpha = 0.3,
                        parse = TRUE,
                        angle = 90,
                        vjust        = 1,
                        segment.size = 0.2,
                        size = 2.2)  +
        
        
        scale_x_continuous(limits = c( -3.5, 3.5)) + 
        theme(legend.position = "none") +
        theme_void()


#### Naive rug

naive_rug <- density.dat %>%  
        filter(model %in% c("naive"),  
               coef == "timew12:setsmultiple") %>%   
        
        
        
        ggplot(aes(log2fc)) + 
        
        
        ## Annotate fold-change scale on the x axis bottom plot 
        
        annotate("segment", x = 0, xend = 0, 
                 y = -0.5, yend = 5.2,
                 lty = 2, 
                 color = "gray20") +
        
        annotate("text", 
                 x = -3.4, 
                 y = 2, 
                 label = naive_text, 
                 size = 2.5, 
                 hjust = 0) +
        
        annotate("segment", 
                 x = c(0, 0, 1), 
                 xend = c(1, 0, 1), 
                 y = c(-1, -0.9, -0.9), 
                 yend = c(-1, -1.1, -1.1), 
                 size = 0.4) +
        
        annotate("text", x = c(0, 1), y = c(-1.5, -1.5), label = c(0, 1), size = 2.2, color = "black") +
        annotate("richtext", x = 0.5, y = -2.2, 
                 hjust = 0.5, 
                 label = "Log<sub>2</sub>(Moderate-/Low-volume)", 
                 size = 2.5, 
                 
                 color = "black", 
                 fill = NA, 
                 label.color = NA,
                 label.padding = unit(0, "pt")) +
        
        
        
        
        geom_density(alpha = 0.2, fill = group.study.color[5]) +
        
        
        # Plot all log2 fold changes
        geom_segment(aes(x = log2fc, xend = log2fc, y =  0, yend = 0 + log2fc), 
                     color = "gray80") +
        
        # Plot log2 fold changes contained in gene set 
        geom_segment(data =  density.dat %>%
                             filter(model %in% c("naive"),  
                                    coef == "timew2pre:setsmultiple", 
                                    SYMBOL %in% go_set), 
                     aes(x = log2fc, xend = log2fc, y =  0, yend = 0 + log2fc), color = "gray20") +
        
        # Plot gene symbols connected to "top leading edge" genes
        
        
        scale_x_continuous(limits = c( -3.5, 3.5)) + 
        theme(legend.position = "none") +
        theme_void()


gene_set_rug <- plot_grid(tissue_offset_rug, 
                          lib_size_rug, 
                          naive_rug, ncol = 1)


go_analysis_week2 <- plot_grid(w2pre_tissue_bubble, 
                               plot_grid(NULL, top_go_week2, NULL, ncol = 1, rel_heights = c(0.25, 1, 0.25)), 
                               ncol = 2, rel_widths = c(0.5, 0.5))

#### Save plots 


saveRDS(gene_set_rug, "./data/derivedData/study1-gene-set-enrichment-week2/gene-set-rug.RDS")





saveRDS(go_analysis_week2, "./data/derivedData/study1-gene-set-enrichment-week2/go_week2.RDS")



############## Complete figure ###################
