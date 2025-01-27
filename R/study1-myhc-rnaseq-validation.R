
##-------------------------------------
## study1-myhc-rnaseq-validation.R
##
## Title: MyHC distributions and RNA-seq validation in Study I
## Purpose: 
## Author:
##
##
##
##-------------------------------------

#### Myosin heavy chain comparison between RNA-seq and Protein family-based normalization. ####

# Notes: 

# To circumvent any issues with normalization when comparing alignment tools, a biological validation
# can be done through correlation analysis of mRNA to protein expression. MyHC are well suited for this 
# purpose. The idea is to normalize to gene-family (Myosins) and compare correlations between alignment tools. 



# Source functions and settings

source("./R/libraries.R")
source("./R/themes.R")



# read DGEList

dge_lists <- readRDS("./data/study-1b/data/derivedData/dge_lists/dge_list.RDS")


## Download gene family information from HUGO

# Myosin heavy chain gene family 

url <- "https://www.genenames.org/cgi-bin/genegroup/download?id=1098&type=node"


myosin_family <- read_delim(url, delim = "\t") %>%
        dplyr::select(symbol = 'Approved symbol', gene = 'Ensembl gene ID') %>%
        print()


results <- list()

for(i in 1:length(dge_lists)) {
        
        
        
        temp  <- dge_lists[[i]][rownames(dge_lists[[i]]) %in% myosin_family$gene,]
        
        myhc_all  <- temp$counts %>%
                data.frame(gene = rownames(.)) %>%
                gather(sample, count, FP11w0L:FP9w2preR) %>%
                inner_join(temp$samples %>%
                                   mutate(sample = rownames(.))) %>%
                inner_join(myosin_family) %>%
                filter(symbol %in% c("MYH1", "MYH2", "MYH7")) %>%
                dplyr::select(subject, time, sets, gene, symbol, count) %>%
                group_by(subject, time, sets) %>%
                mutate(sum_count = sum(count)) %>%
                ungroup() %>%
                mutate(relative_count = (count/sum_count)*100, 
                       method = names(dge_lists[i])) 
        
        
        results[[i]] <- myhc_all
        
        
}

# Results for all lists are compiled
rna_myh <- bind_rows(results) %>%
        mutate(fibertype = if_else(symbol == "MYH1", "type2x", 
                                   if_else(symbol == "MYH2", "type2a", "type1"))) %>%
        mutate(time = if_else(time == "w2pre", "w2", as.character(time))) %>%
        dplyr::select(subject, time, sets, fibertype,method, rna_percentage = relative_count) %>%
        print()




# Load fiber type data (immunohistochemistry)
## Load data ##

fib <- read_excel("./data/study-1/immuno/fibertype_checked.xlsx") 
leg <- read.csv("./data/study-1/oneThreeSetLeg.csv", sep=";")

leg <- leg %>%
        gather(sets, leg, multiple:single) %>%
        filter(include == "incl")


# Type 2 AX counted as 0.5 A, 0.5 X

fib.trans <- fib %>%
        mutate(Type2A = as.integer(type2a + 0.5 * type2ax),
               Type2X = as.integer(type2x + 0.5 * type2ax),
               Type1 = type1,
               type2x = type2x,
               type2ax = type2ax,
               time = ifelse(timepoint == 1, "w0", ifelse(timepoint %in% c(2,3), "w2", "w12"))) %>%
        dplyr::select(subject, time, timepoint, leg, Type1, Type2A, Type2X, type2x, type2ax) %>%
        inner_join(leg) %>%
        gather(fibertype, counts, Type1:type2ax) %>%
        group_by(subject, sex, sets, time, timepoint, fibertype) %>%
        summarise(counts = sum(counts)) %>%
        spread(fibertype, counts) %>%
        rowwise() %>%
        mutate(sum = as.integer(round(Type1 + Type2A + Type2X, 0)),
               type2x_pure = type2x / sum, 
               type2x_hybrids = type2ax / sum,
               type1 = Type1 / sum,
               type2a = Type2A / sum,
               type2x = Type2X / sum) %>%
        ungroup() %>%
        dplyr::select(subject, sex, sets, time, timepoint, sum, type1, type2a, type2x) %>%
        mutate(time = factor(time, levels = c("w0", "w2", "w12")),
               sets = factor(sets, levels = c("single", "multiple"))) %>%
        gather(fibertype, percentage, type1:type2x) %>%
        mutate(sample = paste0(subject, sets, timepoint),
               leg = paste0(subject, sets)) %>%
        group_by(subject, sets, time, fibertype) %>%
        summarise(prot_percentage = mean(percentage, na.rm = TRUE) * 100) %>%
        ungroup() %>%
        print()


# Combine data sets 

myh_combined <- fib.trans %>%
        dplyr::select(subject, sets, time,fibertype, prot_percentage) %>%
        inner_join(rna_myh %>% ungroup()) %>%
        print()



myh_cor_points <- myh_combined %>%
        # filter(sets == "multiple") %>%
        mutate(time = factor(time, levels = c("w0", "w2", "w12"), 
                             labels = c("Week 0", "Week 2", "Week 12"))) %>%
        mutate(method = factor(method,
                               levels = c("star", "hisat", "rsem", "salmon", "kallisto"), 
                               labels = c("STAR", "HISAT2", "RSEM", "Salmon", "kallisto")), 
               fibertype = factor(fibertype, 
                                  levels = c("type1", "type2a", "type2x"), 
                                  labels = c("Type I", "Type IIA", "Type IIX"))) %>%
        ggplot(aes(rna_percentage, prot_percentage, color = fibertype, fill = fibertype)) + 
        geom_abline(intercept = 0, slope = 1, color = "gray40", lty = 2) +
        geom_point(aes(color = NULL), color = "black", shape = 21, alpha = 0.4, size = 0.8) +
        facet_grid(method ~ time) + 
        geom_smooth(method = "lm", color = "black", size = 0.4, se = FALSE) + 
        
        scale_y_continuous(breaks = c(0, 25, 50, 75, 100), 
                           labels = c(0, "", 50, "", 100), 
                           limits = c(-1,100), 
                           expand = c(0,0)) +
        scale_x_continuous(breaks = c(0, 25, 50, 75, 100), 
                           labels = c(0, "", 50, "", 100), 
                           limits = c(-2,102), 
                           expand = c(0,0)) +
        
        scale_fill_manual(values = group.study.color[c(1,5,6)]) +
        
        scale_color_manual(values = group.study.color[c(1,5,6)]) +
        
        
        labs(x = "RNA-seq gene-family percentage", 
             y = "Immunohistochemistry protein-family percentage") +
        
        dissertation_theme()  +
        theme(strip.background = element_rect(fill = "white", 
                                              color = "white"),
              strip.text = element_text(size = 7), 
              legend.position = "none")


myh_cor_coef <- myh_combined %>%
        # filter(sets == "single") %>%
        mutate(time = factor(time, levels = c("w0", "w2", "w12"), 
                             labels = c("Week 0", "Week 2", "Week 12"))) %>%
        group_by(time, method, fibertype) %>%
        summarise(cor = cor.test(rna_percentage, prot_percentage)$estimate, 
                  p.val = cor.test(rna_percentage, prot_percentage)$p.value, 
                  ucl = cor.test(rna_percentage, prot_percentage)$conf.int[1], 
                  lcl = cor.test(rna_percentage, prot_percentage)$conf.int[2]) %>%
        ungroup() %>%
        mutate(method = fct_reorder(method, .x = cor)) %>%
        
        mutate(method = factor(method, labels = c("STAR", "HISAT2", "RSEM", "Salmon", "kallisto")), 
               fibertype = factor(fibertype, 
                                  levels = c("type1", "type2a", "type2x"), 
                                  labels = c("Type I", "Type IIA", "Type IIX"))) %>%
        
        
        ggplot(aes(method, cor, fill = fibertype)) + 
        geom_errorbar(aes(ymin = lcl, ymax = ucl), width = 0, size = 0.4) +
        geom_point(position = position_dodge(width = 0.5), size = 1.5, shape = 21) + 
        facet_grid(fibertype ~ time) +
        
        scale_y_continuous(breaks = c(-0.2, 0, 0.2, 0.4, 0.6, 0.8, 1), 
                           labels = c("",  0,  "",  0.4, "",  0.8, ""), 
                           limits = c(-0.2, 1), 
                           expand = c(0,0)) + 
        
        scale_fill_manual(values = group.study.color[c(1,5,6)]) +
        
        scale_color_manual(values = group.study.color[c(1,5,6)]) +
        
        labs(y = "Correlation coeficient \U00B1 95% CI") +
        
        dissertation_theme()  +
        theme(strip.background = element_rect(fill = "white", 
                                              color = "white"), 
              strip.text = element_text(size = 7),
              legend.position = "none", 
              axis.title.x = element_blank()) +
        scale_x_discrete(guide = guide_axis(n.dodge = 2))





##### FIGURE 2 COMPOSITION #####




myhc_rnaseq_fig <- plot_grid(myh_cor_points, myh_cor_coef, ncol = 1, rel_heights =  c(0.55, 0.45)) +
        
        
        
        draw_plot_label(label=c("a", "b"),
                        x = c(0.02, 0.02), 
                        y = c(0.97, 0.47),
                        hjust=.5, vjust=.5, size = label.size)



saveRDS(myhc_rnaseq_fig, "./data/derivedData/study1-myhc-rnaseq-validation/myhc_rnaseq_fig.RDS")






