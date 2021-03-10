##-------------------------------------
## study1-pairwise-variation-rnaseq.R
##
## Title: Validation of RNA-seq tools using pairwise variations
## Purpose: 
## Author:
##
##
##
##-------------------------------------




# Source functions and settings

source("./R/libraries.R")
source("./R/themes.R")






paired_sd_data <- readRDS(file = "./data/study-1b/data/derivedData/paired_sd_data.RDS")



## Plot results from analysis ##

sd_a <- paired_sd_data %>%
        mutate(method = factor(method, levels = c("star", "hisat", "salmon", "kallisto", "rsem"), 
                               labels = c("STAR", "HISAT2", "Salmon", "kallisto", "RSEM"))) %>%
        filter(hk == TRUE) %>%
        ggplot(aes(A, SD, color = method, fill = method)) + 
        
        geom_smooth(se = FALSE) +
        dissertation_theme() + 
        scale_y_continuous(breaks = c(0, 0.5, 1, 1.5)) +
        scale_x_continuous(breaks = c(0, 5, 10, 15, 20), 
                           labels = c(0, "", 10, "", 20)) + 
        coord_cartesian(ylim=c(0, 1.5),
                        xlim = c(-1, 15), 
                        expand = 0) +
        
        scale_color_manual(values = group.study.color) +
        scale_fill_manual(values = group.study.color) +
        theme(legend.position = c(0.6, 0.8), 
              legend.title = element_blank(), 
              legend.key.height = unit(0.45, "cm")) +
        labs(y = "Average log2-difference\nbetween replicates", 
             x = "Average log2-count (A)")





ncount_medium <- paired_sd_data %>%
        filter(hk == TRUE) %>%
        mutate(abundance.strata = if_else(A <= 1, "low", 
                                          if_else(A>1 & A <6, "medium", 
                                                  if_else(A >= 6, "high", "undetermined")))) %>%
        mutate(abundance.strata = factor(abundance.strata, levels = c("low", "medium", "high"), 
                                         labels = c("Low abundace\n(A<1)",
                                                    "Medium abundace\n(1 \u2264 A \u2265 6)",
                                                    "High abundace\n(A > 6)"))) %>%
        group_by(method, abundance.strata) %>%
        filter(abundance.strata == "Medium abundace\n(1 \u2264 A \u2265 6)") %>%
        summarise(n = n()) %>%
        print()



ann_text <- data.frame(abundance.strata = factor( "Medium count\n(1 \u2264 A \u2265 6)", 
                                                  levels = c("Low abundace\n(A<1)",
                                                             "Medium count\n(1 \u2264 A \u2265 6)",
                                                             "High abundace\n(A > 6)")), 
                       
                       x = c(1, 2, 3, 4, 5), 
                       method = c("rsem", "kallisto", "salmon", "star", "hisat"),
                       n = ncount_medium$n + 30,
                       lab = c("RSEM", 
                               "kallisto", 
                               "Salmon", 
                               "STAR", 
                               "HISAT2"))



strata.labels <- data.frame(x = c(1, 1, 1), 
                            y = c(1, 2, 3), 
                            label =  c("High count\n(A > 6)",
                                       "Medium count\n(1 \u2264 A \u2265 6)", 
                                       "Low count\n(A<1)")) %>%
        ggplot(aes(x, y, label = label)) + geom_text(size = 2.2) + 
        scale_y_continuous(limits = c(0, 3.3)) +
        theme_void() 




n_counts_fig <- paired_sd_data %>%
        filter(hk == TRUE) %>%
        mutate(abundance.strata = if_else(A <= 1, "low", 
                                          if_else(A>1 & A <6, "medium", 
                                                  if_else(A >= 6, "high", "undetermined")))) %>%
        mutate(abundance.strata = factor(abundance.strata, levels = c("low", "medium", "high"), 
                                         labels = c("Low count\n(A<1)",
                                                    "Medium count\n(1 \u2264 A \u2265 6)",
                                                    "High count\n(A > 6)"))) %>%
        group_by(method, abundance.strata) %>%
        summarise(n = n()) %>%
        ungroup() %>%
        ggplot(aes(method, n, fill = method)) + 
        geom_bar(stat = "identity") +
        
        
        scale_fill_manual(values = group.study.color) +
        scale_y_continuous(breaks = c(0, 3000, 6000), 
                           labels = c(0, " ", 6000), 
                           limits = c(0, 6000), 
                           expand = c(0, 0), 
                           name = "n gene counts") +
        facet_grid(abundance.strata ~ .) + 
        geom_text(data = ann_text, aes(label = lab), hjust = 0, size = 2) +
        
        dissertation_theme() +
        theme(strip.text.y = element_blank(), 
              strip.background = element_blank(), 
              legend.position = "none", 
              axis.text.y = element_blank(),
              axis.text.x = element_text(hjust = c(0, 1, 1)), # This produces a warnings but works
              axis.title.y = element_blank(), 
              axis.ticks.y = element_blank()) +
        
        coord_flip()

# Combine plots to get abundance strata labels on the y axis ...
n_counts_fig <- plot_grid(strata.labels, n_counts_fig, nrow = 1, rel_widths = c(0.2, 0.8))


rna_seq_validation_pairwise <- plot_grid(n_counts_fig, sd_a, rel_widths = c(0.5, 0.5), 
          ncol = 2)



saveRDS(sd_a, "./data/derivedData/study1-pairwise-variation-rnaseq/pairwise_fig.RDS")





