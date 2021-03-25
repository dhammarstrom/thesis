##-------------------------------------
## study2-protein-analyses.R
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



library(brms)
source("./R/libraries.R")
source("./R/themes.R")




s6.m1 <- readRDS( "./data/study-2/western-analysis/s6_tx_m1.RDS")
ubf.m1 <- readRDS( "./data/study-2/western-analysis/ubf_tx_m1.RDS")



tx_contrasts <-  c("txS1_int_train -     txpre_int_train =      txS1_con_train",
                   # Interaction effects post-training
                   "txpost_int_train -     txpre_int_train =      txpost_con_train",
                   # Interaction effects detraining (post1w)
                   "txpost_int_detrain -     txpre_int_train =      txpost_con_train",
                   # Whitin group fold change
                   "txS1_int_train =     txpre_int_train",     
                   "txpost_int_train =     txpre_int_train",     
                   "txpost_int_detrain =     txpre_int_train",     
                   "txS1_con_train = 0",
                   "txpost_con_train = 0")



h.s6 <- hypothesis(s6.m1, tx_contrasts)
h.ubf <- hypothesis(ubf.m1, tx_contrasts)

# Create a data frame for the results 
tx_results_western <- data.frame(target = rep(c("rpS6","UBF"), each = 8), 
                                 comparison = rep(c("inter:S1",
                                                    "inter:post", 
                                                    "inter:post1w", 
                                                    
                                                    "int_S1", 
                                                    "int_post", 
                                                    "int_post1w", 
                                                    
                                                    "con_S1",
                                                    "con_post"), 2)) %>%
        cbind(data.frame(rbind(h.s6$hypothesis[,-1], 
                               h.ubf$hypothesis[,-1]) )) %>%
        print()


#### Fold change plot rpS6 and UBF 


prot_fold_change <- tx_results_western %>%
        filter(comparison %in% c("int_S1", "int_post", "int_post1w", "con_S1", "con_post")) %>%
        mutate(Estimate = exp(Estimate), 
               CI.Lower= exp(CI.Lower), 
               CI.Upper = exp(CI.Upper)) %>%
        
        separate(comparison, into = c("group", "time")) %>%
        
        mutate(detrain = if_else(time == "post1w", "detrain", "train"),
               time = if_else(time == "post1w", "post", time), 
               time = factor(time, levels = c("S1", "post"), 
                             labels = c("Session 1", "Post-\ntraining")), 
               group = paste0(group, "_",  detrain), 
               group = factor(group, levels = c("con_train", 
                                                "int_train", 
                                                "int_detrain"), 
                              labels = c("Control", "Training", 
                                         "Training\n+de-training"))) %>%
        
        ggplot(aes(time, Estimate, fill = group)) + 
        
        
        geom_hline(yintercept = 1, lty = 2, color = "gray80") +
        
        geom_errorbar(aes(ymin = CI.Lower, ymax = CI.Upper), 
                      position = position_dodge(width = 0.2),
                      size = 0.3,
                      width = 0) +
        geom_point(position = position_dodge(width = 0.2), shape = 21) + 
        
        facet_wrap(~ target , ncol = 1, strip.position = "top", scales = "free") +
        # scale_y_continuous(limits = c(0.5, 3), breaks = c(1, 2, 3), 
        #                    expand = c(0,0)) +
        
        scale_fill_manual(values = c(group.study.color[1], group.study.color[3],group.study.color[5])) +
        
        
        labs(y = "Fold change from Baseline") +
        dissertation_theme()  +
        theme(strip.placement = "top",
              strip.background = element_rect(fill = "white", color = "white"), 
              strip.text = element_text(size = 7, hjust = 0),
              legend.position = "left",
              legend.title = element_blank(),
              legend.text = element_text(size = 7, margin = margin(t = 0.1, b= 0.1,r = 0.1,l = 0.1, unit = "pt")),
              legend.key.size = unit(0, "cm"),
              legend.margin = margin(-0.5, 0, -0.5, 0, "cm"),
              legend.spacing.x = unit(0.1, 'cm'),
              legend.direction = "vertical",
              axis.title.x = element_blank()) 



### Comparison plot 

prot_interaction <- tx_results_western %>%
        filter(!(comparison %in% c("int_S1", "int_post", "int_post1w", "con_S1", "con_post"))) %>%
        mutate(Estimate = exp(Estimate), 
               CI.Lower= exp(CI.Lower), 
               CI.Upper = exp(CI.Upper)) %>%
        mutate(comparison =  gsub("inter:", "", comparison), 
               comparison = factor(comparison, levels = c("S1", "post", "post1w"), 
                                   labels = c("Session 1", "Post-training", 
                                              "Post-training\n+de-training")), 
               comparison = fct_rev(comparison), 
               robust = if_else(CI.Lower  > 1, "robust", "notrobust")) %>%
        ggplot(aes(Estimate, comparison, alpha = robust)) + 
        labs(x = "Fold change\ncompared to Control") +
        
        geom_vline(xintercept = 1, color = "gray50", lty = 2) +
        geom_errorbarh(aes(xmin = CI.Lower, xmax = CI.Upper), height = 0, size = 0.3) + 
        geom_point(shape = 24, fill = group.study.color[3]) +

        facet_wrap( ~ target, ncol = 1, strip.position = "top") +
        
        scale_alpha_manual(values = c(0.5, 1)) +
        
        dissertation_theme() +
        theme(strip.background = element_rect(color = "white", fill = "white"), 
              
              strip.text = element_text(size = 7, hjust = 0.1),
              axis.title.y = element_blank()  , 
              axis.text.y = element_text(size = 7),
              legend.position = "none")





### Insert image of western blotting  




saveRDS(list(prot_fold_change = prot_fold_change, 
             prot_interaction = prot_interaction), "./data/derivedData/study2-protein-analyses/protein_analyses.RDS")



