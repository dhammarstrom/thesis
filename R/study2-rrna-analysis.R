##-------------------------------------
## study2-rrna-analysis.R
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

########### rRNA vs. control group ############




qpcr_res_int_con <- readRDS("./data/study-2/qpcr-analysis-bayes2/qpcr_res_int_con.RDS")


complete_rrna_comp <- qpcr_res_int_con %>%
        filter(model == "tissue") %>%
        dplyr::select(target, comparison = contrast, estimate, lower.CL , upper.CL) %>%
        
        filter(!(target %in% c("UBTF F4R4",
                               "UBTF F6R6", 
                               "rpS6 F2R2", 
                               "MyHC1 F1R1", 
                               "MyHC2A F5R5", 
                               "MyHC2X F5R5"))) %>%
        print()






# Create an annotation data frame

anno.df <- data.frame(target = unique(complete_rrna_comp$target), 
                      label = c("rRNA 18S",
                                "rRNA 28S",  
                                "rRNA 45S ETS",
                                "rRNA 45S ITS", 
                                "rRNA 47S ETS", 
                                "rRNA 5.8S",
                                "rRNA 5S" ), 
                      time = factor("S1", levels = c("S1", "post")), 
                      estimate = 3.5) %>%
        mutate(target = factor(target, levels = c("rRNA18S F2R2",
                                                  "rRNA5.8S F2R2",
                                                  "rRNA28S F2R2",
                                                  "rRNA5S F3R3", 
                                                  "rRNA45SITS F12R12",
                                                  "rRNA45S F5R5",     
                                                  "rRNA47S F1R1")), 
               comparison = "S1",
               comparison = factor(comparison, levels = c("S1", "post", "post1w"), 
                                   labels = c("Session 1", 
                                              "Post-training", 
                                              "Post-trainin \n+ De-training")),
               target = fct_rev(target))



interaction_effects <- complete_rrna_comp %>%
        
        
        filter(comparison %in% c("inter:S1", "inter:post", "inter:post1w","")) %>%
        
        mutate(target = factor(target, levels = c("rRNA18S F2R2",
                                                  "rRNA5.8S F2R2",
                                                  "rRNA28S F2R2",
                                                  "rRNA5S F3R3", 
                                                  "rRNA45SITS F12R12",
                                                  "rRNA45S F5R5",     
                                                  "rRNA47S F1R1")),
               
               comparison = gsub("inter:", "", comparison), 
               comparison = factor(comparison, levels = c("S1", "post", "post1w"), 
                                   labels = c("Session 1", 
                                              "Post-training", 
                                              "Post-trainin \n+ De-training")),
               estimate = exp(estimate), 
               lower.CL = exp(lower.CL), 
               upper.CL = exp(upper.CL), 
               robust = if_else(lower.CL  > 1, "robust", "notrobust")) %>%
        
        ggplot(aes(estimate, 
                   target, color = robust)) + 
        
    #   geom_text(data = anno.df %>%
    #                     mutate(target = factor(target, levels = c("rRNA18S F2R2",
    #                                                               "rRNA5.8S F2R2",
    #                                                               "rRNA28S F2R2",
    #                                                               "rRNA5S F3R3", 
    #                                                               "rRNA45SITS F12R12",
    #                                                               "rRNA45S F5R5",     
    #                                                               "rRNA47S F1R1"))),
    #             aes(estimate - 2, target,  label = label, color = NULL), 
    #             position = position_nudge(y = 0.25), 
    #             hjust = 0,
    #             size = 2.2) +
    #   
        
        labs(x = "Fold change compared to Control") +
        
        geom_vline(xintercept = 1, color = "gray85", lty = 2) +
        
    geom_errorbarh(aes(xmin = lower.CL, xmax = upper.CL), height = 0, size = 0.3) + 
        
    geom_point(shape = 24, fill = group.study.color[3]) +

        facet_wrap(  comparison ~ .) +
        
        scale_color_manual(values = c("gray50", "gray10")) +
        
        dissertation_theme() +
        theme(strip.background = element_rect(color = "white", fill = "white"), 
              strip.text = element_text(size = 7),
              axis.title.y = element_blank()  , 
              axis.line.y = element_blank(), 
              axis.ticks.y = element_blank(),
              axis.text.y = element_blank(), 
              legend.position = "none")




fold_changes <- complete_rrna_comp %>%
        
        filter(!(comparison %in% c("inter:S1", "inter:post", "inter:post1w",""))) %>%
        separate(comparison, into = c("group", "time"), sep = "_") %>%
        
        
        
        mutate(target = factor(target, levels = c(
                "rRNA18S F2R2",
                "rRNA5.8S F2R2",
                "rRNA28S F2R2",
                "rRNA5S F3R3", 
                "rRNA45SITS F12R12",
                "rRNA45S F5R5",     
                "rRNA47S F1R1")),
               target = fct_rev(target),
               
               estimate = exp(estimate), 
               lower.CL = exp(lower.CL), 
               upper.CL = exp(upper.CL), 
               group = if_else(time == "post1w", "int_detrain", group),
               
               group = factor(group, levels = c("con", "int", "int_detrain"), 
                              labels = c("Control", "Training", "Training\n+De-training")),
               
               time = if_else(time == "post1w", "post", time),
               time = factor(time, levels = c("S1", "post"))) %>%
        
        ggplot(aes(time, estimate, fill = group)) + 
        
        geom_text(data = anno.df, aes(time, estimate, label = label, fill = NULL), 
                  position = position_nudge(x = -0.5), 
                  hjust = 0,
                  size = 2.2) +
        geom_hline(yintercept = 1, lty = 2, color = "gray80") +
        
        # geom_bar(stat = "identity", position = position_dodge(width = 0.3), width = 0.15) + 
        geom_errorbar(aes(ymin = lower.CL, 
                          ymax = upper.CL), 
                      position = position_dodge(width = 0.3), 
                      width = 0, 
                      size = 0.3) + 
        
        geom_point(position = position_dodge(width = 0.3), shape = 21) +
        
        scale_y_continuous(limits = c(0.5, 3.8), breaks = c(1, 2, 3)) +
        scale_x_discrete(limits = c( "S1", "post"), labels = c("Session 1", "Post-\ntraining")) +
        
        scale_fill_manual(values = c(group.study.color[1], group.study.color[2],group.study.color[5])) +
        
        
        labs(y = "Fold change from Baseline") +
       dissertation_theme() +
        
        theme(strip.background = element_blank(), 
              strip.text = element_blank(), 
              legend.position = "top", 
              legend.title = element_blank(),
              legend.text = element_text(size = 7, margin = margin(t = 1, b= 1,r = 0.1,l = 0.1, unit = "pt")),
              legend.key.size = unit(0, "cm"),
              legend.margin = margin(0.5, 0, 0.5, 0, "cm"),
              legend.spacing.x = unit(0.1, 'cm'),
              legend.direction = "horizontal",
              axis.title.x = element_blank()) +
        
        facet_grid(target ~ .)



### Combine 



rrna <- plot_grid(plot_grid(NULL, fold_changes, NULL, ncol = 1, rel_heights = c(0, 1, 0.065)),
          plot_grid(NULL, interaction_effects, NULL, ncol = 1, rel_heights = c(0.1, 1, 0.01)),
          ncol = 2, rel_widths = c(0.5, 0.4))






### Total RNA ctrl vs exp ##############

comp_rna <- readRDS("./data/study-2/total-rna-analysis/comp_rna_bayes.RDS")

tot_rna_interaction <- comp_rna %>%
        mutate(comparison = time_group, 
               target = "totalRNA") %>%
        dplyr::select(target, comparison, estimate:upper.CL) %>%
        filter(comparison %in% c("inter:S1", "inter:post", "inter:post1w","")) %>%
        
        mutate(
                comparison = gsub("inter:", "", comparison), 
                comparison = factor(comparison, levels = c("S1", "post", "post1w"), 
                                    labels = c("Session 1", 
                                               "Post-training", 
                                               "Post-trainin \n+ De-training")),
                estimate = exp(estimate), 
                lower.CL = exp(lower.CL), 
                upper.CL = exp(upper.CL), 
                robust = if_else(lower.CL  > 1, "robust", "notrobust")) %>%
        
        ggplot(aes(estimate, 
                   target, color = robust)) + 
        
        
        labs(x = "Fold change compared to Control") +
        
        geom_vline(xintercept = 1, color = "gray85", lty = 2) +
         geom_errorbarh(aes(xmin = lower.CL, xmax = upper.CL), height = 0, size = 0.3) + 
        geom_point(shape = 24, fill = group.study.color[3]) +
       
        facet_wrap(  comparison ~ .) +
        
        scale_color_manual(values = c("gray50", "gray10")) +
        scale_x_continuous(limits = c(0.8, 2), 
                           expand = c(0, 0), 
                           breaks = c(0.8, 1, 1.2, 1.4, 1.6, 1.8, 2), 
                           labels = c("",  1, "",   1.4, "", 1.8, "")) +
        
        dissertation_theme() +
        theme(strip.background = element_rect(color = "white", fill = "white"), 
              strip.text = element_text(size = 7),
              axis.title.y = element_blank()  , 
              axis.line.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.text.y = element_blank(), 
              legend.position = "none")



tot_rna_fold_change <- comp_rna %>%
        mutate(comparison = time_group, 
               target = "totalRNA") %>%
        dplyr::select(target, comparison, estimate:upper.CL) %>%
        
        filter(!(comparison %in% c("inter:S1", "inter:post", "inter:post1w",""))) %>%
        separate(comparison, into = c("group", "time"), sep = "_") %>%
        
        mutate(  estimate = exp(estimate), 
                 lower.CL = exp(lower.CL), 
                 upper.CL = exp(upper.CL), 
                 group = if_else(time == "post1w", "int_detrain", group),
                 time = if_else(time == "post1w", "post", time),
                 time = factor(time, levels = c("S1", "post"))) %>%
        
        ggplot(aes(time, estimate, fill = group)) + 
        
        geom_hline(yintercept = 1, lty = 2, color = "gray80") +
        
        # geom_bar(stat = "identity", position = position_dodge(width = 0.3), width = 0.15) + 
        geom_errorbar(aes(ymin = lower.CL, 
                          ymax = upper.CL), 
                      position = position_dodge(width = 0.3), 
                          size = 0.3,
                      width = 0) + 
        
        geom_point(position = position_dodge(width = 0.3), shape = 21) +
        
        scale_x_discrete(limits = c( "S1", "post"), labels = c("Session 1", "Post-\ntraining")) +
        scale_y_continuous(limits = c(0.8, 2), breaks = c(1, 1.2, 1.4, 1.6, 1.8, 2), 
                           labels = c(1, 1.2, 1.4, 1.6, 1.8, 2), 
                           expand = c(0, 0)) +
        
        scale_fill_manual(values = c(group.study.color[1], group.study.color[2],group.study.color[5])) +
        labs(y = "Fold change<br>from Baseline") +
        dissertation_theme() +
        
        theme(strip.background = element_blank(), 
              strip.text = element_blank(), 
              axis.title.y = element_markdown(),
              legend.position = "none", 
              axis.title.x = element_blank()) 


tot_rna <- plot_grid(NULL, 
          
          plot_grid(NULL, tot_rna_fold_change, NULL, ncol = 1, rel_heights = c(0.05, 1, 0.05)),
          
          tot_rna_interaction, 
          
          ncol = 3, rel_widths = c(0.05, 0.4, 0.4))


### Combine all plots in training vs. control ###########



rrna_totrna_vs_control <- plot_grid(tot_rna, rrna, 
          ncol = 1, 
          rel_heights = c(0.2, 0.8))



saveRDS(rrna_totrna_vs_control, "./data/derivedData/study2-rrna-analysis/rrna_totrna_vs_control.RDS")



#### Time course plots ############################



### time course rrna ##############

# "trace plots" of each qpcr-estimate


rrna_data <- readRDS("./data/study-2/qpcr-analysis-bayes2/time_course_qpcr.RDS")

# Get Estimates and create a data set for statistical robustness indicators
# Statistical robust is defined as when 0 is not within the 95% credible interval


rrna_timecourse <- rrna_data$estimated_means %>%
        filter(model == "tissue", 
               target %in% c("rRNA18S F2R2",
                             "rRNA5.8S F2R2",
                             "rRNA28S F2R2",
                             "rRNA5S F3R3", 
                             "rRNA45SITS F12R12",
                             "rRNA45S F5R5",     
                             "rRNA47S F1R1")) %>%
        
        mutate(target = factor(target, levels = c("rRNA18S F2R2",
                                                  "rRNA5.8S F2R2",
                                                  "rRNA28S F2R2",
                                                  "rRNA5S F3R3", 
                                                  "rRNA45SITS F12R12",
                                                  "rRNA45S F5R5",     
                                                  "rRNA47S F1R1")),
               target = fct_rev(target), 
               time.c = gsub("S", "", time), 
               time.c = as.numeric(if_else(time.c == "post1w", "12", time.c))) %>%
        print()


diffs <- rrna_data$estimated_diff %>%
        
        filter(model == "tissue") %>%
        
        mutate(target = factor(target, levels = c("rRNA18S F2R2",
                                                  "rRNA5.8S F2R2",
                                                  "rRNA28S F2R2",
                                                  "rRNA5S F3R3", 
                                                  "rRNA45SITS F12R12",
                                                  "rRNA45S F5R5",     
                                                  "rRNA47S F1R1")),
               target = fct_rev(target), 
               time.c = gsub("S", "", time), 
               time.c = as.numeric(if_else(time.c == "post1w", "12", time.c))) %>%
        filter(robust == "robust") %>%
        mutate(cond = "var") %>%
        left_join(rrna_timecourse %>% 
                          group_by(target, time) %>%
                          summarise(emmean = max(emmean))) %>% 
        mutate(emmean = emmean * 1.05, 
               sign = "\U0002A") %>%
        print()


targets <- unique(rrna_timecourse$target)


plots <- list()


for(i in 1:length(targets)) {
        
        dat <- rrna_timecourse %>%
                filter(target == targets[i]) 
        
        
        
        diffs2 <- diffs %>%
                filter(target == targets[i])
        
        plot  <-   dat %>%
                filter(time != "post1w" )  %>%
                ggplot(aes(time.c, emmean, color = cond)) + 
                geom_line() +
                geom_text(data = diffs2, 
                          aes(time.c, emmean, label = sign, color = NULL), 
                          show.legend = FALSE) +
                
                geom_point(data = filter(dat, time == "post1w"), 
                           position = position_nudge(x = 1)) +
                
                scale_color_manual(values = c(group.study.color[1], group.study.color[5])) +
                
                scale_y_continuous(limits = c(floor(min(dat$emmean, na.rm = TRUE)), 
                                              ceiling(max(dat$emmean, na.rm = TRUE))),
                                   breaks = c(floor(min(dat$emmean, na.rm = TRUE)),
                                              ceiling(max(dat$emmean, na.rm = TRUE)) + ((floor(min(dat$emmean, na.rm = TRUE)) - ceiling(max(dat$emmean, na.rm = TRUE)))/2), 
                                              ceiling(max(dat$emmean, na.rm = TRUE)))) +
                
                labs(x = "Session", 
                     y = "Estimated log-abundance per tissue weight (AU)")  +
                dissertation_theme() +
                
                ggtitle(anno.df[anno.df$target == targets[i], 2]) + 
                
                theme(strip.background = element_blank(), 
                      strip.text = element_blank(), 
                      plot.title = element_text(size = 7),
                      
                      legend.position = "none", 
                      axis.title = element_blank(), 
                      axis.text.x = element_blank(), 
                      axis.line.x = element_blank(), 
                      axis.ticks.x = element_blank()) 
        
        
        plots[[i]]  <- plot
        
}


names(plots) <- unique(rrna_timecourse$target)



## Make axis label for subplots
library(grid)
library(gridExtra)

axis_lab <- textGrob("Estimated log-counts per tissue weight (AU)", 
                     gp=gpar( fontsize=7),
                     rot=90)

### create a plot with only x axis

x_axis <- data.frame(Session = c(0, 3, 6, 9, 12), 
                     y       = c(0, 0, 0, 0, 0)) %>%
        ggplot(aes(Session, y)) + geom_blank() + 
        dissertation_theme() + 

        
        scale_x_continuous(breaks = c(0, 3, 6, 9, 12), 
                           limits = c(0, 12.5)) +
        
        labs(x = "Session") +
        
        theme(strip.background = element_blank(), 
              strip.text = element_blank(), 
              legend.position = "none", 
              axis.title.y = element_blank(), 
              axis.text.y = element_blank(), 
              axis.line.y = element_blank(), 
              axis.ticks.y = element_blank()) 



time_traces <-  plot_grid(axis_lab, 
                          plot_grid(plots[[5]], 
                                    plots[[3]], 
                                    plots[[4]], 
                                    plots[[7]], 
                                    plots[[2]], 
                                    plots[[6]], 
                                    plots[[1]], 
                                    x_axis,
                                    
                                    rel_heights = c(rep(1, 7), 0.5),
                                    
                                    align = "v",
                                    
                                    ncol = 1), ncol = 2, rel_widths = c(0.1, 1))






cond_eff_rna_tc <- readRDS("./data/study-2/total-rna-analysis/cond_eff_rna_tc.RDS")






rna_tc_fig <-  cond_eff_rna_tc %>%
        mutate(cond = factor(cond, levels = c("const", "var"), 
                             labels = c("Constant volume", 
                                        "Variable volume")), 
               cond = fct_rev(cond)) %>%
        ggplot(aes(time.c, estimate, color = cond, fill = cond)) +
        geom_line() + 
        geom_ribbon(aes(ymin = lower, ymax = upper, color = NULL), alpha = 0.2) +
        
        labs(x = "Session", 
             y = "RNA     ng mg<sup>-1</sup>") +
        
        scale_x_continuous(limits = c(0, 12), expand = c(0,0), 
                           breaks = c(0, 3, 6, 9, 12)) +
        
        scale_y_continuous(limits = c(250, 650), expand = c(0,0), 
                           breaks = c(250, 300, 350, 400, 450, 500, 550, 600, 650), 
                           labels = c("",  300, "",  400, "",  500, "" , 600, "")) +
        
        scale_color_manual(values = c(group.study.color[5], group.study.color[1]), guide = NULL) +
        scale_fill_manual(values = c(group.study.color[5], group.study.color[1])) +
        
        dissertation_theme() +
        theme(legend.position = c(0.7, 0.15), 
              legend.key.size = unit(0.3, "cm"), 
              legend.title = element_blank(),
              axis.title.y = element_markdown(), 
        ) 



rna_abundance_exp_group <- plot_grid(plot_grid(NULL, rna_tc_fig, NULL, ncol = 1, rel_heights = c(0.7, 1, 0.7)), 
          time_traces, 
          ncol = 2, 
          rel_widths = c(0.7, 0.5))


saveRDS(rna_abundance_exp_group, "./data/derivedData/study2-rrna-analysis/rna_abundance_exp_group.RDS")



