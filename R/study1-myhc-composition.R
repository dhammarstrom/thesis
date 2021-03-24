##-------------------------------------
## study1-myhc-composition.R
##
## Title: MyHC distributions in Study I
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




# Source functions and settings

source("./R/libraries.R")
source("./R/themes.R")



## Load data ##

fib <- read_excel("./data/study-1/immuno/fibertype_checked.xlsx") 
leg <- read.csv("./data/study-1/oneThreeSetLeg.csv", sep=";")

leg <- leg %>%
        gather(sets, leg, multiple:single) %>%
        filter(include == "incl")


## Tranform data ##

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
        dplyr::select(subject, sex, sets, time, timepoint, sum, type1, type2a, type2x, type2x_pure, type2x_hybrids) %>%
        mutate(time = factor(time, levels = c("w0", "w2", "w12")),
               sets = factor(sets, levels = c("single", "multiple"))) %>%
        gather(fibertype, percentage, type1:type2x_hybrids) %>%
        mutate(sample = paste0(subject, sets, timepoint),
               leg = paste0(subject, sets)) 




#### Figure a -- Type IIX fibers distributed among time and participants ###########


overall_type2x <- fib.trans %>%
        
        group_by(subject, time, fibertype) %>%
        summarise(percentage = mean(percentage, na.rm = TRUE)) %>%
        ungroup() %>%

        
        filter(fibertype %in% c("type2x_hybrids", "type2x_pure")) %>%
        pivot_wider(names_from = "time", values_from = "percentage") %>%
        filter(complete.cases(.)) %>%
        pivot_longer(names_to = "time", 
                     values_to = "percentage", 
                     cols = w0:w12) %>%

        mutate(subject = fct_reorder(subject, percentage, max), 
               fibertype = factor(fibertype, levels = c("type2x_pure", 
                                                        "type2x_hybrids"), 
                                  labels = c("Type IIX", 
                                             "Type IIA/X")), 
               time = factor(time, levels = c("w0", "w2", "w12"), 
                             labels = c("Week 0", "Week 2", "Week 12")), 
               percentage = percentage * 100) %>%
        
        
        
        ggplot(aes(time, subject, fill = percentage)) + geom_raster() + 
        
        facet_wrap(~ fibertype, ncol = 2) +
               scale_fill_gradient(breaks = c(0, 15, 30), 
                            limits = c(0, 30), 
                            low = group.study.color[3], 
                            high = group.study.color[1]) +
        labs(fill = "% of\ntotal fibers", 
             y = "Participants") +
        
        guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5)) +
        dissertation_theme() + 
        theme(axis.title.y = element_text(size = 7), 
              
              axis.text.y = element_blank(), 
              axis.ticks = element_blank(), 
              axis.line = element_blank(),
              axis.text.x = element_text(size = 7, angle = 45, hjust = 0.9),
              axis.title.x = element_blank(),
              strip.background = element_blank(), 
              strip.text.x = element_text(size = 7, hjust = 0), 
              legend.title = element_text(size = 7),
              legend.text = element_text(size = 7),
              legend.key.height = unit(0.6, "cm"),
              legend.key.width = unit(0.2, "cm"),
              legend.position="right",
              legend.box="vertical") 
        
 
        
overall_type2x

### mRNA figure #####################



qpcr.model.estimates.full <- read_feather("./data/study-1/derivedData/qpcr_model_estimates.full")

qpcr.means.estimates.full <-  read_feather("./data/study-1/derivedData/qpcr_mean_estimates.full")



diff_2x_mRNA <- qpcr.model.estimates.full %>%
        filter(gene == "MyHC2X", 
               parameter %in% c("timepointw0:setsmultiple", 
                                "timepointw2pre:setsmultiple",
                                "timepointw2post:setsmultiple",
                                "timepointw12:setsmultiple")) %>%
        mutate(timepoint = gsub("timepoint", "", parameter), 
               timepoint = gsub(":setsmultiple", "", timepoint), 
               time = if_else(timepoint == "w0", 0, 
                              if_else(timepoint == "w2pre", 2, 
                                      if_else(timepoint == "w2post", 2.4, 12))), 
               sig = if_else(timepoint == "w0", "ns", "sig")) %>%
        
        ggplot(aes(time, exp(estimate), alpha = sig)) +
        
        geom_hline(yintercept = 1, lty = 2, color = "gray50", size = 0.5) +
        
        geom_errorbar(aes(ymin = exp(lower), ymax = exp(upper)), 
                      width = 0) +
 
        
        geom_point(shape = 21,  size = 2, fill = group.study.color[2]) +
        scale_y_continuous(limits = c(0.25, 1.5), 
                           expand = c(0,0), 
                           breaks = c(0.25, 0.5, 0.75, 1, 1.25, 1.5), 
                           labels = c("", 0.50, "", 1.00, "", 1.50)) +
        scale_x_continuous(breaks = c(0, 2, 12), 
                           limits = c(-0.5, 12.5), 
                           labels = c("Week 0", "Week 2", "Week 12")) +
        
        scale_alpha_manual(values = c(0.5, 1)) +
        
        dissertation_theme() +
        
        ylab("Moderate-/<br>Low-volume") + 
        theme(axis.title.x = element_blank(),
              axis.title.y = element_markdown(size = 7),
              axis.text.y = element_text(size = 7),
              axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
              legend.position = "none")






myhc2x_mrna_means <-  data.frame(qpcr.means.estimates.full) %>%
        mutate(time = if_else(timepoint == "w0", 0, 
                                                 if_else(timepoint == "w2pre", 2, 
                                                         if_else(timepoint == "w2post", 2.4, 12))), 
               timepoint = factor(timepoint, levels = c("w0", "w2pre", "w2post", "w12"),
                                  labels = c("Week 0", "Week 2\nPre-ex", "Week 2\nPost-ex", "Week 12")),
               sets = factor(sets, levels = c("single", "multiple"))) %>%
        filter(gene %in% c("MyHC2X")) %>%
        ggplot(aes(time, emmean, fill = sets, group = paste(sets, gene))) +

        geom_line(position = position_dodge(width = 0.2), lty = 2, color = "gray50") +
        geom_point(shape = 21, position = position_dodge(width = 0.2), size = 2.5) +
        scale_y_continuous(limits = c(-16.5, -12.5), 
                           breaks = seq(from = -16.5, to = -12.5, by = 2), 
                           expand = c(0, 0)) +
        
        scale_x_continuous(breaks = c(0, 2, 12), 
                           limits = c(-0.5, 12.5), 
                           labels = c("Week 0", "Week 2", "Week 12")) +
        
        # Annotations 
        annotate("text", x = c(9, 9), y = c(-13, -15.4), 
                 label = c("Low-volume", "Moderate-volume"), 
                 size = 2.7) +
        # From moderate volume to point
        annotate("curve",
                 x = 11.4,
                 xend = 11.9,
                 y = -15,
                 yend = -14.5,
                 curvature = 0.2,
                 color = "gray50",
                 arrow = arrow(length=unit(0.1,"cm"), type = "closed")) +
        
        # From low volume to point
        annotate("curve",
                 x = 11.25,
                 xend = 11.7,
                 y = -13.3,
                 yend = -13.65,
                 curvature = 0.2,
                 color = "gray50",
                 arrow = arrow(length=unit(0.1,"cm"), type = "closed")) +
        
        
        scale_fill_manual(values = c(group.study.color[5], group.study.color[2])) +
        dissertation_theme() +
        

        ylab("<i>MYH</i><br>Log-abundance") + 
        theme(axis.title.x = element_blank(),
              axis.title.y = element_markdown(size = 7),
              axis.text.y = element_text(size = 7),
              axis.text.x = element_blank(),
              axis.line.x = element_blank(),
              axis.ticks.x = element_blank(),
              legend.position = "none")


######## Percentage plot ############################


stacked_fibertypes_data <- fib.trans %>%
        filter(fibertype %in% c("type1", "type2a", "type2x")) %>%
        group_by(subject, sets, time, fibertype) %>%
        summarise(percentage = (sum((percentage * sum))/sum(sum)), 
                  sum = sum(sum)) %>%
        group_by(sets, time, fibertype) %>%
        summarise(s = sd(percentage), 
                  q10 = quantile(percentage, probs = c(0.1, 0.9))[1], 
                  q90 = quantile(percentage, probs = c(0.1, 0.9))[2],
                  percentage.median = median(percentage), 
                  percentage = (sum((percentage * sum))/sum(sum))) %>%
        
        

        mutate(fibertype = factor(fibertype, levels = c("type1", "type2a", "type2x"), 
                                  labels = c("Type I", 
                                             "Type IIA", 
                                             "Type IIX")), 
               time_cont = if_else(time == "w0", 
                                   1, if_else(time == "w2", 2, 3))) %>%
        print()
        

stacked_fibertypes_fig <- stacked_fibertypes_data %>%        
        ggplot(aes(time_cont, percentage, fill = fibertype)) + 
        
        annotate("rect", xmin = 0.8, xmax = Inf, ymin = 0, ymax = 0.4, 
                 alpha = 0.2) +
        
        # Two different colums as stacked and nudge cannot be combined
        geom_col(data = filter(stacked_fibertypes_data, sets == "single"), 
                 aes(time_cont - 0.15, percentage, fill = fibertype),
                 position = position_stack(),
                 width = 0.2) +

        geom_col(data = filter(stacked_fibertypes_data, sets == "multiple"), 
                aes(time_cont + 0.15, percentage, fill = fibertype),
                 position = position_stack(),
                width = 0.2) +
        # Scales
      scale_y_continuous(limits = c(0, 1),
                         expand = c(0,0), 
                         breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), 
                         labels = c(0, 20, 40, 60, 80, 100)) +
        
        scale_x_continuous(limits = c(-0.2, 3.5), 
                           expand = c(0, 0), 
                           breaks = c(1, 2, 3), 
                           labels = c("Week 0", 
                                      "Week 2", 
                                      "Week 12")) +
        
        scale_fill_manual(values = c(group.study.color[5], 
                                         group.study.color[2],
                                         group.study.color[1])) +
        

        labs(y = "% of total fibers") +
        
        ## Annotations
        annotate("text", 
                 x = 0.3,   
                 y = c(0.8, 0.4, 0.05), 
                 label = c("Type I", 
                           "Type IIA", 
                           "Type IIX"), 
                 color = c(group.study.color[5], 
                           group.study.color[2],
                           group.study.color[1]),
                 size = 2.5) +

        annotate("text", 
                 x = c(1.67, 2.33), 
                 y = 0.65, 
                 label = c("Low-volume", "Moderate-volume"), 
                 size = 2.5,
                 angle = c(90, -90)) +
        
        
        dissertation_theme() +
        
        theme(legend.position = "none", 
              axis.title.x = element_blank()) 



 

fib.sum.avg <-  fib.trans %>%
        filter(fibertype %in% c("type1", "type2a", "type2x")) %>%
        group_by(subject, sets, time, fibertype) %>%
        summarise(percentage = (sum((percentage * sum))/sum(sum)), 
                  sum = sum(sum)) %>%

        group_by(sets, time, fibertype) %>%
        summarise(s = sd(percentage), 
                  q10 = quantile(percentage, probs = c(0.1, 0.9))[1], 
                  q90 = quantile(percentage, probs = c(0.1, 0.9))[2],
                  percentage.median = median(percentage), 
                  percentage = (sum((percentage * sum))/sum(sum))) %>%
        mutate(percentage = percentage * 100) %>%
        print()

       
        
paired_plot_2x <- fib.trans %>%
        filter(fibertype %in% c("type1", "type2a", "type2x")) %>%
        group_by(subject, sets, time, fibertype) %>%
        summarise(percentage = (sum((percentage * sum))/sum(sum)), 
                  sum = sum(sum)) %>%
        filter(fibertype == "type2x") %>%
        mutate(percentage = percentage * 100) %>%   
        ggplot(aes(sets, percentage, group = subject)) + 
        
        geom_line(color = group.study.color[1], size = 0.5) +
        
       
        geom_point(data = filter(fib.sum.avg, sets == "single", fibertype == "type2x"), 
                           position = position_nudge(x = -0.1), 
                           aes(group = NULL), 
                   size = 2.5, 
                   shape = 21, 
                   fill = group.study.color[4]) +
        
        geom_point(data = filter(fib.sum.avg, sets == "multiple",  fibertype == "type2x"), 
                   position = position_nudge(x = 0.1), 
                   aes(group = NULL), 
                   size = 2.5, 
                   shape = 21, 
                   fill = group.study.color[5]) +

        facet_wrap(~ time) +
        
        ## Scales 
        scale_y_continuous(limits = c(0, 40), 
                           expand = c(0, 0), 
                           breaks = c(0, 5, 10, 20, 30, 40), 
                           position = "right") +
        labs(y = "Type IIX % of total fibers") +
        

        dissertation_theme() +

        theme(axis.title.y = element_text(size = 7), 
              axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              strip.background = element_blank(), 
              strip.text.x = element_text(size = 7, hjust = 0), 
              legend.title = element_text(size = 7),
              legend.text = element_text(size = 7),
              legend.key.height = unit(0.25, "cm"),
              panel.spacing.x = unit(0.3, "lines"),
              legend.position="bottom",
              legend.box="horizontal") 
        

### Clean plot to use for background and annotations


back_plot <- data.frame(x = c(0, 1), y = c(0, 1)) %>%
        ggplot(aes(x, y)) + theme_void()



fiber_type_insets <- ggdraw(back_plot) + 
        draw_plot(stacked_fibertypes_fig, 0, 0, 0.5, 1) + 
        annotate("polygon", x = c(0.485, 0.485, 0.58, 0.58), y = c(0.1, 0.45, 0.95, 0.17), 
                 alpha = 0.1) +
          draw_plot(paired_plot_2x, 0.58, 0.15, 0.4, 0.8) +
        annotate("text", 
                 x = c(0.73,0.77), y = c(0.13,0.13), 
                 label = c("Low-\nvolume", 
                           "Moderate-\nvolume"), 
                 color = c("black"),
                 
                 size = 2.5, 
                 hjust = c(1,0 )) 
        
      




## immuno control image #################

myhc_fig_comps <- list(overall_type2x = overall_type2x, 
                       myhc2x_mrna_means = myhc2x_mrna_means, 
                       diff_2x_mRNA = diff_2x_mRNA, 
                       fiber_type_insets = fiber_type_insets)

saveRDS(myhc_fig_comps, "./data/derivedData/study1-myhc-composition/myhc_fig.RDS")


## Model for stats in manuscript
mx <- glmmTMB(percentage ~ time + time:sets + (1|subject) + (1|leg) + (1|sample), 
              data = fib.trans[fib.trans$fibertype == "type2x", ], 
              family = binomial(link = "logit"), 
              weights = sum)



mx.summary <- summary(mx)
mx.summary <- cbind(mx.summary$coefficients$cond, confint(mx)[1:6,])

saveRDS(mx.summary, "./data/derivedData/study1-myhc-composition/myhc2x_stats.RDS")




