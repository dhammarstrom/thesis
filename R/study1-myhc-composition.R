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
              legend.key.height = unit(0.25, "cm"),
              legend.position="bottom",
              legend.box="horizontal") +
        
        scale_fill_gradient(breaks = c(0, 15, 30), 
                            limits = c(0, 30), 
                            low = group.study.color[3], 
                            high = group.study.color[1]) +
        labs(fill = "Percentage of fibers", 
             y = "Participants") +
        
        guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5))
        
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
                                      if_else(timepoint == "w2post", 2.4, 12)))) %>%
        
        ggplot(aes(time, exp(estimate))) +
        geom_errorbar(aes(ymin = exp(lower), ymax = exp(upper)), 
                      width = 0) +
 
        
        geom_point(shape = 21,  size = 2, ) +
        
        scale_x_continuous(breaks = c(0, 2, 12), 
                           limits = c(-0.5, 12.5), 
                           labels = c("Week 0", "Week 2", "Week 12")) +
dissertation_theme() +
        
        ylab("<i>MYH1</i> Moderate- / Low-volume") + 
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
        geom_point(shape = 21, position = position_dodge(width = 0.2), size = 2) +
        scale_y_continuous(limits = c(-16.5, -12.5), 
                           breaks = seq(from = -16.5, to = -12.5, by = 1), 
                           expand = c(0, 0)) +
        
        scale_x_continuous(breaks = c(0, 2, 12), 
                           limits = c(-0.5, 12.5), 
                           labels = c("Week 0", "Week 2", "Week 12")) +
        
        
        dissertation_theme() +

        ylab("<i>MYH1</i> mRNA<br>estimated log-abundance") + 
        theme(axis.title.x = element_blank(),
              axis.title.y = element_markdown(size = 7),
              axis.text.y = element_text(size = 7),
              axis.text.x = element_blank(),
              axis.line.x = element_blank(),
              axis.ticks.x = element_blank(),
              legend.position = "none")


######## Percentage plot ############################


stacked_bar <- fib.trans %>%
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

        filter(fibertype != "type1") %>%
        mutate(fibertype = factor(fibertype, levels = c("type2a", "type2x"), 
                                  labels = c("Type IIA", 
                                             "Type IIX"))) %>%
        ggplot(aes(time, percentage, fill = fibertype, group = sets)) + 
        
        annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = 0.08, 
                 alpha = 0.2) +
        
        
        scale_y_continuous(limits = c(0, 0.6),
                           expand = c(0,0), 
                           breaks = c(0, 0.2, 0.4, 0.6), 
                           labels = c(0, 20, 40, 60)) +

        geom_bar(position = position_dodge(width = 0.5), stat = "identity", width = 0.3) + 
        theme(legend.position = "none") +
        dissertation_theme() +
        
        theme(axis.title.y = element_text(size = 7), 
              legend.position = "none",
              axis.text.y = element_text(7), 
              axis.text.x = element_text(size = 7),
              axis.title.x = element_blank(),
              strip.background = element_blank(), 
              strip.text.x = element_text(size = 7, hjust = 0), 
              legend.title = element_text(size = 7),
              legend.text = element_text(size = 7),
              legend.key.height = unit(0.25, "cm"),
              panel.spacing.x = unit(2, "lines")) 



        print()
 

fib.sum.avg <- fib.sum %>%
        group_by(sets, time, fibertype) %>%
        summarise(s = sd(percentage), 
                  q10 = quantile(percentage, probs = c(0.1, 0.9))[1], 
                  q90 = quantile(percentage, probs = c(0.1, 0.9))[2],
                  percentage.median = median(percentage), 
                  percentage = (sum((percentage * sum))/sum(sum))) %>%
        print()

       
        
paired_plot_2x <- fib.trans %>%
        filter(fibertype %in% c("type1", "type2a", "type2x")) %>%
        group_by(subject, sets, time, fibertype) %>%
        summarise(percentage = (sum((percentage * sum))/sum(sum)), 
                  sum = sum(sum)) %>%
        filter(fibertype == "type2x") %>%
        mutate(percentage = percentage * 100) %>%   
        ggplot(aes(sets, percentage, group = subject)) + 
        
        geom_line(color = "gray50", size = 0.5) +
        
       
        geom_point(data = filter(fib.sum.avg, sets == "single"), 
                           position = position_nudge(x = -0.2), 
                           aes(group = NULL), 
                   size = 2.5, 
                   shape = 21, 
                   fill = group.study.color[4]) +
        
        geom_point(data = filter(fib.sum.avg, sets == "multiple"), 
                   position = position_nudge(x = 0.2), 
                   aes(group = NULL), 
                   size = 2.5, 
                   shape = 21, 
                   fill = group.study.color[1]) +

        facet_wrap(~ time) +

        dissertation_theme() +

        theme(axis.title.y = element_text(size = 7), 
              
              axis.text.y = element_text(7), 
              axis.text.x = element_text(size = 7, angle = 45, hjust = 0.9),
              axis.title.x = element_blank(),
              strip.background = element_blank(), 
              strip.text.x = element_text(size = 7, hjust = 0), 
              legend.title = element_text(size = 7),
              legend.text = element_text(size = 7),
              legend.key.height = unit(0.25, "cm"),
              panel.spacing.x = unit(1, "lines"),
              legend.position="bottom",
              legend.box="horizontal") 
        

### Clean plot to use for background and annotations


back_plot <- data.frame(x = c(0, 1), y = c(0, 1)) %>%
        ggplot(aes(x, y)) + theme_void()



fiber_type_insets <- ggdraw(back_plot) + 
        draw_plot(stacked_bar, 0, 0, 0.5, 1) + 
        annotate("polygon", x = c(0.485, 0.485, 0.58, 0.58), y = c(0.05, 0.18, 0.9, 0.05), 
                 alpha = 0.2) +
        draw_plot(paired_plot_2x, 0.58, 0, 0.4, 1) 




## immuno control image #################

myhc_fig_comps <- list(overall_type2x = overall_type2x, 
                       myhc2x_mrna_means = myhc2x_mrna_means, 
                       diff_2x_mRNA = diff_2x_mRNA, 
                       fiber_type_insets = fiber_type_insets)

saveRDS(myhc_fig_comps, "./data/derivedData/study1-myhc-composition/myhc_fig.RDS")


## Combine figure directly in manuscript as saving the png takes some space 

myhc_fig <- plot_grid(plot_grid(img, overall_type2x, 
                    rel_widths = c(0.72, 0.25), 
                    ncol = 2),
               plot_grid(
                       plot_grid(myhc2x_mrna_means, diff_2x_mRNA, ncol = 1, align = "v"), 
                       fiber_type_insets, 
                       ncol = 2, 
                       rel_widths = c(0.5, 0.5)), 
          ncol = 1, 
          rel_heights = c(0.6, 0.4))
   


























### How many subjects have 0 2X in the start of the training ###

fib.trans %>%
        filter(fibertype == "type2x_pure" & time == "w0") %>%
        filter(percentage == 0) %>%
        group_by() %>%
        summarise(n = n_distinct(subject)) %>%
        print()





type2x <- fib.trans %>%
        filter(fibertype %in% c("type2x")) %>%
        group_by(subject) %>%
        mutate(m = mean(percentage)) %>%
        ungroup() %>%
        mutate(high = if_else(m > median(m), 1, 0)) %>%
        group_by(time, sets, fibertype) %>%
        summarise(s = sd(percentage * 100), 
                  q10 = quantile(percentage * 100, probs = c(0.1, 0.))[1], 
                  q90 = quantile(percentage * 100, probs = c(0.1, 0.85))[2],
                  percentage = (sum((percentage * sum))/sum(sum))*100) %>%
        ggplot(aes(time, percentage, fill = sets)) + 
        geom_errorbar(aes(ymin = q10, ymax = q90), width = 0, position = position_dodge(width = 0.2)) +
        geom_point(position = position_dodge(width = 0.2), shape = 21, size = 3) + 


        scale_y_continuous(breaks = c(0, 7, 14), limits = c(0, 14), expand = c(0,0)) +
        theme(axis.title = element_text(size = 8), 
              axis.text = element_text(size = 8),
              legend.position = "none")







type2_switch_summary <- fib.trans %>%
        filter(fibertype %in% c("type2x_pure", "type2x_hybrids")) %>%
        group_by(subject, sets, time, fibertype) %>%
        summarise(percentage = (sum((percentage * sum))/sum(sum)), 
                  sum = sum(sum)) %>%
        group_by(time, sets, fibertype) %>%
        summarise(s = sd(percentage * 100), 
                  q10 = quantile(percentage * 100, probs = c(0.1, 0.9))[1], 
                  q90 = quantile(percentage * 100, probs = c(0.1, 0.9))[2],
                  percentage = (sum((percentage * sum))/sum(sum)) * 100) %>%
        ungroup() %>%
        mutate(fibertype = factor(fibertype, levels = c("type2x_pure", "type2x_hybrids"), labels = c("Type IIX", "Type IIX/IIA")),
               sets = factor(sets, levels = c("single", "multiple"), labels = c("Single-set", "Multiple-set")), 
               time = factor(time, levels = c("w0", "w2", "w12"), labels = c("Week 0", "Week 2", "Week 12"))) 




type_switch_data <- fib.trans %>%
        filter(fibertype %in% c("type2x_pure", "type2x_hybrids")) %>%
        group_by(subject, sets, time, fibertype) %>%
        summarise(percentage = (sum((percentage * sum))/sum(sum))) %>%
        ungroup() %>%
        mutate(highlight = factor(if_else(subject == "FP46", "FP46", 
                                          if_else(subject == "FP20", "FP20", 
                                                  if_else(subject == "FP9", "FP9", 
                                                          if_else(subject == "FP42", "FP42", "other")))), 
                                  levels = c("FP46", "FP20", "FP9", "FP5", "other")), 
               fibertype = factor(fibertype, levels = c("type2x_pure", "type2x_hybrids"), labels = c("Type IIX", "Type IIX/IIA")),
               sets = factor(sets, levels = c("single", "multiple"), labels = c("Single-set", "Multiple-set")), 
               time = factor(time, levels = c("w0", "w2", "w12"), labels = c("Week 0", "Week 2", "Week 12")), 
               percentage = percentage * 100) %>%
        dplyr::select(subject, time, sets, fibertype, highlight, percentage) 



### Models

type2xpure.m1 <- glmer(percentage ~  time + time:sets  + (1|subject)  + (1|leg) + (1|sample), 
                       data = fib.trans[fib.trans$fibertype == "type2x_pure", ], 
                       family = binomial(link = "logit"), 
                       weights = sum,
                       control = glmerControl(optimizer = "nloptwrap", 
                                              calc.derivs = FALSE))

type2xhybrids.m1 <- glmer(percentage ~  time + time:sets  + (1|subject)  + (1|leg) + (1|sample), 
                          data = fib.trans[fib.trans$fibertype == "type2x_hybrids", ], 
                          family = binomial(link = "logit"), 
                          weights = sum,
                          control = glmerControl(optimizer = "nloptwrap", 
                                                 calc.derivs = FALSE))



text_size_factor <- 1.4

anno_type2switch_plot <- data.frame(time = c("Week 2", "Week 12", "Week 12", "Week 12", "Week 2"), 
                                    y = c(20, 15, 7, 20, 30), 
                                    size.lab = c(time_size , time_size , star_size, time_size , star_size),
                                    fibertype = c("Type IIX", "Type IIX", "Type IIX", "Type IIX/IIA", "Type IIX/IIA"), 
                                    sets = c("Single-set", "Single-set", "Multiple-set", "Single-set", "Multiple-set"), 
                                    label = c(pval(data.frame(summary(type2xpure.m1)$coef)[2,4] , flag = TRUE, sign = time.flag),
                                              pval(data.frame(summary(type2xpure.m1)$coef)[2,4] , flag = TRUE, sign =  time.flag),
                                              pval(data.frame(summary(type2xpure.m1)$coef)[2,4] , flag = TRUE, sign = sets.flag),
                                              pval(data.frame(summary(type2xhybrids.m1)$coef)[3, 4], flag = TRUE, sign =  time.flag),
                                              pval(data.frame(summary(type2xhybrids.m1)$coef)[5, 4], flag = TRUE, sign = sets.flag)))





strip_label_type2switch <- data.frame(time = c("Week 2", "Week 2"),
                                      y = c(37, 37), 
                                      size.lab = c(7.5, 7.5), 
                                      fibertype = c( "Type IIX", "Type IIX/IIA"), 
                                      sets = c("Single-set", "Single-set"), 
                                      label = c("Type IIX", "Type IIX/IIA"))


individual_swith_plot <- type_switch_data %>% 
        ggplot(aes(time, percentage, group = subject)) + 
        geom_line(alpha = 0.2) + 
        #  geom_line(data = filter(type_switch_data, highlight != "other"), 
        #            aes(color = highlight), 
        #            size = 1) +
        facet_grid(sets ~ fibertype) +
        
        
        ### Single-set summary points
        geom_point(data = filter(type2_switch_summary, time == "Week 0", sets == "Single-set"), 
                   aes(group = NULL), position = position_nudge(x = -0.1), 
                   size = 2.5, shape = 21, fill = fill_scale[1]) +
        geom_point(data = filter(type2_switch_summary, time == "Week 2", sets == "Single-set"), 
                   aes(group = NULL), position = position_nudge(x = 0), 
                   size = 2.5, shape = 21, fill = fill_scale[1]) +
        geom_point(data = filter(type2_switch_summary, time == "Week 12", sets == "Single-set"), 
                   aes(group = NULL), position = position_nudge(x = 0.1), 
                   size = 2.5, shape = 21, fill = fill_scale[1]) +
        
        ### Multiple-sets summary points
        geom_point(data = filter(type2_switch_summary, time == "Week 0", sets == "Multiple-set"), 
                   aes(group = NULL), position = position_nudge(x = -0.1), 
                   size = 2.5, shape = 21, fill = fill_scale[2]) +
        geom_point(data = filter(type2_switch_summary, time == "Week 2", sets == "Multiple-set"), 
                   aes(group = NULL), position = position_nudge(x = 0), 
                   size = 2.5, shape = 21, fill = fill_scale[2]) +
        geom_point(data = filter(type2_switch_summary, time == "Week 12", sets == "Multiple-set"), 
                   aes(group = NULL), position = position_nudge(x = 0.1), 
                   size = 2.5, shape = 21, fill = fill_scale[2]) +
        
        geom_text(data = anno_type2switch_plot, aes(time, y, group = "NULL", label = label), size = c(time_size-1.5, 
                                                                                                      time_size-1.5, 
                                                                                                      star_size,
                                                                                                      time_size-1.5,
                                                                                                      star_size)) +
        
        geom_text(data = strip_label_type2switch, aes(time, y, group = "NULL", label = label), size = 2.5) +
        
        scale_color_manual(values = c("#1b9e77", "#d95f02", "#e7298a", "blue", "gray40")) +
        scale_y_continuous(breaks = c(0, 10, 20, 30, 40), limits = c(-2, 40), expand = c(0, 0)) +
        publr_theme() +
        ylab("IHC proportion (%)") +
        theme(legend.position = "none", 
              axis.text.x = element_text(size = figure_text_size, angle = 45, hjust = 1), 
              axis.text.y = element_text(size = figure_text_size),
              axis.title.y = element_text(size = figure_text_size), 
              axis.title.x = element_blank(), 
              strip.text = element_text(size = figure_text_size), 
              strip.background.y = element_rect(fill = "white", color = "white"), 
              strip.background.x = element_blank(), 
              strip.text.x = element_blank(),
              panel.spacing.y = unit(1, "lines")) 




#### May be used in alternative figure ###
type2ax_plot <- plot_grid(ggdraw() + 
                                  annotate("text", x = 0.1, y = 0.95, label = "Type IIX", size = 3) +
                                  annotate("text", x = 0.59, y = 0.95, label = "Type IIX/IIA", size = 3) +
                                  draw_plot(plot = temp, x = 0.05, y = 0, width = 0.95, height = 0.9),
                          plot_grid(pure_group, hybrids_group, ncol = 2), nrow = 2)








### How many subjects have 0 2X in the start of the training ###
# 
# fib.trans %>%
#   filter(fibertype == "type2x" & time == "w0") %>%
#   filter(percentage == 0) %>%
#   print()
# 


#### Alternative figure 5A ########

### Statistics for mhc family figures ####

### Models PROTEIN ####
# A binomial (logit transformation) generalized me-model is choosen to account for the 0-1 bound 
# of proportions data. The number of fibers counted were used to indicate number of "trials" and 
# proportions represents number of "successes" i.e. number of positive fibers.


baseline.2x <- fib.trans %>%
        filter(time == "w0", 
               fibertype == "type2x") %>%
        dplyr::select(subject, leg, sets, percentage) %>%
        group_by(subject, leg) %>%
        summarise(baseline.2x = mean(percentage, na.rm = TRUE)) %>%
        mutate(xbaseline = if_else(baseline.2x > 0.05, "upr", "lwr")) %>%
        print()
        










temp <- fib.trans %>%
        filter(time != "w0", 
               fibertype == "type2x") %>%
        
        inner_join(baseline.2x) %>%
        print()
        







type2a.m1 <- glmer(percentage ~ time + time:sets + (1|subject)  + (1|sample),
                   data = fib.trans[fib.trans$fibertype == "type2a", ], 
                   family = "binomial", 
                   weights = sum,
                   control = glmerControl(optimizer = "nloptwrap", 
                                          calc.derivs = FALSE))






type1.m1 <- glmer(percentage ~ time + time:sets + (1|subject) + (1|leg) + (1|sample),
                  data = fib.trans[fib.trans$fibertype == "type1", ], 
                  family = "binomial", 
                  weights = sum,
                  control = glmerControl(optimizer = "nloptwrap", 
                                         calc.derivs = FALSE))



fibertypes.results <- rbind(data.frame(parameter = rownames(coef(summary(type2x.m1))), 
                                       coef(summary(type2x.m1)), row.names = NULL, 
                                       fibertype = "type2x"),
                            data.frame(parameter = rownames(coef(summary(type2a.m1))),
                                       coef(summary(type2a.m1)), 
                                       row.names = NULL, 
                                       fibertype = "type2a"),
                            data.frame(parameter = rownames(coef(summary(type1.m1))),
                                       coef(summary(type1.m1)), 
                                       row.names = NULL, 
                                       fibertype = "type1")) %>%
        mutate(time = if_else(parameter == "timew0:setsmultiple", "Week 0", 
                              if_else(parameter == "timew2:setsmultiple", "Week 2", 
                                      if_else(parameter == "timew12:setsmultiple", "Week 12", 
                                              if_else(parameter == "timew2", "Week 2",
                                                      if_else(parameter == "timew12", "Week 12", "Week 0")))))) %>% 
        mutate(interaction = if_else(parameter %in% c("(Intercept)", "timew2", "timew12"), "main", "inter")) %>%
        dplyr::select(parameter, Estimate, p.value = Pr...z.., fibertype, time, interaction) %>%
        mutate(percentage = if_else(fibertype == "type2x", 8, 
                                    if_else(fibertype == "type2a", 45, 100)),
               sets = "Single-\nset", 
               time = factor(time, levels = c("Week 0", "Week 2", "Week 12"))) 

### Models mRNA gene family ###

myhc.dat <- read_feather("./data/study-1/qpcr/qpcrdat_replicates") %>%
        filter(target %in% c("MyHC1 F1R1",   "MyHC2A F5R5",  "MyHC2X F5R5")) %>%
        filter(timepoint != "w2post") %>%
        mutate(expression = eff^-cq, 
               count = eff^(37-cq)) %>%
        group_by(subject, leg, timepoint, cdna) %>%
        mutate(total.expression = sum(expression), 
               total.count = sum(count)) %>%
        ungroup() %>%
        mutate(rel.expression = expression/total.expression, 
               rel.count = count / total.count) %>%
        mutate(gene = if_else(target == "MyHC1 F1R1", "MyHC1", 
                              if_else(target == "MyHC2A F5R5", "MyHC2A", "MyHC2X"))) %>%
        dplyr::select(subject, leg, timepoint, gene, rel.expression, rel.count, total.count) %>%
        inner_join(leg) %>%
        mutate(sets = factor(sets, levels = c("single", "multiple")),
               leg = paste(subject, leg), 
               sample = paste0(subject, leg, timepoint)) %>%
        filter(include == "incl")


### As counts in the mRNA-gene family approach is arbitrary (large), the denominator can 
# be said to be unknown. The binomial model therefore not a good idea. Using the beta family here instead.


myhc1x <- glmmTMB(rel.count ~ timepoint + timepoint:sets  + (1|subject) + (1|leg), 
                  data = myhc.dat[myhc.dat$gene == "MyHC1", ], 
                  family = beta_family())

myhc2ax <- glmmTMB(rel.count ~ timepoint + timepoint:sets  + (1|subject/leg), 
                   data = myhc.dat[myhc.dat$gene == "MyHC2A", ], 
                   family = beta_family())

myhc2xx <- glmmTMB(rel.count ~ timepoint + timepoint:sets  + (1|subject/leg), 
                   data = myhc.dat[myhc.dat$gene == "MyHC2X", ], 
                   family = beta_family())

myhc1x.e <- glmmTMB(rel.expression ~ timepoint + timepoint:sets  + (1|subject) + (1|leg), 
                    data = myhc.dat[myhc.dat$gene == "MyHC1", ], 
                    family = beta_family())

myhc2ax.e <- glmmTMB(rel.expression ~ timepoint + timepoint:sets  + (1|subject/leg), 
                     data = myhc.dat[myhc.dat$gene == "MyHC2A", ], 
                     family = beta_family())

myhc2xx.e <- glmmTMB(rel.expression ~ timepoint + timepoint:sets  + (1|subject/leg), 
                     data = myhc.dat[myhc.dat$gene == "MyHC2X", ], 
                     family = beta_family())




### Collect parameters in a summary data frame.

fibertypes.results.mrna <- rbind(data.frame(parameter = rownames(coef(summary(myhc2xx))$cond), 
                                            coef(summary(myhc2xx))$cond, row.names = NULL, 
                                            fibertype = "type2x"),
                                 data.frame(parameter = rownames(coef(summary(myhc2ax))$cond),
                                            coef(summary(myhc2ax))$cond, 
                                            row.names = NULL, 
                                            fibertype = "type2a"),
                                 data.frame(parameter = rownames(coef(summary(myhc1x))$cond),
                                            coef(summary(myhc1x))$cond, 
                                            row.names = NULL, 
                                            fibertype = "type1")) %>%
        mutate(timepoint = if_else(parameter == "timepointw0:setsmultiple", "Week 0", 
                                   if_else(parameter == "timepointw2pre:setsmultiple", "Week 2", 
                                           if_else(parameter == "timepointw12:setsmultiple", "Week 12", 
                                                   if_else(parameter == "timepointw2pre", "Week 2",
                                                           if_else(parameter == "timepointw12", "Week 12", "Week 0")))))) %>% 
        mutate(interaction = if_else(parameter %in% c("(Intercept)", "timepointw2pre", "timepointw12"), "main", "inter")) %>%
        dplyr::select(parameter, Estimate, p.value = Pr...z.., fibertype, timepoint, interaction) %>%
        mutate(percentage = if_else(fibertype == "type2x", 8, 
                                    if_else(fibertype == "type2a", 45, 100)),
               sets = "Single-\nset", 
               timepoint = factor(timepoint, levels = c("Week 0", "Week 2", "Week 12"))) %>%
        rowwise() %>%
        mutate(size = if_else(p.value < 0.05, if_else(Estimate > 0, "<", ">"), ""),
               lab = pval(p.value, flag = TRUE), 
               lab = paste0(size,"\n",lab)) 

#### Fibertype family Figures ######

### Fibertype MHC protein ####
pos <- position_dodge(width = 0.3)
x.fib.text <- 0.6 # Fiber type annotation placement on x axis
errorbar.size <- error_size


mhc_fam_theme <- function() { publr_theme() +
                theme(axis.title.x = element_blank(),
                      axis.text.x = element_text(size = figure_text_size, angle = 45, hjust = 1),
                      axis.text.y = element_text(size = figure_text_size),
                      strip.text = element_text(size = figure_text_size),
                      # axis.line.x = element_blank(),
                      # axis.ticks.x = element_blank(),
                      legend.position = "none",
                      axis.title.y = element_blank()) }


fib2x <- fib.sum %>%
        mutate(fibertype = factor(fibertype, 
                                  levels = c("type1", "type2a", "type2x"))) %>%
        filter(fibertype == "type2x") %>%
        ungroup() %>%
        mutate(time = factor(time, labels = c("Week 0", "Week 2", "Week 12")),
               sets = factor(sets, labels = c("Single-\nset", "Multiple-\nset")), 
               percentage = percentage * 100,
               percentage.median = percentage.median * 100,
               q10 = q10 * 100,
               q90 = q90 * 100) %>%
        ggplot(aes(time, percentage, fill = sets)) + 
        ## Adds a split violin plot to the plot showing the distribution.
        #  geom_split_violin(data = fib.trans %>%
        #                      ungroup() %>%
        #                      mutate(fibertype = factor(fibertype, 
        #                            levels = c("type1", "type2a", "type2x")),
        #                            time = factor(time, labels = c("Week 0", "Week 2", "Week 12")),
        #                            sets = factor(sets, labels = c("Single-\nsets", "Multiple-\nsets")), 
        #                            percentage = percentage * 100),
        #                    position = position_dodge(width = 0),
        #                    trim = TRUE) +
        geom_errorbar(aes(ymin = q10, ymax = q90), width = 0, position = pos, size = errorbar.size) +
        geom_point(shape = 21, position = pos, size = 2) +
        geom_signif(y_position = c(9.2, 5.7), 
                    xmin = c(1.9, 2.9), 
                    xmax = c(2.1, 3.1), 
                    annotation=c(pval(
                            fibertypes.results %>%
                                    filter(time == "Week 2",
                                           fibertype == "type2x",
                                           interaction == "inter") %>%
                                    dplyr::select(p.value) %>%
                                    as.numeric(), 
                            flag = TRUE, 
                            sign = sets.flag),
                            pval(
                                    fibertypes.results %>%
                                            filter(time == "Week 12",
                                                   fibertype == "type2x",
                                                   interaction == "inter") %>%
                                            dplyr::select(p.value) %>%
                                            as.numeric(), 
                                    flag = TRUE, 
                                    sign = sets.flag)), 
                    tip_length = tip_length, 
                    family="sans",
                    textsize = star_size, 
                    size = bar_size,
                    vjust = condition_vjust) +
        geom_signif(y_position = c(6.9), 
                    xmin = c(3), 
                    xmax = c(3), 
                    annotation=c(pval(
                            fibertypes.results %>%
                                    filter(time == "Week 12",
                                           fibertype == "type2x",
                                           interaction == "main") %>%
                                    dplyr::select(p.value) %>%
                                    as.numeric(), 
                            flag = TRUE, 
                            sign = time.flag)), 
                    tip_length = tip_length, 
                    family="sans",
                    textsize = time_size-1.5, 
                    size = bar_size) +
        mhc_fam_theme () +
        scale_fill_manual(values = fill_scale) +
        xlab(" ") +
        scale_y_continuous(expand = c(0,0), breaks = c(0,  5 , 10, 15), limits = c(0, 15)) 



fib2a <- fib.sum %>%
        mutate(fibertype = factor(fibertype, 
                                  levels = c("type1", "type2a", "type2x"))) %>%
        filter(fibertype == "type2a") %>%
        ungroup() %>%
        mutate(time = factor(time, labels = c("Week 0", "Week 2", "Week 12")),
               sets = factor(sets, labels = c("Single-\nset", "Multiple-\nset")), 
               percentage = percentage * 100,
               q10 = q10 * 100,
               q90 = q90 * 100) %>%
        ggplot(aes(time, percentage, fill = sets)) + 
        
        geom_errorbar(aes(ymin = q10, ymax = q90), width = 0, position = pos, size = errorbar.size) +
        geom_point(shape = 21, position = pos, size = 2) +
        
        mhc_fam_theme () +
        scale_fill_manual(values = fill_scale) +
        xlab("Week") +
        scale_y_continuous(expand = c(0,0), breaks = c(30,45, 60, 75), limits = c(30, 75))
#  annotate("text", x = x.fib.text, y = 60 + (75-60)/2, label = "Type 2A", size = 3)


fib1 <- fib.sum %>%
        mutate(fibertype = factor(fibertype, 
                                  levels = c("type1", "type2a", "type2x"))) %>%
        filter(fibertype == "type1") %>%
        ungroup() %>%
        mutate(time = factor(time, labels = c("Week 0", "Week 2", "Week 12")),
               sets = factor(sets, labels = c("Single-\nset", "Multiple-\nset")), 
               percentage = percentage * 100,
               q10 = q10 * 100,
               q90 = q90 * 100) %>%
        ggplot(aes(time, percentage, fill = sets)) + 
        geom_errorbar(aes(ymin = q10, ymax = q90), width = 0, position = pos, size = errorbar.size) +
        geom_point(shape = 21, position = pos, size = 2) +
        mhc_fam_theme () +
        scale_fill_manual(values = fill_scale) +
        xlab(" ") +
        theme(axis.title.y = element_text(size = 8, 
                                          angle = 90,
                                          margin = margin(t = 0, r = 1, b = 0, l = -4))) +
        ylab("IHC proportions") +
        scale_y_continuous(expand = c(0,0), breaks = c(25, 40, 55, 70), limits = c(25, 70)) 
#  annotate("text", x = x.fib.text, y = 55 + (70-55)/2, label = "Type 1", size = 3)




##### Gene family plots ####


gene.family.dat <-  read_feather("./data/qpcrdat_replicates") %>%
        filter(target %in% c("MyHC1 F1R1",   "MyHC2A F5R5",  "MyHC2X F5R5")) %>%
        mutate(eff = if_else(target == "MyHC2A F5R5", 1.75, eff)) %>%
        group_by(subject, leg, timepoint, target) %>%
        summarise(eff = mean(eff, na.rm = TRUE),
                  cq = mean(cq, na.rm = TRUE)) %>%
        mutate(expression = eff ^ -cq, 
               count = eff ^ (37-cq)) %>%
        group_by(subject, leg, timepoint) %>%
        mutate(total.expression = sum(expression), 
               total.count = sum(count)) %>%
        ungroup() %>%
        mutate(rel.expression = count / total.count, 
               rel.count = count / total.count) %>%
        mutate(fibertype = if_else(target == "MyHC1 F1R1", "type1", 
                                   if_else(target == "MyHC2A F5R5", "type2a", "type2x")),
               time = if_else(timepoint == "w0", "w0", 
                              if_else(timepoint == "w12", "w12", "w2"))) %>%
        filter(timepoint != "w2post") %>%
        inner_join(leg) %>%
        filter(include == "incl") %>%
        dplyr::select(subject, timepoint, sets,  fibertype, rel.expression) %>%
        group_by(fibertype, timepoint, sets) %>%
        summarise(q10 = quantile(rel.expression, probs = c(0.1, 0.9))[1], 
                  q90 = quantile(rel.expression, probs = c(0.1, 0.9))[2],
                  percentage = mean(rel.expression)) 



pos <- position_dodge(width = 0.3)
x.fib.text <- 1.5 # Fiber type annotation placement on x axis


### Fibertype gene family ####

gene_fam_theme <- function() { publr_theme() +
                theme(axis.title.x = element_blank(),
                      axis.text.x = element_blank(),
                      axis.text.y = element_text(size = figure_text_size),
                      strip.text = element_text(size = figure_text_size),
                      axis.line.x = element_blank(),
                      axis.ticks.x = element_blank(),
                      legend.position = "none",
                      axis.title.y = element_blank()) }



genefam2x <- gene.family.dat %>%
        ungroup() %>%
        mutate(fibertype = factor(fibertype, 
                                  levels = c("type1", "type2a", "type2x"))) %>%
        filter(fibertype == "type2x") %>%
        ungroup() %>%
        mutate(timepoint = factor(timepoint, levels = c("w0", "w2pre", "w12"),
                                  labels = c("Week 0", "Week 2", "Week 12")),
               sets = factor(sets, levels = c("single", "multiple"),
                             labels = c("Single-\nset", "Multiple-\nset")), 
               percentage = percentage * 100,
               q10 = q10 * 100,
               q90 = q90 * 100) %>%
        ggplot(aes(timepoint, percentage, fill = sets)) + 
        geom_errorbar(aes(ymin = q10, ymax = q90), width = 0, position = pos, size = errorbar.size) +
        geom_point(shape = 21, position = pos, size = 2) +
        geom_signif(y_position = c(19, 31), 
                    xmin = c(1.9, 2.9), 
                    xmax = c(2.1, 3.1), 
                    annotation=c(pval(
                            fibertypes.results.mrna %>%
                                    filter(timepoint == "Week 2",
                                           fibertype == "type2x",
                                           interaction == "inter") %>%
                                    dplyr::select(p.value) %>%
                                    as.numeric(), 
                            flag = TRUE, 
                            sign = sets.flag),
                            pval(
                                    fibertypes.results.mrna %>%
                                            filter(timepoint == "Week 12",
                                                   fibertype == "type2x",
                                                   interaction == "inter") %>%
                                            dplyr::select(p.value) %>%
                                            as.numeric(), 
                                    flag = TRUE, 
                                    sign = sets.flag)), 
                    tip_length = tip_length, 
                    family="sans",
                    textsize = star_size, 
                    size = bar_size,
                    vjust = condition_vjust) +
        geom_signif(y_position = c(22), 
                    xmin = c(2), 
                    xmax = c(2), 
                    annotation=c(pval(
                            fibertypes.results.mrna %>%
                                    filter(timepoint == "Week 2",
                                           fibertype == "type2x",
                                           interaction == "main") %>%
                                    dplyr::select(p.value) %>%
                                    as.numeric(), 
                            flag = TRUE, 
                            sign = time.flag)), 
                    tip_length = tip_length, 
                    family="sans",
                    textsize = time_size-1.5, 
                    size = bar_size) +
        scale_fill_manual(values = fill_scale) +
        gene_fam_theme() +
        xlab(" ") +
        scale_y_continuous(expand = c(0,0), breaks = c(0, 10, 20, 30, 40), limits = c(0, 40)) +
        annotate("text", x = x.fib.text, y = 30 + (40-30)/1.25, label = "Type IIX", size = 2.5)

### Genefam type 2a

genefam2a <- gene.family.dat %>%
        ungroup() %>%
        mutate(fibertype = factor(fibertype, 
                                  levels = c("type1", "type2a", "type2x"))) %>%
        filter(fibertype == "type2a") %>%
        ungroup() %>%
        mutate(timepoint = factor(timepoint, levels = c("w0", "w2pre", "w12"),
                                  labels = c("Week 0", "Week 2", "Week 12")),
               sets = factor(sets, levels = c("single", "multiple"),
                             labels = c("Single-\nset", "Multiple-\nset")), 
               percentage = percentage * 100,
               q10 = q10 * 100,
               q90 = q90 * 100) %>%
        ggplot(aes(timepoint, percentage, fill = sets)) + 
        geom_errorbar(aes(ymin = q10, ymax = q90), width = 0, position = pos, size = errorbar.size) +
        geom_point(shape = 21, position = pos, size = 2) +
        geom_signif(y_position = c(78), 
                    xmin = c( 2.9), 
                    xmax = c( 3.1), 
                    annotation=c(   pval(
                            fibertypes.results.mrna %>%
                                    filter(timepoint == "Week 12",
                                           fibertype == "type2a",
                                           interaction == "inter") %>%
                                    dplyr::select(p.value) %>%
                                    as.numeric(), 
                            flag = TRUE, 
                            sign = sets.flag)), 
                    tip_length = tip_length, 
                    family="sans",
                    textsize = star_size, 
                    size = bar_size,
                    vjust = condition_vjust) +
        geom_signif(y_position = c(82), 
                    xmin = c(2), 
                    xmax = c(2), 
                    annotation=c(   pval(
                            fibertypes.results.mrna %>%
                                    filter(timepoint == "Week 2",
                                           fibertype == "type2a",
                                           interaction == "main") %>%
                                    dplyr::select(p.value) %>%
                                    as.numeric(), 
                            flag = TRUE, 
                            sign = time.flag)), 
                    tip_length = tip_length, 
                    family="sans",
                    textsize = time_size-1.5, 
                    size = bar_size) +
        scale_fill_manual(values = fill_scale) +
        gene_fam_theme() +
        xlab("Week") +
        scale_y_continuous(expand = c(0,0), breaks = c(35, 50, 65, 80, 95), limits = c(35, 95)) +
        annotate("text", x = x.fib.text, y = 80 + (95-80)/1.25, label = "Type IIA", size = 2.5)



genefam1 <- gene.family.dat %>%
        ungroup() %>%
        mutate(fibertype = factor(fibertype, 
                                  levels = c("type1", "type2a", "type2x"))) %>%
        filter(fibertype == "type1") %>%
        ungroup() %>%
        mutate(timepoint = factor(timepoint, levels = c("w0", "w2pre", "w12"),
                                  labels = c("Week 0", "Week 2", "Week 12")),
               sets = factor(sets, levels = c("single", "multiple"),
                             labels = c("Single-\nset", "Multiple-\nset")), 
               percentage = percentage * 100,
               q10 = q10 * 100,
               q90 = q90 * 100) %>%
        ggplot(aes(timepoint, percentage, fill = sets)) + 
        
        geom_errorbar(aes(ymin = q10, ymax = q90), width = 0, position = pos, size = errorbar.size) +
        geom_point(shape = 21, position = pos, size = 2) +
        geom_signif(y_position = c(32, 37), 
                    xmin = c(2, 3), 
                    xmax = c(2, 3), 
                    annotation=c(pval(fibertypes.results.mrna %>%
                                              filter(timepoint == "Week 2",
                                                     fibertype == "type1",
                                                     interaction == "main") %>%
                                              dplyr::select(p.value) %>%
                                              as.numeric(), 
                                      flag = TRUE, 
                                      sign = time.flag),
                                 pval(fibertypes.results.mrna %>%
                                              filter(timepoint == "Week 12",
                                                     fibertype == "type1",
                                                     interaction == "main") %>%
                                              dplyr::select(p.value) %>%
                                              as.numeric(), 
                                      flag = TRUE, 
                                      sign = time.flag)), 
                    tip_length = tip_length, 
                    family="sans",
                    textsize = time_size-1.5, 
                    size = bar_size) +
        gene_fam_theme() + 
        theme(legend.position = "bottom") +
        scale_fill_manual(values = fill_scale) +
        theme(axis.title.y = element_text(size = 8, 
                                          angle = 90,
                                          margin = margin(t = 0, r = 1, b = 0, l = -4))) +
        ylab("GeneFam proportions") +
        xlab(" ") +
        scale_y_continuous(expand = c(0,0), breaks = c(0,10, 20, 30, 40), limits = c(0, 40)) +
        annotate("text", x = x.fib.text, y = 30 + (40-30)/1.25, label = "Type I", size = 2.5)


#### Family plot 

gene_fam <- plot_grid(genefam1 + theme(legend.position = "none"), genefam2a, genefam2x, nrow = 1)
mhc_fam <- plot_grid(fib1, fib2a, fib2x, nrow = 1)

fam_fig <- plot_grid(gene_fam, mhc_fam, ncol = 1, rel_heights = c(1, 1.25))

### Shared legend from 
legend <- get_legend(genefam1 + theme(legend.title = element_blank(),
                                      legend.background = element_rect(color = "white")))



### Fig 5D Myosin heavy chain 2x mRNA #####

pos <- position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8)




img <- ggdraw() + draw_image("./data/immuno/immuno_control_image2.png")


figure3 <- plot_grid(NULL, plot_grid(img, fam_fig, 
                                     plot_grid(myhc2x.mrna.fig, legend.custom.left, nrow = 1, rel_widths = c(1, 0.25)), 
                                     ncol = 1, rel_heights = c(1, 1, 0.6, 0.1)), 
                     rel_widths = c(0.025, 1), ncol = 2) +
        draw_plot_label(label=c("A", "B","C"),
                        x = c(0.03, 0.03, 0.03), 
                        y = c(0.98, 0.61, 0.23),
                        hjust=.5, vjust=.5, size = label.size)



figure3_alt <- plot_grid(NULL, 
                         plot_grid(plot_grid(img, fam_fig, ncol = 1), 
                                   plot_grid(individual_swith_plot,
                                             plot_grid(myhc2x.mrna.fig, legend.custom.left, nrow = 1, rel_widths = c(1, 0.3)), 
                                             ncol = 1, rel_heights = c(1, 0.7)), ncol = 1, rel_heights = c(1, 0.82)), ncol = 2, rel_widths = c(0.1,1)) +
        draw_plot_label(label=c("A", "B","C", "D"),
                        x = c(0.05, 0.05, 0.05, 0.05), 
                        y = c(0.98, 0.72, 0.44, 0.17),
                        hjust=.5, vjust=.5, size = label.size)





# Width of figure = 2x columns 17 cm
# height of figure = 3/4 page = 23 * 0.75 cm

# ggsave("figures/figure5.pdf", plot = figure5, width = 8.5, height = 18, units = "cm", device = "pdf", title = "Hammarstrom et al. 2019 Figure 5")

ggsave("figures/figure3.pdf", plot = figure3_alt, width = 8.5, 
       dpi = 600,
       height = 23, units = "cm", device=cairo_pdf)


#### Save models for summary table 

saveRDS(list(pure2x = type2xpure.m1,
             hybrids = type2xhybrids.m1, 
             prot2x = type2x.m1, 
             prot2a = type2a.m1, 
             prot1 = type1.m1, 
             genefam1 = myhc1x.e,
             genefam2a = myhc2ax.e,
             genefam2x = myhc2xx.e), "./derivedData/models/fibertype.models")






