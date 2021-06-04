##-------------------------------------
## study2-rna-muscle-mass-fig.R
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
library(brms)
library(readxl)


#### Figure 4B Leave one out analysis #########################



loo.results <- readRDS("./data/study-2/ubf-tot-rna-model/leave-one-out.RDS")

leg <- read_excel("./data/study-2/leg_randomization.xlsx")


res <- loo.results$loo.sample %>%
        inner_join(leg) %>%
        
        mutate(Participant = factor(participant, 
                                    levels = paste0("P", c(seq(1:7), 19, 21, 22, 23)), 
                                    labels = paste0("P", seq(1:11))), 
               Participant = fct_rev(Participant)) %>%
        
        
        pivot_longer(names_to = "stat", values_to = "estimates", estimate:upr) %>%
        print()

res2 <- loo.results$loo.participant %>%
        inner_join(leg) %>%
        mutate(Participant = factor(participant, 
                                    levels = paste0("P", c(seq(1:7), 19, 21, 22, 23)), 
                                    labels = paste0("P", seq(1:11))), 
               Participant = fct_rev(Participant)) %>%
        
        
        filter(model == "m3", coef != "sex") %>%
        mutate(coef = factor(coef, levels = c("slope", "intercept"), 
                             labels = c("Total RNA increase\nper session", 
                                        "Average total RNA\n(Session 6)"))) %>%
        print()



loo_panel <- res %>%
        filter(stat != "t_val", model == "m3", coef != "sex") %>%
        
        mutate(coef = factor(coef, levels = c("slope", "intercept"), 
                             labels = c("Total RNA increase\nper session", 
                                        "Average total RNA\n(Session 6)"))) %>%
        
        ggplot(aes(estimates, Participant, fill = stat)) +
        
        
        geom_vline(xintercept = 0, lty = 2, color = "gray50") +
        
        #  geom_density_ridges(color = "white") +
        
        geom_point(position = position_nudge(y = 0.2), shape = 21, alpha = 0.4) +
        
        facet_wrap( ~ coef, scales = "free") +
        
        scale_fill_manual(values = c(group.study.color[5], group.study.color[4], group.study.color[4])) +
        
        geom_point(data = res2, aes(estimate, Participant, fill = NULL)) +
        geom_errorbarh(data = res2, aes(fill = NULL, x = estimate, y = Participant, 
                                        xmax =  upr, xmin = lwr), 
                       height = 0) + 
        
        labs(x = "Estimated change in muscle thickness (mm)\nper one unit change in predictor", 
             y = "Participant") +
        
       dissertation_theme() + 
        theme(strip.background = element_rect(color = "white", fill = "white"),
              strip.text = element_text(size = 8, hjust = 0),
              axis.title.y = element_markdown(size = 7), 
              legend.position = "none")


### Panel of predicted values against observed 



# Read data frame of predicted values
predict_df <- readRDS("./data/study-2/ubf-tot-rna-model/predict_df.RDS")

# get predictions from the model
hyp.pre.m1 <- readRDS("./data/study-2/ubf-tot-rna-model/hyp_pre_m3.RDS")
summary(hyp.pre.m1)
effects_m1 <- conditional_effects(hyp.pre.m1)

mean_diff <- (effects_m1$sex[2, 9] - effects_m1$sex[1, 9]) / 2


prediction_slope <- predict_df %>%
        
        inner_join(leg) %>%
        mutate(cond = factor(cond, levels = c("const", "var"), 
                             labels = c("Constant volume", 
                                        "Variable volume")))  %>%
        
        group_by(participant, leg, cond) %>%
        summarise(slope = mean(slope), 
                  mm_incr = mean(mm_incr)) %>%
        ungroup() %>%
        mutate(Participant = factor(participant, 
                                    levels = paste0("P", c(seq(1:7), 19, 21, 22, 23)), 
                                    labels = paste0("P", seq(1:11))), 
               Participant = fct_rev(Participant)) %>%
        ggplot(aes(slope, mm_incr, fill = cond))  + 
        geom_ribbon(data =   effects_m1$slope, 
                    aes(ymin = lower__ + mean_diff, 
                        ymax = upper__ + mean_diff), 
                    fill = "gray85") +
        geom_line(data =   effects_m1$slope, 
                  aes(slope, estimate__ + mean_diff, fill = NULL)) +
        geom_point(shape = 21, size = 2) + 
        geom_text_repel(aes(label = Participant), 
                        size = 2.5) +
        
        scale_y_continuous(limits = c(-3, 5), 
                           breaks = c(-3, -2, -1, 0, 1, 2, 3, 4, 5), 
                           labels = c("", -2, "", 0, "", 2, "", 4, "")) +
        
        scale_fill_manual(values = c(group.study.color[1], group.study.color[2])) +
        labs(x = "Total RNA increase per session (%)", 
             y = "Muscle thickness change (mm)") +
        dissertation_theme() + 
        theme(axis.title.y = element_markdown(size = 7), 
              legend.position = c(0.25, 0.95),
              legend.title = element_blank(),
              legend.margin = margin(c(0, 1, 0, 0)),
              legend.key.height = unit(0.2, "cm"),
              legend.key.width = unit(0.2, "cm"),
              legend.text = element_text(size = 7, margin = margin(t = 0, b = 0,r = 0,l= 0, unit = "pt")))



prediction_intercept <- predict_df %>%
        
        inner_join(leg) %>%
        mutate(cond = factor(cond, levels = c("const", "var"), 
                             labels = c("Constant volume", 
                                        "Variable volume")))  %>%
        
        
        group_by(participant, leg, cond) %>%
        summarise(intercept = mean(intercept), 
                  mm_incr = mean(mm_incr)) %>%
        ungroup() %>%
        mutate(Participant = factor(participant, 
                                    levels = paste0("P", c(seq(1:7), 19, 21, 22, 23)), 
                                    labels = paste0("P", seq(1:11))), 
               Participant = fct_rev(Participant)) %>%
        ggplot(aes(intercept, mm_incr, fill = cond))  + 
        geom_ribbon(data =   effects_m1$intercept, 
                    aes(ymin = lower__ + mean_diff, 
                        ymax = upper__ + mean_diff), 
                    fill = "gray85") +
        geom_line(data =   effects_m1$intercept, 
                  aes(intercept, estimate__ + mean_diff, fill = NULL)) +
        geom_point(shape = 21, size = 2) + 
        geom_text_repel(aes(label = Participant), 
                        size = 2.5) +
        scale_fill_manual(values = c(group.study.color[1], group.study.color[2])) +
        labs(x = "Average Total RNA at Session 6\n(Standard deviations from the mean)", 
             y = "Muscle thickness change (mm)") +
        
        scale_y_continuous(limits = c(-3, 5), 
                           breaks = c(-3, -2, -1, 0, 1, 2, 3, 4, 5), 
                           labels = c("", -2, "", 0, "", 2, "", 4, "")) +
        dissertation_theme() + 
        theme(axis.title.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.line.y = element_blank(),
              axis.text.y = element_blank(),
              legend.position = "none")


# Combining data to show estimates of each predictor when leavining one out


combined_df <- readRDS("./data/study-2/ubf-tot-rna-model/combined_df.RDS")


# A exploitative plot of linear fits for each participant (RNA to number of sessions).

rna_to_time_estimates <- combined_df %>%
        filter(!(time %in% c("post1w"))) %>% # Removing the de-training estimate as doeas not represent training induced increase.
        
        mutate(Participant = factor(participant, 
                                    levels = paste0("P", c(seq(1:7), 19, 21, 22, 23)), 
                                    labels = paste0("P", seq(1:11)))) %>%
        
        inner_join(leg) %>%
        mutate(cond = factor(cond, levels = c("const", "var"), 
                             labels = c("Constant volume", 
                                        "Variable volume")))  %>%
        
        # filter(!(participant == "P21" & time == "S12" & leg == "R")) %>%
        
        mutate(time.c = time.c) %>% # this centers number of sessions and the intercept becomes
        # the estimate total RNA (on log scale) in the mid of the training intervention.
        
        ggplot(aes(time.c, rna.mg, group = paste(participant, cond), 
                   color = cond, fill = cond)) + 
        geom_point(shape = 21, size = 1.5) + 
        geom_smooth(method = "lm", se = FALSE, size = 0.75) +
        facet_wrap(~ Participant) + 
        
        scale_x_continuous(limits = c(0, 12), 
                           breaks = c(0, 3, 6, 9, 12), 
                           labels = c(0, "", 6, "", 12)) +
        
        scale_fill_manual(values = c(group.study.color[1], group.study.color[2])) +
        scale_color_manual(values = c(group.study.color[1], group.study.color[2])) +
        
        
        labs(y = "Total RNA (ng mg <sup>-1</sup>)", 
             x = "Sessions") +
        
        dissertation_theme() + 
        theme(strip.background = element_rect(color = "white", fill = "white"),
              strip.text = element_text(size = 8, hjust = 0),
              axis.title.y = element_markdown(size = 7), 
              legend.position = "none")




# Combine all elements of plot 



### Whole page width 
tot_rna_muscle_thickness <- plot_grid( 
        plot_grid(NULL, prediction_slope, prediction_intercept, align = "h", rel_widths = c(0.02, 1, 1), ncol = 3), 
        plot_grid(plot_grid(NULL, rna_to_time_estimates, NULL, ncol = 3, rel_widths = c(0.2, 1, 0.2)),
                  
                  loo_panel, ncol = 1),
        ncol = 1, 
        rel_heights = c(1, 1.6)) +
        
        
        draw_plot_label(label=c("a",  "",  "b", "c"),
                        x =   c(0.01, 0.53, 0.1, 0.02), 
                        y =   c(0.98, 0.98, 0.6, 0.25),
                        hjust=.5, vjust=.5, size = label.size)


saveRDS(tot_rna_muscle_thickness, "./data/derivedData/study2-rna-muscle-mass-fig/tot_rna_muscle_thickness.RDS")

