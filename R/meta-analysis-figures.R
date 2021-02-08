##-------------------------------------
## meta-analysis-figures.R
##
## Title: Meta-analysis figures
## Purpose: Create figures (and tables) based on the meta analysis 
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

# Get coefficients ######

# 


library(tidybayes)
library(tidyverse)
library(brms)

# muscle mass model
m1 <- readRDS("./data/meta-volume/models/muscle_size_full_data.RDS")

# strength model
s1 <- readRDS("./data/meta-volume/models/strength_full_data.RDS")

summary(s1)

mods1 <- rbind(s1 %>%
                       spread_draws(b_weekly_sets, r_study[study,weekly_sets]) %>%
                       filter(weekly_sets != "Intercept") %>%
                       # add the grand mean to the group-specific deviations
                       mutate(mu = b_weekly_sets + r_study) %>%
                       ungroup() %>%
                       mutate(Study = paste0(gsub("[0-9]", "", study), 
                                             " et al. (", 
                                             gsub("[A-z]", "", study), 
                                             ")"), 
                              type = "strength"),
               m1 %>%
                       spread_draws(b_weekly_sets, r_study[study,weekly_sets]) %>%
                       filter(weekly_sets != "Intercept") %>%
                       # add the grand mean to the group-specific deviations
                       mutate(mu = b_weekly_sets + r_study) %>%
                       ungroup() %>%
                       mutate(Study = paste0(gsub("[0-9]", "", study), 
                                             " et al. (", 
                                             gsub("[A-ö]", "", study), 
                                             ")"), 
                              type = "mass")) 

saveRDS(mods1, "./data/derivedData/get_coefs/mods1.RDS")

mods1 %>%
        filter(type == "strength") %>%
        ggplot(aes(x = mu, y =reorder(Study, mu, .fun='mean'))) +
        geom_vline(xintercept = fixef(s1)[2, 1], color = "white", size = 1) +
        geom_vline(xintercept = fixef(s1)[2, 3:4], color = "white", linetype = 2) +
        stat_halfeye(.width = .95, size = 1) +
        labs(x = expression(italic("Cohen's d")),
             y = NULL) +
        theme(panel.grid   = element_blank(),
              axis.ticks.y = element_blank(),
              axis.text.y  = element_text(hjust = 0)) 






s2 <- readRDS("./data/derivedData/brms_models/strength_full_data_s2.RDS")


summary(s2)

s2 %>%
        spread_draws(b_weekly_sets, r_study[study,weekly_sets]) %>%
        filter(weekly_sets != "Intercept") %>%
        # add the grand mean to the group-specific deviations
        mutate(mu = b_weekly_sets + r_study) %>%
        ungroup() %>%
        mutate(Study = paste0(gsub("[0-9]", "", study), 
                              " et al. (", 
                              gsub("[A-z]", "", study), 
                              ")"), 
               type = "strength")

m3 <- readRDS("./data/derivedData/brms_models/muscle_size_full_data_m3.RDS")


summary(m3)

get_variables(m3)


# Tendency of smaller effect in upper body compared to the lower body
# the figure below mixes different estimates (interaction and estimated means)

m3 %>%
        recover_types() %>%
        
        spread_draws(`b_weekly_sets`, 
                     `b_weekly_sets:body_halfupper`, 
                     `b_weekly_sets:body_halfwhole`) %>%
        pivot_longer(names_to = "var", values_to = "val", cols =  `b_weekly_sets:body_halfupper`: 
                             `b_weekly_sets:body_halfwhole`) %>%
        
        # add the grand mean to the group-specific deviations
        mutate(mu = val + b_weekly_sets) %>%
        ggplot(aes(x = mu, y = var)) +
        geom_vline(xintercept = fixef(m3)[2, 1], color = "blue", size = 1) +
        geom_vline(xintercept = fixef(m3)[2, 3:4], color = "lightblue", linetype = 2) +
        stat_halfeye(.width = .95, size = 1) +
        labs(x = expression(italic("Cohen's d")),
             y = NULL) +
        
        annotate("text", x = 0.04, y = 1, 
                 label = paste0(sprintf("%.3f", fixef(m3)[5, 1]), " [", 
                                sprintf("%.3f", fixef(m3)[5, 3]), ", ",
                                sprintf("%.3f", fixef(m3)[5, 4]), "]")) +
        
        annotate("text", x = 0.04, y = 2, 
                 label = paste0(sprintf("%.3f", fixef(m3)[4, 1]), " [", 
                                sprintf("%.3f", fixef(m3)[4, 3]), ", ",
                                sprintf("%.3f", fixef(m3)[4, 4]), "]")) +
        
        
        theme(panel.grid   = element_blank(),
              axis.ticks.y = element_blank(),
              axis.text.y  = element_text(hjust = 0)) 







