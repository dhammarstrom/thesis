##-------------------------------------
## study1-strength-mass-gains.R
##
##
## Title: 
## Purpose: Combine strength and muscle mass gains for presentation 
## Author:
##
##
##
##-------------------------------------
## Notes: Combining all data on muscle strength and function from study 1
# to compare effects of multiple and single set training.
# Below all models are fitted and a figure is saved with the estimates.
#       
#
## ------------------------------------


source("./R/libraries.R")
source("./R/themes.R")



### Read data on included/excluded, body composition and cybex data
dxa1 <- read_excel("./data/study-1/body-composition/bodycomp_DXA.xlsx") 
participants <- read.csv("./data/study-1/body-composition/oneThreeSetLeg.csv", sep=";") %>%
        pivot_longer(names_to = "sets", values_to = "leg", cols = multiple:single) %>%
        print()




# Fat free mass 
dxa1 <- dxa1 %>%
        dplyr::select(subject, timepoint, L = lean.left_leg, R = lean.right_leg) %>%
        pivot_longer(names_to = "leg", values_to = "ffm", cols = L:R) %>%
        pivot_wider(names_from = timepoint, values_from = ffm) %>%
        inner_join(participants) %>%
        filter(include == "incl") %>%
        # Calculate change score 
        mutate(change = 100 * ((post/pre-1))) %>%
        mutate(sets = factor(sets, levels = c("single", "multiple"))) %>%
        
        
        print()


# Cybex data 
load("./data/study-1/strength/cybex.rda")


cybex2 <- data.frame(cybex) %>%
        mutate(leg = if_else(leg == "Left", "L", "R")) %>%
        inner_join(participants) %>%
        filter(include == "incl", 
               timepoint %in% c("pre","session1", "post")) %>%
        
        filter(!(subject=="FP6" & timepoint == "post" & speed == 240)) %>% # This removes one record, exclusion due to error in measurement, see protocol FP6
        
        
        
        dplyr::select(subject:leg, sex, sets, speed, torque.max) %>%
        pivot_wider(names_from = timepoint, values_from = torque.max) %>%
        rowwise() %>%
        mutate(pre = max(c(pre, session1), na.rm = TRUE)) %>%
        mutate(change = 100 * ((post/pre-1))) %>%
        mutate(sets = factor(sets, levels = c("single", "multiple")), 
               speed = paste0("s", speed)) %>%
        print()

### Temporary for meta analysis
temp.iso <- data.frame(cybex) %>%
        mutate(leg = if_else(leg == "Left", "L", "R")) %>%
        inner_join(participants) %>%
        filter(include == "incl", 
               timepoint %in% c("pre","session1", "post")) %>%
        
        filter(!(subject=="FP6" & timepoint == "post" & speed == 240)) %>% # This removes one record, exclusion due to error in measurement, see protocol FP6
        
        
        dplyr::select(subject:leg, sex, sets, speed, torque.max) %>%
        pivot_wider(names_from = timepoint, values_from = torque.max) %>%
        rowwise() %>%
        mutate(pre = max(c(pre, session1), na.rm = TRUE)) %>%
        dplyr::select(subject:leg, sex, sets, speed, pre, post) %>%
        pivot_longer(names_to = "timepoint", values_to = "torque.max", cols = pre:post) %>%
        group_by(timepoint, sets,sex, speed) %>%
        summarise(m = mean(torque.max, na.rm = TRUE), 
                  s = sd(torque.max, na.rm = TRUE)) %>%
        pivot_wider(names_from = timepoint, values_from = c(m, s)) %>%
        dplyr::select(sets, sex, speed, m_pre, s_pre, m_post, s_post) %>%
        print()



# mr data 
source("./R/study-1/mr-and-dxa-data.R")

mr <- mr.results %>%
        dplyr::select(subject, leg, timepoint, CSA.avg) %>%
        inner_join(participants) %>%
        filter(include == "incl") %>%
        pivot_wider(names_from = timepoint, values_from = CSA.avg) %>%
        mutate(change = 100 * ((post/pre-1))) %>%
        mutate(sets = factor(sets, levels = c("single", "multiple"))) %>%
        
        print()

# 1RM data 


# Full data sets for use in pre-post analysis
rm <- read_excel("./data/study-1/strength/strengthTests.xlsx", na = "NA") %>%
        filter(!(exercise %in% c("benchpress", "legcurl"))) %>%
        inner_join(participants) %>%
        filter(include == "incl") %>%
        
        dplyr::select(subject,sex,sets, exercise, timepoint, load) %>%
        mutate(load = as.numeric(load),
               sets = factor(sets, levels = c("single", "multiple"))) %>%
        spread(timepoint, load) %>%
        rowwise() %>%
        mutate(pre = max(c(pre, session1))) %>%
        dplyr::select(subject, sex, sets, exercise, pre, post) %>%
        mutate(change = 100 * ((post/pre-1))) %>%
        
        mutate(sets = factor(sets, levels = c("single", "multiple"))) %>%
        
        print()

# For meta                
rm %>%
        dplyr::select(subject:post) %>%
        pivot_longer(names_to = "timepoint", values_to = "kg", cols = pre:post) %>%
        group_by(timepoint, sets,sex, exercise) %>%
        summarise(m = mean(kg, na.rm = TRUE), 
                  s = sd(kg, na.rm = TRUE)) %>%
        pivot_wider(names_from = timepoint, values_from = c(m, s)) %>%
        dplyr::select(sets, sex, exercise, m_pre, s_pre, m_post, s_post) %>%
        print()
# write_csv(path = "./temp.csv")



# Fit models 
# In study 1, there were no interaction effects between gain and sex. 
# Sex is therefore not modelled as an interaction but as a balancing variable
# for the pre values

# muscle size
dxa.m1 <- lmer(change ~ sets + scale(pre) + sex + (1|subject), 
               data = dxa1)
mr.m1 <- lmer(change ~ sets + scale(pre) + sex + (1|subject), 
              data = mr)


# Muscle strength

legext.m1 <- lmer(change ~ sets + scale(pre) + sex + (1|subject), 
                  data = rm[rm$exercise == "legext", ])

legpress.m1 <- lmer(change ~ sets + scale(pre) + sex + (1|subject), 
                    data = rm[rm$exercise == "legpress", ])

# Cybex 
s0.m1 <- lmer(change ~ sets + scale(pre) + sex + (1|subject), 
              data = cybex2[cybex2$speed == "s0", ])

s60.m1 <- lmer(change ~ sets + scale(pre) + sex + (1|subject), 
               data = cybex2[cybex2$speed == "s60", ])
s120.m1 <- lmer(change ~ sets + scale(pre) + sex + (1|subject), 
                data = cybex2[cybex2$speed == "s120", ])
s240.m1 <- lmer(change ~ sets + scale(pre) + sex + (1|subject), 
                data = cybex2[cybex2$speed == "s240", ])



comd.df <- rbind(
        cbind(summary(mr.m1)$coef, data.frame(confint(mr.m1))[3:6,]) %>%
                mutate(variable = "MR", 
                       coef = row.names(.)), 
        cbind(summary(dxa.m1)$coef, data.frame(confint(dxa.m1))[3:6,]) %>%
                mutate(variable = "DXA", 
                       coef = row.names(.)), 
        cbind(summary(legext.m1)$coef, data.frame(confint(legext.m1))[3:6,]) %>%
                mutate(variable = "legext", 
                       coef = row.names(.)), 
        cbind(summary(legpress.m1)$coef, data.frame(confint(legpress.m1))[3:6,]) %>%
                mutate(variable = "legpress", 
                       coef = row.names(.)),  
        cbind(summary(s0.m1)$coef, data.frame(confint(s0.m1))[3:6,]) %>%
                mutate(variable = "s0", 
                       coef = row.names(.)),  
        cbind(summary(s60.m1)$coef, data.frame(confint(s60.m1))[3:6,]) %>%
                mutate(variable = "s60", 
                       coef = row.names(.)),  
        cbind(summary(s120.m1)$coef, data.frame(confint(s120.m1))[3:6,]) %>%
                mutate(variable = "s120", 
                       coef = row.names(.)),  
        cbind(summary(s240.m1)$coef, data.frame(confint(s240.m1))[3:6,]) %>%
                mutate(variable = "s240", 
                       coef = row.names(.))) %>%
        
        
        dplyr::select(variable, 
                      coef, 
                      estimate = Estimate, 
                      SE = `Std. Error`, 
                      t.val = `t value`, 
                      lwr = X2.5.., 
                      upr = X97.5..) %>%
        print()



muscle_strength_fig <- comd.df %>%
        filter(coef == "setsmultiple") %>%
        
        mutate(variable_group = if_else(variable %in% c("MR", "DXA"), "muscle", "strength"),
               variable = factor(variable, 
                                 levels = c("MR", 
                                            "DXA", 
                                            "legext", 
                                            "legpress", 
                                            "s0", 
                                            "s60", 
                                            "s120", 
                                            "s240"), 
                                 labels = c("Cross sectional area (MRI)", 
                                            "Fat free mass (DXA)", 
                                            "Knee-extension 1RM", 
                                            "Leg-press 1RM", 
                                            "Knee extension isometric torque<br>(60&deg;)", 
                                            "Knee extension isokinetic torque<br>(60&deg; sec<sup>-1</sup>)", 
                                            "Knee extension isokinetic torque<br>(120&deg; sec<sup>-1</sup>)",
                                            "Knee extension isokinetic torque<br>(240&deg; sec<sup>-1</sup>)")), 
               variable = fct_rev(variable)) %>%
        
        
        ggplot(aes(estimate, variable, fill = variable_group)) + 
        
        geom_vline(xintercept = 0, lty = 2, color = "gray50") +
        
        geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0, size = line.size) +
        geom_point(shape = 21, size = 2.5) +
        scale_x_continuous(limits = c(-2.5, 15), 
                           expand = c(0, 0),
                           breaks = c(-2.5, 0, 2.5, 5, 7.5, 10, 12.5, 15), 
                           labels = c("", 0, "", 5, "", 10, "", 15)) +
        labs(x = "Difference between volume conditions<br>(%-gain Multiple-sets - %-gain Single-set &#x00B1;95%CI)") +
        
        scale_fill_manual(values = c(group.study.color[5], group.study.color[3])) +
        
        dissertation_theme() +
        theme(axis.text.y = element_markdown(size = 7), 
              legend.position = "none",
              axis.title.y = element_blank(), 
              axis.title.x = element_markdown())




############# Meta analysis figures ############################


library(tidybayes)
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
                                             gsub("[A-ø]", "", study), 
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
                                             gsub("[A-ø]", "", study), 
                                             ")"), 
                              type = "mass")) 

## get fixed effects

s_fx <- s1 %>%
        spread_draws(b_weekly_sets) %>%
        mutate(mu = b_weekly_sets, 
               Study = "Summary", 
               type = "strength") %>%
                dplyr::select(mu, Study, type) %>%
        print()

m_fx <- m1 %>%
        spread_draws(b_weekly_sets) %>%
        mutate(mu = b_weekly_sets, 
               Study = "Summary", 
               type = "mass") %>%
        dplyr::select(mu, Study, type) %>%
        print()



summary_table <- bind_rows(s_fx %>%
        mutate(type = "strength"), 
        m_fx %>%
                mutate(type = "mass")) %>%
        group_by(type) %>%
        summarise(median = median(mu), 
                  cri.lwr = quantile(mu, 0.025), 
                  cri.upr = quantile(mu, 0.975))
        

saveRDS(summary_table, "./data/derivedData/meta/naive_models.RDS")




meta_fig <- mods1 %>%
        dplyr::select(mu, Study, type) %>%
        rbind(m_fx, s_fx) %>%
        
        
        mutate(Study = reorder(Study, mu, .fun = 'mean'), 
               Study = fct_relevel(Study, "Summary"), 
               type = factor(type, levels = c("mass", "strength"), 
                             labels = c("Muscle mass", 
                             "Muscle strength"))) %>%
        
        ggplot(aes(x = mu, y = Study, fill = type)) +
        
        geom_vline(xintercept = 0, lty = 2, color = "gray50") +
        
        stat_cdfinterval(.width = c(.95),
                        
                         size = 0.1) +
        labs(x = "Increased ES per weekly number of sets<br>(Hedges' <i>g</i>, 95% CrI)",
             y = NULL) +
        scale_color_manual(values = c(group.study.color[5], group.study.color[2])) +
        scale_fill_manual(values = c(group.study.color[5], group.study.color[2])) +
 
        dissertation_theme() +
        
        theme(panel.grid   = element_blank(),
              legend.position = "none",
              axis.title.x = element_markdown(size = 7),
              axis.ticks.y = element_blank(),
              axis.text.y  = element_text(hjust = 0), 
              strip.text = element_text(size = 8, hjust = 0),
              strip.background = element_blank()) +
        facet_grid(. ~ type, scales = "free")



mass_fig_1 <- mods1 %>%
        filter(type == "mass") %>%
        dplyr::select(mu, Study, type) %>%
        rbind(m_fx) %>%
        
        mutate(Study = reorder(Study, mu, .fun = 'median'), 
               Study = fct_relevel(Study, "Summary")) %>%
        
        ggplot(aes(x = mu, y = Study)) +
        
        geom_vline(xintercept = 0, lty = 2, color = "gray50") +
        
        stat_halfeye(.width = .95, 
                     size = 0.5, 
                     shape = 21,
                     fill = group.study.color[5]) +
        
        labs(x = "Increased ES per weekly number<br>of sets (Hedges' <i>g</i>, 95% CrI)",
             y = NULL) +
        
        dissertation_theme() +
        
        theme(panel.grid   = element_blank(),
              axis.title.x = element_markdown(size = 7),
              axis.ticks.y = element_blank(),
              axis.text.y  = element_text(hjust = 0)) 



### Strength model 1


strength_fig_1 <- mods1 %>%
        filter(type == "strength") %>%
        dplyr::select(mu, Study, type) %>%
        rbind(s_fx) %>%
        
        mutate(Study = reorder(Study, mu, .fun = 'median'), 
               Study = fct_relevel(Study, "Summary")) %>%
        
        ggplot(aes(x = mu, y = Study)) +
        
        geom_vline(xintercept = 0, lty = 2, color = "gray50") +
        
        stat_halfeye(.width = .95, 
                     shape = 21,
                     size = 0.5, 
                     fill = group.study.color[3]) +
        labs(x = "Increased ES per weekly number<br>of sets (Hedges' <i>g</i>, 95% CrI)",
             y = NULL) +
        
        dissertation_theme() +
        
        theme(panel.grid   = element_blank(),
              axis.title.x = element_markdown(size = 7),
              axis.ticks.y = element_blank(),
              axis.text.y  = element_text(hjust = 0)) 






## Combine figure 

muscle_strength_fig_comb  <- plot_grid(muscle_strength_fig, 
                                       NULL,
                                       meta_fig, 
          ncol = 1, 
          rel_heights = c(0.4,0.02, 0.6)) + 
       
        draw_plot_label(label=c("a", "b"),
                        x = c(0.02, 0.02), 
                        y = c(0.99, 0.59),
                        hjust=.5, vjust=.5, size = 12)












saveRDS(muscle_strength_fig_comb, "./data/derivedData/study1_muscle-strength-size/muscle-strength-fig.RDS")




