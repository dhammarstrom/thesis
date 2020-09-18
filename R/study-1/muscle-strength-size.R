##-------------------------------------
## study-1/muscle-size-strength-gains.R
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
dxa <- dxa1 %>%
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
rm <- read_excel("./data/study-1/strength/strengthTests.xlsx") %>%
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
                


# Fit models 
# In study 1, there were no interaction effects between gain and sex. 
# Sex is therefore not modelled as an interaction but as a balancing variable
# for the pre values

# muscle size
dxa.m1 <- lmer(change ~ sets + scale(pre) + sex + (1|subject), 
               data = dxa)
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
        
        mutate(variable = factor(variable, 
                                 levels = c("MR", 
                                            "DXA", 
                                            "legext", 
                                            "legpress", 
                                            "s0", 
                                            "s60", 
                                            "s120", 
                                            "s240"), 
                                 labels = c("Cross sectional area", 
                                            "Fat free mass", 
                                            "Knee-extension 1RM", 
                                            "Leg-press 1RM", 
                                            "Knee extension isometric torque<br>(60&deg;)", 
                                            "Knee extension isokinetic torque<br>(60&deg; sec<sup>-1</sup>)", 
                                            "Knee extension isokinetic torque<br>(120&deg; sec<sup>-1</sup>)",
                                            "Knee extension isokinetic torque<br>(240&deg; sec<sup>-1</sup>)")), 
               variable = fct_rev(variable)) %>%

        
        ggplot(aes(estimate, variable)) + 
        
        geom_vline(xintercept = 0, lty = 2, color = "gray50") +
        
        geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.2) +
        geom_point() +
        scale_x_continuous(limits = c(-2.5, 15), 
                           breaks = c(-2.5, 0, 2.5, 5, 7.5, 10, 12.5, 15), 
                           labels = c("", 0, "", 5, "", 10, "", 15)) +
        labs(x = "Difference between volume conditions (%MOD - %LOW &#x00B1;95%CI)") +
        
        
        
                dissertation_theme() +
        theme(axis.text.y = element_markdown(size = 7), 
              
              axis.title.y = element_blank(), 
              axis.title.x = element_markdown())


saveRDS(muscle_strength_fig, "./data/derivedData/study1_muscle-strength-size/muscle-strength-fig.RDS")




