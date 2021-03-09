##-------------------------------------
## study2-muscle-mass-strength.R
##
## Title:
## Purpose: 
## Author:
##
##
##
##-------------------------------------
## Notes:
# Muscle mass and strength in Study II.
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




comp_strength <-  readRDS("./data/study-2/strength-bayes/comp_strength.RDS")



######### Strength plot ###############################

strength <- read_excel("./data/study-2/tr010_humac.xlsx") %>%
        inner_join(read_excel("./data/study-2/leg_randomization.xlsx")) %>%
        mutate(time = if_else(timepoint %in% c("B1", "B2", "fam"), "baseline", timepoint), 
               time = if_else(time == "post_ctrl", "post", time),
               cond = if_else(cond == "ctrl_leg", "ctrl", cond)) %>%
        dplyr::select(participant, sex, time, leg, group,cond, isokinetic_torque, isometric_torque) %>%
        group_by(participant,sex,leg, time, group, cond) %>%
        
        # For plotting purposes, the baseline averaged over all attempts
        
        summarise(isok = mean(isokinetic_torque, na.rm = TRUE), 
                  isom = mean(isometric_torque, na.rm = TRUE)) %>%
        pivot_longer(names_to = "type", values_to= "torque", cols = isok:isom) %>%
        pivot_wider(names_from = time, values_from = torque) %>%
        mutate(post = post - baseline, 
               post1w = post1w - baseline) %>%
        pivot_longer(names_to = "time", values_to = "change", cols = post:post1w) %>%
        filter(!(is.na(change))) %>%
        group_by(type) %>%
        mutate(baseline = baseline - mean(baseline, na.rm = TRUE)) %>%
        mutate(time_cond = paste0(time, "_", cond), 
               time_group = paste0(time, "_", group), 
               time_group = factor(time_group, levels = c("post_con", 
                                                          "post_int", 
                                                          "post1w_int"))) %>%
        print()



raw_scores <- strength %>%
        ggplot(aes(time_group, change, fill = type)) + 
        
        geom_hline(yintercept = 0, color = "gray50", lty = 2) +
        
        geom_point(position = position_jitterdodge(jitter.width = 0.05, dodge.width = 0.25), 
                   shape = 21, size = 1.2, alpha = 0.4) + 
        
        # Isokinetic strength
        
        geom_point(data = comp_strength[1:3,], 
                   aes(time_group, estimate, fill = NULL), 
                   position = position_nudge(x = -0.25), color = "white", size = 0.5, shape = 15) + 
        # Isometric strength
        
        scale_fill_manual(values = c(group.study.color[2], group.study.color[5])) +
        
        geom_point(data = comp_strength[6:8,], 
                   aes(time_group, estimate, fill = NULL), 
                   position = position_nudge(x = 0.25), color = "white", size = 0.5, shape = 15) + 
        
        labs(y = expression(Delta*" Torque (Nm)")) +
        
        dissertation_theme() +
        theme(axis.line.x = element_blank(), 
              axis.text.x = element_blank(), 
              axis.title.x = element_blank(), 
              axis.ticks.x = element_blank(), 
              legend.position = "none")



# Strength compare plot




comp_panel <- comp_strength %>%
        
        filter(time_group %in% c("inter:post_int", "inter:post1w_int")) %>%
        
        mutate(time_group = gsub("inter:", "", time_group)) %>%
        dplyr::select(estimate, lower.CL, upper.CL, type, time_group) %>%
        
        
        rbind(data.frame(estimate = NA, lower.CL = NA, upper.CL = NA, type = c("isom", "isok"),
                         time_group = "post_con")) %>%
        
        mutate(time_group = factor(time_group,  
                                   levels = c("post_con", 
                                              "post_int", 
                                              "post1w_int"), 
                                   labels =c("Control-group", "Training-group\nS12", "Training-group\nDe-train")), 
               type = factor(type, levels = c("isok", "isom"), 
                             labels = c("Isokinetic",
                                        "Isometric"))) %>%
        
        ggplot(aes(time_group, estimate, fill = type)) + 
        
        geom_hline(yintercept = 0, color = "gray90") +
        
        geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL, color = type), width = 0, 
                      position = position_dodge(width = 0.25)) + 
        geom_point(size = 1.5, 
                   shape = 24,
                   position = position_dodge(width = 0.25)) + 
        
        scale_y_continuous(limits = c(-20, 40), 
                           breaks = c(-20, -10, 0, 10, 20, 30, 40), 
                           labels = c(-20, "",   0,  "", 20,  "", 40), 
                           expand = c(0, 0)) +
        
        
        scale_fill_manual(values = c(group.study.color[2], group.study.color[5])) +
        scale_color_manual(values = c(group.study.color[2], group.study.color[5])) +
        
        annotate("text", x = 1.45, y = 23, label = "Isokinetic", 
                 size = 2.5) +
        
        annotate("curve",
                 x = 1.55,
                 xend = 1.85,
                 y = 18,
                 yend = 14,
                 curvature = 0.2,
                 color = "gray50",
                 arrow = arrow(length=unit(0.1,"cm"), type = "closed")) +
        
        
        annotate("text", x = 2.5, y = 22.5, label = "Isometric", 
                 size = 2.5) +
        
        annotate("curve",
                 x = 2.3,
                 xend = 2.15,
                 y = 18,
                 yend = 15,
                 curvature = -0.2,
                 color = "gray50",
                 arrow = arrow(length=unit(0.1,"cm"), type = "closed")) +

        labs(y = expression(Delta*"Training - "*Delta*"Control")) +
        
       dissertation_theme() +
        
        theme(axis.title.x = element_blank(), 
              legend.position = "none", 
              legend.spacing.y = unit(-0.5, 'cm'),
              legend.key.size = unit(0.2, "cm"),
              legend.text = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0.3), size = 7))



### Save UL results as graph
strength_results <- plot_grid(raw_scores, comp_panel, nrow = 2, align = "v")






## Muscle thickness ############################################




comp_us <-  readRDS("./data/study-2/us-freq-bayes/comp_us_bayes.RDS") %>%
        dplyr::select(time_group, estimate, lower.CL, upper.CL)





us_data <- rbind(read_excel("./data/study-2/ultrasound/ultrasound_data.xlsx") %>%
                         inner_join(read_csv("./data/study-2/ultrasound/ultrasound_codekey.csv")) %>%
                         mutate(leg = gsub("VL", "", leg)) %>%
                         inner_join(read_excel("./data/study-2/leg_randomization.xlsx")) %>%
                         dplyr::select(participant, time, leg, sex, cond, code, length) %>%
                         mutate(group = if_else(participant %in% paste("P", 1:7, sep = ""), "experiment", "control")) %>%
                         group_by(participant, time, leg, sex, cond, group) %>%
                         summarise(thickness = mean(length, na.rm = TRUE)) %>%
                         ungroup(), 
                 read_excel("./data/study-2/ultrasound/ultrasound_data_2019.xlsx") %>%
                         inner_join(read_csv("./data/study-2/ultrasound/ultrasound_codekey_2019.csv")) %>%
                         mutate(leg = gsub("VL", "", leg)) %>%
                         inner_join(read_excel("./data/study-2/leg_randomization.xlsx")) %>%
                         dplyr::select(participant, time, leg, sex, cond, code, length) %>%
                         mutate(group = if_else(participant %in% paste("P", c(1:7, 19:23), sep = ""), "experiment", "control")) %>%
                         group_by(participant, time, leg, sex, cond, group) %>%
                         summarise(thickness = mean(length, na.rm = TRUE)) %>%
                         ungroup()) %>%
        
        mutate(time_pp = if_else(time == "post1w", "post", time),
               # The de-training period get its own coefficient
               detrain = if_else(time == "post1w" & group == "experiment", "detrain", "train"),
               # The effect of de training will be added to the model --within-- the intervention group 
               detrain = factor(detrain, levels = c("train", "detrain")), 
               time = factor(time, levels = c("pre", "post", "post1w")), 
               time_pp = factor(time_pp, levels = c("pre", "post"))) %>%
        
        print()

# 





us_temp <- us_data %>%
        dplyr::select(participant:sex, group, thickness) %>%
        pivot_wider(names_from = time, values_from = thickness) %>%
        mutate(post = post - pre, 
               post1w = post1w - pre) %>%
        pivot_longer(names_to = "time", values_to = "change", cols = post:post1w) %>%
        filter(!is.na(change)) %>%
        mutate(group = if_else(group == "experiment", "int", "con"), 
               time_group = paste0(time, "_", group)) %>%
        
        print()





# 
# 
raw_scores <- us_temp %>%
        ggplot(aes(time_group, change)) + 
        
        geom_hline(yintercept = 0, color = "gray50", lty = 2) +
        
        geom_jitter(width = 0.05, size = 1.5, alpha = 0.4, 
                    shape = 21,
                    fill = group.study.color[5]) + 

        
        labs(y = expression(Delta*" Muscle thickness (mm)")) +
        
        dissertation_theme() +
        theme(axis.line.x = element_blank(), 
              axis.text.x = element_blank(), 
              axis.title.x = element_blank(), 
              axis.ticks.x = element_blank())

# scale_y_continuous(limits = c(-8, 8), 
#                    breaks = c(-2.5, 0, 2.5, 5), 
#                    labels = c(-2.5,   0,  2.5, 5), 
#                    expand = c(0, 0))



comp_panel <- comp_us %>%
        
        filter(time_group %in% c("inter:post_int", "inter:post1w_int")) %>%
        
        mutate(time_group = gsub("inter:", "", time_group)) %>%
        
        rbind(data.frame(time_group = "post_con", estimate = NA, lower.CL = NA, upper.CL = NA )) %>%
        
        mutate(time_group = factor(time_group,  
                                   levels = c("post_con", 
                                              "post_int", 
                                              "post1w_int"), 
                                   labels = c("Control-group", "Training-group\nS12", "Training-group\nDe-train"))) %>% 
        
        ggplot(aes(time_group, estimate)) + 
        
        
        geom_hline(yintercept = 0, color = "gray90") +
        
        geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0) + 
        geom_point(size = 2, shape = 24, fill = group.study.color[1]) + 
        
        scale_y_continuous(limits = c(-0.5, 2.5), 
                           breaks = c(-0.5, 0, 0.5, 1, 1.5, 2, 2.5), 
                           labels = c("",   0,  "", 1,  "", 2,  ""), 
                           expand = c(0, 0)) +
        
        labs(y = expression(Delta*"Training - "*Delta*"Control")) +
        
        dissertation_theme() +
        
        theme(axis.title.x = element_blank())


### Save UL results as graph
ul_results <- plot_grid(raw_scores, comp_panel, nrow = 2, align = "v")





### Combine plots #############################

ultra_sound_strength <- plot_grid(ul_results, 
          strength_results, 
          ncol = 2, rel_widths = c(1, 1)) 



saveRDS(ultra_sound_strength, "./data/derivedData/study2-muscle-mass-strength/ultra_sound_strength.RDS")



