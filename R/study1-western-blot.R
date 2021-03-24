##-------------------------------------
## study1-western-blot.R
##
## Title: Western blot analysis in study I
## Purpose: Display results from western blot analysis in study I
## Author:
##
##
##
##-------------------------------------
## Notes:
# Western blot data are presented in the thesis as diffs between conditions with 95% CI
#
#
#
#
#
#
#
## ------------------------------------


# Load libraries and functions

source("./R/libraries.R")
source("./R/themes.R")

#### Read data and combine data frames ####
# read western data
west <- read_excel("./data/study-1/westernResults.xlsx", sheet="run1")
west2 <- read_excel("./data/study-1/westernResults.xlsx", sheet="run2")

# sample setups
sample <- read_excel("./data/study-1/westernResults.xlsx", sheet="samples")%>%
        filter(include=="incl")

# leg training load-partitioning
leg <- read.csv("./data/study-1/oneThreeSetLeg.csv", sep=";")

leg <- leg %>%
        gather(sets, leg, multiple:single)%>%
        filter(include=="incl")

# combine sample and leg sets
sl <- inner_join(sample, leg)


#### Mean center expression data and collect data from all runs ####

# West acute only contains w2pre and w2post

### Extract pan p70 values for normalisation 
west <- west %>%
        mutate(totalprotein=(mean.gray1 + mean.gray2)/2 - (mean.graybg1+ mean.graybg2)/2) %>%
        inner_join(sl) %>%
        dplyr::select(run, gel, well, subject, leg, sets, timepoint, date, 
                      totalprotein, p.p70, t.p70, p.mTOR, t.mTOR, p.s6, t.s6, p.4EBP1, t.4EBP1) %>%
        group_by(date, gel) %>%
        mutate(totalprotein = totalprotein / max(totalprotein, na.rm = TRUE)) %>%
        gather(target, expression, p.p70:t.4EBP1) %>%
        mutate(expression=expression,
               run = "west1") %>%
        ungroup()



west2 <- west2 %>%
        mutate(totalprotein=(mean.gray1+mean.gray2)/2 - (mean.graybg1+ mean.graybg2)/2) %>%
        inner_join(sl) %>%
        dplyr::select(run, gel, well, subject, leg, sets, timepoint, date, 
                      totalprotein, t.mTOR, p.s6, t.s6) %>%
        group_by(date, gel) %>%
        mutate(totalprotein = totalprotein / max(totalprotein, na.rm = TRUE)) %>%
        gather(target, expression, t.mTOR:t.s6) %>%
        mutate(expression=expression,
               run = "west2") %>%
        # This removes bad gels from analysis
        filter(!(subject %in% c("FP21", "FP22", "FP24", "FP25", "FP27", "FP28", "FP29", "FP30", "FP31") & target == "t.mTOR")) %>%
        ungroup()


## Comnbine data from full runs
wb.raw <- rbind(west, west2) 


wb.raw.acute <- wb.raw %>% 
        filter(timepoint %in% c("w2pre", "w2post")) %>%
        spread(target, expression) %>%
        mutate(p.mTOR = p.mTOR/t.mTOR, 
               p.mTOR.tp = p.mTOR/totalprotein,
               t.mTOR = t.mTOR/totalprotein,
               p.p85.tp = p.p70/totalprotein, 
               p.p85 = p.p70/t.p70,
               p.s6.tp = p.s6/totalprotein, 
               p.s6 = p.s6/t.s6,
               t.s6 = t.s6/totalprotein,
               t.p70 = t.p70/totalprotein) %>%
        dplyr::select(run, gel, well, subject, sets, timepoint, 
                      p.mTOR, 
                      p.mTOR.tp,
                      t.mTOR,
                      p.p85.tp, 
                      p.p85,
                      p.s6.tp,
                      p.s6,
                      t.p70,
                      t.s6) %>%
        gather(target, expression, p.mTOR:t.s6) %>%
        mutate(sample = paste0(subject, sets, timepoint, gel, well), 
               timepoint = factor(timepoint, levels = c("w2pre", "w2post")), 
               sets = factor(sets, levels = c("single", "multiple"))) %>%
        dplyr::select(subject, target, sets, timepoint, expression) %>%
        ungroup() %>%
        print()


west.acute <-   read_excel("./data/study-1/westernResults.xlsx", sheet="runAcute") %>%
        mutate(totalprotein=(mean.gray1+mean.gray2)/2 - (mean.graybg1+ mean.graybg2)/2)%>%
        inner_join(sl)%>%
        dplyr::select(run, gel, well, subject, leg, sets, timepoint, date, 
                      totalprotein, p.p70 = p.p70.new, t.p70, t.mTOR, p.mTOR)%>%
        group_by(date, gel)%>%
        mutate(totalprotein = totalprotein / max(totalprotein, na.rm = TRUE)) %>%
        ungroup() %>%
        mutate(p.p70 = p.p70/totalprotein,
               p.mTOR = p.mTOR/t.mTOR, 
               p.mTOR.tp = p.mTOR/totalprotein,
               t.mTOR = t.mTOR/totalprotein) %>%
        gather(target, expression, p.p70:p.mTOR.tp)%>%
        mutate(run = "westacute")%>%
        ungroup() %>%
        filter(target %in% c("p.p70", "t.mTOR", "p.mTOR")) %>%
        mutate(timepoint = factor(timepoint, levels=c("w2pre", "w2post")),
               sets = factor(sets, levels = c("single", "multiple")),
               leg = factor(paste0(subject, sets))) %>%
        data.frame() %>%
        dplyr::select(subject, target, sets, timepoint, expression) %>%
        ungroup() %>%
        print()


wb.raw.acute <- rbind(wb.raw.acute, west.acute) %>%
        group_by(subject, sets, timepoint, target) %>%
        summarise(expression = log(mean(expression, na.rm = TRUE))) %>%
        print()



temp <- wb.raw.acute %>%
        pivot_wider(names_from = timepoint, values_from = expression) %>%
        mutate(change = w2post - w2pre) %>%
        print()



########### Raw data plot for each target


wb.raw.acute %>%
        filter(target %in% c("p.p70", 
                             "p.p85", 
                             "p.mTOR", 
                             "p.s6.tp")) %>%
        group_by(target, timepoint) %>%
        summarise(m = mean(expression), 
                  )





wb.raw.acute %>%
        filter(target %in% c("p.p70", 
                             "p.p85", 
                             "p.mTOR", 
                             "p.s6.tp")) %>%
        ggplot(aes(timepoint, expression, color = sets)) + 
        
        geom_boxplot(position = position_dodge(width = 0.2), 
                     width = 0.2) +
        
        geom_point(position = position_jitterdodge(dodge.width = 0.2, 
                                                   jitter.width = 0.1), 
                   shape = 21, alpha = 0.3) +
        
        facet_grid(target ~ ., scales = "free")






p.p70 <- lmer(change ~ w2pre + sets + (1|subject), 
              data = temp[temp$target == "p.p70", ]) 

p.p85 <- lmer(change ~ w2pre + sets + (1|subject), 
              data = temp[temp$target == "p.p85", ]) 

p.mTOR <- lmer(change ~ w2pre + sets + (1|subject), 
               data = temp[temp$target == "p.mTOR", ]) 

p.s6.tp <- lmer(change ~ w2pre + sets + (1|subject), 
                data = temp[temp$target == "p.s6.tp", ]) 


p.p70.unadjusted <- lmer(change ~  sets + (1|subject), 
              data = temp[temp$target == "p.p70", ]) 

p.p85.unadjusted  <- lmer(change ~  sets + (1|subject), 
              data = temp[temp$target == "p.p85", ]) 

p.mTOR.unadjusted  <- lmer(change ~  sets + (1|subject), 
               data = temp[temp$target == "p.mTOR", ]) 

p.s6.tp.unadjusted  <- lmer(change ~  sets + (1|subject), 
                data = temp[temp$target == "p.s6.tp", ]) 




######### Small figures showing fold change in each group ################


sets.fc <- rbind(data.frame(emmeans(p.p70, specs = ~ sets)) %>%
        mutate(target = "p.p70"), 
      data.frame(emmeans(p.p85, specs = ~ sets)) %>%
              mutate(target = "p.p85"), 
      data.frame(emmeans(p.mTOR, specs = ~ sets)) %>%
              mutate(target = "p.mTOR"), 
      data.frame(emmeans(p.s6.tp, specs = ~ sets)) %>%
              mutate(target = "p.s6.tp")) %>%
        mutate(fold.change = exp(emmean)) %>%
        print()
  


# Settings for annotate in indiviual plots
percent_from_bottom <- 90
text_size_label <- 2.5

      
        
 p.p70.fcfig <- sets.fc %>%
         filter(target == "p.p70") %>%
        ggplot(aes(sets, fold.change, fill = sets)) + 
        geom_bar(stat = "identity", 
                 width = 0.2) +
         
         annotate("richtext", label = "p70-S6K1<sup>Thr389</sup>", 
                  y = (10/100) * percent_from_bottom, x = 0.5, hjust = 0, 
                  size = text_size_label,
                  label.color = NA, 
                  fill = NA) +
         
         scale_y_continuous(limits = c(0,10), 
                            expand = c(0,0), 
                            breaks = c(0, 5, 10)) +
         
         scale_fill_manual(values = c(group.study.color[3], group.study.color[4])) +
         
         theme_bw() +
         theme(panel.grid = element_blank(), 
               panel.border = element_blank(), 
               axis.line.x  = element_line(size = line.size), 
               axis.line.y = element_line(size = line.size), 
               axis.ticks.x = element_line(size = line.size), 
               axis.ticks.y = element_line(size = line.size), 
               axis.text.y =  element_markdown(color = "black", size = 7), 
               
               # Individual plots settings
               axis.title.y = element_blank(),
               axis.text.x =  element_blank(), 
               axis.title.x = element_blank(),
               legend.title = element_blank(), 
               legend.background = element_rect(fill = "white"),
               legend.margin = margin(t = 0, r = 1, b = 1, l = 1, unit = "pt"),
               legend.key = element_rect(fill = "white"),
               legend.position = "none")

 p.p85.fcfig <- sets.fc %>%
         filter(target == "p.p85") %>%
         ggplot(aes(sets, fold.change, fill = sets)) + 
         geom_bar(stat = "identity", 
                  width = 0.2) +
         
         annotate("richtext", label = "p85-S6K1<sup>Thr412</sup>", 
                  y = (2/100) * percent_from_bottom, x = 0.5, hjust = 0, 
                  size = text_size_label,
                  label.color = NA, 
                  fill = NA) +
         
         scale_y_continuous(limits = c(0,2), 
                            expand = c(0,0), 
                            breaks = c(0, 1, 2)) +
         
         scale_fill_manual(values = c(group.study.color[3], group.study.color[4])) +
         
         theme_bw() +
         theme(panel.grid = element_blank(), 
               panel.border = element_blank(), 
               axis.line.x  = element_line(size = line.size), 
               axis.line.y = element_line(size = line.size), 
               axis.ticks.x = element_line(size = line.size), 
               axis.ticks.y = element_line(size = line.size), 
               axis.text.y =  element_markdown(color = "black", size = 7), 
               
               # Individual plots settings
               axis.title.y = element_blank(),
               axis.text.x =  element_blank(), 
               axis.title.x = element_blank(),
               legend.title = element_blank(), 
               legend.background = element_rect(fill = "white"),
               legend.margin = margin(t = 0, r = 1, b = 1, l = 1, unit = "pt"),
               legend.key = element_rect(fill = "white"),
               legend.position = "none")

 p.mtor.fcfig <- sets.fc %>%
         filter(target == "p.mTOR") %>%
         ggplot(aes(sets, fold.change, fill = sets)) + 
         geom_bar(stat = "identity", 
                  width = 0.2) +
         
         
         annotate("richtext", label = "mTOR<sup>Ser2448</sup>", 
                  y = (2/100) * percent_from_bottom, x = 0.5, hjust = 0, 
                  size = text_size_label,
                  label.color = NA, 
                  fill = NA) +
         scale_y_continuous(limits = c(0,2), 
                            expand = c(0,0), 
                            breaks = c(0, 1, 2)) +
         
         scale_fill_manual(values = c(group.study.color[3], group.study.color[4])) +
         
         theme_bw() +
         theme(panel.grid = element_blank(), 
               panel.border = element_blank(), 
               axis.line.x  = element_line(size = line.size), 
               axis.line.y = element_line(size = line.size), 
               axis.ticks.x = element_line(size = line.size), 
               axis.ticks.y = element_line(size = line.size), 
               axis.text.y =  element_markdown(color = "black", size = 7), 
               
               # Individual plots settings
               axis.title.y = element_blank(),
               axis.text.x =  element_blank(), 
               axis.title.x = element_blank(),
               legend.title = element_blank(), 
               legend.background = element_rect(fill = "white"),
               legend.margin = margin(t = 0, r = 1, b = 1, l = 1, unit = "pt"),
               legend.key = element_rect(fill = "white"),
               legend.position = "none")
 

p.s6.fcfig <-  sets.fc %>%
         filter(target == "p.s6.tp") %>%
        
        mutate(sets = factor(sets,  levels = c("single", "multiple"), 
                             labels = c("Low-volume", "Moderate-volume"))) %>%
        
         ggplot(aes(sets, fold.change, fill = sets)) + 
         geom_bar(stat = "identity", 
                  width = 0.2) +
        
        
        annotate("richtext", label = "rpS6<sup>Ser235/236</sup>", 
                 y = (6/100) * percent_from_bottom, x = 0.5, hjust = 0, 
                 size = text_size_label, 
                 label.color = NA, 
                 fill = NA) +
         
         scale_y_continuous(limits = c(0,6), 
                            expand = c(0,0), 
                            breaks = c(0, 3, 6)) +
         
         scale_fill_manual(values = c(group.study.color[3], group.study.color[4])) +
         
        labs(x = "") +
        
        
         theme_bw() +
         theme(panel.grid = element_blank(), 
               panel.border = element_blank(), 
               axis.line.x  = element_line(size = line.size), 
               axis.line.y = element_line(size = line.size), 
               axis.ticks.x = element_line(size = line.size), 
               axis.ticks.y = element_line(size = line.size), 
               axis.text.y =  element_markdown(color = "black", size = 7), 
               
               # Individual plots settings
               axis.title.y = element_blank(),
               axis.text.x =  element_text(color = "black", size = 7, angle = 45, 
                                           hjust = 1), 
               axis.title.x = element_text(color = "black", size = 7),
               legend.title = element_blank(), 
               legend.background = element_rect(fill = "white"),
               legend.margin = margin(t = 0, r = 1, b = 1, l = 1, unit = "pt"),
               legend.key = element_rect(fill = "white"),
               legend.position = "none")
 


y_axis_label <- data.frame(x = 1, y = 1) %>%
        ggplot(aes(x, y)) + 
        theme_void() + 
        annotate("text", x = 1, y = 0.5,  label = "Fold-change from Pre-Ex", 
                 angle = 90, size = 2.4)



fold_change_plot <- plot_grid(y_axis_label, 
        
        plot_grid( 
          p.p85.fcfig, 
          p.p70.fcfig,
          p.mtor.fcfig,
          p.s6.fcfig, 
          nrow = 4, 
          rel_heights = c(1, 1, 1, 1.7)), 
        rel_widths = c(0.2, 1), 
        nrow = 1)









fold_changes_diffs <- rbind( data.frame(exp(confint(p.p70))) %>%
                                     mutate(coef = rownames(.), 
                                            target = "p.p70", 
                                            est = exp(coef(summary(p.p70))[3, 1]), 
                                            adjust = "baseline") %>%
                                     filter(coef == "setsmultiple"), 
                             
                             data.frame(exp(confint(p.p85))) %>%
                                     mutate(coef = rownames(.), 
                                            target = "p.p85", 
                                            est = exp(coef(summary(p.p85))[3, 1]), 
                                            adjust = "baseline") %>%
                                     filter(coef == "setsmultiple"), 
                             
                             data.frame(exp(confint(p.mTOR))) %>%
                                     mutate(coef = rownames(.), 
                                            target = "p.mTOR", 
                                            est = exp(coef(summary(p.mTOR))[3, 1]), 
                                            adjust = "baseline") %>%
                                     filter(coef == "setsmultiple"), 
                             
                             data.frame(exp(confint(p.s6.tp))) %>%
                                     mutate(coef = rownames(.), 
                                            target = "p.s6.tp", 
                                            est = exp(coef(summary(p.s6.tp))[3, 1]), 
                                            adjust = "baseline") %>%
                                     filter(coef == "setsmultiple"), 
                       ##### Unadjusted data ##########################################  
                             data.frame(exp(confint(p.p70.unadjusted))) %>%
                                     mutate(coef = rownames(.), 
                                            target = "p.p70", 
                                            est = exp(coef(summary(p.p70.unadjusted))[2, 1]), 
                                            adjust = "non") %>%
                                     filter(coef == "setsmultiple"), 
                             
                             data.frame(exp(confint(p.p85.unadjusted))) %>%
                                     mutate(coef = rownames(.), 
                                            target = "p.p85", 
                                            est = exp(coef(summary(p.p85.unadjusted))[2, 1]), 
                                            adjust = "non") %>%
                                     filter(coef == "setsmultiple"), 
                             
                             data.frame(exp(confint(p.mTOR.unadjusted))) %>%
                                     mutate(coef = rownames(.), 
                                            target = "p.mTOR", 
                                            est = exp(coef(summary(p.mTOR.unadjusted))[2, 1]), 
                                            adjust = "non") %>%
                                     filter(coef == "setsmultiple"), 
                             
                             data.frame(exp(confint(p.s6.tp.unadjusted))) %>%
                                     mutate(coef = rownames(.), 
                                            target = "p.s6.tp", 
                                            est = exp(coef(summary(p.s6.tp.unadjusted))[2, 1]), 
                                            adjust = "non") %>%
                                     filter(coef == "setsmultiple") ) %>%
        print()


fold_change_fig <- fold_changes_diffs %>%
        filter(adjust == "baseline") %>%
        mutate(target = factor(target, levels = c("p.p85", 
                                                  "p.p70", 
                                                  "p.mTOR", 
                                                  "p.s6.tp"), 
                               labels = c("p85-S6K1<sup>Thr412</sup>", 
                                          "p70-S6K1<sup>Thr389</sup>", 
                                          "mTOR<sup>Ser2448</sup>", 
                                          "rpS6<sup>Ser235/236</sup>")), 
               target = fct_rev(target)) %>%
        
        ggplot(aes(target, est)) + 
        geom_errorbar(aes(ymin = X2.5.., ymax = X97.5..), 
                       width = 0, 
                       position = position_dodge(width = 0.2 )) +
        geom_point(shape = 21, 
                   fill = group.study.color[5], 
                   size = 2.5, 
                   position = position_dodge(width = 0.2 )) +
        
        scale_y_continuous(limits = c(0.8, 2.4), 
                           expand = c(0,0), 
                           breaks = c(0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4), 
                           labels = c("",  1,  "", 1.4,  "", 1.8, "", 2.2, "")) +
       
         geom_hline(yintercept = 1, lty = 2, color = "grey85") +
        

        

        labs(y = expression(paste(Delta, "Moderate / ", Delta,  "Low"))) + 
        theme_bw() +
        theme(panel.grid = element_blank(), 
              panel.border = element_blank(), 
              axis.line.x  = element_line(size = line.size), 
              axis.line.y = element_blank(), 
              axis.ticks.x = element_line(size = line.size), 
              axis.ticks.y = element_blank(), 
              axis.text.y =  element_blank(),
              axis.title.y = element_blank(),
              axis.text.x =  element_text(color = "black", size = 7), 
              axis.title.x = element_text(color = "black", size = 7), 
              legend.title = element_blank(), 
              legend.background = element_rect(fill = "white"),
              legend.margin = margin(t = 0, r = 1, b = 1, l = 1, unit = "pt"),
              legend.key = element_rect(fill = "white"),
              legend.position = "none") + 
        coord_flip()





fc_plot <- plot_grid(fold_change_plot, 
          plot_grid(fold_change_fig, NULL, nrow = 2, rel_heights = c(1, 0.07)),  
          nrow = 1, rel_widths = c(0.7, 1))
          







west.img <- ggdraw() + draw_image("./data/study-1/westernImages/signalling_context.png")


western_fig <- plot_grid(west.img, fc_plot, 
                         labels = "auto", 
                         rel_widths = c(0.6, 1))



saveRDS(western_fig, "./figures/results/study1-western-blot.RDS")



####################################### Correlation with RT induced growth #########

# Using S6K1 p70

p.85 <- wb.raw.acute %>%
  filter(target == "p.p85") %>%
  group_by(subject, sets, timepoint) %>%
  summarise(expression = mean(expression, na.rm = TRUE)) %>%
  pivot_wider(names_from = timepoint, 
              values_from = expression) %>%
  mutate(fc.p85 = w2post / w2pre, 
         w2post.p85 = w2post, 
         w2pre.p85 = w2pre) %>%
  dplyr::select(subject:sets, fc.p85:w2pre.p85) %>%
  print()



p70 <- west.acute %>%
  filter(target == "p.p70") %>%
  group_by(subject, sets, timepoint) %>%
  summarise(expression = mean(expression)) %>%
  pivot_wider(names_from = timepoint, 
              values_from = expression) %>%
  mutate(fc = w2post / w2pre) %>%
  print()

# getting mr data 

# Total RNA analysis #
# loads include and condition data
include <- read.csv("./data/study-1/oneThreeSetLeg.csv", sep=";") %>%
  gather(sets,leg, multiple:single) %>%
  mutate(leg) %>%
  dplyr::select(subject, sex, include, leg, sets) %>%
  filter(include == "incl")


# mr data 
source("./R/study-1/mr-and-dxa-data.R")

p70_csa <-  mr.results %>%
  inner_join(include) %>%
  dplyr::select(subject,  sets, timepoint,  CSA.avg) %>%
  pivot_wider(names_from = timepoint, values_from = CSA.avg) %>%
  mutate(change = 100 * ((post/pre-1))) %>%
  mutate(sets = factor(sets, levels = c("single", "multiple"))) %>%

  inner_join(p70) %>%
  inner_join(p.85) %>%
  print()
  

saveRDS(p70_csa, "./data/derivedData/study1-western-blot/csa_p70.RDS")








