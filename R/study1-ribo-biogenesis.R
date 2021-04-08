##-------------------------------------
## study1-ribo-biogenesis.R
##
## Title: Ribosome biogenesis markers from study 1
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



# Total RNA and ribo-analysis #######

source("./R/libraries.R")
source("./R/themes.R")




# Total RNA analysis #
# loads include and condition data
include <- read.csv("./data/study-1/oneThreeSetLeg.csv", sep=";") %>%
        gather(sets,leg, multiple:single) %>%
        mutate(condition=leg) %>%
        dplyr::select(subject, sex, include, condition, sets)

# loads sample weight and total RNA concentration in trizol extracts 
total.rna <- read_excel("./data/study-1/sample_weight.xlsx") %>%
        dplyr::select(subject, timepoint, leg, weight, conc1.5, elution.volume)%>% # total RNA measured in 1:5 dilution
        mutate(condition=leg)%>%
        inner_join(include)%>%
        filter(include=="incl")%>%
        mutate(rna = (conc1.5*5*elution.volume)/weight)%>% # estimates of total RNA
        dplyr::select(subject, sex, sets, timepoint, rna, weight, conc1.5, elution.volume)%>%
        mutate(timepoint = factor(timepoint, levels = c("w0", "w2pre", "w2post", "w12")),
               conc = conc1.5*5*elution.volume)%>% # scale RNA concentration based on dilution
        filter(timepoint !="w2post")


# Some samples were underestimated because of sample loss during 
# sample prep. This was due to loss of RNA in washing steps were piece of 
# RNA pellet was occasionally lost. To estimate deviance from weight to RNA relationship
# a linear model was constructed log(concentration)~log(weight) and
# deviance of more than 4 units in standardized residuals from zero meant exclusion
# from further modeling/plotting

# model for RNA concentration to weight relationship
l.mod <- lm(log(conc) ~ log(weight), data = total.rna) 

# Calculates standardized residuals
total.rna$resid<- resid(l.mod)/sd(resid(l.mod))
# store predicted values for diagnostic plotting
total.rna$pred<- predict(l.mod)

# plot predicted vs. standardized residuals

# creates a outlier variable
total.rna<-total.rna%>%
        mutate(outlier = if_else(resid < -4 | resid > 4 , "out", "in")) # 4 units deviance from zero

tot.rna<-total.rna%>%
        filter(outlier != "out")%>% 
        group_by(subject, timepoint, sets)%>%
        summarise(rna = mean(rna, na.rm=T))%>%
        ungroup()

## Data for statistics ##
stats.total.rna <- tot.rna %>%
        mutate(log.rna = log(rna),
               leg = paste0(subject, sets), 
               sets = factor(sets, levels = c("single", "multiple")))

# Model using nlme

rna.m0 <- lme(log(rna) ~ timepoint + timepoint:sets,
              random = list(subject = ~1), 
              control=list(msMaxIter=120,
                           opt = "nloptwrap"), method = "REML",
              stats.total.rna)

rna.m1 <- lme(log(rna) ~ timepoint + timepoint:sets,
              random = list(subject = ~1, leg = ~1), 
              control=list(msMaxIter=120,
                           opt = "nloptwrap"), method = "REML",
              stats.total.rna)

rna.m2 <- lme(log(rna) ~ timepoint + timepoint:sets,
              random = list(subject = ~1,
                            leg = pdDiag(~1 + timepoint)), 
              control=list(msMaxIter=120,
                           opt = "nloptwrap",
                           msVerbose=FALSE), method = "REML",
              stats.total.rna)

# anova(rna.m0, rna.m1, rna.m2) 
# Evidence for including a random slope for time but not random intercept for sets.





### Diagnostic plots not run ###
# qqnorm(resid(rna.m0), main = "total.rna.model log-values"); qqline(resid(rna.m0))
# plot(rna.m0, resid(.)~fitted(.)|timepoint)


# Patterns in the residual plot, the simplest model is the best model. 
# m0 is prefered.


############ Descriptive statistics ###################

# total RNA % change




# Gene estimates 
qpcr.means.estimates <- read_feather("./data/study-1/derivedData/qpcr_mean_estimates_rrna")

descriptive_ribo_markers <- rbind(qpcr.means.estimates %>%
        filter(gene %in% c("rRNA18s", 
                           "rRNA28s", 
                           "rRNA58s", 
                           "rRNA45s")) %>%
        dplyr::select(sets:emmean) %>%
        pivot_wider(names_from = timepoint, 
                    values_from = emmean) %>%
        mutate(w2pre = exp(w2pre - w0), 
               w12 = exp(w12 - w0)) %>%
        group_by(sets, gene) %>%
        summarise(w2pre = (mean(w2pre, na.rm = TRUE) - 1) * 100, 
                  w12 = (mean(w12, na.rm = TRUE) - 1) * 100) %>%
        dplyr::select(sets, w2pre, w12, target = gene), 
      
      stats.total.rna %>%
              dplyr::select(subject, timepoint, sets, rna) %>%
              pivot_wider(names_from = timepoint, 
                          values_from = rna) %>%
              
              mutate(w2pre = w2pre/w0, 
                     w12 = w12/w0) %>%
              
              group_by(sets) %>%
              summarise(w2pre = (mean(w2pre, na.rm = TRUE) - 1) * 100, 
                        w12 = (mean(w12, na.rm = TRUE) - 1) * 100)  %>%
              mutate(target = "totalrna"))


saveRDS(descriptive_ribo_markers, 
        "./data/derivedData/study1-ribo-biogenesis/descriptive_ribo_markers.RDS")





### Summary plot of model estimates

qpcr.model.estimates <- read_feather("./data/study-1/derivedData/qpcr_model_estimates_rrna")




training_ribo_markers  <- rbind(data.frame(intervals(rna.m0)$fixed) %>%
        mutate(timepoint = rownames(.), 
               target = "totalrna") %>%
        filter(timepoint %in% c("timepointw0:setsmultiple", 
                                "timepointw2pre:setsmultiple", 
                                "timepointw12:setsmultiple")) %>%
        dplyr::select(target, 
                      timepoint, 
                      estimate = est., 
                      lower, 
                      upper), 
      
      qpcr.model.estimates %>%
        filter(parameter %in% c("timepointw0:setsmultiple", 
                                "timepointw2pre:setsmultiple", 
                                "timepointw12:setsmultiple"), 
               gene != "generRNA45S") %>%
        dplyr::select(target = gene, 
                      timepoint = parameter,
                      estimate, 
                      lower, 
                      upper)) %>%
        
        mutate(target = gsub("gene", "", target), 
               timepoint = gsub(":setsmultiple", "", timepoint), 
               timepoint = gsub("timepoint", "", timepoint), 
               
               target = factor(target, 
                               levels = c("rRNA28s",
                                          "rRNA58s",
                                          
                                          "rRNA18s",
                                          "rRNA45s",  
                                          "totalrna"), 
                               labels = c("rRNA 28S",
                                          "rRNA 5.8S",
                                          "rRNA 18S",
                                          "pre-rRNA 45S", 
                                          "Total RNA")), 
               timepoint = factor(timepoint, 
                                  levels = c("w12", 
                                             "w2pre", 
                                             "w0"), 
                                  labels = c("Week 12", 
                                             "Week 2", 
                                             "Week 0")), 
               sig = if_else(estimate > 0 & lower > 0 | estimate < 0 & upper < 0, "sig", "ns")) %>%
   
        
        ggplot(aes(target, exp(estimate), fill = timepoint, alpha = sig)) + 
        
        geom_hline(yintercept = 1, lty = 2, color = "grey85") +
        
        geom_errorbar(aes(ymin = exp(lower), ymax = exp(upper)), 
                      position = position_dodge(width = 0.5), 
                      width = 0) +
        
        geom_point(position = position_dodge(width = 0.5), 
                   shape = 21, size = 2.5) + 
        
        scale_alpha_manual(values = c(0.4, 1)) +
        
        scale_fill_manual(values = c(group.study.color[1], 
                                     group.study.color[5], 
                                     group.study.color[6])) +
       
    
        
         scale_y_continuous(limits = c(0.5, 1.5), 
                            expand = c(0, 0),
                           breaks = c(0.5, 0.75, 1, 1.25, 1.5), 
                           labels = c(0.5, "", 1, "", 1.5)) +

        labs(y = "Moderate-volume / Low-volume") + 
        theme_bw() +
        theme(panel.grid = element_blank(), 
              panel.border = element_blank(), 
              axis.line.x  = element_line(size = line.size), 
              axis.line.y = element_blank(), 
              axis.ticks.x = element_line(size = line.size), 
              axis.ticks.y = element_blank(), 
              axis.text.y =  element_text(color = "black", size = 7), 
              axis.title.y = element_blank(),
              axis.text.x =  element_text(color = "black", size = 7), 
              axis.title.x = element_text(color = "black", size = 7), 
              legend.title = element_blank(), 
              legend.background = element_rect(fill = "white"),
              legend.margin = margin(t = 0, r = 1, b = 1, l = 1, unit = "pt"),
              legend.key = element_rect(fill = "white"),
              legend.position = "none") + 
        coord_flip() +
        # Annotate time points
        annotate("text", 
                 x = c(5.2, 
                       5.05,
                       4.7), 
                       
                 y =  c(0.7, 
                        1.35, 
                        1.32),
                 label = c("Week 0",
                         "Week 2",
                         "Week 12"), 
                 color = "gray50",        
                 size = 2.2) +
        
        # Legend arrow week 2
        annotate("curve",
                 x = 5.1,
                 xend = 5.08,
                 y = 1.25,
                 yend = 1.17,
                 curvature = 0.2,
                 color = "gray50",
                 arrow = arrow(length=unit(0.1,"cm"), type = "closed")) +
        # legend arrow week 0
        annotate("curve",
                 x = 5.1,
                 xend =  5.12,
                 y =  0.74,
                 yend =  0.87,
                 curvature = 0.2,
                 color = "gray50",
                 arrow = arrow(length=unit(0.1,"cm"), type = "closed")) +
        # legend arraw week 12
        annotate("curve",
                 x =  4.6,
                 xend =  4.7,
                 y = 1.2,
                 yend = 1.07,
                 curvature = -0.2,
                 color = "gray50",
                 arrow = arrow(length=unit(0.1,"cm"), type = "closed"))




training_ribo_markers






### Acute fold changes c-myc ###


qpcr.model.estimates.acute <- read_feather("./data/study-1/derivedData/qpcr_model_estimates.acute")
qpcr.means.estimates.acute <- read_feather("./data/study-1/derivedData/qpcr_mean_estimates.acute")



acute_qpcr <- qpcr.means.estimates.acute %>%
        filter(gene %in% c("Myc", "rRNA45S")) %>%
        dplyr::select(sets, timepoint, gene, emmean) %>%
        pivot_wider(names_from = timepoint, 
                    values_from = emmean) %>%
        mutate(w2post = exp(w2post - w2pre)) %>%
        print()



myc_acute <- acute_qpcr %>%
        filter(gene == "Myc") %>%

        ggplot(aes(sets, w2post, fill = sets)) + 
        geom_bar(stat = "identity", 
                 width = 0.2) +
        
        annotate("richtext", label = "c-Myc", 
                 y = (10/100) * 90, x = 0.5, hjust = 0, 
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

rrna45_acute <- acute_qpcr %>%
        filter(gene == "rRNA45S") %>%
        mutate(sets = factor(sets, levels = c("single", "multiple"), 
                             labels = c("Low-\nvolume", 
                                        "Moderate-\nvolume"))) %>%
        
        ggplot(aes(sets, w2post, fill = sets)) + 
        geom_bar(stat = "identity", 
                 width = 0.2) +
        
        annotate("richtext", label = "rRNA 45S", 
                 y = (2/100) * 90, x = 0.5, hjust = 0, 
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
              axis.text.x =  element_text(color = "black", size = 7, angle = 45, 
                                          hjust = 1, lineheight = 0.6), 
              axis.title.x = element_blank(),
              legend.title = element_blank(), 
              legend.background = element_rect(fill = "white"),
              legend.margin = margin(t = 0, r = 1, b = 1, l = 1, unit = "pt"),
              legend.key = element_rect(fill = "white"),
              legend.position = "none")







acute_estimates <- qpcr.model.estimates.acute %>%
        filter(gene %in% c("Myc","rRNA45S"), 
               parameter == "timepointw2post:setsmultiple") %>%
        
        mutate(gene = factor(gene, levels = c("rRNA45S", 
                                              "Myc"), 
                             labels = c( "rRNA45S", 
                                         "c-Myc"))) %>%
        
        
        ggplot(aes(exp(estimate), gene, alpha = gene)) +
        geom_vline(xintercept = 1, lty = 2, color = "gray80") +
        geom_errorbarh(aes(xmin = exp(lower), xmax = exp(upper)), height = 0) +
        geom_point(fill = group.study.color[2], size = 3, shape = 21) +
        labs(x = expression(paste(Delta, "Moderate / ", Delta, "Low"))) +
        scale_x_continuous(limits = c(0.5, 2.5), 
                           breaks = c(0.5, 1, 1.5, 2, 2.5), 
                           labels = c(0.5, 1, 1.5, 2, 2.5), 
                           expand = c(0, 0)) +
        
        scale_alpha_manual(values = c(0.4, 1)) +
        
 
        theme_bw() +
        theme(panel.grid = element_blank(), 
              panel.border = element_blank(), 
              axis.line.x  = element_line(size = line.size), 
              axis.line.y = element_blank(), 
              axis.ticks.x = element_line(size = line.size), 
              axis.ticks.y = element_blank(), 
              axis.text.y =  element_blank(),  
              
              # Individual plots settings
              axis.title.y = element_blank(),
              axis.text.x =  element_text(color = "black", size = 7), 
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


ribo_markers_fig <- plot_grid(training_ribo_markers, 
        plot_grid(plot_grid(y_axis_label, 
          
                plot_grid(myc_acute, 
                        rrna45_acute, nrow = 2, 
                        rel_heights = c(0.8, 1)), 
                
                 nrow = 1, rel_widths = c(0.1, 1)), 
           
           acute_estimates, 
          rel_widths = c(1, 1)), 
        rel_widths = c(1, 0.8), 
        nrow = 1, 
        labels = "auto")








saveRDS(ribo_markers_fig, "./figures/results/study1-ribo-markers.RDS")




















