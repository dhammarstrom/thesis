##-------------------------------------
## study-1-determinants.R
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



### Figure 1. Strength and muscle mass differences in 


include <- read.csv("./data/study-1/oneThreeSetLeg.csv", sep=";") %>%
        gather(sets,leg, multiple:single) %>%
        mutate(condition=leg) %>%
        dplyr::select(subject, sex, include, condition, sets)




# mr data 
source("./R/study-1/mr-and-dxa-data.R")





### Quantify avergae percentage error ###
csa.error <- mr.results %>%
        mutate(condition=leg) %>%
        inner_join(read.csv("./data/study-1/oneThreeSetLeg.csv", sep=";") %>%
                           gather(sets,leg, multiple:single) %>%
                           mutate(condition = leg) %>%
                           dplyr::select(subject, sex, include, condition, sets)) %>%
        filter(include=="incl") %>%
        dplyr::select(subject, timepoint, sex, sets, CSA.avg) %>%
        filter(timepoint == "pre") %>%
        group_by(subject, sex) %>%
        summarise(m = mean(CSA.avg)) %>%
        group_by(sex) %>%
        mutate(m.centered = m - mean(m, na.rm = TRUE)) %>%
        group_by() %>%
        mutate(swc = sd(m.centered) * 0.2) %>%
        group_by(sex) %>%
        summarise(swc = 100 * (mean(swc)/mean(m)), 
                  n = n()) %>%
        group_by() %>%
        summarise(swc = sum(swc * n) / sum(n)) %>%
        print()





benefit <- mr.results %>%
        mutate(condition=leg) %>%
        inner_join(read.csv("./data/study-1/oneThreeSetLeg.csv", sep=";") %>%
                           gather(sets,leg, multiple:single) %>%
                           mutate(condition = leg) %>%
                           dplyr::select(subject, sex, include, condition, sets)) %>%
        filter(include=="incl") %>%
        dplyr::select(subject, timepoint, sex, sets, CSA.avg) %>%
        spread(timepoint, CSA.avg) %>%
        mutate(error = as.numeric(csa.error), 
               change = ((post/pre)-1) * 100) %>%
        inner_join(read_excel("./data/study-1/body-composition/bodycomp_DXA.xlsx") %>%
                           inner_join(read.csv("./data/study-1/oneThreeSetLeg.csv", sep=";")) %>%
                           filter(include =="incl") %>% 
                           group_by(subject) %>%
                           summarise(weight = mean(weight),
                                     height = mean(height))) %>% 
        group_by(subject) %>%
        mutate(mean.pre = mean(pre)) %>%
        ungroup() %>%
        dplyr::select(subject, sex, sets, change, mean.pre, error) %>%
        spread(sets, change) %>%
        mutate(diff = multiple - single, 
               benefit = if_else(multiple - single > error, 1, 0),
               benefit.single = if_else(single - multiple > error, 1, 0)) %>%
        dplyr::select(subject, sex, benefit, benefit.single, diff) 





temp.stats <- mr.results%>%
        mutate(condition=leg)%>%
        inner_join(include)%>%
        filter(include=="incl")%>%
        dplyr::select(subject, timepoint, sex, sets, CSA.avg)%>%
        spread(timepoint, CSA.avg)%>%
        mutate(csa.change=((post / pre)-1)*100)%>%
        dplyr::select(subject, csa.change, sets)%>%
        spread(sets, csa.change)







buildPoly <- function(limits.x = c(-3, 14), limits.y = c(-3, 14), slope = 1, intercept = 0, above = TRUE){
        #Assumes ggplot expand = c(0.05,0)
        xrTru <- limits.x
        yrTru <- limits.y
        
        #Find where the line crosses the plot edges
        yCross <- (yrTru - intercept) / slope
        xCross <- (slope * xrTru) + intercept
        
        #Build polygon by cases
        if (above & (slope >= 0)){
                rs <- data.frame(x=-Inf,y=Inf)
                if (xCross[1] < yrTru[1]){
                        rs <- rbind(rs,c(-Inf,-Inf),c(yCross[1],-Inf))
                }
                else{
                        rs <- rbind(rs,c(-Inf,xCross[1]))
                }
                if (xCross[2] < yrTru[2]){
                        rs <- rbind(rs,c(Inf,xCross[2]),c(Inf,Inf))
                }
                else{
                        rs <- rbind(rs,c(yCross[2],Inf))
                }
        }
        if (!above & (slope >= 0)){
                rs <- data.frame(x= Inf,y= -Inf)
                if (xCross[1] > yrTru[1]){
                        rs <- rbind(rs,c(-Inf,-Inf),c(-Inf,xCross[1]))
                }
                else{
                        rs <- rbind(rs,c(yCross[1],-Inf))
                }
                if (xCross[2] > yrTru[2]){
                        rs <- rbind(rs,c(yCross[2],Inf),c(Inf,Inf))
                }
                else{
                        rs <- rbind(rs,c(Inf,xCross[2]))
                }
        }
        if (above & (slope < 0)){
                rs <- data.frame(x=Inf,y=Inf)
                if (xCross[1] < yrTru[2]){
                        rs <- rbind(rs,c(-Inf,Inf),c(-Inf,xCross[1]))
                }
                else{
                        rs <- rbind(rs,c(yCross[2],Inf))
                }
                if (xCross[2] < yrTru[1]){
                        rs <- rbind(rs,c(yCross[1],-Inf),c(Inf,-Inf))
                }
                else{
                        rs <- rbind(rs,c(Inf,xCross[2]))
                }
        }
        if (!above & (slope < 0)){
                rs <- data.frame(x= -Inf,y= -Inf)
                if (xCross[1] > yrTru[2]){
                        rs <- rbind(rs,c(-Inf,Inf),c(yCross[2],Inf))
                }
                else{
                        rs <- rbind(rs,c(-Inf,xCross[1]))
                }
                if (xCross[2] > yrTru[1]){
                        rs <- rbind(rs,c(Inf,xCross[2]),c(Inf,-Inf))
                }
                else{
                        rs <- rbind(rs,c(yCross[1],-Inf))
                }
        }
        
        return(rs)
}


mr.cor.data <- mr.results %>%
        mutate(condition=leg) %>%
        inner_join(include) %>%
        inner_join(benefit) %>%
        filter(include=="incl") %>%
        dplyr::select(subject, timepoint, sex, sets, CSA.avg, benefit) %>%
        spread(timepoint, CSA.avg) %>%
        mutate(csa.change=((post / pre)-1) * 100) %>%
        dplyr::select(subject, sex, sets,benefit, csa.change) %>%
        spread(sets, csa.change)%>%
        mutate(sex = factor(sex, levels=c("female", "male"), 
                            labels=c("Female", "Male")))

triangle <- buildPoly(limits.x = c(-3, 14), limits.y = c(-3, 14),
                      slope=1,intercept=2.41,above=FALSE)

### CSA correlation figure #### 

mr.cor.fig <- mr.cor.data %>%
        ggplot(aes(single, multiple, fill = as.factor(benefit)))+
        geom_abline(intercept = 2.41, slope = 1, color = "gray40") +
        geom_abline(intercept = -2.41, slope = 1, color = "gray40") +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
        geom_point(size = 2.5, shape = 21, stroke = 0.2)+
        geom_polygon(data = triangle, aes(x = x,
                                          y = y, 
                                          shape = NULL, 
                                          color = NULL), 
                     fill = "gray80", alpha = 0.2) +
        scale_x_continuous(limits=c(-3, 14), 
                           expand = c(0,0), 
                           breaks = c(-2, 0, 2, 4, 6, 8, 10, 12, 14), 
                           labels = c("", 0, "", 4, "", 8, "", 12, ""))+
        scale_y_continuous(limits=c(-3, 14), 
                           expand = c(0,0), 
                           breaks = c(-2, 0, 2, 4, 6, 8, 10, 12, 14), 
                           labels = c("", 0, "", 4, "", 8, "", 12, ""))+

        dissertation_theme() +
        
        scale_fill_manual(values = c(group.study.color[1], group.study.color[5])) +
        xlab("Single-set (%-change)")  +
        ylab("Multiple-set (%-change)") +
        
        
        theme(legend.position = "none",
              legend.title = element_blank(),
              legend.text = element_text(size = 8), 
              axis.text = element_text(size = 8),
              axis.title = element_text(size = 8),
              panel.border = element_blank()) +
        annotate("text", -Inf, Inf, label = "CSA", hjust = -0.1, vjust = 1.5, size = 2.8)


#### Strength fig 

rm_str <- read_excel("./data/study-1/strength/strengthTests.xlsx")%>%
        filter(!(exercise %in% c("benchpress", "legcurl")))%>%
        mutate(condition = leg) %>%
        inner_join(include)%>%
        filter(include=="incl")%>%
        dplyr::select(subject, sex, sets, exercise, timepoint, load)%>%
        mutate(load = as.numeric(load),
               sets = factor(sets, levels = c("single", "multiple"))) %>%
        spread(timepoint, load) %>%
        mutate(baseline = pmax(pre, session1)) %>%
        dplyr::select(subject, sex, sets, exercise, baseline, week2, week5, week9, post) %>%
        gather(week, load, baseline:post) %>%
        mutate(week = factor(week, levels = c("baseline", "week2", "week5", "week9", "week12", "post"))) %>%
        group_by(exercise) %>%
        mutate(str = load / max(load, na.rm = TRUE)) %>%
        dplyr::select(-load) %>%
        filter(week %in% c("baseline", "post")) %>%
        print()



# load cybex data 
load("./data/study-1/strength/cybex.rda")


cybex <- data.frame(cybex)

# clean data sets
cybex$condition<-NA
cybex[cybex$leg=="Left",]$condition<-"L"
cybex[cybex$leg=="Right",]$condition<-"R"

cybex <- cybex %>%
        inner_join(include) %>%
        filter(include=="incl") %>%
        dplyr::select(subject, timepoint, sets, sex, type, speed, torque.max, position, actual.speed)%>%
        filter(!(subject=="FP6" & timepoint == "post" & speed == 240)) # This removes one record, exclusion due to error in measurement, see protocol FP6


# Clean cybex data set and combine with 1RM 
total_strength <- cybex %>%
        ungroup() %>%
        dplyr::select(subject, timepoint, sets,speed, sex, torque.max) %>%
        spread(timepoint, torque.max) %>%
        rowwise()%>%
        mutate(baseline = max(pre, session1, na.rm=T)) %>%
        ungroup() %>%
        dplyr::select(subject, sex, sets, exercise = speed, baseline, week2, week5, week9, post) %>%
        gather(week, load, baseline:post) %>%
        group_by(exercise) %>%
        mutate(str = load/max(load, na.rm = TRUE)) %>%
        dplyr::select(-load) %>%
        filter(week %in% c("baseline", "post")) %>%
        rbind(rm_str) %>% # This includes all strength tests
        mutate(ex_fam = if_else(exercise %in% c("legext", "legpress"), 
                                "rm", if_else(exercise == 0, "isom", "isok"))) %>%
        print()



### Smallest worthwile change calculation ###
strength.cv <- total_strength %>%
        filter(week == "baseline") %>%
        group_by(subject, sex, ex_fam) %>%
        summarise(m = mean(str, na.rm = TRUE)) %>%
        group_by(sex, ex_fam) %>%
        mutate(m.centered = m - mean(m, na.rm = TRUE)) %>%
        group_by(ex_fam) %>%
        mutate(swc = (sd(m.centered) * 0.2)) %>%
        group_by(sex, ex_fam) %>%
        summarise(swc = 100 * (mean(swc)/mean(m)), 
                  n = n()) %>%
        group_by() %>%
        summarise(swc = sum(swc * n) / sum(n) ) %>%
        print()


strength_data <- total_strength %>%
        group_by(subject, sex, sets, week, ex_fam) %>%
        summarise(str = mean(str, na.rm = TRUE)) %>%
        group_by(subject, sex, sets, week) %>%
        summarise(str = mean(str, na.rm = TRUE)) %>%
        spread(week, str) %>%
        mutate(change = ((post/baseline)-1)*100) %>%
        
        dplyr::select(subject, sex, sets, change) %>%
        spread(sets, change) %>%
        mutate(cv = data.frame(strength.cv)[1,1]) %>%
        mutate(benefit.strength = if_else(multiple - single > cv, 1, 0)) %>%
        print()







triangle.rm <- buildPoly(limits.x = c(-5, 70), limits.y = c(-5, 70),
                         slope = 1, intercept = as.numeric(strength.cv[1,1]), above = FALSE) 

strength_data %>%
        group_by(benefit.strength) %>%
        summarise(n = n())

#### Strength correlations figure ####

strength.benefit <- strength_data %>%
        ungroup() %>%
        
        ggplot(aes(single, multiple, fill = as.factor(benefit.strength))) + 
        geom_abline(mapping = aes(intercept = cv, slope = 1),  color = "gray40") +
        geom_abline(mapping = aes(intercept = -cv,  slope = 1), color = "gray40") +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
        geom_point(size = 2.5, shape = 21, stroke = 0.2) +
        geom_polygon(data = triangle.rm, aes(x = x,
                                             y = y, 
                                             shape = NULL, 
                                             color = NULL), 
                     fill = "gray80", alpha = 0.2) +
        scale_x_continuous(limits=c(-5, 70), 
                           expand = c(0,0), 
                           breaks = c(0, 10, 20, 30, 40, 50, 60, 70),
                           labels = c(0, "", 20, "", 40, "", 60, ""))+
        scale_y_continuous(limits=c(-5, 70), 
                           expand = c(0,0), 
                           breaks = c(0, 10, 20, 30, 40, 50, 60, 70),
                           labels = c(0, "", 20, "", 40, "", 60, ""))+
        
        dissertation_theme() +
        
        scale_fill_manual(values = c(group.study.color[1], group.study.color[5])) +
        
        xlab(expression(paste("Single-set (%-change)")))  +
        ylab(expression(paste("Multiple-set (%-change)"))) +
        theme(legend.position = "none",
              legend.title = element_blank(),
              axis.text = element_text(size = 8),
              axis.title = element_text(size = 8),
              strip.background = element_blank(),
              strip.text = element_text(size = 8),
              panel.border = element_blank()) +
        annotate("text", -Inf, Inf, label = "Strength", hjust = -0.1, vjust = 1.5, size = 2.8)



#### Combine and save figure 


csa_str_cor <- plot_grid(mr.cor.fig, strength.benefit, ncol = 2, rel_widths = c(.5, .5))


saveRDS(csa_str_cor, "./data/derivedData/study-1-determinants/csa_str_cor.RDS")


### Save correlations between sets in csa and strength



corr <- data.frame(variable = c("strength", "csa"), 
           estimate = c(cor.test(strength_data$single, strength_data$multiple)$estimate, 
                        cor.test(temp.stats$single, temp.stats$multiple)$estimate), 
           ci.lwr = c(cor.test(temp.stats$single, temp.stats$multiple)$conf.int[1], 
                      cor.test(temp.stats$single, temp.stats$multiple)$conf.int[1]),
           ci.upr = c(cor.test(temp.stats$single, temp.stats$multiple)$conf.int[2], 
                      cor.test(temp.stats$single, temp.stats$multiple)$conf.int[2]))

saveRDS(corr, "./data/derivedData/study-1-determinants/csa_str_cor_df.RDS")


############ Modeling #############################



# Get univariate coefficients

results_univariate <- readRDS("./data/derivedData/study-1-benefit-analysis/univariate_analysis.Rda")


# Maka a transformation for log-transformation of negative values


univariate_regression <- results_univariate %>%
        dplyr::select(variable,  
                      estimate.csa, 
                      csa.ci.lwr,   csa.ci.upr, csa.ci.lwr.80, csa.ci.upr.80, 
                      estimate.str,
                      str.ci.lwr,  str.ci.upr, str.ci.lwr.80, str.ci.upr.80) %>%

        unite("csa", estimate.csa:csa.ci.upr.80, sep = "_") %>%
        unite("str", estimate.str:str.ci.upr.80, sep = "_") %>%

        dplyr::select(variable, csa, str) %>%


        pivot_longer(names_to = "type", 
                     values_to = "estimate", 
                     cols = csa:str) %>%
    
        
        separate(estimate, into = c("est", "lwr", "upr", "lwr.80", "upr.80"), sep = "_", convert = TRUE) %>%

       # filter(variable != "cortisol.ba") %>%
        
        mutate(sig = if_else(est > 0 & lwr.80 > 0 | est < 0 & upr.80 < 0, 
                             "sig", "ns"), 
               variable = fct_reorder(variable, est, mean), 
               type = factor(type, levels = c("csa", "str"), 
                             labels = c("CSA", "Strength"))) %>%
 
        
     #   filter(variable != "cortisol.ba") %>%
        
        ggplot(aes(est, variable, alpha = sig)) + 
        
        
        geom_vline(xintercept = 0, lty = 2, color = "gray50") +
        
        ## 95% CI
        geom_errorbarh(aes(xmin = lwr, 
                           xmax = upr), 
                       height = 0) +
        
        ## 95% CI
        geom_errorbarh(aes(xmin = lwr.80, 
                           xmax = upr.80), 
                       height = 0, size = 1.5) +
        
        geom_point(size = 2.5, shape = 21, fill = group.study.color[5])   +
        
        dissertation_theme() +
        
        scale_alpha_manual(values = c(0.4, 1)) +
        
        facet_grid(. ~ type) +
        
        labs(y = "Variable", 
             x = "Benefit vs. No benefit of moderate-volume training (Difference SD-units)") +
        
        theme(legend.position = "none",
              legend.title = element_blank(),
              axis.text = element_text(size = 7),
              axis.title.x = element_text(size = 7),
              axis.title.y = element_blank(),
              strip.background = element_blank(),
              panel.border = element_blank()) 

        
        
         
     

saveRDS(univariate_regression, "./data/derivedData/study-1-benefit-analysis/univariate_regression_results_fig.Rds")




### Model selection process ##########



csa_models <- readRDS("./data/derivedData/study-1-benefit-analysis/csa_logit_models.Rds")
str_models <- readRDS("./data/derivedData/study-1-benefit-analysis/str_logit_models.Rds")


csa_model_reduction <- csa_models %>%
        filter(!(term %in% c("sexmale", "(Intercept)")), 
                 !(model %in% c("m5", "reduced"))) %>%
        
        mutate(ci = paste0("[", 
                           sprintf("%.2f", exp(X2.5..)), 
                           ", ", 
                           sprintf("%.2f", exp(X97.5..)),
                           "]"), 
               model = factor(model, 
                              levels = c("m1", "m2", "m3", "m4", "m5"), 
                              labels = c("Model 1", 
                                         "Model 2", 
                                         "Model 3", 
                                         "Model 4", 
                                         "Model 5"))) %>% 
        
        mutate(term = factor(term, 
                             levels = c("rna_w2pre",
                                        "lean.massb_low",
                                        "testo.ba",
                                        "gh.tr", 
                                        "type2xlo" ))) %>%

        ggplot(aes(exp(estimate), term, fill = term)) + 
        
       geom_vline(xintercept = 1, lty = 2, color = "gray50", 
                  size = 0.4) +
        
        geom_text(aes(x = exp(estimate), y = term, label = ci), 
                  size = 2.5, 
                  position = position_nudge(y = -0.22, x = 0.25)) +

        geom_point(shape = 21, size = 2) +
        
        scale_x_continuous(limits = c(0.1, 1.8), 
                           breaks = c(0.5, 1, 1.5), 
                           expand = c(0,0)) +
        
        facet_grid(. ~ model) + 
        dissertation_theme() +
        
        scale_fill_manual(values = group.study.color) +
        
         theme(legend.position = "none",
              legend.title = element_blank(),
              axis.text = element_text(size = 7),
              axis.title = element_blank(),
              
              axis.line.y = element_blank(), 
              axis.ticks.y = element_blank(),
              
              panel.background = element_rect(fill = "gray95"),
              strip.background = element_blank(),
              strip.text = element_text(size = 7, hjust = 0),
              panel.border = element_blank()) 

        
        
str_model_reduction <- str_models %>%
        filter(!(term %in% c("sexmale", "(Intercept)")), 
               !(model %in% c("reduced"))) %>%
        
        mutate(ci = paste0("[", 
                           sprintf("%.2f", exp(X2.5..)), 
                           ", ", 
                           sprintf("%.2f", exp(X97.5..)),
                           "]"), 
               model = factor(model, 
                              levels = c("m1", "m2", "m3"), 
                              labels = c("Model 1", 
                                         "Model 2", 
                                         "Model 3"))) %>% 
        
       mutate(term = factor(term, 
                            levels = c("rna_w2pre",
                                       "p.p70_w2post",
                                       "cortisol.ba"))) %>%
       
        ggplot(aes(exp(estimate), term, fill = term)) + 
        
        geom_vline(xintercept = 1, lty = 2, color = "gray50", 
                   size = 0.4) +
        
        geom_text(aes(x = exp(estimate), y = term, label = ci), 
                  size = 2.5, 
                  position = position_nudge(y = -0.22, x = 0.25)) +
        
        geom_point(shape = 21, size = 2) +
        
        scale_x_continuous(limits = c(0.1, 1.8), 
                           breaks = c(0.5, 1, 1.5), 
                           expand = c(0,0)) +
        
        facet_grid(. ~ model) + 
        dissertation_theme() +
        
        scale_fill_manual(values = group.study.color) +
        
        labs(y = "Variable", 
             x = "Odds-ratio, benefit of moderate-volume training") +
        
        theme(legend.position = "none",
              legend.title = element_blank(),
              axis.text = element_text(size = 7),
              axis.title.x = element_text(size = 7, hjust = 0.2),
              
              axis.line.y = element_blank(), 
              axis.ticks.y = element_blank(),
              axis.title.y = element_blank(),
              
              panel.background = element_rect(fill = "gray95"),
              strip.background = element_blank(),
              strip.text = element_text(size = 7, hjust = 0),
              panel.border = element_blank()) 

## Combine plots 



model_reduction_plot <- plot_grid(csa_model_reduction, 
          plot_grid(NULL, str_model_reduction, rel_widths = c(0.20, 0.75), ncol = 2), 
          ncol = 1, rel_heights = c(0.4, 0.3))

saveRDS(model_reduction_plot, "./data/derivedData/study-1-determinants/model_reduction_plot.RDS")                    






