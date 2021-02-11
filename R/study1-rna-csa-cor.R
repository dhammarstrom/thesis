##-------------------------------------
## study1-rna-csa-cor.R
##
## Title: Correlations in study 1
## Purpose: Compile data for cerrelation analysis
## Author:
##
##
##
##-------------------------------------
## Notes: RNA content is likely related to protein synthesis rates (ref Millard 1973). 
# To explore RNA measured at different timepoint influence muscle growth, naive correlation
# and modeling controlling for sex, baseline measurements in muscle size and other values for RNA
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



#### Setting up the data ##########################

temp <- tot.rna %>%
        inner_join(mr) %>%
        pivot_wider(names_from = timepoint, 
                    values_from = rna) %>%
        
        
        print()
        
        
tot.rna %>%
        inner_join(mr) %>%      
        ggplot(aes(rna, change, color = sex, size = pre)) + geom_point() + 
        facet_wrap(~ timepoint) + 
        geom_smooth(method = "lm")

m.scaled0 <- lme(change ~  sets + scale(w0) +  scale(w2pre) + scale(w12), 
                random = list(subject = ~ 1), 
           
                data = temp, 
                na.action = na.omit)

plot(m.scaled0) # signs of heterogeneity, check if varExp helps

m.scaled <- lme(change ~  sets + scale(w0) +  scale(w2pre) + scale(w12), 
                    random = list(subject = ~ 1), 
                    weights = varExp(form = ~ fitted(.)), 
                    data = temp, 
                    na.action = na.omit)

anova(m.scaled0, m.scaled) 
plot(m.scaled)
# this removes the bad pattern in the residual plot and improves the fit

## Check if there is colinearity (rna measures could be correlated)

car::vif(m.scaled)

# Not an issue...


# Unscaled estimates for the regression table
m.unscaled <- lme(change ~  sets + w0 +  w2pre + w12, 
                         random = list(subject = ~ 1), 
                         weights = varExp(form = ~ fitted(.)), 
                 data = temp, 
                 na.action = na.omit)



summary(m.scaled)
summary(m.unscaled)

# Combine all values for regression table ############

rna_csa_change_models <- cbind(data.frame(coef(summary(m.unscaled))),
      data.frame(intervals(m.unscaled)$fixed)[,c(1,3)], 
      data.frame(intervals(m.scaled)$fixed)[,c(2,1,3)] %>%
              dplyr::select(scaled.est =  est.,
                            scaled.lower = lower, 
                            scaled.upper = upper)) %>%
        mutate(term = rownames(.)) %>%
        dplyr::select(term, 
                      Estimate = Value, 
                      SE = Std.Error, 
                      df = DF,
                      t.value, 
                      lower, 
                      upper, 
                      scaled.est, 
                      scaled.lower, 
                      scaled.upper) %>%
        
        print()
      
saveRDS(rna_csa_change_models, "./data/derivedData/study1-rna-csa-cor/regression.RDS")


## For plotting, individual time/RNA must be unscaled, but other scaled (to simplify plotting) 
# Refitting...
w0.unscaled <- lme(change ~  sets + w0 +  scale(w2pre) + scale(w12), 
                random = list(subject = ~ 1), 
                weights = varExp(form = ~ fitted(.)), 
                data = temp, 
                na.action = na.omit)

w2pre.unscaled <- lme(change ~  sets + scale(w0) +  w2pre + scale(w12), 
                   random = list(subject = ~ 1), 
                   weights = varExp(form = ~ fitted(.)), 
                   data = temp, 
                   na.action = na.omit)

w12.unscaled <- lme(change ~  sets + scale(w0) +  scale(w2pre) + w12, 
                      random = list(subject = ~ 1), 
                      weights = varExp(form = ~ fitted(.)), 
                      data = temp, 
                      na.action = na.omit)

## Create a plot of each RNA/timepoint againts change score in CSA


## Week 0


week0 <- tot.rna %>%
        inner_join(mr) %>%  
        filter(timepoint == "w0") %>%
        ggplot(aes(rna, change)) + geom_point(shape = 21, fill = group.study.color[5], alpha = 0.5) + 
        geom_smooth(method = "lm", 
                    se = FALSE, 
                    color = "gray50", 
                    lty = 2) +
        
        geom_segment(aes(x = 200, 
                         xend = 425, 
                         y = coef(summary(w0.unscaled))[1,1] + (coef(summary(w0.unscaled))[2,1])/2 + coef(summary(w0.unscaled))[3,1] * 200  ,
                         yend = coef(summary(w0.unscaled))[1,1] + (coef(summary(w0.unscaled))[2,1])/2 + coef(summary(w0.unscaled))[3,1] * 425), 
                     color = "gray30") +
        
        scale_x_continuous(limits = c(190, 700), 
                          expand = c(0,0),
                          breaks = c(200, 300, 400, 500, 600, 700)) +
        labs(x = "", 
             y = "CSA change (%)", 
             title = "Week 2") +
        
        
        theme_bw() +
        theme(panel.grid = element_blank(), 
              panel.border = element_rect(), 
              axis.line.x  = element_line(size = line.size), 
              axis.line.y = element_blank(), 
              axis.ticks.x = element_line(size = line.size), 
              axis.ticks.y = element_line(size = line.size), 
              axis.text.y =  element_text(color = "black", size = 7), 
              axis.title.y = element_text(color = "black", size = 7), 
              axis.text.x =  element_text(color = "black", size = 7), 
              axis.title.x = element_text(color = "black", size = 7), 
              title = element_text(color = "black", size = 7), 
              legend.title = element_blank(), 
              legend.background = element_rect(fill = "white"),
              legend.margin = margin(t = 0, r = 1, b = 1, l = 1, unit = "pt"),
              legend.key = element_rect(fill = "white"),
              legend.position = "none")
        
week2 <- tot.rna %>%
        inner_join(mr) %>%  
        filter(timepoint == "w2pre") %>%
        ggplot(aes(rna, change)) + geom_point(shape = 21, fill = group.study.color[5], alpha = 0.5) + 
        geom_smooth(method = "lm", 
                    se = FALSE, 
                    color = "gray50", 
                    lty = 2) +
        
        geom_segment(aes(x = 300, 
                         xend = 690, 
                         y = coef(summary(w2pre.unscaled))[1,1] + (coef(summary(w2pre.unscaled))[2,1])/2 + coef(summary(w2pre.unscaled))[4,1] * 300  ,
                         yend = coef(summary(w2pre.unscaled))[1,1] + (coef(summary(w2pre.unscaled))[2,1])/2 + coef(summary(w2pre.unscaled))[4,1] * 690), 
                     color = "gray30") +
        
        scale_x_continuous(limits = c(190, 700), 
                           expand = c(0,0),
                           breaks = c(200, 300, 400, 500, 600, 700)) +
        labs(x = "RNA per muscle weight (ng mg<sup>-1</sup>", 
             title = "Week 2") +
        
        theme_bw() +
        theme(panel.grid = element_blank(), 
              panel.border = element_rect(), 
              axis.line.x  = element_line(size = line.size), 
              axis.line.y = element_blank(), 
              axis.ticks.x = element_line(size = line.size), 
              axis.ticks.y = element_blank(), 
              axis.text.y =  element_blank(),
              axis.title.y = element_blank(),
              axis.text.x =  element_markdown(color = "black", size = 7), 
              axis.title.x = element_markdown(color = "black", size = 7), 
              title = element_text(color = "black", size = 7), 
              legend.title = element_blank(), 
              legend.background = element_rect(fill = "white"),
              legend.margin = margin(t = 0, r = 1, b = 1, l = 1, unit = "pt"),
              legend.key = element_rect(fill = "white"),
              legend.position = "none") 

week12 <- tot.rna %>%
        inner_join(mr) %>%  
        filter(timepoint == "w12") %>%
        ggplot(aes(rna, change)) + geom_point(shape = 21, fill = group.study.color[5], alpha = 0.5) + 
        geom_smooth(method = "lm", 
                    se = FALSE, 
                    color = "gray50", 
                    lty = 2) +
        
        geom_segment(aes(x = 200, 
                         xend = 590, 
                         y = coef(summary(w12.unscaled))[1,1] + (coef(summary(w12.unscaled))[2,1])/2 + coef(summary(w12.unscaled))[5,1] * 200  ,
                         yend = coef(summary(w12.unscaled))[1,1] + (coef(summary(w12.unscaled))[2,1])/2 + coef(summary(w12.unscaled))[5,1] * 590), 
                     color = "gray30") +
        
        scale_x_continuous(limits = c(190, 700), 
                           expand = c(0,0),
                           breaks = c(200, 300, 400, 500, 600, 700)) +
        labs(x = "", 
             title = "Week 12") +
        theme_bw() +
        theme(panel.grid = element_blank(), 
              panel.border = element_rect(), 
              axis.line.x  = element_line(size = line.size), 
              axis.line.y = element_blank(), 
              axis.ticks.x = element_line(size = line.size), 
              axis.ticks.y = element_blank(), 
              axis.text.y =  element_blank(),
              axis.title.y = element_blank(),
              axis.text.x =  element_text(color = "black", size = 7), 
              axis.title.x = element_text(color = "black", size = 7), 
              title = element_text(color = "black", size = 7), 
              legend.title = element_blank(), 
              legend.background = element_rect(fill = "white"),
              legend.margin = margin(t = 0, r = 1, b = 1, l = 1, unit = "pt"),
              legend.key = element_rect(fill = "white"),
              legend.position = "none") 

  


rna_csa_fig  <- plot_grid(week0, week2, week12, nrow = 1, 
                          align = "h",
          rel_widths = c(1, 0.95, 0.95))
        
    

saveRDS(rna_csa_fig, "./figures/study1-rna-csa.RDS")

