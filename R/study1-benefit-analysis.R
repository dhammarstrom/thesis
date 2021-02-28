##-------------------------------------
## study1-benefit-analysis.R
##
## Title: Benefit analysis in study 1
## Purpose: 
## Author:
##
##
##
##-------------------------------------


##############################################
#
#
# Analysis of Determinants of additional benefit of multiple sets 
#
# 
#
#
###############################################

## Libraries ##
source("./R/libraries.R")


### Content: 
## 1. Quantification of smallest worthwhile change and classification of benefit of multiple sets for muscle CSA
## 1.b and Strength data
## 2. RNA content 
## 3. mTOR pathway activation
## 4. Blood parameters
## 5. Myosin heavy chain distribution
## 6. Body composition
## 7. Training parameters

### Notes on model building strategies ###

# Data reduction is employed to "concentrate" the data set in the spirit  
# of Harrell 2015 Regression Modeling Strategies. All variables except training variables are
# screened for relation to SWC-grouping using a linear model controlling for sex. 
# This would be equialent to two-sample t-tests, as suggested by Hosmer, Lemeshow and Sturdivant (2013), with 
# the addition of the sex effect. In an earlier version, univariate (with sex a covariate) logistic regression was used: 

# benefit ~ sex + variable_i

# to increase interpretability, differences between SWC-groups was instead tested in a linear model using 

# variable_i ~ sex + benefit

# Differences associated with p-values < 0.2 was considered for multivariate logistic regression. The choice of
# t-test (or equialent regression) is unproblematic as described in Hosmer et al.

# The question in the logistic modelling is how different predictor affect the odds of benefit of multiple sets. 
# We choose this approach as we are not primarily interested in the magnitude of the difference between sets-conditions 
# but rather, on an individual level, benefit or not. This would in some sense "scale" the problem to each individual.

# Dichotomization of training-induced differences between sets-conditions are made on 
# the smalles worthwile change (SWC) calculated as:

# swc_j = 100 * (SD_mc * 0.2)/M_j

# where the swc is calculated per sex (j) based on standard deviation from both sexes using the standard deviation of
# mean centered pre-training between-subject data (SD_mc) times 0.2 and divided with sex-specific mean values (M_j). 
# The swc is calculated as the weighted mean of men and women.

# The sex specific calculation is done as both baseline CSA and strength has strong bimodal characteristics (sex effect).
# To get an representatiove SD pooled data is used. To get unbiased swc, sex specific values are calculated prior to avering.

# if the difference in hypertrophy between sets is greater than the average swc
# benefit of multiple-sets is declared. The same basic principle applies to molecular measurements
# (mTOR and riobsome biogenesis) when dichomization is done. See below for details



##### Functions for variable and model checks #####


### Linearity check plots linear predictors to check for linearity in the logit

# A key assumption of binomial regression with a continuous predictor is that the predictor
# is linear in the logit. On a small data-set this is not easily check, however it seems that 
# design variables (creating categorical variables from quartiles) and checking linearity in estimated
# coefficients works OK. This was is the primary linearity-check perfomed for each variable that 
# was identified in the screening process.

# In bimodal data (with biological confirmed sex differences, e.g. endocrine variables), design variables are constructed 
# per sex. 


linearity_check <-
        function(data, variable, formula, per.sex = TRUE) {
                temp <- data.frame(data)
                
                
                temp$fctr <- rep(NA, nrow(temp))
                
                temp$cont <- temp[, colnames(temp) == variable]
                
                
                
                if (per.sex) {
                        temp[temp$sex == "female",]$fctr <-
                                as.character(cut(
                                        temp[temp$sex == "female", colnames(temp) == variable],
                                        quantile(
                                                temp[temp$sex == "female", colnames(temp) == variable],
                                                na.rm = TRUE,
                                                prob = seq(0, 1, by = 0.25),
                                                type = 7
                                        ),
                                        include.lowest = TRUE,
                                        labels = c("f1", "f2", "f3", "f4")
                                ))
                        
                        temp[temp$sex == "male",]$fctr <-
                                as.character(cut(
                                        temp[temp$sex == "male", colnames(temp) == variable],
                                        quantile(
                                                temp[, colnames(temp) == variable],
                                                na.rm = TRUE,
                                                prob = seq(0, 1, by = 0.25),
                                                type = 7
                                        ),
                                        include.lowest = TRUE,
                                        labels = c("f1", "f2", "f3", "f4")
                                ))
                        
                        temp2 <-
                                aggregate(temp[, colnames(temp) == variable], list(temp$fctr, temp$sex), median)
                        
                        
                        
                        colnames(temp2) <- c("factor", "sex", "m")
                        
                        form <- as.formula("benefit ~ sex + fctr")
                        
                        model <- temp %>%
                                glm(form, data = ., method = "brglmFit")
                        
                        
                        
                        temp3 <-  data.frame(coef(model)) %>%
                                mutate(factor = row.names(.)) %>%
                                mutate(
                                        factor = if_else(factor == "(Intercept)", "f1", gsub("fctr", "", factor)),
                                        coef = if_else(factor == "f1", 0, coef.model.)
                                ) %>%
                                inner_join(temp2)
                        
                        plot1 <- temp3 %>%
                                ggplot(aes(m, coef)) + geom_line() + ggtitle(variable) + facet_grid(sex ~ .)
                        
                        
                        
                        plot2 <-
                                ggplot() + geom_density(data = temp, aes(cont, color = sex)) +
                                geom_vline(data = temp3, aes(
                                        xintercept = m,
                                        color = sex,
                                        lty = factor
                                )) + facet_grid(sex ~ .)
                        
                        
                        return(cowplot::plot_grid(plot1, plot2, ncol = 2))
                        
                }
                
                if (!per.sex) {
                        temp$fctr <-
                                cut(
                                        temp[, colnames(temp) == variable],
                                        quantile(
                                                temp[, colnames(temp) == variable],
                                                na.rm = TRUE,
                                                prob = seq(0, 1, by = 0.25),
                                                type = 7
                                        ),
                                        include.lowest = TRUE,
                                        labels = c("f1", "f2", "f3", "f4")
                                )
                        
                        
                        temp$cont <- temp[, colnames(temp) == variable]
                        
                        
                        temp2 <-
                                aggregate(temp[, colnames(temp) == variable], list(temp$fctr), mean)
                        
                        colnames(temp2) <- c("factor", "m")
                        
                        
                        form <- as.formula("benefit ~ sex + fctr")
                        
                        model <- temp %>%
                                glm(form, data = ., method = "brglmFit")
                        
                        
                        temp3 <-  data.frame(coef(model)) %>%
                                mutate(factor = row.names(.)) %>%
                                mutate(
                                        factor = if_else(factor == "(Intercept)", "f1", gsub("fctr", "", factor)),
                                        coef = if_else(factor == "f1", 0, coef.model.)
                                ) %>%
                                inner_join(temp2)
                        
                        plot1 <- temp3 %>%
                                ggplot(aes(m, coef)) + geom_line() + ggtitle(variable)
                        
                        
                        
                        plot2 <- ggplot() + geom_density(data = temp, aes(cont))
                        
                        
                        return(cowplot::plot_grid(plot1, plot2, ncol = 2))
                        
                        
                }
                
        }

### Logit loess plot ### 
# Not recommended for use in small samples (n < 50). As an extra check to the design variables it identifies potential 
# influential data points. May thus be used as complementary method.

# Code from: https://thestatsgeek.com/2014/09/13/checking-functional-form-in-logistic-regression-using-loess/

logitloess <- function(x, y, s, cf) {
        logit <- function(pr) {
                log(pr / (1 - pr))
        }
        
        if (missing(s)) {
                locspan <- 0.7
        } else {
                locspan <- s
        }
        
        loessfit <- predict(loess(y ~ x, span = locspan))
        pi <- pmax(pmin(loessfit, 0.9999), 0.0001)
        logitfitted <- logit(pi)
        
        plot(x, logitfitted, ylab = "logit", col = cf)
        
}



### Check for confunding effect of variable
# delta_beta_fun calculates % change of predictors when removing other predictors.
# In variable selection, a delta-beta of less than ~20% suggests that the the removed 
# variable did not affect other variables that much (see Hosmer et al.). The function 
# calculates delta-beta for a full and reduced model.

delta_beta_fun <- function(m1, m2) {
        coef(summary({
                m2
        })) %>%
                data.frame() %>%
                mutate(parameter = rownames(.),
                       estimate_reduced = Estimate) %>%
                dplyr::select(parameter, estimate_reduced) %>%
                filter(parameter != "(Intercept)") %>%
                inner_join(
                        coef(summary({
                                m1
                        })) %>%
                                data.frame() %>%
                                mutate(parameter = rownames(.),
                                       estimate_full = Estimate) %>%
                                dplyr::select(parameter, estimate_full)
                ) %>%
                mutate(delta.beta = abs((estimate_reduced - estimate_full) / estimate_full) * 100) %>%
                return()
        
        
}

##### Quantification of measurement error and classification of benefit of multiple-sets #####
# 
# The average CSA variation (%CV) between legs at baseline, expressed as a mean percentage was regarded as the 
# measurement error in an earlier iteration of the manuscript. In order to classify as benefit of multiple-sets, the difference between gains 
# must be larger than the between legs measurement error. 

# In peer-review, the smallest worthwile change (SWC) was brought up as an alternative to the between leg variation.
# We have previously calculated the measurement error as the average whitin subject coefficient of variation. This 
# is a conservative measure of measurement error. The SWC is more intuitive and comparative to other studies.

# SWC is calculated as 0.2 * SD and expressed as a percentage to the mean. To account for sex-differences
# (CSA and strength has bimodal distribution), the standard deviation is calculated on data mean-centered per
# sex. The swc = 0.2 * sd_sex is divided by the sex specific mean and the weighted average of the two are
# used as the estimate of the SWC.


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
        dplyr::select(subject, sex, benefit, benefit.single, diff) %>%
        print()



##### Strength data #####

# In the review process, we were challenged to change the measurement error to SWC to potentially find other predictors
# of multiple-set benefit. In the process we decided to combined all strength measures into one using a weighted average of 
# transformed strength values x_i/max(bar(x)). Weights are calculated per test-family (1RM, isometric and isokinetic tests). The
# combined measure should in theory reflect generic muscular strength more accurate.


### Combine all strength data into one data set 

include <- read.csv("./data/study-1/oneThreeSetLeg.csv", sep=";") %>%
        gather(sets,leg, multiple:single) %>%
        mutate(condition = leg) %>%
        dplyr::select(subject, sex, include, condition, sets)


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
        group_by(subject, sex, sets, week, ex_fam) %>%
        summarise(str = mean(str, na.rm = TRUE)) %>%
        group_by(subject, sex, sets, week) %>%
        summarise(str = mean(str, na.rm = TRUE)) %>%
        print()

# Smallest wortwhile change in strength exercises is calculated as the average of the different
# exercise groups. Failing to calculate per exercise group would underestimate the SD as averaging would 
# lead to "regression towards the mean effect".

strength.error <- cybex %>%
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
        rbind(rm_str) %>% # This includes all strength tests %>%
        mutate(ex_fam = if_else(exercise %in% c("legext", "legpress"), 
                                "rm", if_else(exercise == 0, "isom", "isok"))) %>%
        dplyr::select(subject, week, sex, sets, str, ex_fam) %>%
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


## Save total_strength for descriptive analysis ##
# saveRDS(total_strength, "./derivedData/total_strength.Rda")



str.benefit <- total_strength %>%
        spread(week, str) %>%
        mutate(rel_increase = ((post/baseline)-1)*100) %>%
        ungroup() %>%
        group_by(subject) %>%
        mutate(strength.baseline = mean(baseline, na.rm = TRUE)) %>%
        group_by(sex) %>%
        mutate(strength.baseline = strength.baseline - mean(strength.baseline, na.rm = TRUE)) %>% # center strength baseline per sex
        ungroup() %>%
        dplyr::select(subject, sex, sets, strength.baseline, rel_increase) %>%
        spread(sets, rel_increase) %>%
        mutate(cv = as.numeric(strength.error)) %>%
        mutate(strength.diff = multiple - single, 
               benefit.strength = if_else(multiple - single > cv, 1, 0),
               benefit.strength.single = if_else(single - multiple > cv, 1, 0)) %>%
        dplyr::select(subject, sex, strength.baseline, benefit = benefit.strength, benefit.single = benefit.strength.single, diff = strength.diff) %>%
        print()




strength <- total_strength %>%
        inner_join(read_excel("./data/study-1/body-composition/bodycomp_DXA.xlsx") %>%
                           inner_join(read.csv("./data/study-1/oneThreeSetLeg.csv", sep=";")) %>%
                           filter(include =="incl" & timepoint == "pre") %>% 
                           dplyr::select(subject, weight) %>%
                           mutate(body.mass = weight)) %>% 
        filter(week == "baseline") %>%
        mutate(rel.strength = (str/body.mass) * 1000) %>%
        group_by(subject, sex) %>%
        summarise(rel.strength = mean(rel.strength, na.rm = TRUE)) %>%
        group_by(sex) %>%
        mutate(rel.strength = rel.strength - mean(rel.strength), 
               rel.strength.factor =   as.character(cut(rel.strength, quantile(rel.strength, na.rm = TRUE, 
                                                                               prob = seq(0, 1, by = 0.25),
                                                                               type = 7), include.lowest = TRUE, 
                                                        labels = c("f1", "f2", "f3", "f4")))) %>%
        print()

### Save data sets in two versions for variable selection

strength.cont <- strength  %>% dplyr::select(subject, sex, rel.strength)
strength.fact <- strength  %>% dplyr::select(subject, sex, rel.strength.factor)



#### RNA content #####

include <- read.csv("./data/study-1/oneThreeSetLeg.csv", sep = ";") %>%
        gather(sets, leg, multiple:single) %>%
        mutate(condition = leg) %>%
        dplyr::select(subject, sex, include, condition, sets)

# loads sample weight and total RNA concentration in trizol extracts 
total.rna <- read_excel("./data/study-1/sample_weight.xlsx") %>%
        dplyr::select(subject, timepoint, leg, weight, conc1.5, elution.volume) %>% # total RNA measured in 1:5 dilution
        mutate(condition = leg) %>%
        inner_join(include) %>%
        filter(include == "incl") %>%
        mutate(rna = (conc1.5 * 5 * elution.volume) / weight) %>% # estimates of total RNA
        dplyr::select(subject,
                      sex,
                      sets,
                      timepoint,
                      rna,
                      weight,
                      conc1.5,
                      elution.volume) %>%
        mutate(timepoint = factor(timepoint, levels = c("w0", "w2pre", "w2post", "w12")),
               conc = conc1.5 * 5 * elution.volume)# %>% # scale RNA concentration based on dilution
#  filter(timepoint !="w2post")



# Some samples may be underestimated because of sample loss during 
# sample prep. This was due to loss of RNA in washing steps were piece of 
# RNA pellet was lost. To estimate deviance from weight to RNA relationship
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
total.rna <- total.rna%>%
        mutate(outlier = if_else(resid < -4 | resid > 4 , "out", "in")) # 4 units deviance from zero

### RNA variables 

# Calculate measurement error in RNA measurement

rna.cv <- total.rna %>%
        filter(outlier != "out") %>%
        group_by(subject, timepoint, sets) %>%
        summarise(rna = mean(rna, na.rm = T)) %>%
        filter(timepoint == "w0") %>%
        spread(sets, rna) %>%
        rowwise() %>%
        mutate(s2 = ((multiple - single) ^ 2) / 2,
               m2 = (mean(c(multiple, single), na.rm = TRUE)) ^ 2,
               s2m2 = s2 / m2) %>%
        ungroup() %>%
        summarise(cv = 100 * sqrt(mean(s2m2))) %>%
        data.frame() %>%
        print()




tot.rna <- total.rna %>%
        filter(outlier != "out") %>%
        group_by(subject, timepoint, sets) %>%
        summarise(rna = mean(rna, na.rm = T)) %>%
        spread(sets, rna) %>%
        mutate(
                rna_sets_ratio = ((multiple / single) - 1) * 100,
                rna_dic = if_else(rna_sets_ratio > rna.cv, 1, 0)
        ) %>%
        dplyr::select(subject, timepoint, rna_sets_ratio, rna_dic) %>%
        unite(rna, c(rna_sets_ratio, rna_dic)) %>%
        spread(timepoint, rna) %>%
        separate(w0,
                 c("rna_w0", "rna.dic_w0"),
                 sep = "_",
                 convert = TRUE) %>%
        separate(w2pre,
                 c("rna_w2pre", "rna.dic_w2pre"),
                 sep = "_",
                 convert = TRUE) %>%
        separate(w2post,
                 c("rna_w2post", "rna.dic_w2post"),
                 sep = "_",
                 convert = TRUE) %>%
        separate(w12,
                 c("rna_w12", "rna.dic_w12"),
                 sep = "_",
                 convert = TRUE) %>%
        dplyr::select(subject,
                      # rna_w0,
                      # rna.dic_w0,
                      rna_w2pre,
                      rna_w2post,
                      # rna.dic_w2pre,
                      rna_w12, #,
                      rna.dic_w12) %>%
        #  mutate(rna_w2pre = if_else(subject == "FP15", 5.682378, rna_w2pre)) %>% # For comparisons imputed data (see below)
        dplyr::select(-rna_w2post) %>%
        ungroup() %>%
        print()

tot.rna.cont <-
        tot.rna %>% dplyr::select(subject, rna_w2pre, rna_w12)



### Western blot data: mTOR pathway activation ####

# The mTOR pathway is activated acutely in response to the fifth session. The fold-change activation between sets are
# used as a ratio of acute activation. 

#### Read data and combine data frames ####
# read western data
west <-
        read_excel("./data/study-1/westernResults.xlsx", sheet = "run1")
west2 <-
        read_excel("./data/study-1/westernResults.xlsx", sheet = "run2")
west.acute <-
        read_excel("./data/study-1/westernResults.xlsx", sheet = "runAcute")

# sample setups
sl <-
        read_excel("./data/study-1/westernResults.xlsx", sheet = "samples") %>%
        filter(include == "incl") %>%
        inner_join(
                read.csv("./data/study-1/oneThreeSetLeg.csv", sep = ";") %>%
                        gather(sets, leg, multiple:single) %>%
                        filter(include == "incl")
        )


#### Mean center expression data and collect data from all runs ####

# West acute only contains w2pre and w2post, used to validate p70 expression
# possibly supplementary data
west.acute <- west.acute %>%
        mutate(totalprotein = (mean.gray1 + mean.gray2) / 2 - (mean.graybg1 + mean.graybg2) /
                       2) %>%
        inner_join(sl) %>%
        dplyr::select(run,
                      gel,
                      well,
                      subject,
                      leg,
                      sets,
                      timepoint,
                      date,
                      totalprotein,
                      p.p70) %>%
        gather(target, expression, totalprotein:p.p70) %>%
        group_by(date, gel, target) %>%
        mutate(expression = expression / mean(expression, na.rm = T),
               run = "westacute") %>%
        ungroup()

west <- west %>%
        mutate(totalprotein = (mean.gray1 + mean.gray2) / 2 - (mean.graybg1 + mean.graybg2) /
                       2) %>%
        inner_join(sl) %>%
        dplyr::select(
                run,
                gel,
                well,
                subject,
                leg,
                sets,
                timepoint,
                date,
                totalprotein,
                p.p70
        ) %>%
        gather(target, expression, totalprotein:p.p70) %>%
        group_by(date, gel, target) %>%
        mutate(expression = expression / mean(expression, na.rm = T),
               run = "west1") %>%
        ungroup()

west2 <- west2 %>%
        mutate(totalprotein = (mean.gray1 + mean.gray2) / 2 - (mean.graybg1 + mean.graybg2) /
                       2) %>%
        inner_join(sl) %>%
        dplyr::select(run,
                      gel,
                      well,
                      subject,
                      leg,
                      sets,
                      timepoint,
                      date,
                      totalprotein,
                      t.mTOR,
                      p.s6,
                      t.s6) %>%
        gather(target, expression, totalprotein:t.s6) %>%
        group_by(date, gel, target) %>%
        mutate(expression = expression / mean(expression, na.rm = T),
               run = "west2") %>%
        ungroup()

## Comnbine data from full runs
wb.raw <- rbind(west, west.acute)


wb.data <- wb.raw %>%
        spread(target, expression) %>%
        mutate(
                norm.expr = p.p70 / totalprotein,
                timepoint = factor(timepoint, levels = c("w0", "w2pre", "w2post", "w12")),
                sets = factor(sets, levels = c("single", "multiple")),
                leg = factor(paste0(subject, sets))
        ) %>%
        group_by(run, subject, sets, timepoint) %>%
        summarise(norm.expr = mean(norm.expr, na.rm = TRUE)) %>%
        filter(timepoint == "w2post") %>%
        group_by(run) %>%
        mutate(norm.expr = norm.expr / max(norm.expr)) %>%
        group_by(subject, sets) %>%
        summarise(norm.expr = mean(norm.expr, na.rm = TRUE)) %>%
        spread(sets, norm.expr) %>%
        mutate(p.p70_w2post = multiple / single) %>%
        print()



### Saves continuous variables for primary screening
wb.data.cont <- wb.data %>% dplyr::select(subject, p.p70_w2post)

## Saves full data set (with dichotomized variables)
wb.data <- wb.data %>%
        dplyr::select(subject, p.p70_w2post) %>%
        data.frame() 




#### Blood parameters #####


blood_data <- read_excel("./data/study-1/blood/blood_parameters.xlsx") %>%
        gather(analyte, concentration, andro:tsh) %>%
        filter(analyte %in% c("cortisol", "gh", "igf1", "testo")) %>%
        inner_join(read.csv("./data/study-1/oneThreeSetLeg.csv", sep = ";")) %>%
        filter(include == "incl") %>%
        dplyr::select(subject, timepoint, analyte, concentration, sex) %>%
        filter(timepoint != "w12") %>%
        mutate(
                timepoint = factor(timepoint,
                                   levels = c("w0", "w2pre", "w2post", "w2post60")),
                time_group = if_else(
                        analyte %in% c("testo", "cortisol"),
                        if_else(timepoint %in% c("w0", "w2pre", "w2post", "w2post60"), "ba", "na"),
                        if_else(timepoint %in% c("w2post", "w2post60"), "tr", "ba")
                )
        ) %>%
        
        # "tr" for training, "ba" for baseline, week2 for mean week 0-2
        group_by(subject, sex, analyte, time_group) %>%
        summarise(m.conc = mean(concentration, na.rm = TRUE)) %>%
        ungroup() %>%
        mutate(analyte_time = paste(analyte, time_group, sep = ".")) %>%
        dplyr::select(subject, sex, analyte_time, m.conc) %>%
        spread(analyte_time, m.conc) %>%
        dplyr::select(subject, sex, cortisol.ba, testo.ba, gh.tr, igf1.ba, igf1.tr) %>%
        # Female values below detection are coded as 0
        mutate(testo.ba = if_else(
                sex == "female" &
                        !is.na(cortisol.ba) & is.na(testo.ba),
                0,
                testo.ba
        )) %>%
        print()


vitd <- read_excel("./data/study-1/blood/blood_parameters.xlsx") %>%
        gather(analyte, concentration, andro:tsh) %>%
        filter(analyte %in% c("vitd")) %>%
        group_by(subject) %>%
        summarise(vitd = mean(concentration, na.rm = TRUE)) %>%
        print()



# Save continuous data set
blood_data.cont <- blood_data %>%
        inner_join(vitd)



#### Myosin heavy chain distribution #####

## Tranform data ##

# Type 2 AX counted as 0.5 A, 0.5 X
fibertype_pre <- read_excel("./data/study-1/immuno/fibertype_checked.xlsx") %>%
        mutate(
                Type2A = type2a + 0.5 * type2ax,
                Type2X = type2x + 0.5 * type2ax,
                Type1 = type1,
                time = ifelse(timepoint == 1, "w0", ifelse(timepoint %in% c(2, 3), "w2", "w12"))
        ) %>%
        dplyr::select(subject, time, leg, Type1, Type2A, Type2X) %>%
        inner_join(
                read.csv("./data/study-1/oneThreeSetLeg.csv", sep = ";") %>%
                        gather(sets, leg, multiple:single) %>%
                        filter(include == "incl")
        ) %>%
        gather(fibertype, counts, Type1:Type2X) %>%
        group_by(subject, sex, sets, time, fibertype) %>%
        summarise(counts = sum(counts)) %>%
        spread(fibertype, counts) %>%
        rowwise() %>%
        mutate(
                sum = as.integer(round(Type1 + Type2A + Type2X, 0)),
                type1 = (Type1 / sum) * 100,
                type2a = (Type2A / sum) * 100,
                type2x = (Type2X / sum) * 100
        ) %>%
        ungroup() %>%
        dplyr::select(subject, sex, sets, time, sum, type1:type2x) %>%
        mutate(time = factor(time, levels = c("w0", "w2", "w12")),
               sets = factor(sets, levels = c("single", "multiple"))) %>%
        gather(fibertype, percentage, type1:type2x) %>%
        dplyr::select(subject, sex, sets, time, fibertype, percentage) %>%
        
        spread(sets, percentage) %>%
        filter(time %in%  c("w0", "w2")) %>%
        rowwise() %>%
        mutate(avg = mean(c(single, multiple), na.rm = TRUE)) %>%
        ungroup() %>%
        dplyr::select(subject, sex, fibertype, time, avg) %>%
        group_by(fibertype, subject, sex) %>%
        summarise(avg = mean(avg, na.rm = TRUE)) %>%
        spread(fibertype, avg) %>%
        #  mutate(type2x.dic = if_else(type2x > 5, 1, 0)) %>% # used for additional preliminary models
        dplyr::select(subject, sex, type2a, type2x,  type1) %>%
        print()

### Save continous data set

fibertype.cont <- fibertype_pre



#### Body composition ####

body.comp <- read_excel("./data/study-1/body-composition/bodycomp_DXA.xlsx") %>%
        inner_join(read.csv("./data/study-1/oneThreeSetLeg.csv", sep = ";")) %>%
        filter(include == "incl" & timepoint == "pre") %>%
        dplyr::select(subject, sex, BMD.whole, fat.whole, lean.whole) %>%
        rowwise() %>%
        mutate(whole = sum(c(BMD.whole, fat.whole, lean.whole))) %>%
        ungroup() %>%
        mutate(lean.mass = (lean.whole / whole) * 100) %>%
        group_by(sex) %>%
        mutate(lean.mass = lean.mass - mean(lean.mass)) %>%
        dplyr::select(subject, sex, lean.mass) %>%
        print()

lean.mass.cont <- body.comp

#### Training adherence  ######

train1 <- read_excel("./data/study-1/trainingdata.xlsx", sheet = 2)

train1$performed.sessions <- (train1$total.sessions/31) *100
train1$supervised.sessions <- (1 - (train1$nonsupervised.sessions / train1$total.sessions))*100

train1$subject <- paste("FP", train1$subject, sep = "")

train1 <- train1 %>%
        dplyr::select(subject, total.sessions, supervised.sessions) %>%
        mutate(total.sessions.dic = if_else(total.sessions == 31, 1, 0),
               supervised.sessions.dic = if_else(supervised.sessions == 100, 1, 0)) %>%
        dplyr::select(subject, total.sessions.dic, supervised.sessions.dic) %>%
        print()



## Save data sets for descriptive statistics

saveRDS(benefit, "./data/derivedData/study-1-benefit-analysis/benefit.Rda") 
saveRDS(str.benefit, "./data/derivedData/study-1-benefit-analysis/strbenefit.Rda")
saveRDS(list(csa.error = csa.error, strength.cv = strength.error), "./data/derivedData/study-1-benefit-analysis/cvBenefit.Rda")

#### Combine data-sets ######

# Continuous data sets #

benefit_csa_complete_cont <- benefit[,c(1,2,3)] %>%
        inner_join(tot.rna.cont) %>%
        inner_join(wb.data.cont) %>%
        inner_join(blood_data.cont) %>%
        inner_join(strength.cont) %>%
        inner_join(str.benefit %>%
                           dplyr::select(-benefit, -benefit.single, -diff)) %>%
        
        inner_join(body.comp) %>%
        inner_join(fibertype.cont) %>%
        # inner_join(train1) %>%
        filter(!(subject %in% c("FP15","FP23"))) %>% # Filters complete cases
        print()

benefit_strength_complete_cont <- str.benefit[,c(1,2,4)] %>%
        inner_join(tot.rna.cont) %>%
        inner_join(wb.data.cont) %>%
        inner_join(blood_data.cont) %>%
        inner_join(strength.cont) %>%
        inner_join(str.benefit %>%
                           dplyr::select(-benefit, -benefit.single, -diff)) %>%
        inner_join(body.comp) %>%
        inner_join(fibertype.cont) %>%
        # inner_join(train1) %>%
        filter(!(subject %in% c("FP15", "FP23"))) %>% # Filters complete cases
        print()


### Save data for descriptive statistics

saveRDS(benefit_csa_complete_cont, "./data/derivedData/study-1-benefit-analysis/benefit_csa_rawdata.Rda")
saveRDS(benefit_strength_complete_cont, "./data/derivedData/study-1-benefit-analysis/benefit_str_rawdata.Rda")


univar_csa <- colnames(benefit_csa_complete_cont)[-c(1, 2, 3)]

univar_str <- colnames(benefit_strength_complete_cont)[-c(1, 2, 3)]


### Creates a result-df for the screening loop

results <- data.frame(variable = rep(NA, length(univar_csa)),
                      estimate.csa = rep(NA, length(univar_csa)),
                      se.csa = rep(NA, length(univar_csa)), 
                      tval.csa = rep(NA, length(univar_csa)),
                      pval.csa = rep(NA, length(univar_csa)),
                      ## Added ci for diff between benefits
                      csa.ci.lwr = rep(NA, length(univar_csa)),
                      csa.ci.upr = rep(NA, length(univar_csa)),
                      csa.ci.lwr.80 = rep(NA, length(univar_csa)),
                      csa.ci.upr.80 = rep(NA, length(univar_csa)),
                      
                      
                      
                      sex.pval.csa = rep(NA, length(univar_csa)),
                      sex.inter.pval.csa = rep(NA, length(univar_csa)),
                      m.b.csa = rep(NA, length(univar_csa)),
                      cil.b.csa = rep(NA, length(univar_csa)),
                      ciu.b.csa = rep(NA, length(univar_csa)),
                      m.nb.csa = rep(NA, length(univar_csa)),
                      cil.nb.csa = rep(NA, length(univar_csa)),
                      ciu.nb.csa = rep(NA, length(univar_csa)),
                      m.b.male.csa = rep(NA, length(univar_csa)),
                      cil.b.male.csa = rep(NA, length(univar_csa)),
                      ciu.b.male.csa = rep(NA, length(univar_csa)),
                      m.nb.male.csa = rep(NA, length(univar_csa)),
                      cil.nb.male.csa = rep(NA, length(univar_csa)),
                      ciu.nb.male.csa = rep(NA, length(univar_csa)),
                      m.b.female.csa = rep(NA, length(univar_csa)),
                      cil.b.female.csa = rep(NA, length(univar_csa)),
                      ciu.b.female.csa = rep(NA, length(univar_csa)),
                      m.nb.female.csa = rep(NA, length(univar_csa)),
                      cil.nb.female.csa = rep(NA, length(univar_csa)),
                      ciu.nb.female.csa = rep(NA, length(univar_csa)),
                      estimate.str = rep(NA, length(univar_csa)),
                      se.str = rep(NA, length(univar_csa)), 
                      tval.str = rep(NA, length(univar_csa)),
                      pval.str = rep(NA, length(univar_csa)),
                      ## Added ci for diff between benefits
                      str.ci.lwr = rep(NA, length(univar_csa)),
                      str.ci.upr = rep(NA, length(univar_csa)),
                      str.ci.lwr.80 = rep(NA, length(univar_csa)),
                      str.ci.upr.80 = rep(NA, length(univar_csa)),
                      
                      sex.pval.str = rep(NA, length(univar_csa)),
                      sex.inter.pval.str = rep(NA, length(univar_csa)),
                      m.b.str = rep(NA, length(univar_csa)),
                      cil.b.str = rep(NA, length(univar_csa)),
                      ciu.b.str = rep(NA, length(univar_csa)),
                      m.nb.str = rep(NA, length(univar_csa)),
                      cil.nb.str = rep(NA, length(univar_csa)),
                      ciu.nb.str = rep(NA, length(univar_csa)),
                      m.b.male.str = rep(NA, length(univar_csa)),
                      cil.b.male.str = rep(NA, length(univar_csa)),
                      ciu.b.male.str = rep(NA, length(univar_csa)),
                      m.nb.male.str = rep(NA, length(univar_csa)),
                      cil.nb.male.str = rep(NA, length(univar_csa)),
                      ciu.nb.male.str = rep(NA, length(univar_csa)),
                      m.b.female.str = rep(NA, length(univar_csa)),
                      cil.b.female.str = rep(NA, length(univar_csa)),
                      ciu.b.female.str = rep(NA, length(univar_csa)),
                      m.nb.female.str = rep(NA, length(univar_csa)),
                      cil.nb.female.str = rep(NA, length(univar_csa)),
                      ciu.nb.female.str = rep(NA, length(univar_csa)), 
                      n = rep(NA, length(univar_csa)))



# This loop takes ever single predictor and tests for group differences in csa and strength 
# swc-groups (multiple-benefit).

for(i in 1:length(univar_csa)) {
        
        # An interaction between variables a sex 
        formula.csa.inter <- as.formula(paste("scale(",univar_csa[i], ") ~ sex * as.factor(benefit)"))
        formula.str.inter <- as.formula(paste("scale(",univar_str[i], ")~ sex * as.factor(benefit)"))
        
        # The linear model controlling for sex
        formula.csa <- as.formula(paste("scale(",univar_csa[i], ")~ sex + as.factor(benefit)"))
        formula.str <- as.formula(paste("scale(",univar_str[i], ")~ sex + as.factor(benefit)"))
        
        
        model.csa <- lm(formula.csa, data = benefit_csa_complete_cont)
        model.str <- lm(formula.str, data = benefit_strength_complete_cont)
        
        
        model.csa.inter <- lm(formula.csa.inter, data = benefit_csa_complete_cont)
        model.str.inter <- lm(formula.str.inter, data = benefit_strength_complete_cont)
        
        
        if(univar_csa[i] != univar_str[i]) stop() # Adds an check to make sure the variables are aligned in the two df
        
        # Extract estimates
        csa.em <- data.frame(emmeans(model.csa, specs = ~"benefit"))
        csa.sex.em <- data.frame(emmeans(model.csa.inter, specs = ~"benefit|sex"))
        
        str.em <- data.frame(emmeans(model.str, specs = ~"benefit"))
        str.sex.em <- data.frame(emmeans(model.str.inter, specs = ~"benefit|sex"))
        
        # Saves the results
        results[i,1] <- univar_csa[i]
        
        results[i, 2]  <- coef(summary(model.csa))[3, 1]
        results[i, 3]  <- coef(summary(model.csa))[3, 2]
        results[i, 4]  <- coef(summary(model.csa))[3, 3]
        results[i, 5]  <- coef(summary(model.csa))[3, 4]
        results[i, 6] <- confint(model.csa)[3,1]
        results[i, 7] <- confint(model.csa)[3,2]
        results[i, 8] <- confint(model.csa, level = 0.8)[3,1]
        results[i, 9] <- confint(model.csa, level = 0.8)[3,2]
        
        
        results[i, 10]  <- coef(summary(model.csa))[2, 4]
        results[i, 11]  <- coef(summary(model.csa.inter))[4, 4]
        results[i, 12]  <- csa.em[2, 2]
        results[i, 13]  <- csa.em[2, 5]
        results[i, 14] <- csa.em[2, 6]
        results[i, 15] <- csa.em[1, 2]
        results[i, 16] <- csa.em[1, 5]
        results[i, 17] <- csa.em[1, 6]
        results[i, 18] <- csa.sex.em[4, 3]
        results[i, 19] <- csa.sex.em[4, 6]
        results[i, 20] <- csa.sex.em[4, 7]
        results[i, 21] <- csa.sex.em[3, 3]
        results[i, 22] <- csa.sex.em[3, 6]
        results[i, 23] <- csa.sex.em[3, 7]
        results[i, 24] <- csa.sex.em[2, 3]
        results[i, 25] <- csa.sex.em[2, 6]
        results[i, 26] <- csa.sex.em[2, 7]
        results[i, 27] <- csa.sex.em[1, 3]
        results[i, 28] <- csa.sex.em[1, 6]
        results[i, 29] <- csa.sex.em[1, 7]
        
        
        results[i, 30] <- coef(summary(model.str))[3, 1]
        results[i, 31] <- coef(summary(model.str))[3, 2]
        results[i, 32] <- coef(summary(model.str))[3, 3]
        results[i, 33] <- coef(summary(model.str))[3, 4]
        results[i, 34] <- confint(model.str)[3,1]
        results[i, 35] <- confint(model.str)[3,2]
        results[i, 36] <- confint(model.str, level = 0.8)[3,1]
        results[i, 37] <- confint(model.str, level = 0.8)[3,2]
        
        results[i, 38] <- coef(summary(model.str))[2, 4]
        results[i, 39] <- coef(summary(model.str.inter))[4, 4]
        results[i, 40] <- str.em[2, 2]
        results[i, 41] <- str.em[2, 5]
        results[i, 42] <- str.em[2, 6]
        results[i, 43] <- str.em[1, 2]
        results[i, 44] <- str.em[1, 5]
        results[i, 45] <- str.em[1, 6]
        results[i, 46] <- str.sex.em[4, 3]
        results[i, 47] <- str.sex.em[4, 6]
        results[i, 48] <- str.sex.em[4, 7]
        results[i, 49] <- str.sex.em[3, 3]
        results[i, 50] <- str.sex.em[3, 6]
        results[i, 51] <- str.sex.em[3, 7]
        results[i, 52] <- str.sex.em[2, 3]
        results[i, 53] <- str.sex.em[2, 6]
        results[i, 54] <- str.sex.em[2, 7]
        results[i, 55] <- str.sex.em[1, 3]
        results[i, 56] <- str.sex.em[1, 6]
        results[i, 57] <- str.sex.em[1, 7]
        results[i, 58] <- nrow(model.str$model)
        
        
}

# Saves for tables (see below)
results_univariate <- results

saveRDS(results_univariate, "./data/derivedData/study-1-benefit-analysis/univariate_analysis.Rda")






## Inclusion of factor variables ##
library(brglm2)
method.fit <- "brglmFit"

# Training variables were factorized due to funky numbers... 


csa.dic <- benefit[,c(1,2,3)] %>%
        inner_join(train1) %>%
        print()

str.dic <- str.benefit[,c(1,2,4)] %>%
        inner_join(train1) %>%
        print()


total.session.csa <- glm(benefit ~ total.sessions.dic, data = csa.dic, method = method.fit, family = binomial("logit"))
supervised.session.csa <- glm(benefit ~ supervised.sessions.dic, data = csa.dic, method = method.fit, family = binomial("logit"))

total.session.str <- glm(benefit ~ total.sessions.dic, data = str.dic, method = method.fit, family = binomial("logit"))
supervised.session.str <- glm(benefit ~ supervised.sessions.dic, data = str.dic, method = method.fit, family = binomial("logit"))


training_data_for_table <- list(csa.dic = csa.dic, str.dic = str.dic, 
                                total.session.csa = total.session.csa, 
                                supervised.session.csa = supervised.session.csa, 
                                total.session.str = total.session.str,
                                supervised.session.str = supervised.session.str)


#### Baseline training 

baseline_training <- read_excel("./data/study-1/baseline_training.xlsx") %>%
        mutate(activity = if_else(as.numeric(amount) > 0, 1, 0)) %>%
        dplyr::select(subject, activity, any_strength_type) %>%
        print()

csa.dic.baseline <- benefit[,c(1,2,3)] %>%
        inner_join(baseline_training) %>%
        print()

str.dic.baseline <- str.benefit[,c(1,2,4)] %>%
        inner_join(baseline_training) %>%
        print()



training.csa <- glm(benefit ~ activity, data = csa.dic.baseline, method = method.fit, family = binomial("logit"))
str.training.csa <- glm(benefit ~ any_strength_type, data = csa.dic.baseline, method = method.fit, family = binomial("logit"))

training.str <- glm(benefit ~ activity, data = str.dic.baseline, method = method.fit, family = binomial("logit"))
str.training.str <- glm(benefit ~ any_strength_type, data = str.dic.baseline, method = method.fit, family = binomial("logit"))


training_data_for_table <- list(csa.dic = csa.dic, 
                                str.dic = str.dic, 
                                csa.dic.baseline = csa.dic.baseline,
                                str.dic.baseline = str.dic.baseline,
                                training.csa = training.csa,
                                str.training.csa = str.training.csa,
                                training.str = training.str,
                                str.training.str = str.training.str,
                                total.session.csa = total.session.csa, 
                                supervised.session.csa = supervised.session.csa, 
                                total.session.str = total.session.str,
                                supervised.session.str = supervised.session.str)


## Saves for tables
saveRDS(training_data_for_table, "./data/derivedData/study-1-benefit-analysis/training_data_for_table.Rds")

# No training variables had relationship with benefit #


## Classify variables 

variables.csa <- results %>%
        dplyr::select(variable, estimate.csa, pval.csa, sex.pval.csa, sex.inter.pval.csa) %>%
        mutate(pval.csa = round(pval.csa, 3), 
               sex.pval.csa = round(sex.pval.csa, 3), 
               sex.inter.pval.csa = round(sex.inter.pval.csa, 3)) %>%
        mutate(selection = if_else(pval.csa < 0.2, "include", "exclude")) %>%
        data.frame() %>%
        print()


variables.str <- results %>%
        dplyr::select(variable, estimate.str, pval.str, sex.pval.str, sex.inter.pval.str) %>%
        mutate(pval.str = round(pval.str, 3), 
               sex.pval.str = round(sex.pval.str, 3), 
               sex.inter.pval.str = round(sex.inter.pval.str, 3)) %>%
        mutate(selection = if_else(pval.str < 0.2, "include", "exclude")) %>%
        print()



#### Check continuous variables for linearity in the logit ####


inluded.csa <- variables.csa[variables.csa$selection == "include", 1]

csa.df <- benefit_csa_complete_cont %>%
        dplyr::select(subject, sex, benefit, one_of(inluded.csa)) %>%
        filter(complete.cases(.)) %>%
        data.frame() %>%
        print()


inluded.csa_testo <- inluded.csa[-2]


lin_check_results_csa <- list()

for(i in 1:length(inluded.csa_testo)) {
        
        p1 <- linearity_check(csa.df, inluded.csa_testo[i])
        
        
        lin_check_results_csa[[i]] <- p1
        
}

# with(csa.df, logitloess(type2x,benefit, cf = sex))



### Strength 

included.str <- variables.str[variables.str$selection == "include", 1]

str.df <- benefit_strength_complete_cont %>%
        dplyr::select(subject, sex, benefit, one_of(included.str)) %>%
        filter(complete.cases(.)) %>%
        print()



lin_check_results_str <- list()

for(i in 1:length(included.str)) {
        
        p1 <- linearity_check(str.df, included.str[i])
        
        
        lin_check_results_str[[i]] <- p1
        
}





# with(benefit_csa_complete_cont, logitloess(lean.mass, benefit, cf = sex))



##### Create dichotomization/factorization of variables #####


# Selection and categorization of variables are based training-induced differences 
# between sets-conditions (mTOR activation, ribosome biogenesis) and 
# baseline characteristics (muscle mass, fiber-type composition, strength and blood parameters). 
# Blood variables with training-induced changes (IGF-1 and GH) are analyzed as mean baseline- or
# training induced values. Blood variables with no apperent training induction (testosterone and cortisol)
# are analyzed as mean baseline values. Vit-D is used as the mean of pre and post values.

# All variables were tested for linearity in the logit using design variables (data divided in quartiles). 
# If a variable was not linear in the logit one of X alternatives were used (1) biological relevant thresholds were used for
# dichotomization (e.g. Vit-D insufficiency) (2) smallest wortwhile change based on 
# baseline variation (volume dependent, training induced variables e.g. mTOR), (3) dichotomization on the median and/or detection limit 
# (e.g. testosterone, median for males and detection-limit for females). Bimodal variables (e.g. % lean mass) are mean-centered 
# per sex. Ref on centering: 
# https://afni.nimh.nih.gov/pub/dist/doc/htmldoc/STATISTICS/center.html


# Vit D at baseline was not linear in the logit and therefore factorized.
# Average Vit d were factorized to Normal values "a_norm" 75-100 nmol/L (30-40 ng/ml). 
# or "b_insufdefic" levels less than 75 nmol/L (20 ng/ml), including insufficiency and deficiency




lin_check_results_csa[[2]]

csa_data <- csa.df %>%
        mutate(type2x = if_else(type2x > median(type2x), "hi", "lo")) %>%
        group_by(sex) %>%
        mutate(lean.mass = if_else(lean.mass > median(lean.mass), "a_high", "b_low"),
               testo.ba = if_else(sex == "female", if_else(testo.ba > 0, 1, 0), 
                                  if_else(testo.ba > median(testo.ba), 1, 0))) %>%
        ungroup() %>%
        print()


lin_check_results_str[[2]]

str_data <- str.df %>%
        group_by(sex) %>%
        mutate(cortisol.ba = cortisol.ba - mean(cortisol.ba)) %>%
        ungroup() %>%
        print()



##### CSA models #######

library(brglm2)
method.fit <- "brglmFit"



form.csa <- formula(paste0("benefit ~ sex + ", 
                           paste0(inluded.csa, collapse = "+"), 
                           collapse = "+" ))



m1 <- glm(form.csa, data = csa_data, 
          family = binomial(link = "logit"), 
          method = method.fit)


summary(m1)
car::vif(m1)


form.csa <- formula(paste0("benefit ~ sex + ", 
                           paste0(inluded.csa[-c(5)], collapse = "+"), 
                           collapse = "+" ))

# ... and fitted
m2 <- glm(form.csa, data = csa_data, 
          family = binomial(link = "logit"), 
          method = method.fit)


summary(m2)
delta_beta_fun(m1, m2)
car::vif(m2)



form.csa <- formula(paste0("benefit ~ sex + ", 
                           paste0(inluded.csa[-c(5, 3)], collapse = "+"), 
                           collapse = "+" ))

# ... and fitted
m3 <- glm(form.csa, data = csa_data, 
          family = binomial(link = "logit"), 
          method = method.fit)


summary(m3)
delta_beta_fun(m2, m3)
car::vif(m3)



form.csa <- formula(paste0("benefit ~ sex + ", 
                           paste0(inluded.csa[-c(2, 3, 5)], collapse = "+"), 
                           collapse = "+" ))

# ... and fitted
m4 <- glm(form.csa, data = csa_data, 
          family = binomial(link = "logit"), 
          method = method.fit)


summary(m4)
delta_beta_fun(m3, m4)
car::vif(m4)



form.csa <- formula(paste0("benefit ~ sex + ", 
                           paste0(inluded.csa[-c(2, 3, 4, 5)], collapse = "+"), 
                           collapse = "+" ))

# ... and fitted
m5 <- glm(form.csa, data = csa_data, 
          family = binomial(link = "logit"), 
          method = method.fit)


summary(m5)
delta_beta_fun(m4, m5)
car::vif(m5)
anova(m4, m5, test = "LRT")



### Excluded analysis CSA ####


excluded.csa <- variables.csa[variables.csa$selection == "exclude", 1]

excluded.results <- list()

for(i in 1:length(excluded.csa)) {
        
        benefit.csa <- benefit_csa_complete_cont %>%
                dplyr::select(one_of(c("subject","benefit", "sex",  "rna_w2pre" , excluded.csa[i]))) %>%
                filter(complete.cases(.)) 
        
        
        form <- formula(paste0(paste0("benefit ~ sex +  rna_w2pre", collapse = "+" ),"+", excluded.csa[i]))
        
        
        mx <- glm(form, data = benefit.csa, 
                  family = binomial(link = "logit"), 
                  method = method.fit)
        
        n <- length(coef(mx))
        
        
        results <- data.frame(variable = rep(NA, length(coef(mx)) - 3),
                              coef = rep(NA, length(coef(mx)) - 3),
                              estimate = rep(NA, length(coef(mx)) - 3),
                              ci.lower = rep(NA, length(coef(mx)) - 3),
                              ci.upper = rep(NA, length(coef(mx)) - 3),
                              deviance = rep(NA, length(coef(mx)) - 3),
                              p.value = rep(NA,  length(coef(mx)) - 3), 
                              p.value.wald = rep(NA, length(coef(mx)) - 3))
        
        results[c(1:(n-3)), 1] <- excluded.csa[i]
        results[c(1:(n-3)), 2] <- names(coef(mx))[c(4:n)]
        results[c(1:(n-3)), 3] <- round(summary(mx)$coef[c(4:n),1],3)
        results[c(1:(n-3)), 4] <- confint(mx)[c(4:n),1]
        results[c(1:(n-3)), 5] <- confint(mx)[c(4:n),2]
        results[c(1:(n-3)), 6] <- data.frame(anova(mx, test = "Chisq"))[4,2]
        results[c(1:(n-3)), 7] <- data.frame(anova(mx, test = "Chisq"))[4,5]
        results[c(1:(n-3)), 8] <- summary(mx)$coef[c(4:n),4]
        
        excluded.results[[i]] <- results 
        
}


excluded.analysis <- bind_rows(excluded.results)


excl.csa_data <- csa_data %>%
        ungroup() %>%
        inner_join(benefit_csa_complete_cont %>%  
                           dplyr::select(subject, "igf1.ba"), by = "subject") %>%
        filter(complete.cases(.)) %>%
        mutate(igf1.ba = igf1.ba) %>%
        print()

linearity_check(excl.csa_data, "igf1.ba")

# The only variable with significant relation from excluded analysis is not linear in the logit.
# Transforming to tertier factors


excl.csa_data <- csa_data %>%
        ungroup() %>%
        inner_join(benefit_csa_complete_cont %>%  
                           dplyr::select(subject, "igf1.ba"), by = "subject") %>%
        filter(complete.cases(.)) %>%
        mutate(igf1.ba = cut(igf1.ba, quantile(igf1.ba, 
                                               na.rm = TRUE, 
                                               prob = seq(0, 1, by = 0.33), type = 7), include.lowest = TRUE, 
                             labels = c("f1", "f2", "f3")),
               igf1.ba = factor(igf1.ba, levels = c("f3", "f1", "f2"))) %>%
        print()


form.csa <- formula("benefit ~ sex + rna_w2pre + lean.mass + igf1.ba")


# ... and fitted
mx <- glm(form.csa, data = excl.csa_data, 
          family = binomial(link = "logit"), 
          method = method.fit)


summary(mx)



# Transforming the IGF1 variable, no evidence or significant relationship exists... 
# It is thus concluded that variable selection seleceted w2pre and lean.mass (at p < 0.1). Both seems relevant and 
# should be further discussed/examined.

##### Check model assumptions CSA ###############


## Influential data points ##

library(car)


full <- glm(benefit ~ sex + rna_w2pre + lean.mass, data = csa_data, 
            family = binomial(link = "logit"), 
            method = method.fit)

influencePlot(full)



## Identified as influential
reduced7 <- update(full, subset = c(-7))
reduced9 <- update(full, subset = c(-9))
reduced19 <- update(full, subset = c(-19))
reduced25 <- update(full, subset = c(-25))

## Not identified as influetial
reduced21 <- update(full, subset = c(-21))
reduced22 <- update(full, subset = c(-22))


compareCoefs(full, reduced7, reduced9,  reduced19, reduced25, reduced21, reduced22)

# Influential data points are identified, but do not change the main conclusion (increased OR with increased w2pre).

full <- glm(benefit ~ sex + rna_w2pre + lean.mass, data = csa_data, 
            family = binomial(link = "logit"), 
            method = method.fit)

reduced <- glm(benefit ~ rna_w2pre + lean.mass, data = csa_data, 
               family = binomial(link = "logit"), 
               method = method.fit)

null <- glm(benefit ~  sex, data = csa_data, 
            family = binomial(link = "logit"), 
            method = method.fit)

anova(reduced, full,null, test = "LRT")

# Comparing the final model with and without sex does not improve the fit. Sex not important for the model but opotentially important
# in model fitting. 


##### Strength models #######




form.str <- formula(paste0("benefit ~ sex + ", 
                           paste0(included.str, collapse = "+"), 
                           collapse = "+" ))



m1.str <- glm(form.str, data = str_data, 
              family = binomial(link = "logit"), 
              method = method.fit)


summary(m1.str)
car::vif(m1.str)


form.str <- formula(paste0("benefit ~ sex + ", 
                           paste0(included.str[-c(3)], collapse = "+"), 
                           collapse = "+" ))

# ... and fitted
m2.str <- glm(form.str, data = str_data, 
              family = binomial(link = "logit"), 
              method = method.fit)


summary(m2.str)
delta_beta_fun(m1.str, m2.str)
car::vif(m2.str)



form.str <- formula(paste0("benefit ~ sex + ", 
                           paste0(included.str[-c(2, 3)], collapse = "+"), 
                           collapse = "+" ))

# ... and fitted
m3.str <- glm(form.str, data = str_data, 
              family = binomial(link = "logit"), 
              method = method.fit)


summary(m3.str)
delta_beta_fun(m2.str, m3.str)
car::vif(m3.str)
anova(m2.str, m3.str, test = "LRT")



#### Excluded analysis strength #####

### Excluded analysis CSA ####


excluded.str <- variables.str[variables.str$selection == "exclude", 1]

excluded.results <- list()

for(i in 1:length(excluded.str)) {
        
        benefit.str <- benefit_strength_complete_cont %>%
                dplyr::select(one_of(c("subject","benefit", "sex",  "rna_w2pre" , excluded.str[i]))) %>%
                filter(complete.cases(.)) 
        
        
        form <- formula(paste0(paste0("benefit ~ sex +  rna_w2pre", collapse = "+" ),"+", excluded.str[i]))
        
        
        mx <- glm(form, data = benefit.str, 
                  family = binomial(link = "logit"), 
                  method = method.fit)
        
        n <- length(coef(mx))
        
        
        results <- data.frame(variable = rep(NA, length(coef(mx)) - 3),
                              coef = rep(NA, length(coef(mx)) - 3),
                              estimate = rep(NA, length(coef(mx)) - 3),
                              ci.lower = rep(NA, length(coef(mx)) - 3),
                              ci.upper = rep(NA, length(coef(mx)) - 3),
                              deviance = rep(NA, length(coef(mx)) - 3),
                              p.value = rep(NA,  length(coef(mx)) - 3), 
                              p.value.wald = rep(NA, length(coef(mx)) - 3))
        
        results[c(1:(n-3)), 1] <- excluded.str[i]
        results[c(1:(n-3)), 2] <- names(coef(mx))[c(4:n)]
        results[c(1:(n-3)), 3] <- round(summary(mx)$coef[c(4:n),1],3)
        results[c(1:(n-3)), 4] <- confint(mx)[c(4:n),1]
        results[c(1:(n-3)), 5] <- confint(mx)[c(4:n),2]
        results[c(1:(n-3)), 6] <- data.frame(anova(mx, test = "Chisq"))[4,2]
        results[c(1:(n-3)), 7] <- data.frame(anova(mx, test = "Chisq"))[4,5]
        results[c(1:(n-3)), 8] <- summary(mx)$coef[c(4:n),4]
        
        excluded.results[[i]] <- results 
        
}


excluded.analysis.str <- bind_rows(excluded.results)

excluded_reanalysis <- excluded.analysis.str[excluded.analysis.str$p.value.wald < 0.2, 1]


### Check linearity
benefit_strength_complete_cont %>%
        dplyr::select(one_of(c("subject","benefit", "sex",  "rna_w2pre" , excluded_reanalysis[3]))) %>%
        filter(complete.cases(.)) %>%
        linearity_check(excluded_reanalysis[3])



excl.str_data <- str_data %>%
        ungroup() %>%
        inner_join(benefit_strength_complete_cont %>%  
                           dplyr::select(subject, "testo.ba", "rel.strength", "lean.mass"), by = "subject") %>%
        mutate(lean.mass = if_else(lean.mass > median(lean.mass), "a_high", "b_low")) %>%
        group_by(sex) %>%
        mutate(testo.ba = as.factor(if_else(sex == "female", if_else(testo.ba > 0, 1, 0),
                                            if_else(testo.ba > median(testo.ba), 1, 0)))) %>%
        ungroup() %>%
        filter(complete.cases(.)) %>%
        print()


form.str <- formula("benefit ~ sex + rna_w2pre  + lean.mass  + rel.strength + testo.ba")


mx.str <- glm(form.str, data = excl.str_data, 
              family = binomial(link = "logit"), 
              method = method.fit)


summary(mx.str)


form.str <- formula("benefit ~ sex + rna_w2pre +  rel.strength + lean.mass")

mx1.str <- glm(form.str, data = excl.str_data, 
               family = binomial(link = "logit"), 
               method = method.fit)


summary(mx1.str)
delta_beta_fun(mx.str, mx1.str)


form.str <- formula("benefit ~ sex + rna_w2pre +  rel.strength")

mx2.str <- glm(form.str, data = excl.str_data, 
               family = binomial(link = "logit"), 
               method = method.fit)


summary(mx2.str)
delta_beta_fun(mx1.str, mx2.str)

form.str <- formula("benefit ~ sex + rna_w2pre")

mx3.str <- glm(form.str, data = excl.str_data, 
               family = binomial(link = "logit"), 
               method = method.fit)


summary(mx3.str)
delta_beta_fun(mx2.str, mx3.str)

##### Check model assumptions CSA ###############


## Influential data points ##

library(car)


full.str <- glm(benefit ~ sex + rna_w2pre, data = str_data, 
                family = binomial(link = "logit"), 
                method = method.fit)

influencePlot(full.str)



## Identified as influential
reduced7 <- update(full.str, subset = c(-7))
reduced13 <- update(full.str, subset = c(-13))
reduced14 <- update(full.str, subset = c(-14))


## Not identified as influetial
reduced21 <- update(full.str, subset = c(-21))
reduced22 <- update(full.str, subset = c(-22))


compareCoefs(full.str, reduced7, reduced13,  reduced14, reduced21, reduced22)

# Influential data points are identified, but do not change the main conclusion (increased OR with increased w2pre).

full.str <- glm(benefit ~ sex + rna_w2pre, data = str_data, 
                family = binomial(link = "logit"), 
                method = method.fit)

reduced.str <- glm(benefit ~ rna_w2pre, data = str_data, 
                   family = binomial(link = "logit"), 
                   method = method.fit)



null.str <- glm(benefit ~  sex, data = str_data, 
                family = binomial(link = "logit"), 
                method = method.fit)

anova(reduced.str, full.str, null.str, test = "LRT")



# Conclusion on strength models #
# Stepwise removal of predictors selected in first and second screening did not reveal any other predictors than
# rna_w2pre. This predictor is considered further.


#### Saving models for table ######


confint(m1)

csa_models <- rbind(cbind(broom::tidy(m1) %>%
                            mutate(model = "m1", 
                                   dependent = "csa", 
                                   mod.comp = NA),data.frame(confint(m1))),
                   cbind(broom::tidy(m2) %>%
                            mutate(model = "m2", 
                                   dependent = "csa", 
                                   mod.comp = data.frame(anova(m1, test = "LRT"))[7, 5]),data.frame(confint(m2))),
                   cbind(broom::tidy(m3) %>%
                            mutate(model = "m3", 
                                   dependent = "csa",
                                   mod.comp = data.frame(anova(m2, m3, test = "LRT"))[2, 5]),data.frame(confint(m3))),
                         cbind(broom::tidy(m4) %>%
                            mutate(model = "m4", 
                                   dependent = "csa",
                                   mod.comp = data.frame(anova(m3, m4, test = "LRT"))[2, 5]),data.frame(confint(m4))),
                            cbind(broom::tidy(m5) %>%
                            mutate(model = "m5", 
                                   dependent = "csa",
                                   mod.comp = data.frame(anova(m4, m5, test = "LRT"))[2, 5]),data.frame(confint(m5))),
                            cbind(broom::tidy(reduced) %>%
                            mutate(model = "reduced", 
                                   dependent = "csa",
                                   mod.comp = data.frame(anova(m4, reduced, test = "LRT"))[2, 5]),data.frame(confint(reduced))))

str_models <- rbind(cbind(broom::tidy(m1.str) %>%
                            mutate(model = "m1", 
                                   dependent = "str",
                                   mod.comp = NA),data.frame(confint(m1.str))),
                          cbind(broom::tidy(m2.str) %>%
                            mutate(model = "m2", 
                                   dependent = "str",
                                   mod.comp = data.frame(anova(m1.str, test = "LRT"))[5,5]),data.frame(confint(m2.str))),
                            cbind(broom::tidy(m3.str) %>%
                            mutate(model = "m3", 
                                   dependent = "str",
                                   mod.comp = data.frame(anova(m2.str, m3.str, test = "LRT"))[2, 5]),data.frame(confint(m3.str))),
                    
                            cbind(broom::tidy(reduced.str) %>%
                            mutate(model = "reduced", 
                                   dependent = "str",
                                   mod.comp = data.frame(anova(m2.str, reduced.str, test = "LRT"))[2, 5]),data.frame(confint(reduced.str))))



# saveRDS(list(reduced.str = reduced.str, reduced.csa = reduced), "./derivedData/logit_models.Rds")

saveRDS(csa_models, "./data/derivedData/study-1-benefit-analysis/csa_logit_models.Rds")
saveRDS(str_models, "./data/derivedData/study-1-benefit-analysis/str_logit_models.Rds")



