##-------------------------------------
## meta-analysis.R
##
## Title: Meta analysis effect of inter-session volume on muscle strength and mass 
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



# References for calculations
#
# 
# https://onlinelibrary.wiley.com/doi/pdf/10.1002/9780470743386.ch4
# http://handbook-5-1.cochrane.org/
# https://mvuorre.github.io/post/2016/09/29/meta-analysis-is-a-special-case-of-bayesian-multilevel-modeling/
#
#


library(readxl)
library(tidyverse)
# library(metafor)
# library(nlme)
# library(MAd)
library(brms); library(rstan)
library(tidybayes)

# settings for rstan
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Read data
dat <- read_excel("./data/meta-volume/raw-data.xlsx")



# Calculate effect sizes
es <- dat %>%
        # calculate pre sd and post m in Wittke 2017
        mutate(pre_sd = if_else(study == "Wittke2017", sqrt(n_participants) * (pre_ciu - pre_cil)/(2 * qt(0.975, 40-1)), pre_sd), 
               post_m = if_else(study == "Wittke2017", pre_m + change_m, post_m), 
               post_sd = if_else(study == "Wittke2017", pre_sd, post_sd)) %>%
        
        # Create variables for modeling 
        
        mutate(age_cat = if_else(age_m < 30, "young", 
                                 if_else(age_m > 50, "old", "middle")), 
               age_cat = if_else(study == "Cannon2010", "mixed", age_cat), 
               age_cat = factor(age_cat, levels = c("young", 
                                                    "middle", 
                                                    "old", 
                                                    "mixed")), 
               training_status = factor(training_status, levels = c("untrained", 
                                                                    "trained")), 
               sex = factor(sex, levels = c("male", 
                                            "female", 
                                            "mixed")), 
               body_half = factor(body_half, levels = c("lower", 
                                                        "upper", 
                                                        "whole"))) %>%
        # other variables:
        # sex, training_status, body_half, type_technique, outcome_group
        
        # Filter studies
        filter(group_type == "INTERVENTION") %>%
        
        filter( !(study %in% c("Radaelli2013", "Radaelli2014"))) %>%
        
        # Calculate weekly sets
        mutate(weekly_sets = n_sessions_week * n_sets * n_exercises) %>%
        data.frame() %>%
        
        group_by(study, outcome) %>%
        mutate(sd = (n_participants - 1) * pre_sd^2, 
               n = sum((n_participants - 1)),
               sd = sqrt(sum(sd) / n), 
               change = post_m - pre_m, 
               es = change / sd) %>%
        mutate(j = (1 - (3/((4*(n_participants-1) - 1)))),
               es.ub = es * j) %>% # unbiased es
        group_by(muscle_size) %>% # groups by strength or muscle size
        mutate(pop.es = mean(es.ub)) %>% # average es per outcome
        ungroup() %>%
        
        # Estimate r from previous studies or by pre-post SDs?
        # This does not works so well, as pre post SD leads to biased r (> 1)
        mutate(sd.diff = sqrt(((pre_sd^2)/(n_participants)) + ((post_sd^2)/(n_participants))),  # for use in correlation estimate
               r = (pre_sd^2 + post_sd^2 - sd.diff^2) / (2*pre_sd*post_sd)) %>%
        # This is a better approach, estimating from previous studies. 
        mutate(r = if_else(outcome_group %in% c("dxa_whole", "dxa_extr"), 0.98,  # Based on Hammarström 2020
                           if_else(outcome_group %in% c("MRI"), 0.97, # Based on Hammarstrom 2020
                                   if_else(outcome_group %in% c("MT"), 0.93,  # Hammarstrom Ofsteng 2020
                                           if_else(outcome_group %in% c("isometric","isokinetic"), 0.85, # Based on Hammarstrom 2020
                                                   if_else(outcome_group %in% c("1RM", "5RM"), 0.85, 0.8)))))) %>% # Based on Hammarstrom 2020, air displacement low estimate
        
        group_by(study,group, outcome) %>%
        mutate(vd = ((2 * (1 - r)) / n_participants) * ((n_participants - 1) /(n_participants -3)) * (1 + (n_participants/(2*(1-r)))*pop.es^2) - ((pop.es^2)/(j^2)),  
               vd2 = ((1/n_participants) + ((pop.es^2) / (2*n_participants))) * 2 * (1-r), # https://onlinelibrary.wiley.com/doi/pdf/10.1002/9780470743386.ch4
               vd2 = vd2 * j^2,
               vd = vd * j^2,
               se = sqrt(vd2),
               
               HI_LO = factor(HI_LO, levels = c("LO", "HI")), 
               muscle_size = if_else(muscle_size == TRUE, "muscle_size", "muscle_strength")) %>%
        print()




### Save data set for descriptive stats ################

saveRDS(es, "./data/meta-volume/models/raw_data_es.RDS")




### brms models muscle mass ##################

# Models:

# m1: no interactions
# m2: age:sets interaction
# m3: training status:sets interaction
# m4: sex:sets interaction
# m5: body half:sets interaction
# m6: type measuremet:sets interaction
# m7: n weeks:sets interaction
# m8: all interactions




m1 <- brm(data = es[es$muscle_size == "muscle_size",], 
          family = gaussian,
          es.ub | se(se) ~ weekly_sets + 
                  age_cat +
                  type_technique +
                  training_status +
                  sex +
                  n_weeks +
                  body_half +
                  (weekly_sets||study),
          prior = c(prior(normal(0, 1), class = Intercept),
                    prior(cauchy(0, 1), class = sd)),
          iter = 8000, warmup = 1000, cores = 4, chains = 4,
          seed = 14)

saveRDS(m1, "./data/meta-volume/models/muscle_size_full_data.RDS")
m1 <- readRDS("./data/meta-volume/models/muscle_size_full_data.RDS")
summary(m1)

pp_check(m1)


m2 <-  brm(data = es[es$muscle_size == "muscle_size",], 
           family = gaussian,
           es.ub | se(se) ~ weekly_sets + 
                   age_cat + age_cat:weekly_sets +
                   type_technique +
                   training_status +
                   sex +
                   n_weeks +
                   body_half +
                   (weekly_sets||study),
           prior = c(prior(normal(0, 1), class = Intercept),
                     prior(cauchy(0, 1), class = sd)),
           iter = 12000, warmup = 2000, cores = 4, chains = 4,
           control = list(adapt_delta = 0.95),
           seed = 14)






saveRDS(m2, "./data/meta-volume/models/muscle_size_full_data_m2.RDS")
m2 <- readRDS("./data/meta-volume/models/muscle_size_full_data_m2.RDS")
pp_check(m2)
summary(m2)


m3 <-  brm(data = es[es$muscle_size == "muscle_size",], 
           family = gaussian,
           es.ub | se(se) ~ weekly_sets + 
                   age_cat + 
                   type_technique +
                   training_status + training_status:weekly_sets +
                   sex +
                   n_weeks +
                   body_half +
                   (weekly_sets||study),
           prior = c(prior(normal(0, 1), class = Intercept),
                     prior(cauchy(0, 1), class = sd)),
           iter = 12000, warmup = 2000, cores = 4, chains = 4,
           control = list(adapt_delta = 0.95),
           seed = 14)



saveRDS(m3, "./data/meta-volume/models/muscle_size_full_data_m3.RDS")
m3 <- readRDS("./data/meta-volume/models/muscle_size_full_data_m3.RDS")
pp_check(m3)
summary(m3)


######### 

m4 <-  brm(data = es[es$muscle_size == "muscle_size",], 
           family = gaussian,
           es.ub | se(se) ~ weekly_sets + 
                   age_cat + 
                   type_technique +
                   training_status +
                   sex +  sex:weekly_sets +
                   n_weeks +
                   body_half +
                   (weekly_sets||study),
           prior = c(prior(normal(0, 1), class = Intercept),
                     prior(cauchy(0, 1), class = sd)),
           iter = 12000, warmup = 2000, cores = 4, chains = 4,
           control = list(adapt_delta = 0.95),
           seed = 14)



saveRDS(m4, "./data/meta-volume/models/muscle_size_full_data_m4.RDS")
m4 <- readRDS("./data/meta-volume/models/muscle_size_full_data_m4.RDS")
pp_check(m4)
summary(m4)



m5 <-  brm(data = es[es$muscle_size == "muscle_size",], 
           family = gaussian,
           es.ub | se(se) ~ weekly_sets + 
                   age_cat + 
                   type_technique +
                   training_status +
                   sex +  
                   n_weeks +
                   body_half + body_half:weekly_sets +
                   (weekly_sets||study),
           prior = c(prior(normal(0, 1), class = Intercept),
                     prior(cauchy(0, 1), class = sd)),
           iter = 12000, warmup = 2000, cores = 4, chains = 4,
           control = list(adapt_delta = 0.95),
           seed = 14)



saveRDS(m5, "./data/meta-volume/models/muscle_size_full_data_m5.RDS")
m5 <- readRDS("./data/meta-volume/models/muscle_size_full_data_m5.RDS")
pp_check(m5)
summary(m5)


m6 <-  brm(data = es[es$muscle_size == "muscle_size",], 
           family = gaussian,
           es.ub | se(se) ~ weekly_sets + 
                   age_cat + 
                   type_technique + type_technique:weekly_sets +
                   training_status +
                   sex +  
                   n_weeks +
                   body_half + 
                   (weekly_sets||study),
           prior = c(prior(normal(0, 1), class = Intercept),
                     prior(cauchy(0, 1), class = sd)),
           iter = 12000, warmup = 2000, cores = 4, chains = 4,
           control = list(adapt_delta = 0.95),
           seed = 14)



saveRDS(m6, "./data/meta-volume/models/muscle_size_full_data_m6.RDS")
m6 <- readRDS("./data/meta-volume/models/muscle_size_full_data_m6.RDS")
pp_check(m6)
summary(m6)


###

m7 <-  brm(data = es[es$muscle_size == "muscle_size",], 
           family = gaussian,
           es.ub | se(se) ~ weekly_sets + 
                   age_cat + 
                   type_technique + 
                   training_status +
                   sex +  
                   n_weeks + n_weeks:weekly_sets +
                   body_half + 
                   (weekly_sets||study),
           prior = c(prior(normal(0, 1), class = Intercept),
                     prior(cauchy(0, 1), class = sd)),
           iter = 12000, warmup = 2000, cores = 4, chains = 4,
           control = list(adapt_delta = 0.95),
           seed = 14)



saveRDS(m7, "./data/meta-volume/models/muscle_size_full_data_m7.RDS")
m7 <- readRDS("./data/meta-volume/models/muscle_size_full_data_m7.RDS")
pp_check(m7)
summary(m7)




#### Full model muscle mass ####################

m8 <-  brm(data = es[es$muscle_size == "muscle_size",], 
           family = gaussian,
           es.ub | se(se) ~ weekly_sets + 
                   age_cat + age_cat:weekly_sets + 
                   type_technique + type_technique:weekly_sets +
                   training_status + training_status:weekly_sets+
                   sex +  sex:weekly_sets+
                   n_weeks + n_weeks:weekly_sets +
                   body_half + body_half:weekly_sets+
                   (weekly_sets||study),
           prior = c(prior(normal(0, 1), class = Intercept),
                     prior(cauchy(0, 1), class = sd)),
           iter = 12000, warmup = 2000, cores = 4, chains = 4,
           control = list(adapt_delta = 0.95),
           seed = 14)



saveRDS(m8, "./data/meta-volume/models/muscle_size_full_data_m8.RDS")
m8 <- readRDS("./data/meta-volume/models/muscle_size_full_data_m8.RDS")
pp_check(m8)
summary(m8)



##### Strength models #############

# Models:

# s1: no interactions
# s2: age:sets interaction
# s3: training status:sets interaction
# s4: sex:sets interaction
# s5: body half:sets interaction
# s6: type measuremet:sets interaction
# s7: n weeks:sets interaction
# s8: all interactions




s1 <- brm(data = es[es$muscle_size == "muscle_strength",], 
          family = gaussian,
          es.ub | se(se) ~ weekly_sets + 
                  age_cat +
                  type_technique +
                  training_status +
                  sex +
                  n_weeks +
                  body_half +
                  (weekly_sets||study),
          prior = c(prior(normal(0, 1), class = Intercept),
                    prior(cauchy(0, 1), class = sd)),
          iter = 8000, warmup = 1000, cores = 4, chains = 4,
          seed = 14)

saveRDS(s1, "./data/meta-volume/models/strength_full_data.RDS")
s1 <- readRDS("./data/meta-volume/models/strength_full_data.RDS")
pp_check(s1)
summary(s1)


s2 <-  brm(data = es[es$muscle_size == "muscle_strength",], 
           family = gaussian,
           es.ub | se(se) ~ weekly_sets + 
                   age_cat + age_cat:weekly_sets +
                   type_technique +
                   training_status +
                   sex +
                   n_weeks +
                   body_half +
                   (weekly_sets||study),
           prior = c(prior(normal(0, 1), class = Intercept),
                     prior(cauchy(0, 1), class = sd)),
           iter = 8000, warmup = 2000, cores = 4, chains = 4,
           control = list(adapt_delta = 0.9),
           seed = 14)

saveRDS(s2, "./data/meta-volume/models/strength_full_data_s2.RDS")
s2 <- readRDS("./data/meta-volume/models/strength_full_data_s2.RDS")
pp_check(s2)
summary(s2)





s3 <-  brm(data = es[es$muscle_size == "muscle_strength",], 
           family = gaussian,
           es.ub | se(se) ~ weekly_sets + 
                   age_cat + 
                   type_technique +
                   training_status + training_status:weekly_sets +
                   sex +
                   n_weeks +
                   body_half +
                   (weekly_sets||study),
           prior = c(prior(normal(0, 1), class = Intercept),
                     prior(cauchy(0, 1), class = sd)),
           iter = 8000, warmup = 2000, cores = 4, chains = 4,
           control = list(adapt_delta = 0.9),
           seed = 14)

saveRDS(s3, "./data/meta-volume/models/strength_full_data_s3.RDS")
s3 <- readRDS("./data/meta-volume/models/strength_full_data_s3.RDS")
pp_check(s3)
summary(s3)



s4 <-  brm(data = es[es$muscle_size == "muscle_strength",], 
           family = gaussian,
           es.ub | se(se) ~ weekly_sets + 
                   age_cat + 
                   type_technique +
                   training_status +
                   sex +  sex:weekly_sets +
                   n_weeks +
                   body_half +
                   (weekly_sets||study),
           prior = c(prior(normal(0, 1), class = Intercept),
                     prior(cauchy(0, 1), class = sd)),
           iter = 8000, warmup = 2000, cores = 4, chains = 4,
           control = list(adapt_delta = 0.9),
           seed = 14)

saveRDS(s4, "./data/meta-volume/models/strength_full_data_s4.RDS")
s4 <- readRDS("./data/meta-volume/models/strength_full_data_s4.RDS")
pp_check(s4)
summary(s4)


s5 <-  brm(data = es[es$muscle_size == "muscle_strength",], 
           family = gaussian,
           es.ub | se(se) ~ weekly_sets + 
                   age_cat + 
                   type_technique +
                   training_status +
                   sex +  
                   n_weeks +
                   body_half + body_half:weekly_sets +
                   (weekly_sets||study),
           prior = c(prior(normal(0, 1), class = Intercept),
                     prior(cauchy(0, 1), class = sd)),
           iter = 8000, warmup = 2000, cores = 4, chains = 4,
           control = list(adapt_delta = 0.9),
           seed = 14)

saveRDS(s5, "./data/meta-volume/models/strength_full_data_s5.RDS")
s5 <- readRDS("./data/meta-volume/models/strength_full_data_s5.RDS")
pp_check(s5)
summary(s5)



s6 <-  brm(data = es[es$muscle_size == "muscle_strength",], 
           family = gaussian,
           es.ub | se(se) ~ weekly_sets + 
                   age_cat + 
                   type_technique + type_technique:weekly_sets +
                   training_status +
                   sex +  
                   n_weeks +
                   body_half +
                   (type_technique * weekly_sets||study),
           prior = c(prior(normal(0, 1), class = Intercept),
                     prior(cauchy(0, 1), class = sd)),
           iter = 8000, warmup = 2000, cores = 4, chains = 4,
           control = list(adapt_delta = 0.9),
           seed = 14)

saveRDS(s6, "./data/meta-volume/models/strength_full_data_s6.RDS")
s6 <- readRDS("./data/meta-volume/models/strength_full_data_s6.RDS")
pp_check(s6)
summary(s6)



s7 <-  brm(data = es[es$muscle_size == "muscle_strength",], 
           family = gaussian,
           es.ub | se(se) ~ weekly_sets + 
                   age_cat + 
                   type_technique + weekly_sets +
                   training_status +
                   sex +  
                   n_weeks +  n_weeks:weekly_sets +
                   body_half +
                   (type_technique * weekly_sets||study),
           prior = c(prior(normal(0, 1), class = Intercept),
                     prior(cauchy(0, 1), class = sd)),
           iter = 8000, warmup = 2000, cores = 4, chains = 4,
           control = list(adapt_delta = 0.9),
           seed = 14)

saveRDS(s7, "./data/meta-volume/models/strength_full_data_s7.RDS")
s7 <- readRDS("./data/meta-volume/models/strength_full_data_s7.RDS")
pp_check(s7)
summary(s7)




## Full model strength 

s8 <-  brm(data = es[es$muscle_size == "muscle_strength",], 
           family = gaussian,
           es.ub | se(se) ~ weekly_sets + 
                   age_cat + age_cat:weekly_sets + 
                   type_technique + type_technique:weekly_sets +
                   training_status + training_status:weekly_sets+
                   sex +  sex:weekly_sets+
                   n_weeks + n_weeks:weekly_sets +
                   body_half + body_half:weekly_sets+
                   (weekly_sets||study),
           prior = c(prior(normal(0, 1), class = Intercept),
                     prior(cauchy(0, 1), class = sd)),
           iter = 12000, warmup = 2000, cores = 4, chains = 4,
           control = list(adapt_delta = 0.9),
           seed = 14)

saveRDS(s8, "./data/meta-volume/models/strength_full_data_s8.RDS")
s8 <- readRDS("./data/meta-volume/models/strength_full_data_s8.RDS")
pp_check(s8)
summary(s8)


########################### Combine estimates from models


# m1: no interactions
# m2: age:sets interaction
# m3: training status:sets interaction
# m4: sex:sets interaction
# m5: body half:sets interaction
# m6: type measuremet:sets interaction
# m7: n weeks:sets interaction
# m8: all interactions


m1 <- readRDS("./data/meta-volume/models/muscle_size_full_data.RDS")
m2 <- readRDS("./data/meta-volume/models/muscle_size_full_data_m2.RDS")
m3 <- readRDS("./data/meta-volume/models/muscle_size_full_data_m3.RDS")
m4 <- readRDS("./data/meta-volume/models/muscle_size_full_data_m4.RDS")
m5 <- readRDS("./data/meta-volume/models/muscle_size_full_data_m5.RDS")
m6 <- readRDS("./data/meta-volume/models/muscle_size_full_data_m6.RDS")
m7 <- readRDS("./data/meta-volume/models/muscle_size_full_data_m7.RDS")


# Strength models

s2 <- readRDS("./data/meta-volume/models/strength_full_data_s2.RDS")
s3 <- readRDS("./data/meta-volume/models/strength_full_data_s3.RDS")
s4 <- readRDS("./data/meta-volume/models/strength_full_data_s4.RDS")
s5 <- readRDS("./data/meta-volume/models/strength_full_data_s5.RDS")
s6 <- readRDS("./data/meta-volume/models/strength_full_data_s6.RDS")
s7 <- readRDS("./data/meta-volume/models/strength_full_data_s7.RDS")



# Extract interactions for table


muscle_interactions <- cbind(data.frame(variable = c("Age: Middle", 
                              "Age: Older", 
                              "Age: Mixed", 
                              "Training status: Trained",
                              "Sex: Female", 
                              "Sex: Mixed", 
                              "Body portion: Upper-body", 
                              "Body portion: Whole body", 
                              "Measuremet technique: Indirect/Unspecific", 
                              "Study length (Weeks)")),
        rbind(summary(m2)$fixed[c(13:15),], # Age category
        summary(m3)$fixed[13,], # Training status
        summary(m4)$fixed[c(13:14),], # sex
        summary(m5)$fixed[c(13:14),], # body half
        summary(m6)$fixed[c(13),], # measuremet type
        summary(m7)$fixed[13, ]) %>% # n weeks
        data.frame() %>%
        dplyr::select(m.est = Estimate, 
                      m.err = Est.Error, 
                      m.lwr = l.95..CI,
                      m.upr = u.95..CI)) 
      
rownames(muscle_interactions ) <- NULL
        

strength_interactions <- cbind(data.frame(variable = c("Age: Middle", 
                                                     "Age: Older", 
                                                     "Age: Mixed", 
                                                     "Training status: Trained",
                                                     "Sex: Female", 
                                                     "Sex: Mixed", 
                                                     "Body portion: Upper-body", 
                                                     "Measuremet technique: Indirect/Unspecific", 
                                                     "Study length (Weeks)")),
                             rbind(summary(s2)$fixed[c(12:14),], # Age category
                                   summary(s3)$fixed[12,], # Training status
                                   summary(s4)$fixed[c(12:13),], # sex
                                   summary(s5)$fixed[c(12),], # body half
                                 
                                   summary(s6)$fixed[c(12),], # measuremet type
                                   summary(s7)$fixed[12, ]) %>% # n weeks
                                     data.frame() %>%
                                     dplyr::select(s.est = Estimate, 
                                                   s.err = Est.Error, 
                                                   s.lwr = l.95..CI,
                                                   s.upr = u.95..CI))  %>%
  rbind(data.frame(variable = "Body portion: Whole body",
                   s.est = NA,
                     s.err = NA,
                     s.lwr = NA,
                     s.upr = NA)) %>%
  print()

rownames(strength_interactions ) <- NULL



  

summary_table_interactions <- muscle_interactions %>%
  inner_join(strength_interactions)  %>%
  
  mutate(m.est = round(m.est, digits = 3), 
         m.lwr = round(m.lwr, digits = 3), 
         m.upr = round(m.upr, digits = 3), 
         
         s.est = round(s.est, digits = 3), 
         s.lwr = round(s.lwr, digits = 3), 
         s.upr = round(s.upr, digits = 3)) %>%
  dplyr::select(variable, m.est, m.lwr, m.upr, s.est, s.lwr, s.upr) %>%
  
  print()
  
  
saveRDS(summary_table_interactions, 
        "./data/meta-volume/models/summary_table_interactions.RDS")




### Calculate effects from interaction analyses


meta_hypotheses <- list(strength_upper_body_weeklysets = hypothesis(s5, "weekly_sets + weekly_sets:body_halfupper = 0"), 
                        muscle_upper_body_weeklysets = hypothesis(m5, "weekly_sets + weekly_sets:body_halfupper = 0"), 
                        strength_unspecific = hypothesis(s6, "weekly_sets + weekly_sets:type_techniqueunspecific = 0"))
                        

saveRDS(meta_hypotheses, 
        "./data/meta-volume/models/meta_hypotheses.RDS")




### Combined models effect of strength vs muscle_size measure and weekly sets


c1 <-  brm(data = es, 
           family = gaussian,
           es.ub | se(se) ~ weekly_sets + 
                   muscle_size + muscle_size:weekly_sets +
                   age_cat + 
                   training_status +
                   sex +  
                   n_weeks +
                   body_half +
                   (muscle_size * weekly_sets||study),
           prior = c(prior(normal(0, 1), class = Intercept),
                     prior(cauchy(0, 1), class = sd)),
           iter = 8000, warmup = 2000, cores = 4, chains = 4,
           control = list(adapt_delta = 0.9),
           seed = 14)

saveRDS(c1, "./data/meta-volume/models/combined_full_data_c1.RDS")
c1 <- readRDS("./data/meta-volume/models/combined_full_data_c1.RDS")
pp_check(c1)
summary(c1)







###### Explorative plots
es %>%
        ggplot(aes(weekly_sets, es.ub, color = study)) + geom_point() + 
        facet_wrap(~muscle_size)


es %>%
        ggplot(aes(log(pre_m), log(pre_sd^2), color = study)) + geom_point() + 
        facet_wrap( ~ muscle_size)


