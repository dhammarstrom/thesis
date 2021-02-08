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
library(metafor)
library(nlme)
library(MAd)
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
               
               HI_LO = factor(HI_LO, levels = c("LO", "HI"))) %>%
        print()



es %>%
        filter(study == "Hammarström2020") %>%
        dplyr::select(study, group, r, es,es.ub, vd, vd2, se)

es %>%
        ggplot(aes(r, vd2, color = study)) + geom_point()


### brms model 1 ##################


m1 <- brm(data = es[es$muscle_size == "TRUE",], 
          family = gaussian,
          es.ub | se(se) ~ weekly_sets + (weekly_sets||study),
          prior = c(prior(normal(0, 1), class = Intercept),
                    prior(cauchy(0, 1), class = sd)),
          iter = 8000, warmup = 1000, cores = 4, chains = 4,
          seed = 14)

saveRDS(m1, "./data/meta-volume/models/muscle_size_full_data.RDS")
m1 <- readRDS("./data/meta-volume/models/muscle_size_full_data.RDS")


pp_check(m1)


m2 <-  brm(data = es[es$muscle_size == "TRUE",], 
           family = gaussian,
           es.ub | se(se) ~ weekly_sets + 
                   age_cat + age_cat:weekly_sets +
                   training_status + training_status:weekly_sets +
                   sex + sex:weekly_sets +
                   body_half + body_half:weekly_sets +
                   (weekly_sets||study),
           prior = c(prior(normal(0, 1), class = Intercept),
                     prior(cauchy(0, 1), class = sd)),
           iter = 12000, warmup = 2000, cores = 4, chains = 4,
           control = list(adapt_delta = 0.95),
           seed = 14)

saveRDS(m2, "./data/derivedData/brms_models/muscle_size_full_data_m2.RDS")
m2 <- readRDS("./data/derivedData/brms_models/muscle_size_full_data_m2.RDS")

summary(m2)

## Reduce the model keeping body half ############

m3 <-  brm(data = es[es$muscle_size == "TRUE",], 
           family = gaussian,
           es.ub | se(se) ~ weekly_sets + 
                   body_half + body_half:weekly_sets +
                   (weekly_sets||study),
           prior = c(prior(normal(0, 1), class = Intercept),
                     prior(cauchy(0, 1), class = sd)),
           iter = 8000, warmup = 2000, cores = 4, chains = 4,
           control = list(adapt_delta = 0.9),
           seed = 14)

saveRDS(m3, "./data/derivedData/brms_models/muscle_size_full_data_m3.RDS")
m3 <- readRDS("./data/derivedData/brms_models/muscle_size_full_data_m3.RDS")

summary(m3)

##### Strength models #############
# Start with a weekly_sets only model


s1 <- brm(data = es[es$muscle_size == "FALSE",], 
          family = gaussian,
          es.ub | se(se) ~ weekly_sets + (weekly_sets||study),
          prior = c(prior(normal(0, 1), class = Intercept),
                    prior(cauchy(0, 1), class = sd)),
          iter = 8000, warmup = 1000, cores = 4, chains = 4,
          seed = 14)

saveRDS(s1, "./data/meta-volume/models/strength_full_data.RDS")
s1 <- readRDS("./data/meta-volume/models/strength_full_data.RDS")


s2 <-  brm(data = es[es$muscle_size == "FALSE",], 
           family = gaussian,
           es.ub | se(se) ~ weekly_sets + 
                   age_cat + age_cat:weekly_sets +
                   training_status + training_status:weekly_sets +
                   sex + sex:weekly_sets +
                   body_half + body_half:weekly_sets +
                   (weekly_sets||study),
           prior = c(prior(normal(0, 1), class = Intercept),
                     prior(cauchy(0, 1), class = sd)),
           iter = 8000, warmup = 2000, cores = 4, chains = 4,
           control = list(adapt_delta = 0.9),
           seed = 14)

saveRDS(s2, "./data/derivedData/brms_models/strength_full_data_s2.RDS")
s2 <- readRDS("./data/derivedData/brms_models/strength_full_data_s2.RDS")

summary(s2)



















es %>%
        filter(study == "Rønnestad2007") %>%
        dplyr::select(study, es, es.ub, vd) %>%
        print()



###### Explorative plots
es %>%
        ggplot(aes(weekly_sets, es.ub, color = study)) + geom_point() + 
        facet_wrap(~muscle_size)


es %>%
        ggplot(aes(log(pre_m), log(pre_sd^2), color = study)) + geom_point() + 
        facet_wrap( ~ muscle_size)



### Average effect sizes per study/group in Hi_LO comparison

es.agg <- es %>%
        mutate(group_id = paste0(study, ":", HI_LO, ":", muscle_size)) %>%
        agg(id = group_id, es = es.ub, var = vd, cor = 0, data = .) %>%
        separate(id, into = c("study", "group", "muscle_size"), sep = ":") %>%
        mutate(group = factor(group, levels = c("LO", "HI")), 
               se.d = sqrt(var)) %>%
        print()

es %>%
        filter(study == "Martorelli2017") %>%
        dplyr::select(study,group, es.ub, es, vd) %>%
        print()

es.agg %>%
        filter(study == "Hammarström2020") %>%
        print()



## Model with data aggregated
m1 <- brm(data = es.agg[es.agg$muscle_size == "TRUE",], 
          family = gaussian,
          es | se(se.d) ~ group + (group || study),
          prior = c(prior(normal(0, 1), class = Intercept),
                    prior(cauchy(0, 1), class = sd)),
          iter = 8000, warmup = 2000, cores = 4, chains = 4,
          seed = 14)

summary(m1)


saveRDS(m1, "./data/derivedData/brms_models/muscle_size_aggregated_data.RDS")
m1 <- readRDS("./data/derivedData/brms_models/muscle_size_aggregated_data.RDS")

summary(m1)

get_variables(m2)

m1 %>%
        spread_draws(b_groupHI, r_study[study,groupHI]) %>%
        filter(groupHI != "Intercept") %>%
        # add the grand mean to the group-specific deviations
        mutate(mu = b_groupHI + r_study) %>%
        ungroup() %>%
        mutate(study = str_replace_all(study, "[.]", " ")) %>% 
        
        # plot
        ggplot(aes(x = mu, y = reorder(study, mu))) +
        geom_vline(xintercept = fixef(m1)[2, 1], color = "white", size = 1) +
        geom_vline(xintercept = fixef(m1)[2, 3:4], color = "white", linetype = 2) +
        stat_halfeye(.width = .95, size = 1) +
        labs(x = expression(italic("Cohen's d")),
             y = NULL) +
        theme(panel.grid   = element_blank(),
              axis.ticks.y = element_blank(),
              axis.text.y  = element_text(hjust = 0))






## Model without aggregated data 

m2 <- brm(data = es[es$muscle_size == "TRUE",], 
          family = gaussian,
          es.ub | se(se) ~ HI_LO + (HI_LO || study),
          prior = c(prior(normal(0, 1), class = Intercept),
                    prior(cauchy(0, 1), class = sd)),
          iter = 8000, warmup = 2000, cores = 4, chains = 4,
          seed = 14)

saveRDS(m2, "./data/derivedData/brms_models/muscle_size_full_data.RDS")
m2 <- readRDS("./data/derivedData/brms_models/muscle_size_full_data.RDS")

summary(m2)


m2 %>%
        spread_draws(b_HI_LOHI, r_study[study,HI_LOHI]) %>%
        filter(HI_LOHI != "Intercept") %>%
        # add the grand mean to the group-specific deviations
        mutate(mu = b_HI_LOHI + r_study) %>%
        ungroup() %>%
        mutate(study = str_replace_all(study, "[.]", " ")) %>% 
        
        # plot
        ggplot(aes(x = mu, y = reorder(study, mu))) +
        geom_vline(xintercept = fixef(m1)[2, 1], color = "white", size = 1) +
        geom_vline(xintercept = fixef(m1)[2, 3:4], color = "white", linetype = 2) +
        stat_halfeye(.width = .95, size = 1) +
        labs(x = expression(italic("Cohen's d")),
             y = NULL) +
        theme(panel.grid   = element_blank(),
              axis.ticks.y = element_blank(),
              axis.text.y  = element_text(hjust = 0))


posterior_samples(m2) %>% 
        select(starts_with("sd")) %>% 
        gather(key, tau) %>% 
        mutate(key = str_remove(key, "sd_") %>% str_remove(., "__Intercept")) %>% 
        
        ggplot(aes(x = tau, fill = key)) +
        geom_density(color = "transparent", alpha = 2/3) +
        scale_fill_viridis_d(NULL, end = .85) +
        scale_y_continuous(NULL, breaks = NULL) +
        xlab(expression(tau)) +
        theme(panel.grid = element_blank())


m2 %>%
        spread_draws(b_Intercept, r_outcome_group[outcome_group,]) %>%
        # add the grand mean to the group-specific deviations
        mutate(mu = b_Intercept + r_outcome_group) %>%
        ungroup() %>%
        mutate(outcome = str_replace_all(outcome_group, "[.]", " ")) %>% 
        
        # plot
        ggplot(aes(x = mu, y = reorder(outcome_group, mu))) +
        geom_vline(xintercept = fixef(m2)[1, 1], color = "white", size = 1) +
        geom_vline(xintercept = fixef(m2)[1, 3:4], color = "white", linetype = 2) +
        geom_halfeyeh(.width = .95, size = 2/3) +
        labs(x = expression(italic("Cohen's d")),
             y = NULL) +
        theme(panel.grid   = element_blank(),
              axis.ticks.y = element_blank(),
              axis.text.y  = element_text(hjust = 0))


### First model with high vs low 
get_prior(es.ub | se(vd) ~ HI_LO + (HI_LO| study) + (1|outcome_group), 
          data = es[es$muscle_size == "TRUE",])


m2 <- brm(data = es[es$muscle_size == "TRUE",], 
          family = gaussian,
          es.ub | se(vd) ~ HI_LO + (HI_LO || study),
          prior = c(prior(normal(0, 1), class = Intercept),
                    prior(cauchy(0, 1), class = sd)),
          iter = 10000, warmup = 5000, cores = 4, chains = 4,
          control = list(adapt_delta = 0.98),
          seed = 14)

summary(m2)


### Plot of outcomes 
m2 %>%
        spread_draws(b_Intercept, r_outcome_group[outcome_group,]) %>%
        # add the grand mean to the group-specific deviations
        mutate(mu = b_Intercept + r_outcome_group) %>%
        ungroup() %>%
        mutate(outcome = str_replace_all(outcome_group, "[.]", " ")) %>% 
        
        # plot
        ggplot(aes(x = mu, y = reorder(outcome_group, mu))) +
        geom_vline(xintercept = fixef(m2)[1, 1], color = "white", size = 1) +
        geom_vline(xintercept = fixef(m2)[1, 3:4], color = "white", linetype = 2) +
        geom_halfeyeh(.width = .95, size = 2/3) +
        labs(x = expression(italic("Cohen's d")),
             y = NULL) +
        theme(panel.grid   = element_blank(),
              axis.ticks.y = element_blank(),
              axis.text.y  = element_text(hjust = 0))

### Plot of studies
m2 %>%
        spread_draws(b_Intercept, r_study[study,]) %>%
        # add the grand mean to the group-specific deviations
        mutate(mu = b_Intercept + r_study) %>%
        ungroup() %>%
        mutate(outcome = str_replace_all(study, "[.]", " ")) %>% 
        
        # plot
        ggplot(aes(x = mu, y = reorder(study, mu))) +
        geom_vline(xintercept = fixef(m2)[1, 1], color = "white", size = 1) +
        geom_vline(xintercept = fixef(m2)[1, 3:4], color = "white", linetype = 2) +
        geom_halfeyeh(.width = .95, size = 2/3) +
        labs(x = expression(italic("Cohen's d")),
             y = NULL) +
        theme(panel.grid   = element_blank(),
              axis.ticks.y = element_blank(),
              axis.text.y  = element_text(hjust = 0))





m_str <- rma.mv(yi = es, V = var, mods = ~ group , random = ~ 1|study,  
                method = "REML",
                data = es.agg[es.agg$muscle_size == "FALSE",])

m_mass <- rma.mv(yi = es, V =  var, mods = ~ group , random = ~ 1|study,  
                 method = "REML",
                 data = es[es.agg$muscle_size == "TRUE",])


?rma.mv


summary(m_mass)

forest(m_mass)
funnel(m_mass)




?rma
m_mass <- rma.mv(yi = es.ub, V =  vd, mods = ~ weekly_sets , random = ~ 1|study/outcome,  
                 method = "REML",
                 data = es[es$muscle_size == "TRUE",])

summary(m_str)
summary(m_mass)


m <- lme(es.ub ~ weekly_sets, random = list(study = ~ 1, outcome = ~ 1), 
         weights = varFixed(~ vd),
         control = lmeControl(sigma = 1),
         data = es)
intervals(m)

emmeans(m, specs = ~ HI_LO)



ungroup() %>%
        mutate(es = (post_m - pre_m)/pooled.sd,
               weekly_sets = n_sets * n_exercises) %>%
        dplyr::select(study, group, sex, HI_LO,weekly_sets,n_weeks,es, pre_m, pre_sd, post_m, post_sd, n_participants) %>%
        mutate(r = 0.5) %>%
        
        print()








print()

temp2 <- escalc(data = temp, measure = "SMCR", m2i = pre_m, m1i = post_m, sd2i = pre_sd, sd1i = post_sd, 
                ni = n_participants, ri = r) %>%
        data.frame()



temp2 %>%
        ggplot(aes(weekly_sets, yi)) + geom_point()




m <- rma.mv(yi = yi,V =  vi, mods = ~ weekly_sets + n_weeks, random = ~ 1|study/group,  data = temp2)

summary(m)
profile(m)


