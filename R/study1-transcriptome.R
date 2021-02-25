##-------------------------------------
## study1-transcriptome.R
##
## Title: Transcriptome analyses in results/discussion
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


# Note that some analyses are commented out to improve rendering speed. 
# E.g. comparison of models for different random effects formulations. 

# Load libraries and functions

source("./R/libraries.R")
source("./R/themes.R")



source("./R/study1b/figure_source.R")
source("./R/study1b/modified_upset.R")



firstup <- function(x) {
        substr(x, 1, 1) <- toupper(substr(x, 1, 1))
        x
}


new_line_fun <- function(x){
        
        l <- sapply(strsplit(x, " "), length)
        
        if(l > 2){
                
                paste(paste(strsplit(x, " ")[[1]][1:2], collapse = " "), 
                      "\n", 
                      paste(strsplit(x, " ")[[1]][3:l], collapse = " "), 
                      collaspe = "")
                
        } else {
                x
        }
        
}


######## Fig A: Muscle weight figure ############



# Muscle weights in cDNA synthesis are loaded ##
mw <- read_excel("./data/study-1b/data/RNAamount.xlsx") %>%
        mutate(RNA.per.mg = (conc * elution.volume) / prot.mrna1) %>% 
        filter(include == "incl") %>%
        inner_join(read_csv2("./data/study-1b/data/oneThreeSetLeg.csv") %>%
                           gather("sets", "leg", multiple:single)) %>%
        filter(rnaseq_include == "incl") %>%
        filter(timepoint != "w2post") %>%
        dplyr::select(subject, sets, time = timepoint, RNA.per.mg) %>%
        mutate(tissue = 1000 / RNA.per.mg, 
               tl = log(tissue), 
               sets = factor(sets, levels = c("single", "multiple")))  %>%
        filter(!(subject == "FP15" & time == "w2pre")) 
# Filters out one subject with known discrepancy of RNA to muscle weight. 

# m0 <- lme(log(tissue) ~ time + time:sets, random = list(subject = ~1), data = mw, 
#           control = lmeControl(msMaxIter = 100, msVerbose = TRUE, 
#                                opt = "optim"))

m1 <- lme(log(tissue) ~ time + time:sets, random = list(subject = ~1 + time), data = mw, 
          control = lmeControl(msMaxIter = 100, msVerbose = TRUE, 
                               
                               
                               
                               opt = "optim"))

# anova(m0, m1)
# plot(m1)



muscle_weight_fig <- emmeans(m1, specs = ~"time|sets") %>%
        data.frame() %>%
        mutate(time = factor(time, levels = c("w0", "w2pre", "w12"), labels = c("Week 0", "Week 2", "Week 12")), 
               sets = factor(sets, levels = c("single", "multiple"), labels = c("LOW", "MOD"))) %>%
        ggplot(aes(time, exp(emmean), fill = sets)) + 
        geom_errorbar(aes(ymin = exp(lower.CL), ymax = exp(upper.CL)), 
                      position = position_dodge(width = 0.35), 
                      width = 0, 
                      size = line.size) +
        geom_point(position = position_dodge(width = 0.35), 
                   shape = 21, 
                   size = 2.5) +
        
        scale_fill_manual(values = c(group.study.color[1], group.study.color[2])) +
        scale_y_continuous(limits = c(2, 4), 
                           breaks = c(2, 2.5, 3,3.5, 4), 
                           labels = c(2,  "", 3, "", 4),
                           expand = c(0, 0)) +
        
        labs(fill = "", 
             x = "Time-point", 
             y = "Muscle biopsy mass<br>in cDNA synthesis (mg)") +
        
        theme_bw() +
        theme(panel.grid = element_blank(), 
              panel.border = element_blank(), 
              axis.line.y  = element_line(size = line.size), 
              axis.line.x = element_blank(), 
              axis.ticks.y = element_line(size = line.size), 
              axis.ticks.x = element_blank(), 
              axis.text.y = element_markdown(color = "black", size = 12), 
              axis.title.y = element_markdown(color = "black", size = 12),
              axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              legend.title = element_blank(), 
              legend.background = element_rect(fill = "gray99", color = "gray75"),
              legend.margin = margin(t = 1, r = 3, b = 2, l = 1, unit = "pt"),
              legend.key = element_rect(fill = "white"),
              legend.position = c(0.20, 0.9))




### ################### Total mRNA counts (library size) per time-point ##########



# read DGEList

dge_lists <- readRDS("./data/study-1b/data/derivedData/dge_lists/dge_list.RDS")



# using rsem for this example #

selected.method <- "rsem"


# Muscle weight in cDNA synthesis

# Muscle weights in cDNA synthesis are loaded ##
lib.size.data <- read_excel("./data/study-1b/data/RNAamount.xlsx") %>%
        mutate(RNA.per.mg = (conc * elution.volume) / prot.mrna1) %>% 
        filter(include == "incl") %>%
        inner_join(read_csv2("./data/study-1b/data/oneThreeSetLeg.csv") %>%
                           gather("sets", "leg", multiple:single)) %>%
        filter(rnaseq_include == "incl") %>%
        filter(timepoint != "w2post") %>%
        dplyr::select(subject, sets, time = timepoint, RNA.per.mg) %>%
        mutate(tissue = 1000 / RNA.per.mg, 
               tl = log(tissue))  %>%
        filter(!(subject == "FP15" & time == "w2pre")) %>% 
        # This filters out one subject with bad estimation of RNA to mg tissue estimate
        # See script in training study for details
        inner_join(dge_lists[[selected.method]]$samples %>% 
                           mutate(eff.lib = lib.size * norm.factors)) %>%
        mutate(eff.lib = as.integer(round(eff.lib, 0)), 
               eff.lib.scaled = eff.lib/10^6) %>%
        mutate(time = factor(time, levels = c("w0", "w2pre", "w12")), 
               sets = factor(sets, levels = c("single", 
                                              "multiple"))) 



# Modeling effective library size 
# lib.m0 <- lme(log(eff.lib.scaled) ~ time + time:sets, random = list(subject = ~ 1), 
#               data = lib.size.data)
# 

lib.m1 <- lme(log(eff.lib.scaled) ~ time + time:sets, random = list(subject = ~ 1),
              weights = varIdent(form = ~ 1|time), 
              data = lib.size.data)


# Modeling library size per tissue weight


# library per tissue
# lt.m0 <- lme(log(eff.lib.scaled/tissue) ~ time + time:sets, random = list(subject = ~1), 
#              data = lib.size.data)
# 
# 

lt.m1 <- lme(log(eff.lib.scaled/tissue) ~ time + time:sets, random = list(subject = ~1), 
             data = lib.size.data, weights = varIdent(form = ~ 1|time))



# Accounting for heteroscedasticity? 
#plot(lt.m0, main = "Library size per tissue, m0"); plot(lt.m1, main = "Library size per tissue, m1") 

# anova(lt.m0, lt.m1) 

# Clear evidence of better models using weights

# Random slopes? 

# lt.m2 <- lme(log(eff.lib.scaled/tissue) ~ time + time:sets, random = list(subject = ~1 + time), 
#              data = lib.size.data, weights = varIdent(form = ~ 1|time),
#              control = lmeControl(msMaxIter = 100, msVerbose = TRUE, 
#                                   opt = "optim"))
# 
# anova(lt.m1, lt.m2) # not called for



# Extract statistics




# Plotting the data 
per.rna <- emmeans(lib.m1, specs = ~ "time|sets") %>%
        data.frame() %>%
        mutate(time = factor(time, levels = c("w0", "w2pre", "w12"), labels = c("Week 0", "Week 2", "Week 12")), 
               sets = factor(sets, levels = c("single", "multiple"), labels = c("Single-set", "Multiple-set"))) %>%
        
        ggplot(aes(time, exp(emmean), fill = sets )) +
        geom_errorbar(aes(ymin = exp(lower.CL), ymax = exp(upper.CL)), 
                      position = position_dodge(width = 0.35), 
                      width = 0, 
                      size = line.size) +
        geom_point(position = position_dodge(width = 0.35), 
                   shape = 21, 
                   size = 2.5) +
        scale_fill_manual(values = c(group.study.color[1], group.study.color[2])) +
        scale_y_continuous(limits = c(8, 16), 
                           breaks = c(8, 10, 12, 14, 16), 
                           labels = c(8, "", 12, "", 16),
                           expand = c(0, 0)) +
        
        
        labs(fill = "", 
             x = "Time-point", 
             y = bquote("Effective library size<br>(million counts)")) +
        theme_bw() +
        theme(panel.grid = element_blank(), 
              panel.border = element_blank(), 
              axis.line.y  = element_line(size = line.size), 
              axis.line.x = element_blank(), 
              axis.ticks.y = element_line(size = line.size), 
              axis.ticks.x = element_blank(), 
              axis.text.y = element_markdown(color = "black", size = 12), 
              axis.title.y = element_markdown(color = "black", size = 12),
              axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              legend.title = element_blank(), 
              legend.background = element_rect(fill = "white"),
              legend.margin = margin(t = 0, r = 1, b = 1, l = 1, unit = "pt"),
              legend.key = element_rect(fill = "white"),
              legend.position = "none")



per.tissue <- emmeans(lt.m1, specs = ~ "time|sets") %>%
        data.frame() %>%
        mutate(time = factor(time, levels = c("w0", "w2pre", "w12"), labels = c("Week 0", "Week 2", "Week 12")), 
               sets = factor(sets, levels = c("single", "multiple"), labels = c("Single-set", "Multiple-set"))) %>%
        
        ggplot(aes(time, exp(emmean), fill = sets )) +
        geom_errorbar(aes(ymin = exp(lower.CL), ymax = exp(upper.CL)), 
                      position = position_dodge(width = 0.35), 
                      width = 0, 
                      size = line_size) +
        geom_point(position = position_dodge(width = 0.35), 
                   shape = 21, 
                   size = 2.5) +
        
        scale_fill_manual(values = c(color.scale[1], color.scale[2])) +
        scale_y_continuous(limits = c(2, 6), 
                           breaks = c(2, 3, 4, 5, 6), 
                           labels = c(2, "", 4, "", 6),
                           expand = c(0, 0)) +
        
        # Sets effects stats segment
        
        # geom_segment(data = data.frame(x1 = c(1.9, 2.9), 
        #                                x2 = c(2.1, 3.1), 
        #                                y1 = c(2.92, 3.07), 
        #                                y2 = c(2.92, 3.07)), 
        #              aes(x = x1, xend = x2, y = y1, yend = y2), 
        #              inherit.aes = FALSE) +
        
        labs(fill = "", 
             x = "Time-point", 
             y = "Effective library size<br>(million counts mg^-1)") +
        theme_bw() +
        theme(panel.grid = element_blank(), 
              panel.border = element_blank(), 
              axis.line.y  = element_line(size = line_size), 
              axis.line.x = element_blank(), 
              axis.ticks.y = element_line(size = line_size), 
              axis.ticks.x = element_blank(), 
              axis.text.y = element_markdown(color = "black", size = 12), 
              axis.title.y = element_markdown(color = "black", size = 12),
              axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              legend.title = element_blank(), 
              legend.background = element_rect(fill = "white"),
              legend.margin = margin(t = 0, r = 1, b = 1, l = 1, unit = "pt"),
              legend.key = element_rect(fill = "white"),
              legend.position = "none")

# intervals  ##############


rna.diff <- intervals(m1)$fixed %>%
        data.frame() %>%
        mutate(coef = rownames(.)) %>%
        filter(coef %in% c("timew0:setsmultiple", 
                           "timew12:setsmultiple", 
                           "timew2pre:setsmultiple")) %>%
        mutate(coef = factor(coef, levels = c("timew0:setsmultiple", 
                                              
                                              "timew2pre:setsmultiple",
                                              "timew12:setsmultiple"), 
                             labels = c("Week 0", 
                                        "Week 2", 
                                        "Week 12"))) %>%
        
        ggplot(aes(coef, exp(est.))) + 
        
        geom_hline(yintercept = 1, color = "gray50", lty = 2) +
        
        geom_errorbar(aes(ymin = exp(lower), 
                          ymax = exp(upper)), 
                      width = 0.1) +
        geom_point(shape = 24, size = 2.5, fill = "red" ) +
        labs(x = "Time-point", 
             y = "MOD / LOW") +
        
        scale_y_continuous(limits = c(0.75, 1.25), 
                           expand = c(0,0), 
                           breaks = c(0.75, 0.875, 1, 1.125, 1.25), 
                           labels = c(0.75, "", 1, "", 1.25)) +
        
        
        theme_bw() +
        theme(panel.grid = element_blank(), 
              panel.border = element_blank(), 
              axis.line.y  = element_line(size = line_size), 
              axis.line.x = element_blank(), 
              axis.ticks.y = element_line(size = line_size), 
              axis.ticks.x = element_blank(), 
              axis.text.y =  element_markdown(color = "black", size = 12), 
              axis.title.y = element_markdown(color = "black", size = 12),
              axis.text.x =  element_markdown(color = "black", size = 12), 
              axis.title.x = element_blank(),
              legend.title = element_blank(), 
              legend.background = element_rect(fill = "white"),
              legend.margin = margin(t = 0, r = 1, b = 1, l = 1, unit = "pt"),
              legend.key = element_rect(fill = "white"),
              legend.position = "none")




lib.diff <- intervals(lib.m1)$fixed %>%
        data.frame() %>%
        mutate(coef = rownames(.)) %>%
        filter(coef %in% c("timew0:setsmultiple", 
                           "timew12:setsmultiple", 
                           "timew2pre:setsmultiple")) %>%
        mutate(coef = factor(coef, levels = c("timew0:setsmultiple", 
                                              
                                              "timew2pre:setsmultiple",
                                              "timew12:setsmultiple"), 
                             labels = c("Week 0", 
                                        "Week 2", 
                                        "Week 12"))) %>%
        
        ggplot(aes(coef, exp(est.))) + 
        
        geom_hline(yintercept = 1, color = "gray50", lty = 2) +
        
        geom_errorbar(aes(ymin = exp(lower), 
                          ymax = exp(upper)), 
                      width = 0.1) +
        geom_point(shape = 24, size = 2.5, fill = "red" ) +
        labs(x = "Time-point", 
             y = "") +
        
        scale_y_continuous(limits = c(0.75, 1.25), 
                           expand = c(0,0), 
                           breaks = c(0.75, 0.875, 1, 1.125, 1.25), 
                           labels = c(0.75, "", 1, "", 1.25)) +
        
        theme_bw() +
        theme(panel.grid = element_blank(), 
              panel.border = element_blank(), 
              axis.line.y  = element_line(size = line_size), 
              axis.line.x = element_blank(), 
              axis.ticks.y = element_line(size = line_size), 
              axis.ticks.x = element_blank(), 
              axis.text.y =  element_markdown(color = "black", size = 12), 
              axis.title.y = element_markdown(color = "black", size = 12),
              axis.text.x =  element_markdown(color = "black", size = 12), 
              axis.title.x = element_blank(),
              legend.title = element_blank(), 
              legend.background = element_rect(fill = "white"),
              legend.margin = margin(t = 0, r = 1, b = 1, l = 1, unit = "pt"),
              legend.key = element_rect(fill = "white"),
              legend.position = "none")


tissue.diff <- intervals(lt.m1)$fixed %>%
        data.frame() %>%
        mutate(coef = rownames(.)) %>%
        filter(coef %in% c("timew0:setsmultiple", 
                           "timew12:setsmultiple", 
                           "timew2pre:setsmultiple")) %>%
        mutate(coef = factor(coef, levels = c("timew0:setsmultiple", 
                                              
                                              "timew2pre:setsmultiple",
                                              "timew12:setsmultiple"), 
                             labels = c("Week 0", 
                                        "Week 2", 
                                        "Week 12"))) %>%
        
        ggplot(aes(coef, exp(est.))) + 
        
        geom_hline(yintercept = 1, color = "gray50", lty = 2) +
        
        geom_errorbar(aes(ymin = exp(lower), 
                          ymax = exp(upper)), 
                      width = 0.1) +
        geom_point(shape = 24, size = 2.5, fill = "red" ) +
        labs(x = "Time-point", 
             y = "") +
        
        scale_y_continuous(limits = c(0.75, 1.25), 
                           expand = c(0,0), 
                           breaks = c(0.75, 0.875, 1, 1.125, 1.25), 
                           labels = c(0.75, "", 1, "", 1.25)) +
        
        
        theme_bw() +
        theme(panel.grid = element_blank(), 
              panel.border = element_blank(), 
              axis.line.y  = element_line(size = line_size), 
              axis.line.x = element_blank(), 
              axis.ticks.y = element_line(size = line_size), 
              axis.ticks.x = element_blank(), 
              axis.text.y =  element_markdown(color = "black", size = 12), 
              axis.title.y = element_markdown(color = "black", size = 12),
              axis.text.x =  element_markdown(color = "black", size = 12), 
              axis.title.x = element_blank(),
              legend.title = element_blank(), 
              legend.background = element_rect(fill = "white"),
              legend.margin = margin(t = 0, r = 1, b = 1, l = 1, unit = "pt"),
              legend.key = element_rect(fill = "white"),
              legend.position = "none")



abc <- plot_grid(
        plot_grid(muscle_weight_fig, NULL, rna.diff, ncol = 1, 
                  rel_heights = c(1, 0.05, 0.6), 
                  align = "v"),
        
        plot_grid(per.rna, NULL, lib.diff, ncol = 1, 
                  rel_heights = c(1, 0.05, 0.6), 
                  align = "v"),
        
        plot_grid(per.tissue,NULL, tissue.diff, ncol = 1, 
                  rel_heights = c(1, 0.05, 0.6), 
                  align = "v"),
        ncol = 3)

## Save library size comparisons ###################


saveRDS(abc, "./data/figures/rna_libsize_norm.RDS")

