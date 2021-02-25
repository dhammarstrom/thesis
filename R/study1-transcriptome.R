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


source("./R/study1-modified_upset.R")



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
               sets = factor(sets, levels = c("single", "multiple"), labels = c("Low-volume", "Moderate-volume"))) %>%
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
              axis.text.y = element_markdown(color = "black", size = 7), 
              axis.title.y = element_markdown(color = "black", size = 7),
              axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              legend.title = element_blank(), 
              legend.background = element_blank(),
              legend.text = element_text(size = 7, 
                                         margin = margin(t = -1, r = 0, b = 0, l = -9), 
                                         lineheight = 0.5),
              legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
              legend.spacing.y = unit(0.03, "cm"),
              legend.key = element_rect(fill = "white"),
              legend.key.height = unit(0.2, "cm"),
              legend.position = c(0.5, 0.85))




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
              axis.text.y = element_markdown(color = "black", size = 7), 
              axis.title.y = element_markdown(color = "black", size = 7),
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
                      size = line.size) +
        geom_point(position = position_dodge(width = 0.35), 
                   shape = 21, 
                   size = 2.5) +
        
        scale_fill_manual(values = c(group.study.color[1], group.study.color[2])) +
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
              axis.line.y  = element_line(size = line.size), 
              axis.line.x = element_blank(), 
              axis.ticks.y = element_line(size = line.size), 
              axis.ticks.x = element_blank(), 
              axis.text.y = element_markdown(color = "black", size = 7), 
              axis.title.y = element_markdown(color = "black", size = 7),
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
                      width = 0) +
        geom_point(shape = 24, size = 2, fill = group.study.color[2] ) +
        labs(x = "Time-point", 
             y = "Moderate-volume/<br>Low-volume") +
        
        scale_y_continuous(limits = c(0.75, 1.25), 
                           expand = c(0,0), 
                           breaks = c(0.75, 0.875, 1, 1.125, 1.25), 
                           labels = c(0.75, "", 1, "", 1.25)) +
        
        
        theme_bw() +
        theme(panel.grid = element_blank(), 
              panel.border = element_blank(), 
              axis.line.y  = element_line(size = line.size), 
              axis.line.x = element_blank(), 
              axis.ticks.y = element_line(size = line.size), 
              axis.ticks.x = element_blank(), 
              axis.text.y =  element_markdown(color = "black", size = 7), 
              axis.title.y = element_markdown(color = "black", size = 7),
              axis.text.x =  element_text(color = "black", size = 7, angle = 45, hjust = 1), 
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
                      width = 0) +
        geom_point(shape = 24, size = 2, fill = group.study.color[2]) +
        labs(x = "Time-point", 
             y = "") +
        
        scale_y_continuous(limits = c(0.75, 1.25), 
                           expand = c(0,0), 
                           breaks = c(0.75, 0.875, 1, 1.125, 1.25), 
                           labels = c(0.75, "", 1, "", 1.25)) +
        
        theme_bw() +
        theme(panel.grid = element_blank(), 
              panel.border = element_blank(), 
              axis.line.y  = element_line(size = line.size), 
              axis.line.x = element_blank(), 
              axis.ticks.y = element_line(size = line.size), 
              axis.ticks.x = element_blank(), 
              axis.text.y =  element_markdown(color = "black", size = 7), 
              axis.title.y = element_markdown(color = "black", size = 7),
              axis.text.x =  element_text(color = "black", size = 7, angle = 45, hjust = 1), 
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
                      width = 0) +
        geom_point(shape = 24, size = 2, fill = group.study.color[2]) +
        labs(x = "Time-point", 
             y = "") +
        
        scale_y_continuous(limits = c(0.75, 1.25), 
                           expand = c(0,0), 
                           breaks = c(0.75, 0.875, 1, 1.125, 1.25), 
                           labels = c(0.75, "", 1, "", 1.25)) +
        
        
        theme_bw() +
        theme(panel.grid = element_blank(), 
              panel.border = element_blank(), 
              axis.line.y  = element_line(size = line.size), 
              axis.line.x = element_blank(), 
              axis.ticks.y = element_line(size = line.size), 
              axis.ticks.x = element_blank(),
              axis.text.y =  element_markdown(color = "black", size = 7), 
              axis.title.y = element_markdown(color = "black", size = 7),
              axis.text.x =  element_text(color = "black", size = 7, angle = 45, hjust = 1), 
              axis.title.x = element_blank(),
              legend.title = element_blank(), 
              legend.background = element_rect(fill = "white"),
              legend.margin = margin(t = 0, r = 1, b = 1, l = 1, unit = "pt"),
              legend.key = element_rect(fill = "white"),
              legend.position = "none")



abc <- plot_grid(
        plot_grid(muscle_weight_fig, NULL, rna.diff, ncol = 1, 
                  rel_heights = c(0.65, 0.05, 0.6), 
                  align = "v"),
        
        plot_grid(per.rna, NULL, lib.diff, ncol = 1, 
                  rel_heights = c(0.65, 0.05, 0.6), 
                  align = "v"),
        
        plot_grid(per.tissue,NULL, tissue.diff, ncol = 1, 
                  rel_heights = c(0.65, 0.05, 0.6), 
                  align = "v"),
        ncol = 3)

## Save library size comparisons ###################


saveRDS(abc, "./data/derivedData/study1-transcriptome/muscle-weight-lib-size.RDS")



#### Volcano plots 


#### Volcano plots #################

# Read DE data 

mm2_2 <- readRDS(file = "./data/study-1b/data/derivedData/DE/mixedmodel2_results2.RDS")


mm2_1 <- readRDS(file = "./data/study-1b/data/derivedData/DE/mixedmodel2_results1.RDS")

### Sets comparisons #####

volcano_data <- mm2_1 %>%
        filter(coef %in% c("timew0:setsmultiple",  "timew2pre:setsmultiple", "timew2post:setsmultiple", 
                           "timew12:setsmultiple")) %>%
        mutate(coef = factor(coef, levels = c("timew0:setsmultiple", "timew2pre:setsmultiple", "timew12:setsmultiple"),
                             labels = c("Week 0: Multiple - Single-set", 
                                        "Week 2: Multiple - Single-set", 
                                        "Week 12: Multiple - Single-set")),
               model = factor(model, levels = c("naive", "lib_size_normalized", "tissue_offset_lib_size_normalized"), 
                              labels = c("No normalization", 
                                         "Effective library size", 
                                         "Per tissue mass +\nEffective library size"))) %>%
        
        group_by(coef, model) %>%
        mutate(adj.p = p.adjust(p.val, method = "fdr"), 
               pthreshold = if_else(adj.p > 0.05, "ns", "s"), 
               log2fc = estimate/log(2), 
               fcthreshold = if_else(abs(log2fc) > 0.5, "s", "ns"), 
               sig = if_else(fcthreshold == "s" & pthreshold == "s", "s", "ns"), 
               regulation = if_else(sig == "s" & log2fc > 0.5, "up", 
                                    if_else(sig == "s" & log2fc < -0.5, "down", "ns")), 
               regulation = factor(regulation, levels = c("ns", "up", "down"))) 



##### Settings specifically for figure 4 ##################

## Color themes in figure source 
# using volcano.regulation.color
volcano.regulation.color <- c("gray30", group.study.color[2], group.study.color[4])
## Sizes volcano plots

volcano_size <- c(1.2, 1.9, 1.9)
volcano_alpha <- c(0.2, 0.8, 0.8)


############## Week 2 Volcano plots #####################################


week2_naive <- volcano_data %>%
        filter(model == "No normalization" & coef == "Week 2: Multiple - Single-set") %>%
        
        ggplot(aes(log2fc, -log10(p.val), fill = regulation, alpha = regulation, size = regulation)) + 
        geom_point(shape = 21) +
        
        # Annotate plot
        annotate("text", x = -3.5, y = 7.7, label = "Naive", hjust = 0, size = 2.5) +
        
        # Scales
        scale_color_manual(values = volcano.regulation.color[c(1,3)]) + 
        scale_fill_manual(values = volcano.regulation.color[c(1,3)]) +
        scale_size_manual(values = volcano_size) +
        scale_alpha_manual(values = volcano_alpha) +
        
        scale_x_continuous(breaks = c(-4, -2, 0, 2, 4), limits = c(-4, 4), expand = c(0, 0)) +
        scale_y_continuous(breaks = c(0, 2, 4, 6, 8), limits = c(0, 8), expand = c(0,0)) +
        
        labs(x = " ", 
             y = "") +
        dissertation_theme() +
        
        theme(panel.grid = element_blank(), 
              panel.border = element_blank(), 
              axis.line.y  = element_line(size = line.size), 
              axis.line.x = element_line(size = line.size), 
              axis.ticks.y = element_line(size = line.size), 
              axis.ticks.x = element_line(size = line.size), 
              axis.text.y =  element_markdown(color = "black", size = 7), 
              axis.title.y = element_markdown(color = "black", size = 7),
              axis.text.x =  element_text(color = "black", size = 7), 
              axis.title.x = element_blank(),
              legend.title = element_blank(), 
              legend.background = element_rect(fill = "white"),
              legend.margin = margin(t = 0, r = 1, b = 1, l = 1, unit = "pt"),
              legend.key = element_rect(fill = "white"),
              legend.position = "none")

week2_libsize <- volcano_data %>%
        filter(model == "Effective library size" & coef == "Week 2: Multiple - Single-set") %>%
        
        ggplot(aes(log2fc, -log10(p.val), fill = regulation, alpha = regulation, size = regulation)) + 
        geom_point(shape = 21) +
        
        annotate("text", x = -3.5, y = 7.7, label = "Effective library-\nsize", hjust = 0,vjust = 1, size = 2.5) +
        
        # Scales
        scale_color_manual(values = volcano.regulation.color) + 
        scale_fill_manual(values = volcano.regulation.color) +
        scale_size_manual(values = volcano_size) +
        scale_alpha_manual(values = volcano_alpha) +
        
        scale_x_continuous(breaks = c(-4, -2, 0, 2, 4), limits = c(-4, 4), expand = c(0, 0)) +
        scale_y_continuous(breaks = c(0, 2, 4, 6, 8), limits = c(0, 8), expand = c(0,0)) +
        
        
        labs(x = "Log<sub>2</sub>(MOD / LOW)", 
             y = "") +
        
        dissertation_theme() +
        
        theme(panel.grid = element_blank(), 
              panel.border = element_blank(), 
              axis.line.y  = element_line(size = line.size), 
              axis.line.x = element_line(size = line.size), 
              axis.ticks.y = element_line(size = line.size), 
              axis.ticks.x = element_line(size = line.size), 
              axis.text.y =  element_markdown(color = "black", size = 7), 
              axis.title.y = element_markdown(color = "black", size = 7),
              axis.text.x =  element_text(color = "black", size = 7), 
              axis.title.x = element_blank(),
              legend.title = element_blank(), 
              legend.background = element_rect(fill = "white"),
              legend.margin = margin(t = 0, r = 1, b = 1, l = 1, unit = "pt"),
              legend.key = element_rect(fill = "white"),
              legend.position = "none")


week2_tissue <- volcano_data %>%
        filter(model == "Per tissue mass +\nEffective library size" & coef == "Week 2: Multiple - Single-set") %>%
        
        ggplot(aes(log2fc, -log10(p.val), fill = regulation, alpha = regulation, size = regulation)) + 
        geom_point(shape = 21) +
        
        annotate("text", x = -3.5, y = 7.7, label = "Tissue offset", hjust = 0, size = 2.5) +
        
        # Scales
        scale_color_manual(values = volcano.regulation.color) + 
        scale_fill_manual(values = volcano.regulation.color) +
        scale_size_manual(values = volcano_size) +
        scale_alpha_manual(values = volcano_alpha) +
        
        scale_x_continuous(breaks = c(-4, -2, 0, 2, 4), limits = c(-4, 4), expand = c(0, 0)) +
        scale_y_continuous(breaks = c(0, 2, 4, 6, 8), limits = c(0, 8), expand = c(0,0)) +
        
        labs(x = " ", 
             y = "-Log<sub>10</sub>(<em>P</em>-value)") +
        
        dissertation_theme() +
        
        theme(panel.grid = element_blank(), 
              panel.border = element_blank(), 
              axis.line.y  = element_line(size = line.size), 
              axis.line.x = element_line(size = line.size), 
              axis.ticks.y = element_line(size = line.size), 
              axis.ticks.x = element_line(size = line.size), 
              axis.text.y =  element_markdown(color = "black", size = 7), 
              axis.title.y = element_markdown(color = "black", size = 7),
              axis.text.x =  element_text(color = "black", size = 7), 
              axis.title.x = element_blank(),
              legend.title = element_blank(), 
              legend.background = element_rect(fill = "white"),
              legend.margin = margin(t = 0, r = 1, b = 1, l = 1, unit = "pt"),
              legend.key = element_rect(fill = "white"),
              legend.position = "none")

#### Upset plots for sets-effects ###############


venn_list_up <- list() 
venn_list_down <- list()

Coef <- c("timew2pre:setsmultiple", "timew2pre:setsmultiple", "timew2pre:setsmultiple",
          "timew12:setsmultiple", "timew12:setsmultiple", "timew12:setsmultiple")
Model <- c("naive", "lib_size_normalized", "tissue_offset_lib_size_normalized", 
           "naive", "lib_size_normalized", "tissue_offset_lib_size_normalized")

i <- j <-  5

for(i in 1:6){
        
        temp_up <- mm2_1 %>%
                group_by(model, coef) %>%
                mutate(adj.p = p.adjust(p.val, method = "fdr"), 
                       pthreshold = if_else(adj.p > 0.05, "ns", "s"), 
                       log2fc = estimate/log(2), 
                       fcthreshold = if_else(abs(log2fc) > 0.5, "s", "ns"), 
                       sig = if_else(fcthreshold == "s" & pthreshold == "s", "s", "ns"), 
                       regulation = if_else(sig == "s" & log2fc > 0.5, "up", 
                                            if_else(sig == "s" & log2fc < -0.5, "down", "ns")), 
                       regulation = factor(regulation, levels = c("ns", "up", "down"))) %>%
                
                filter(sig != "ns", 
                       regulation == "up") %>%
                filter(coef == Coef[i] & model == Model[i]) %>%
                
                data.frame()
        
        
        temp_down <- mm2_1 %>%
                group_by(model, coef) %>%
                mutate(adj.p = p.adjust(p.val, method = "fdr"), 
                       pthreshold = if_else(adj.p > 0.05, "ns", "s"), 
                       log2fc = estimate/log(2), 
                       fcthreshold = if_else(abs(log2fc) > 0.5, "s", "ns"), 
                       sig = if_else(fcthreshold == "s" & pthreshold == "s", "s", "ns"), 
                       regulation = if_else(sig == "s" & log2fc > 0.5, "up", 
                                            if_else(sig == "s" & log2fc < -0.5, "down", "ns")), 
                       regulation = factor(regulation, levels = c("ns", "up", "down"))) %>%
                
                filter(sig != "ns", 
                       regulation == "down") %>%
                filter(coef == Coef[i] & model == Model[i]) %>%
                
                data.frame()
        
        
        venn_list_up[[i]] <- temp_up$gene
        venn_list_down[[i]] <- temp_down$gene
        
}


names(venn_list_up) <- c("naive_w2", "els_w2", "to_w2", "naive_w12", "els_w12", "to_w12")
names(venn_list_down) <- c("naive_w2", "els_w2", "to_w2", "naive_w12", "els_w12", "to_w12")



upset.list.w2 <- modified_upset(list.up = venn_list_up[1:3], 
                                list.down = venn_list_down[1:3], 
                                set.order = c("to_w2", "els_w2", "naive_w2"),
                                set.names = c("Tissue-offset", 
                                              "Effective\nlibrary size", 
                                              "Naive"), # Set names, order of display in set size fig, 
                                # should correspond to order
                                legend.labels = c("Moderate-volume > Low-volume", "Moderate-volume < Low-volume"),
                                legend.title = "Gene regulation Week 2",
                                line.size = 0.4, # size of lines in figures (may be changed later)
                                regulation.col = c(group.study.color[2], group.study.color[4]), # Regulation coloring (up vs. down) 
                                bar.width = 0.4, # Width of bars in bar plots
                                text.size = 7, # text size in all figures
                                point.size = 3)  # point size in intersection matrix )


upset.list.w12 <- modified_upset(list.up = venn_list_up[4:6], 
                                 list.down = venn_list_down[4:6], 
                                 set.order = c("to_w12", "els_w12", "naive_w12"),
                                 set.names = c("Tissue-offset", 
                                               "Effective\nlibrary size", 
                                               "Naive"), # Set names, order of display in set size fig, 
                                 # should correspond to order
                                 legend.labels = c("Moderate-volume > Low-volume", "Moderate-volume < Low-volume"),
                                 legend.title = "Gene regulation Week 12",
                                 line.size = 0.4, # size of lines in figures (may be changed later)
                                 regulation.col = c(group.study.color[2], group.study.color[4]), # Regulation coloring (up vs. down) 
                                 bar.width = 0.4, # Width of bars in bar plots
                                 text.size = 7, # text size in all figures
                                 point.size = 3)  # point size in intersection matrix )





upset_w2pre <- ggdraw(plot_grid(NULL,
                                plot_grid(upset.list.w2$set.size.fig, upset.list.w2$intersection.matrix, 
                                          ncol = 2, align = "h"), nrow = 2)) +
        draw_plot(upset.list.w2$intersection.fig, x = 0.47, y = 0.47, width = 0.54, height = 0.55) +
        draw_plot(upset.list.w2$legend, x = 0.05, y = 0.55, width = 0.3, height = 0.3)





upset_w12 <- ggdraw(plot_grid(NULL,
                              plot_grid(upset.list.w12$set.size.fig, upset.list.w12$intersection.matrix, 
                                        ncol = 2, align = "h"), nrow = 2)) +
        draw_plot(upset.list.w12$intersection.fig, x = 0.45, y = 0.47, width = 0.56, height = 0.55) +
        draw_plot(upset.list.w12$legend, x = 0.05, y = 0.55, width = 0.3, height = 0.3)







####### Week 12 volcano plots ##################################

week12_naive <- volcano_data %>%
        filter(model == "No normalization" & coef == "Week 12: Multiple - Single-set") %>%
        
        ggplot(aes(log2fc, -log10(p.val), fill = regulation, alpha = regulation, size = regulation)) + 
        geom_point(shape = 21) +
        
        # Scales
        scale_color_manual(values = volcano.regulation.color[c(1,3)]) + 
        scale_fill_manual(values = volcano.regulation.color[c(1,3)]) +
        scale_size_manual(values = volcano_size) +
        scale_alpha_manual(values = volcano_alpha) +
        
        scale_x_continuous(breaks = c(-4, -2, 0, 2, 4), limits = c(-4, 4), expand = c(0, 0)) +
        scale_y_continuous(breaks = c(0, 2, 4, 6, 8), limits = c(0, 8), expand = c(0,0)) +
        
        annotate("text", x = -3.5, y = 7.7, label = "Naive", hjust = 0, size = 2.2) +
        
        labs(x = "", 
             y = "") +
        
        dissertation_theme() +
        
        theme(panel.grid = element_blank(), 
              panel.border = element_blank(), 
              axis.line.y  = element_line(size = line.size), 
              axis.line.x = element_line(size = line.size), 
              axis.ticks.y = element_line(size = line.size), 
              axis.ticks.x = element_line(size = line.size), 
              axis.text.y =  element_markdown(color = "black", size = 7), 
              axis.title.y = element_markdown(color = "black", size = 7),
              axis.text.x =  element_text(color = "black", size = 7), 
              axis.title.x = element_blank(),
              legend.title = element_blank(), 
              legend.background = element_rect(fill = "white"),
              legend.margin = margin(t = 0, r = 1, b = 1, l = 1, unit = "pt"),
              legend.key = element_rect(fill = "white"),
              legend.position = "none")

week12_libsize <- volcano_data %>%
        filter(model == "Effective library size" & coef == "Week 12: Multiple - Single-set") %>%
        
        ggplot(aes(log2fc, -log10(p.val), fill = regulation, alpha = regulation, size = regulation)) + 
        geom_point(shape = 21) +
        
        annotate("text", x = -3.5, y = 7.7, label = "Effective library-\nsize", hjust = 0,vjust = 1, size = 2.2) +
        
        # Scales
        scale_color_manual(values = volcano.regulation.color[c(1,3)]) + 
        scale_fill_manual(values = volcano.regulation.color[c(1,3)]) +
        scale_size_manual(values = volcano_size) +
        scale_alpha_manual(values = volcano_alpha) +
        
        scale_x_continuous(breaks = c(-4, -2, 0, 2, 4), limits = c(-4, 4), expand = c(0, 0)) +
        scale_y_continuous(breaks = c(0, 2, 4, 6, 8), limits = c(0, 8), expand = c(0,0)) +
        
        
        labs(x = "Log<sub>2</sub>(Moderate- / Low-volume)", 
             y = "") +
        
        dissertation_theme() +
        
        theme(panel.grid = element_blank(), 
              panel.border = element_blank(), 
              axis.line.y  = element_line(size = line.size), 
              axis.line.x = element_line(size = line.size), 
              axis.ticks.y = element_line(size = line.size), 
              axis.ticks.x = element_line(size = line.size), 
              axis.text.y =  element_markdown(color = "black", size = 7), 
              axis.title.y = element_markdown(color = "black", size = 7),
              axis.text.x =  element_text(color = "black", size = 7), 
              axis.title.x = element_blank(),
              legend.title = element_blank(), 
              legend.background = element_rect(fill = "white"),
              legend.margin = margin(t = 0, r = 1, b = 1, l = 1, unit = "pt"),
              legend.key = element_rect(fill = "white"),
              legend.position = "none")


week12_tissue <- volcano_data %>%
        filter(model == "Per tissue mass +\nEffective library size" & coef == "Week 12: Multiple - Single-set") %>%
        
        ggplot(aes(log2fc, -log10(p.val), fill = regulation, alpha = regulation, size = regulation)) + 
        geom_point(shape = 21) +
        
        annotate("text", x = -3.5, y = 7.7, label = "Tissue-offset", hjust = 0, size = 2.2) +
        
        # Scales
        scale_color_manual(values = volcano.regulation.color) + 
        scale_fill_manual(values = volcano.regulation.color) +
        scale_size_manual(values = volcano_size) +
        scale_alpha_manual(values = volcano_alpha) +
        
        scale_x_continuous(breaks = c(-4, -2, 0, 2, 4), limits = c(-4, 4), expand = c(0, 0)) +
        scale_y_continuous(breaks = c(0, 2, 4, 6, 8), limits = c(0, 8), expand = c(0,0)) +
        
        
        labs(x = "", 
             y = "-Log<sub>10</sub>(<em>P</em>-value)") +
        
        dissertation_theme() +
        
        theme(panel.grid = element_blank(), 
              panel.border = element_blank(), 
              axis.line.y  = element_line(size = line.size), 
              axis.line.x = element_line(size = line.size), 
              axis.ticks.y = element_line(size = line.size), 
              axis.ticks.x = element_line(size = line.size), 
              axis.text.y =  element_markdown(color = "black", size = 7), 
              axis.title.y = element_markdown(color = "black", size = 7),
              axis.text.x =  element_text(color = "black", size = 7), 
              axis.title.x = element_blank(),
              legend.title = element_blank(), 
              legend.background = element_rect(fill = "white"),
              legend.margin = margin(t = 0, r = 1, b = 1, l = 1, unit = "pt"),
              legend.key = element_rect(fill = "white"),
              legend.position = "none")




######### Volcano and upset configuration 

week2_grid <-  plot_grid(plot_grid(week2_tissue, week2_libsize, week2_naive, ncol = 3, align = "h"),
                         NULL,
                         plot_grid(NULL, upset_w2pre, NULL, ncol = 3, rel_widths = c(0.01, 0.6, 0.01)), 
                         rel_heights  = c(0.4,0.05, 0.6), nrow = 3)


week12_grid <-  plot_grid(plot_grid(week12_tissue, week12_libsize, week12_naive,   ncol = 3, align = "h"),
                          
                          plot_grid(NULL, upset_w12, NULL,  ncol = 3, rel_widths = c(0.01, 0.6, 0.01)), 
                          rel_heights = c(0.4, 0.6), nrow = 2)



## Save week 2 grid #################

saveRDS(week2_grid, "./data/derivedData/study1-transcriptome/week2_volcano_upset.RDS")

saveRDS(week12_grid, "./data/derivedData/study1-transcriptome/week12_volcano_upset.RDS")














