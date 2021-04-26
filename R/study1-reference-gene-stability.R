##-------------------------------------
## study1-reference-gene-stability.R
##
## Title: Reference gene stability in Study I
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


#library(devtools)
#install_github("dhammarstrom/generefer")

library(generefer)

# Load data, average expression values (Cq) for each target, sample, timepoint and condition 
qpcrdat <- read_feather("./data/study-1/qpcr/qpcrdat_avg")

# load data with subject/study info
qpcrdat <- read.csv("./data/study-1/oneThreeSetLeg.csv", sep = ";") %>%
        gather(sets, leg, multiple:single) %>%
        inner_join(qpcrdat) %>%
        filter(include == "incl") %>%
        filter(!(target %in%  c("Lambda MIX", "Lambda Kit", "PPIA F1R1"))) %>%
        dplyr::select(subject, sex, sets, timepoint, target, cq, eff)














qpcrdat$target.symbol <- str_split_fixed(qpcrdat$target, " ", 2)[,1] # make a gene symbol without primer id


### Descriptive figure 

gene_descriptive <- qpcrdat %>%
        group_by(target) %>%
        mutate(sample = paste(subject, sets, timepoint, sep="_"),
               expression = log(eff^-cq)) %>%
        filter(!(target %in% c("IGF1 F11R11", "MSTN F2R2"))) %>%
        
        separate(target, into = c("target", "primer"), sep = " ") %>%
        mutate(target = fct_reorder(target, expression, .fun = mean)) %>%

        
        ggplot(aes(target, expression, fill = target)) + geom_violin() +
        ylab("Raw abundance<br>(Efficiency<sup>-Cq</sup>)") +
        scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
        scale_fill_discrete() +
        dissertation_theme() + 
        theme(axis.title.x = element_blank(), 
              axis.title.y = element_markdown(), 
              legend.position = "none")



saveRDS(gene_descriptive, "./data/derivedData/study1-reference-gene-stability/gene_descriptive.RDS")







# Determine stable reference genes using genorm 
genorm.results <- qpcrdat %>%
        group_by(target) %>%
        mutate(sample = paste(subject, sets, timepoint, sep="_"),
               expression = (eff^-cq)/max((eff^-cq))) %>%
        filter(!(target %in% c("IGF1 F11R11", "MSTN F2R2"))) %>%
        data.frame() %>%
        tidy_genorm(sample = "sample", gene = "target.symbol", expression = "expression", log = FALSE)

# 3 selected genes in genorm
gn.results <- genorm.results$M[order(genorm.results$M$M),]
selected.genes.gn <- gn.results[c(1:2), 1]


# geNorm plots
genorm.results$M$gene <- str_split(genorm.results$M$gene, " ", simplify = TRUE)[,1]

labels <- c(genorm.results$M$gene[1:9], paste(genorm.results$M$gene[10:11], collapse = "\n"))

m.plot <- genorm.results$M %>%
        ggplot(aes(n.genes, M.avg)) + geom_point() + geom_line() +
        scale_x_reverse(lim=c(11,2), 
                        breaks = c(11,10,9,8,7,6,5,4,3,2), 
                        labels = labels, 
                        guide = guide_axis(n.dodge = 2)) +

        ylab("geNorm Average\nexpression-stability") +

        dissertation_theme() + 
        theme(axis.title.x = element_blank())

        
        

v.plot <- genorm.results$pairwise.variations %>%
        mutate(NFpair = factor(NFpair, levels = c("V2/3","V3/4","V4/5","V5/6","V6/7","V7/8","V8/9","V9/10","V10/11")))%>%
        ggplot(aes(NFpair, pairwise.variation)) + geom_bar(stat = "identity") +
        scale_y_continuous(lim = c(0,0.2)) +
        geom_hline(yintercept = 0.15, color = "gray40", linetype = 2) +
        xlab("Normalization-factor pairs") + ylab("geNorm pairwise\nvariation (V)") +
        scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
        dissertation_theme() + 
        theme(axis.title.x = element_blank())


normfinder.results <- qpcrdat %>%
        mutate(sample = paste(subject, sets, timepoint, sep="_"), 
               group = paste(sets, timepoint, sep="_"), 
               expression = (eff^-cq) * 10^8) %>%
        filter(!(target %in% c("IGF1 F11R11", "MSTN F2R2", "Lambda MIX", "Lambda Kit"))) %>%
        data.frame() %>%
        normfinder(sample = "sample", gene = "target.symbol", cq = "cq", ctVal = TRUE, pStabLim = 0.15)

nf.results <- normfinder.results[[3]][order(normfinder.results[[3]][,3]),]

selected.genes.nf <- c(as.character(nf.results[1,1]),  as.character(nf.results[1,2]))

normfinder.ranking <- normfinder.results[[1]]

normfinder.ranking$Gene <- rownames(normfinder.ranking)

rownames(normfinder.ranking) <- NULL

normfinder.ranking <- normfinder.ranking %>%
        dplyr::select(Gene, GroupDif, GroupSD, Stability)


## Normfinder plots
nf.gene.plot <- normfinder.results[[1]] %>%
        mutate(Gene = rownames(.)) %>%
        arrange(-Stability) %>%
        mutate(Gene = factor(Gene, Gene)) %>%
        ggplot(aes(Gene, Stability)) + geom_bar(stat = "identity") +
        ylab("Normfinder\nexpression-stability") +
        scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
        dissertation_theme() + 
        theme(axis.title.x = element_blank())


nf.pair.plot <- normfinder.results[[3]] %>%
        mutate(Gene_pair = paste(Gene1,"\n", Gene2)) %>%
        arrange(-Stability) %>%
        mutate(Gene_pair = factor(Gene_pair, Gene_pair)) %>%
        ggplot(aes(Gene_pair, Stability)) + geom_bar(stat = "identity") +
        xlab("Gene-pair") + 
        ylab("Normfinder pairwise\nexpression stability") +
        scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
        dissertation_theme() + 
        theme(axis.title.x = element_blank())




selection_process <- plot_grid(m.plot, v.plot, nf.gene.plot, nf.pair.plot, ncol = 1, align = "h")

saveRDS(selection_process, "./data/derivedData/study1-reference-gene-stability/selection_process.RDS")





# Normalization factors variation over design variables

normalization.factors <- read_feather("./data/study-1/qpcr/reference-gene-selections/normalization_factors.feather") %>%
        group_by(nf.id) %>%
        mutate(centered.nf = nf-mean(nf, na.rm = TRUE), 
               sample = paste0(subject, timepoint, sets), 
               sets = factor(sets, levels = c("single", "multiple"))) %>%
        
        ungroup() 


nf.summary <- normalization.factors %>%
        group_by(normalization) %>%
        mutate(centered.expr = avg.expr-mean(avg.expr)) %>%
        group_by(normalization, timepoint, sets) %>%
        summarise(centered.expr = mean(centered.expr), 
                  centered.nf = mean(centered.nf), 
                  nf = mean(nf)) 


#### Models normalization factors ##########
gapdh <- lme(nf ~ timepoint + timepoint:sets, 
             random = list(subject = ~ 1, sample = ~ 1),
              data = normalization.factors[normalization.factors$nf.id == "gapdh",])

gn <- lme(nf ~ timepoint + timepoint:sets, 
          random = list(subject = ~ 1, sample = ~ 1),
           data = normalization.factors[normalization.factors$nf.id == "gn",])

mm <- lme(nf ~ timepoint + timepoint:sets, 
          random = list(subject = ~ 1, sample = ~ 1),
           data = normalization.factors[normalization.factors$nf.id == "mm",])

nf <- lme(nf ~ timepoint + timepoint:sets, 
          random = list(subject = ~ 1, sample = ~ 1),
           data = normalization.factors[normalization.factors$nf.id == "nf",])


coefs <- rbind(intervals(nf)$fixed %>%
        data.frame() %>%
        mutate(coef = rownames(.), 
               method = "nf"), 
      intervals(gn)$fixed %>%
              data.frame() %>%
              mutate(coef = rownames(.), 
                     method = "gn"), 
      intervals(mm)$fixed %>%
              data.frame() %>%
              mutate(coef = rownames(.), 
                     method = "mm"), 
      intervals(gapdh)$fixed %>%
              data.frame() %>%
              mutate(coef = rownames(.), 
                     method = "gapdh"))  %>%
        
        filter(coef != "(Intercept)") %>%
        mutate(coef = factor(coef, 
                             levels = c("timepointw2pre", 
                                        "timepointw2post",
                                        "timepointw12", 
                                        "timepointw0:setsmultiple", 
                                        "timepointw2pre:setsmultiple", 
                                        "timepointw2post:setsmultiple", 
                                        "timepointw12:setsmultiple"), 
                             labels = c("Week 2 Pre-ex\nvs. Week 0", 
                                        "Week 2 Post-ex\nvs. Week 0", 
                                        "Week 12\nvs. Week 0", 
                                        "Mod vs. Low\nat Week 0", 
                                        "Mod vs. Low\nat Week 2 Pre-ex", 
                                        "Mod vs. Low\nat Week 2 Post-ex", 
                                        "Mod vs. Low\nat Week 12")), 
               coef = fct_rev(coef),
               method = factor(method, levels = c("mm", "nf", "gn", "gapdh"), 
                               labels = c("Mixed-model", 
                                          "Normfinder", 
                                          "geNorm", 
                                          "GAPDH"))) %>%
        
        mutate(est. = exp(est.), 
               lower = exp(lower), 
               upper = exp(upper)) %>%

        ggplot(aes(coef, est., shape = method, fill = method, color = method)) + 
        geom_hline(yintercept = 1, color = "gray50", lty = 2) +
        geom_point(position = position_dodge(width = 0.5)) +
        geom_errorbar(aes(ymin = lower, ymax = upper), 
                       width = 0, 
                      position = position_dodge(width = 0.5)) +
        
        scale_color_manual(values = group.study.color[c(1, 6, 2, 5)]) + 
        scale_fill_manual(values = group.study.color[c(1, 6, 2, 5)]) + 
        scale_shape_manual(values = c(21, 22, 23, 24)) +
        scale_y_continuous(limits = c(0.85, 1.15), 
                           breaks = c(0.9, 0.95, 1, 1.05, 1.1)) +
        dissertation_theme() +
        
        annotate("text", 
                         x = c(7.36, 7.24, 7.1, 6.7), 
                 y = c(0.92, 1.01, 1.03, 1.02),
                 color = group.study.color[c(5, 2, 6, 1)],
                 hjust = 0,
                 label = c("GAPDH", "geNorm", "Normfinder", "Mixed-model"), 
                 size = 2.5) +
        
      #  annotate("curve",
      #           x = 1.55,
      #           xend = 1.85,
      #           y = 18,
      #           yend = 14,
      #           curvature = 0.2,
      #           color = "gray50",
      #           arrow = arrow(length=unit(0.1,"cm"), type = "closed")) +
      #  
        labs(y = "Fold change") +
        theme(legend.position = "none", 
              axis.title.y = element_blank()) +
        
        coord_flip()
        
    

nf.summary <- normalization.factors %>%
        mutate(timepoint = factor(timepoint, levels=c("w0", "w2pre", "w2post", "w12"), 
                                  labels = c("Week 0", "Week 2 \n Pre", "Week 2 \n Post", "Week 12")),
               sets = factor(sets, levels = c("single", "multiple"), 
                             labels = c("Low-<br>volume", "Moderate-<br>volume")),
               normalization = factor(nf.id, levels = c("mm", "nf", "gn", "gapdh"), 
                                      labels = c("Mixed-model selection:<br>B2m + CHMP2A", 
                                                 "Normfinder:<br>B2m + RPL32", 
                                                 "geNorm:<br>RPL32 + TBP", 
                                                 "Single gene:<br>GAPDH"))) %>%
        group_by(normalization) %>%
        mutate(centered.expr = avg.expr-mean(avg.expr)) %>%
        group_by(normalization, timepoint, sets) %>%
        summarise(centered.expr = mean(centered.expr), 
                  centered.nf = mean(centered.nf), 
                  nf = mean(nf)) 


    


reference_genes <- normalization.factors %>%
        mutate(timepoint = factor(timepoint, levels=c("w0", "w2pre", "w2post", "w12"), 
                                  labels = c("Week 0", "Week 2 \n Pre", "Week 2 \n Post", "Week 12")),
               sets = factor(sets, levels = c("single", "multiple"), 
                             labels = c("Low-<br>volume", "Moderate-<br>volume")),
               normalization = factor(nf.id, levels = c("mm", "nf", "gn", "gapdh"), 
                                      labels = c("Mixed-model selection:<br>B2m + CHMP2A", 
                                                 "Normfinder:<br>B2m + RPL32", 
                                                 "geNorm:<br>RPL32 + TBP", 
                                                 "Single gene:<br>GAPDH"))) %>%
        group_by(timepoint, subject, sets, normalization) %>%
        summarise(avg.expr = mean(avg.expr, na.rm = TRUE),
                  nf = mean(nf, na.rm = TRUE)) %>%
        group_by(normalization) %>%
        mutate(centered.expr = avg.expr-mean(avg.expr),
               centered.nf = nf - mean(nf)) %>%
        ungroup() %>%
        
        ggplot(aes(timepoint, nf, group = paste0(subject, sets))) + 
        
        geom_line(alpha = 0.1) + 
        
        geom_line(data = nf.summary, aes(timepoint, nf, group = sets), 
                  position = position_dodge(width = 0.15)) +
        
        geom_point(data = nf.summary, size = 3, 
                   alpha = 0.8,
                   aes(fill = sets, group = sets, shape = sets), 
                   position = position_dodge(width = 0.15)) + 
        
        scale_fill_manual(values = c(group.study.color[1], group.study.color[5])) +
        scale_shape_manual(values = c(24, 23)) +
        
        facet_wrap(~normalization, ncol = 1) +
        
        dissertation_theme() + 
        
        theme(strip.background = element_blank(), 
              strip.text = element_markdown(size = 8, hjust = 0), 
              legend.text = element_markdown(size = 7),
              legend.key.height = unit(0.18, "cm"),
              legend.title = element_blank(),
              axis.title.x = element_blank(), 
              axis.title.y = element_blank(),
              legend.position = "bottom")
        
        
        
        
        
ref_set_comp <- plot_grid(reference_genes,  coefs, ncol = 2, rel_widths  = c(1, 1.2))

saveRDS(ref_set_comp, "./data/derivedData/study1-reference-gene-stability/ref_set_comp.RDS")



