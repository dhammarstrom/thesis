##-------------------------------------
## study2-ubf-rps6-rna.R
##
## Title: UBF an d
## Purpose: 
## Author:
##
##
##
##-------------------------------------
## Notes: This script makes a figure supporting information from Table 
# Effect of UBF and rpS6 levels, sessions and de-training on RNA-levels
#
#
#
#
#
#
#
#
## ------------------------------------



comb.s6.m1 <- readRDS("./data/study-2/ubf-tot-rna-model/comb_m1_s6.RDS")


comb.m1 <- readRDS("./data/study-2/ubf-tot-rna-model/comb_m1.RDS")

ubf_cond <- conditional_effects(comb.m1)
rps6_cond <- conditional_effects(comb.s6.m1)



ubf_fig <- ubf_cond$scaled_ubf %>%
        ggplot(aes(scaled_ubf, exp(estimate__))) + 
        
        geom_ribbon(aes(ymin = exp(lower__), ymax = exp(upper__)), 
                    fill = group.study.color[4]) +
        
        
        geom_point(data = comb.m1$data, 
                   aes(scaled_ubf, exp(log.rna)),
                   shape = 21, fill = group.study.color[5]) +
        
        scale_y_continuous(limits = c(195, 900), 
                           breaks = c(300, 500, 700, 900), 
                           expand = c(0,0)) +
        
        labs(x = "UBF protein (SD units)", 
             y = "Total RNA (ng mg<sup>-1</sup>)") +
        
        geom_line() +
        dissertation_theme() +
        theme(axis.title.y = element_markdown())


rps6_fig <- rps6_cond$scaled_rps6 %>%
        ggplot(aes(scaled_rps6, exp(estimate__))) + 
               
        
        geom_ribbon(aes(ymin = exp(lower__), ymax = exp(upper__)), 
                    fill = group.study.color[4]) +
        geom_line() +
        scale_y_continuous(limits = c(195, 900), 
                           breaks = c(300, 500, 700, 900), 
                           expand = c(0,0)) +
        
        geom_point(data = comb.s6.m1$data, 
                   aes(scaled_rps6, exp(log.rna)), 
                   shape = 21, fill = group.study.color[5]) +
        

        labs(x = "rpS6 protein (SD units)", 
             y = "") +
        
        
        dissertation_theme() +
        theme(axis.ticks.y = element_blank(), 
              axis.title.y = element_blank(), 
              axis.text.y = element_blank(),
              axis.line.y = element_blank())
        


ubf_rps6_fig <- plot_grid(ubf_fig, rps6_fig, 
                          rel_widths = c(1, 0.8)) + 
        draw_plot_label(x = c(0, 0.52), 
                        y = c(1, 1), 
                        label = c("a", "b"))




saveRDS(ubf_rps6_fig, "./data/derivedData/study2-ubf-rps6-rna/ubf-rps6-rna-fig.RDS")







   