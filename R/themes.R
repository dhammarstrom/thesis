

#### Figure theme/color/settings for dissertation  ####

library(ggplot2)


color_palette <-  



line.size <- 0.3


dissertation_theme <- function() {
        theme_bw() +
        theme(axis.title = element_text(size = 7),
              axis.text = element_text(size = 7, color = "black"),
              axis.line = element_line(color = "black", size = line.size),
              axis.ticks = element_line(color = "black", size = line.size),
              legend.background = element_blank(),
              panel.border = element_blank(),
              panel.grid = element_blank())
        
}




