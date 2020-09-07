######### Figure comparison between DXA and MR and example image of MR

# Themes contain all figure options
source("./R/themes.R")


# Load data

source("./R/study-1/mr-and-dxa-data.R")

legs <- read_csv2("./data/study-1/body-composition/oneThreeSetLeg.csv") %>%
        pivot_longer(names_to = "set", values_to = "leg", cols = multiple:single) %>%
        print()


# Agreement between DXA and MR 
dxa %>%
        pivot_wider(names_from = timepoint, values_from = g) %>%
        mutate(dxa_rel_change = post/pre) %>%
        
        print()

comb_dat <- mr.results %>%
        
        dplyr::select(subject, leg, timepoint, CSA.avg) %>%
        pivot_wider(names_from = timepoint, values_from = CSA.avg) %>%
        mutate(mr_rel_change = 100 * ((post/pre)-1)) %>%
        dplyr::select(subject, leg, mr_rel_change) %>%
        inner_join(dxa %>%
                           pivot_wider(names_from = timepoint, values_from = g) %>%
                           mutate(dxa_rel_change = 100 * ((post/pre)-1)) %>%
                           dplyr::select(subject, leg, dxa_rel_change)) %>%
        inner_join(legs) %>%
        filter(include == "incl") %>%
        print()


rbind(mr.results %>%

        dplyr::select(subject, leg, timepoint, CSA.avg) %>%
        pivot_wider(names_from = leg, values_from = CSA.avg) %>%
              mutate(method = "MRI"), 
        dxa %>%
                           pivot_wider(names_from = leg, values_from = g) %>%
              mutate(method = "DXA")) %>%
        filter(timepoint == "pre") %>%
        inner_join(legs %>%
                           filter(leg == "R") %>%
                           dplyr::select(subject, include)) %>%
        filter(include == "incl") %>%
        rowwise() %>%
        print()
        
        
        mutate(diff = R - L, 
               m = mean(c(R, L)), 
               rel.diff = abs(diff)/m) %>%
        dplyr::select(subject, method, rel.diff) %>%
        pivot_wider(names_from = method, values_from = rel.diff) %>%
         ggplot(aes(MRI, DXA)) + geom_point()
       print()

       cor.test(comb_dat$mr_rel_change, comb_dat$dxa_rel_change)
               


corr_plot <- comb_dat %>%
        
        ggplot(aes(mr_rel_change, dxa_rel_change)) + 
        geom_abline(intercept = 0, slope = 1) +
        geom_point() + 
        scale_x_continuous(limits = c(-10, 15), expand = c(0, 0)) +
        scale_y_continuous(limits = c(-10, 15), expand = c(0, 0)) + 
        
        annotate("text", x = -5, y = -3, label = "y = x", 
                 size = 2.5, angle = 45) +
        
        labs(x = "Percentage change MRI", 
             y = "Percentage change DXA") +

        dissertation_theme()
        

# Calculate limits of agreement by crude method
# Resource for calculations: https://www-users.york.ac.uk/~mb55/meas/glucose.htm
# and @RN2540.

stats <- comb_dat %>%
        rowwise() %>%
        mutate(m = mean(c(mr_rel_change, dxa_rel_change)), 
               diff = mr_rel_change - dxa_rel_change) %>%
        ungroup() %>%
        summarise(m.diff = mean(diff), 
                  sd.diff = sd(diff)) %>%
        mutate(upper.loa = m.diff + 1.96 * sd.diff, 
               lower.loa = m.diff - 1.96 * sd.diff) %>%
        print()

# Regression methods for the LOA

diff_data <- comb_dat %>%
        rowwise() %>%
        mutate(m = mean(c(mr_rel_change, dxa_rel_change)), 
               diff = mr_rel_change - dxa_rel_change)


m1 <- lm(diff ~ m, diff_data)

summary(m1)

# Calculate the absolute 

diff_data$fitted <- fitted(m1)

m2 <- diff_data %>%
        mutate(resid = abs(diff - fitted)) %>%
        lm(resid ~ m, data = .)


# Take the coefficients and multiply by sqrt(pi/2)
# this is the error term (multiply by 1.96)

m1.coefs <- coef(m1)
m2.coefs <- coef(m2) * sqrt(pi/2)
 
model.loa <- data.frame(lower.loa = m1.coefs[1] + m1.coefs[2] * seq(from = -5, to = 15, by = 5) - 
        1.96 * (m2.coefs[1] + m2.coefs[2] * seq(from = -5, to = 15, by = 5) ), 
        upper.loa = m1.coefs[1] + m1.coefs[2] * seq(from = -5, to = 15, by = 5)  + 
                1.96 * (m2.coefs[1] + m2.coefs[2] * seq(from = -5, to = 15, by = 5) ), 
        m = seq(from = -5, to = 15, by = 5) , 
        bias = m1.coefs[1] + m1.coefs[2] * seq(from = -5, to = 15, by = 5) )


loa_plot <- comb_dat %>%
        rowwise() %>%
        mutate(m = mean(c(mr_rel_change, dxa_rel_change)), 
               diff = mr_rel_change - dxa_rel_change) %>%
        ungroup() %>%
        ggplot(aes(m, diff)) + 
        geom_point() +
        
        # Model based limits of agreement
        geom_line(data = model.loa, aes(x = m, y = lower.loa), color = "gray30") +
        geom_line(data = model.loa, aes(x = m, y = upper.loa), color = "gray30") +
        geom_line(data = model.loa, aes(x = m, y = bias), color = "gray30") +
        
        # Crude limits of agreement
        geom_hline(yintercept = stats$upper.loa, lty = 2, color = "gray70") +
        geom_hline(yintercept = stats$lower.loa, lty = 2, color = "gray70") +
        geom_hline(yintercept = stats$m.diff, lty = 2, color = "gray70") +
        

        
        labs(x = "Mean percentage change", 
             y  = "MRI - DXA") +
        
        annotate("text", x = -2, y = 10.2, label = "Crude 95% LOA", 
                 size = 2.5, color = "gray70") +
        annotate("text", x = 1, y = 13, label = "Regression method, 95% LOA", 
                 size = 2.5, color = "gray30") +
        
        annotate("text", x = 12, y = 2.5, label = "Crude bias", 
                 size = 2.5, color = "gray70") +
        annotate("text", x = 13, y = -1, label = "Bias (regression method)", 
                 size = 2.5, color = "gray30") +

        dissertation_theme() 
        


### Draw example images 

pre.img <- ggdraw() + 
        draw_image("./data/study-1/body-composition/pre.png")
post.img <- ggdraw() + 
        draw_image("./data/study-1/body-composition/post.png")
comp.img <- ggdraw() + 
        draw_image("./data/study-1/body-composition/comparison.png")


pre_title <- ggdraw() + draw_label("Pre-\ntraining", 
                                   angle = 0, 
                                   size = 8, vjust = 1.5, hjust = 0.4)
post_title <- ggdraw() + draw_label("Post-\ntraining", 
                                    angle = 0, 
                                    size = 8, vjust = 1.5, hjust = 0.4)
comp_title <- ggdraw() + draw_label("\nPost - Pre", 
                                    angle = 0, 
                                    size = 8, vjust = 1.5, hjust = 0.4)




mr <- plot_grid(NULL, 
                plot_grid(pre_title, pre.img, ncol = 1,   rel_heights = c(0.2, 1)),
                plot_grid(post_title, post.img, ncol = 1, rel_heights = c(0.2, 1)),
                plot_grid(comp_title, comp.img, ncol = 1, rel_heights = c(0.2, 1)),
                NULL, ncol = 5, rel_widths = c(0.01, 0.33, 0.33, 0.33, 0.01))





#### Combine plot 

mr_dxa_fig <- plot_grid(mr, 
          plot_grid(corr_plot, loa_plot, ncol = 2), 
          rel_heights = c(1, 1), ncol = 1) + 
        draw_plot_label(label=c("A", "B","C"),
                        x = c(0.02, 0.02, 0.51), 
                        y = c(0.98, 0.49, 0.49),
                        hjust=.5, vjust=.5, 
                        size = label.size)


saveRDS(mr_dxa_fig, "./figures/methods/mr-dxa.RDS")







