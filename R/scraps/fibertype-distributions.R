


source("./R/libraries.R")
source("./R/themes.R")

# Sample set up (included samples)

legs <- read_csv2("./data/study-1/body-composition/oneThreeSetLeg.csv") %>%
        pivot_longer(names_to = "sets", values_to = "leg", cols = multiple:single) %>%
        print()


# Total number of fibers per sample

low_fiber_counts <- read_excel("./data/study-1/fiber-type/fibertype.xlsx") %>%
        inner_join(legs) %>%
        filter(include == "incl") %>%
        mutate(total = type1 + type2a + type2ax + type2x) %>%
        filter(total < 200) %>%
        print()
        
total_histogram <- read_excel("./data/study-1/fiber-type/fibertype.xlsx") %>%
        inner_join(legs) %>%
        filter(include == "incl") %>%
        mutate(total = type1 + type2a + type2ax + type2x) %>%        
        
        ggplot(aes(total)) + geom_histogram(binwidth = 30, 
                                            fill = "blue", alpha = 0.4) + 
        dissertation_theme() + 
        scale_x_continuous(limits = c(0, 1600),
                           expand = c(0,0),
                           breaks = c(0, 200, 400, 600, 800, 1000, 1200, 1400, 1600), 
                           labels = c(0, "", 400, "", 800, "", 1200, "", 1600)) + 
        scale_y_continuous(limits = c(0, 20), 
                           breaks = c(0, 5, 10, 15, 20), 
                           expand = c(0,0)) +
        
        labs(x = "Total number of fibers per sample", 
             y = "Counts") + 
        theme(axis.text.x = element_text(hjust = c(0, rep(0.5, 7), 1)))
        
### Comparison between "duplicate samples"

corr_plot_fiber_types <- read_excel("./data/study-1/fiber-type/fibertype.xlsx") %>%
        inner_join(legs) %>%
        filter(include == "incl") %>%
        mutate(total  = type1 + type2a + type2ax + type2x, 
               
               typeI = (type1 / total)*100,
               typeII = ((type2a + type2ax +type2x)/total) *100,
               
               type1  =(type1 / total)*100, 
               type2a =(type2a / total)*100, 
               type2ax=(type2ax / total)*100, 
               type2x =(type2x / total)*100 ) %>% 
    
        filter(timepoint %in% c(2, 3)) %>%
        dplyr::select(subject, timepoint, leg,total, type1:typeII) %>%
        pivot_longer(names_to = "type", values_to = "percentage", cols = c(type1:type2x, typeI:typeII)) %>%
        mutate(timepoint = paste0("t", timepoint)) %>%
        pivot_wider(names_from = timepoint, values_from = c(percentage, total)) %>%

        mutate(type = factor(type, levels = c("type1", 
                                              "type2a", 
                                              "type2ax", 
                                              "type2x", 
                                              "typeI", 
                                              "typeII"), 
                             labels = c("Type I", 
                                        "Type IIA", 
                                        "Type IIAX", 
                                        "Type IIX", 
                                        "Type I_tot", 
                                        "Type II_tot"))) %>%
        
        
        ggplot(aes(percentage_t2, percentage_t3)) +  geom_point(alpha = 0.4, 
                                                                 color = "blue") + 
        facet_wrap(~type, ncol = 2, strip.position = "top")  +

        scale_x_continuous(limits = c(0, 100)) + 
        scale_y_continuous(limits = c(0, 100)) + 
        
        geom_abline(intercept = 0, slope = 1) + 
        dissertation_theme() + 
        theme(strip.background = element_blank(), 
              panel.border = element_rect(color = "black", fill = NA))
        


plot_grid(total_histogram, corr_plot_fiber_types, ncol = 1, rel_heights = c(1, 2)) + 
        draw_plot_label(label=c("A", "B"),
                        x = c(0.02, 0.02), 
                        y = c(0.98, 0.6),
                        hjust=.5, vjust=.5, 
                        size = label.size)






read_excel("./data/study-1/fiber-type/fibertype.xlsx") %>%
    inner_join(legs) %>%
    filter(include == "incl") %>%
    mutate(total  = type1 + type2a + type2ax + type2x, 
           
           typeI = (type1 / total)*100,
           typeII = ((type2a + type2ax +type2x)/total) *100,
           
           type1  =(type1 / total)*100, 
           type2a =(type2a / total)*100, 
           type2ax=(type2ax / total)*100, 
           type2x =(type2x / total)*100 ) %>% 
    
    filter(timepoint %in% c(1)) %>%
    dplyr::select(subject, timepoint, leg,total, type1:type2x, typeI:typeII) %>%
    pivot_longer(names_to = "type", values_to = "percentage", cols = c(type1:type2x, typeI:typeII)) %>%
    mutate(timepoint = paste0("t", timepoint)) %>%
    pivot_wider(names_from = leg, values_from = c(percentage, total)) %>%
    
    mutate(type = factor(type, levels = c("type1", 
                                          "type2a", 
                                          "type2ax", 
                                          "type2x", 
                                          "typeI", 
                                          "typeII"), 
                         labels = c("Type I", 
                                    "Type IIA", 
                                    "Type IIAX", 
                                    "Type IIX", 
                                    "Type I_tot", 
                                    "Type II_tot"))) %>%

        rowwise() %>%
        mutate(mean.perc = mean(c(percentage_L, percentage_R)), 
               diff.perc = percentage_L - percentage_R, 
               diff.sq = (percentage_L - percentage_R)^2,
               diff.mean.sq = (diff.perc/mean.perc)^2,
               min.count = min(c(total_L, total_R))) %>%
        filter(mean.perc > 0) %>%
        
        group_by(type) %>%
        
        summarise(m = mean(mean.perc), 
                  S2 = sum(diff.sq)/(2*n()), 
                  S = sqrt(S2), 
                  RMS = 100 * sqrt(sum(diff.mean.sq) / (2*n())), 
                  ws = 100 * (sqrt(sum(diff.perc^2) / (2*n()) ))/m) %>%
        mutate(cv  = 100 * (S/m)) %>%
        
        print()
        
        
        
 
 read_excel("./data/study-1/fiber-type/fibertype.xlsx") %>%
         inner_join(legs) %>%
         filter(include == "incl") %>%
         mutate(total  = type1 + type2a + type2ax + type2x, 
                type1  =(type1 / total)*100, 
                type2a =(type2a / total)*100, 
                type2ax=(type2ax / total)*100, 
                type2x =(type2x / total)*100) %>%  
         filter(timepoint %in% c(2, 3)) %>%
         dplyr::select(subject, timepoint, leg,total, type1:type2x) %>%
         pivot_longer(names_to = "type", values_to = "percentage", cols = type1:type2x) %>%
         mutate(timepoint = paste0("t", timepoint)) %>%
         pivot_wider(names_from = timepoint, values_from = c(percentage, total)) %>%
         
         mutate(type = factor(type, levels = c("type1", 
                                               "type2a", 
                                               "type2ax", 
                                               "type2x"), 
                              labels = c("Type I", 
                                         "Type IIA", 
                                         "Type IIAX", 
                                         "Type IIX"))) %>%
         
         rowwise() %>%
         mutate(mean.perc = mean(c(percentage_t2, percentage_t3)), 
                diff.perc = percentage_t2 - percentage_t3, 
                diff.sq = (percentage_t2 - percentage_t3)^2,
                diff.mean.sq = (diff.perc/mean.perc)^2,
                min.count = min(c(total_t2, total_t3))) %>%
         filter(mean.perc > 0, 
                type %in% c("Type I", "Type IIA")) %>%  
    
        ggplot(aes(mean.perc * , abs(diff.perc))) + geom_point(alpha = 0.4, 
                                                       color = "blue") + 
        facet_wrap(~type, ncol = 2, strip.position = "top")  +
        


        
        
        dissertation_theme() + 
        labs(x = "Minimum count of each pair", 
            y = "Absolute paired percentage-point difference") +
        theme(strip.background = element_blank()) 
         
        





        
        
print()



