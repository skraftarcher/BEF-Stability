#visualizing ecomem persistence temporally
source("scripts/install_packages_function.R")
lp("tidyverse")

#this function's input is one data set; the visual survey
empers.heat<-function(vis){
  #convert wide to long format
  vislong <- vis %>%
    pivot_longer(cols = -c(1:4), #skip metadata columns
                 names_to = "taxa",
                 values_to = "abund") %>%
    mutate(pa = ifelse(abund == 0, 0, 1)) #create a presence/absence column
  
  #basic heatmap
  persistence_plot <- ggplot(vislong, aes(x = sampling, y = taxa,
                                          fill = factor(pa))) +
    geom_tile(color = "white", linewidth = 0.5) +
    facet_wrap(~blockID, ncol = 4) +
    scale_fill_manual(values = c("0" = "gray90", "1" = "#2c7bb6"),
                      labels = c("Absent", "Present"),
                      name = "Status") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 2.5)) +
    labs(x = "Sampling Period", 
         y = "Species", 
         title = "Species Temporal Persistence by Block")
  print(persistence_plot)
  
  #save the plot
  ggsave("ecomem_species_persistence_heatmap.png",
         plot = persistence_plot,
         width - 12,
         height = 8,
         dpi = 300)
}
