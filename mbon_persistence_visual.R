#visualizing mbon persistence temporally
source("scripts/install_packages_function.R")
lp("tidyverse")

#this function's input is two data sets; the environment dataset and
#the abundance dataset
empers.heat<-function(ds.env,ds.abund){
  #add row numbers so we can join them
  ds.env <- ds.env %>%
    mutate(row_id = row_number())
  
  ds.abund <- ds.abund %>%
    mutate(row_id = row_number())
  
  #join the datasets
  mbon_combined <- ds.env %>%
    left_join(ds.abund, by = "row_id")
  
  #get the names of metadata columns
  metadata_cols <- names(ds.env)
  
  #transform to long format
  mbon_long <- mbon_combined %>%
    pivot_longer(cols = -all_of(metadata_cols),
                 names_to = "taxa",
                 values_to = "abund") %>%
    mutate(pa = ifelse(abund == 0, 0, 1)) #create presence/absence
  
  #see what sites exist and save to a list
  site_list <- mbon_long %>%
    distinct(site) %>%
    pull(site)
  
  #iterate through sites
  for(site_name in site_list) {
    site_data <- mbon_long %>%
      filter(site == site_name)
    
    #Filter to one site
    mbon_site <- mbon_long %>%
      filter(site == site_name)
    
    print(paste("Generating heatmap for site:", site_name))
    
    #basic heatmap
    mbon_plot <- ggplot(mbon_site, aes(x = sampling, y = taxa,
                                       fill = factor(pa))) +
      geom_tile(color = "white", linewidth = 0.5) +
      facet_wrap(~tray, ncol = 3) +
      scale_fill_manual(values = c("0" = "gray90", "1" = "#2c7bb6"),
                        labels = c("Absent", "Present"),
                        name = "Status") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            axis.text.y = element_text(size = 2.5)) +
      labs(x = "Sampling Period", 
           y = "Species", 
           title = "Taxa at", site_name)
    print(mbon_plot)
    
    #save the plot
    ggsave("mbon_species_persistence_heatmap_", site_name, ".png",
           plot = mbon_plot,
           width - 12,
           height = 8,
           dpi = 300)
  }
}
