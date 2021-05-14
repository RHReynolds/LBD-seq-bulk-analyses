#---Load Libraries ----------------------------------------------------------------------------------------------------------------------####
library(ggsci)
library(scales)
library(tidyverse)

#---Theme--------------------------------------------------------------------------------------------------------------------------------####

theme_rhr <-  theme_bw(base_family = "Helvetica") + 
  theme(panel.grid.major.x = element_blank(),
        legend.position = "top",
        strip.text = element_text(size = 8),
        axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(vjust = 0.6),
        axis.title = element_text(size = 10),
        panel.spacing = unit(0.1, "lines"))

#---Colour schemes-----------------------------------------------------------------------------------------------------------------------####

# Function for import
colour_schemes <- function(){
  
  library(ggsci)
  
  ## Cell class colours
  cell_class <- data.frame(group = c("Astrocyte", "Immune", "Oligo", "Vascular", "Neuron", "Unclassified"),
                           colour_hex = pal_npg("nrc", alpha = 1)(10)[c(4,2,1,6,3,10)])
  
  ## Disease colours
  # Green theme
  disease_greens <- data.frame(group = c("Control", "PD", "PDD", "DLB"),
                               colour_hex = c("#868686", "#c4ebe3", "#73cab9", "#19967d"))
  
  # Discrete theme
  disease_discrete <- data.frame(group = c("Control", "PD", "PDD", "DLB"),
                                 colour_hex = pal_jco("default", alpha = 0.9)(10)[c(3,9,1,2)])
  
  colour_schemes_list <- list(cell_class, disease_greens, disease_discrete)
  names(colour_schemes_list) <- c("cell_class", "disease_greens", "disease_discrete")
  
  return(colour_schemes_list)
  
  
}

#---Plot Colour schemes------------------------------------------------------------------------------------------------------------------####

colour_schemes_list <- colour_schemes()

colour_schemes_list %>% 
  lapply(., function(df){
    
    df$colour_hex %>% 
      as.character() %>% 
      scales::show_col()
    
  })
