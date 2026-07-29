library(tidyverse)
library(dplyr)
library(knitr)
library(cowplot)
library(grid)
#helper functions to specify values for specific drugs 

lower_limit <- function(drug_name){
  if (drug_name == "ABZ"){
    return(c(-200,50))
  }
  if (drug_name == "CLO"){
    return(c(-600,50))
  }
  if(drug_name == "DEC"){
    return(c(-300, 0))
  }
  if (drug_name == "EMO"){
    return(c(-700,50))
  }
  if (drug_name == "FBZ"){
    return(c(-300,50))
  }
  if (drug_name == "IVM"){
    return(c(-700, 100))
  }
  if (drug_name == "LEV"){
    return(c(-600,100))
  }
  if (drug_name == "MOX"){
    return(c(-800, 100))
  }
  if (drug_name == "PYR"){
    return(c(-300, 100))
  }
  if(drug_name == "MILB"){
    return(c(-700, 200))
  }
  if (drug_name == ""){
    return(c())
  }
}

drug_map <- function(drug_name){
  if (drug_name == "ABZ"){
    return("Albendazole")
  }
  if (drug_name == "CLO"){
    return("Closantel")
  }
  if(drug_name == "DEC"){
    return("DEC")
  }
  if (drug_name == "EMO"){
    return("Emodepside")
  }
  if (drug_name == "FBZ"){
    return("Fenbendazole")
  }
  if (drug_name == "IVM"){
    return("Ivermectin")
  }
  if (drug_name == "LEV"){
    return("Levamisole")
  }
  if (drug_name == "MOX"){
    return("Moxidectin")
  }
  if (drug_name == "PYR"){
    return("Pyrantel")
  }
  if (drug_name == "MILB"){
    return ("Milbemycin")
  }
}

axis_breaks <- function(drug_name){
  if (drug_name == "ABZ"){
    return(c(0.1170, 1.8750, 15.0000, 120.0000))
  }
  if (drug_name == "CLO"){
    return(c(0.097, 1.560, 12.500, 100))
  }
  if (drug_name == "DEC"){
    return(c(0.645, 10.300, 82.500, 660.000))
  }
  if (drug_name == "EMO"){
    return(c(0.00097, 0.01560, 0.12500, 1.00000))
  }
  if (drug_name == "FBZ"){
    return(c(0.117, 1.875, 15.000, 120.000))
  }
  if (drug_name == "IVM"){
    return(c(0.000029, 0.00046, 0.0037, 0.030))
  }
  if (drug_name == "LEV"){
    return(c(0.1360, 2.1875, 17.5000, 140.0000))
  }
  if (drug_name == "MOX"){
    return(c(0.00097, 0.01560, 0.12500, 1.00000))
  }
  if (drug_name == "PYR"){
    return(c(0.292, 4.680, 37.500, 300.000))
  }
  if (drug_name == "MILB"){
    return(c(0.00097, 0.01560, 0.12500, 1.00000))
  }
}

axis_labels <- function(drug_name){
  if (drug_name == "ABZ"){
    return(c("0.117", "1.875", "15", "120"))
  }
  if (drug_name == "CLO"){
    return(c("0.097", "1.56", "12.5", "100"))
  }
  if(drug_name == "DEC"){
    return(c("0.645", "10.3", "82.5", "660"))
  }
  if (drug_name == "EMO"){
    return(c("0.00097", "0.0156", "0.125", "1.0"))
  }
  if (drug_name == "FBZ"){
    return(c("0.117", "1.875", "15", "120"))
  }
  if (drug_name == "IVM"){
    return(c("0.000029", "0.00046", "0.0037", "0.03"))
  }
  if (drug_name == "LEV"){
    return(c("0.136", "2.1875", "17.5", "140"))
  }
  if (drug_name == "MOX"){
    return(c("0.00097", "0.0156", "0.125", "1"))
  }
  if (drug_name == "PYR"){
    return(c("0.292", "4.68", "37.5" , "300"))
  }
  if (drug_name == "MILB"){
    return(c("0.00097", "0.0156", "0.125", "1"))
  }
}


#read in the csv files
make_plot <- function(drug_name){
  
  directory <- "projects/Olivia/cbriggsae_drc/data/"
  drc3 <- "20250602_Cbriggsae_DRC3_"
  drc4 <- "20250619_Cbriggsae_DRC4"
  drugs_drc3 <- c("ABZ", "CLO", "DEC", "EMO", "FBZ", "IVM", "LEV", "PYR")
  drugs_drc4 <- c("MOX","MILB")
  file_name <- "/20260617_length_delta_nemasize.csv"
  
  if (drug_name %in% drugs_drc3){
    pathway <- paste0(directory,drc3,drug_name, file_name)
  }
  
  if (drug_name %in% drugs_drc4){
    pathway <- paste0(directory,drc4,drug_name, file_name)
  }
  
  length_delta <- read.csv(pathway)
  length_delta <- as.data.frame(length_delta)
  
  dose_average_NS <- length_delta %>%
    mutate(concentration_um = as.numeric(concentration_um)) %>%
    group_by(strain, concentration_um) %>%
    summarise(mean_median_wormlength_um_delta = mean(median_wormlength_um_delta),
              se_median_wormlength_um_delta = sd(median_wormlength_um_delta)) %>%
    dplyr::mutate(log.concentration = round(log(concentration_um), 2))
  
  strain_colors <- c(
    "PB420" = "#53886C",
    "ECA2666" = "#466EB4",
    "ED3102" = "#00A0E1",
    "JU3237" = "#F7D42A",
    "NIC1667" = "#7C1D6F",
    "VX34" = "#CC799D"
  )
  
  dose_average_NS$strain <- factor(dose_average_NS$strain, levels = c("PB420", "ECA2666", "ED3102", "JU3237", "NIC1667", "VX34"))
  if(drug_name == "PYR"){
    SE.DRC_NS <- dose_average_NS[which(dose_average_NS$log.concentration != -Inf),] %>% ggplot(., aes(x = concentration_um, y = mean_median_wormlength_um_delta, color = strain)) +
      theme_bw(base_size = 12) +
      geom_line(aes(group = strain), position = position_dodge2(width = 0.03), size = 0.75) +
      geom_pointrange(aes(ymin = mean_median_wormlength_um_delta - se_median_wormlength_um_delta,
                          ymax = mean_median_wormlength_um_delta + se_median_wormlength_um_delta),
                      position = position_dodge2(width = 0.03), size = 0.4) +
      scale_color_manual(values = strain_colors, labels = c("PB420" = "CGC2"), name = "Strain") +
      scale_y_continuous("Normalized Animal Length (µm)", limits = lower_limit(drug_name), breaks = seq(-800, 0, 100)) +
      ggtitle(drug_map(drug_name)) + 
      scale_x_log10(breaks = axis_breaks(drug_name),
                    labels = axis_labels(drug_name)) +
      theme(panel.grid = element_blank(),
            plot.title = element_text(face = "bold", family = "Arial", size = 10, hjust = 0),  # hjust = 0 for left align
            axis.title.x = element_blank(),
            axis.title.y = element_blank())
  }else{
    SE.DRC_NS <- dose_average_NS[which(dose_average_NS$log.concentration != -Inf),] %>% ggplot(., aes(x = concentration_um, y = mean_median_wormlength_um_delta, color = strain)) +
    theme_bw(base_size = 12) +
    geom_line(aes(group = strain), position = position_dodge2(width = 0.03), size = 0.75) +
    geom_pointrange(aes(ymin = mean_median_wormlength_um_delta - se_median_wormlength_um_delta,
                        ymax = mean_median_wormlength_um_delta + se_median_wormlength_um_delta),
                    position = position_dodge2(width = 0.03), size = 0.4) +
    scale_color_manual(values = strain_colors, labels = c("PB420" = "CGC2", "ECA2666, ED3102, JU3237, NIC1667, VX34"), name = "Strain") +
    scale_y_continuous("Normalized Animal Length (µm)", limits = lower_limit(drug_name), breaks = seq(-800, 0, 100)) +
    ggtitle(drug_map(drug_name))+ 
    scale_x_log10(breaks = axis_breaks(drug_name),
                  labels = axis_labels(drug_name)) +
    theme(panel.grid = element_blank(),
          legend.position = "none",
          plot.title = element_text(face = "bold", family = "Arial", size = 10, hjust = 0),  # hjust = 0 for left align
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
  }
  
  
  
  
  print(SE.DRC_NS)
  
  return(SE.DRC_NS)
}

ABZ_plot <- make_plot("ABZ")  
CLO_plot <- make_plot("CLO")
EMO_plot <- make_plot("EMO")
IVM_plot <- make_plot("IVM")
FBZ_plot <- make_plot("FBZ")
PYR_plot <- make_plot("PYR")  + theme(legend.position = "none")
LEV_plot <- make_plot("LEV")
MILB_plot <- make_plot("MILB")
MOX_plot <- make_plot("MOX")
legend <- get_legend(PYR_plot + theme(legend.position = "right"))

combine_plot <- plot_grid(ABZ_plot,FBZ_plot, NULL, IVM_plot, MILB_plot,
          MOX_plot, LEV_plot, PYR_plot, NULL, EMO_plot, CLO_plot, legend, nrow = 4, ncol = 3) +
  theme(plot.background = element_rect(fill = "white", colour = "white"))
final_plot <- ggdraw(combine_plot)+  draw_label("Concentration (μM)", x = 0.5, y = -0.045, vjust = -1, angle = 0, size = 30) +
  draw_label("Normalized Worm Length (µm)", x = -0.03, y = 0.5, vjust = 1.5, angle = 90, size = 30) +
  theme(plot.margin = margin(50, 50, 50, 50), plot.background = element_rect(fill = "white", color = "white"))
  
ggsave(filename= "final_plot.png", device = png, plot = final_plot, 
         path = "projects/Olivia/cbriggsae_drc/figures", dpi = 300, width = 20, height = 15)
  

rect <- rectGrob(
  gp = gpar(fill = "grey", col = "black", alpha = 0.75)
)


ggdraw(xlim = c(0,15), ylim = c(0,25)) +
  draw_grob(rect, x = 13, y = 20.2, width = 2, height = 4.8) + 
  draw_plot(ABZ_plot, x = 1, y = 20.2, width = 5, height = 4.8) + 
  draw_plot(FBZ_plot, x = 7, y = 20.2, width = 5, height = 4.8) +
  draw_grob(rect, x = 13, y = 10.6, width = 2, height = 9.6) +
  draw_plot(IVM_plot, x = 1, y = 15.4, width = 5, height = 4.8) +
  draw_plot (MILB_plot, x = 7, y = 15.4, width = 5, height = 4.8) +
  draw_plot(MOX_plot, x = 1, y = 10.6, width = 5, height = 4.8) +
  draw_plot(legend, x = 7, y = 11.6, width = 2.9, height = 2.9) +
  draw_grob(rect, x = 13, y = 5.8, width = 2, height = 4.8) + #
  draw_plot(LEV_plot, x = 1, y = 5.8, width = 5, height = 4.8) +
  draw_plot(PYR_plot, x = 7,  y = 5.8, width = 5, height = 4.8) +
  draw_grob(rect, x = 13, y = 0, width = 2, height = 5.8) +
  draw_plot(CLO_plot, x = 1, y = 1, width = 5, height = 4.8) +
  draw_plot(EMO_plot, x = 7, y = 1, width = 5, height = 4.8) +
  draw_label("Benzimidazoles", x = 14, y = 22.5, fontface = "bold", size = 10, angle = -90, colour = "black") +
  draw_label("Macrocyclic Lactones", x = 14, y = 15.5, fontface = "bold", size = 10, angle = -90, colour = "black") +
  draw_label("Nicotinic", x = 14.5, y = 8, fontface = "bold", size = 10, angle = -90, colour = "black") +
  draw_label("Acetylcholine", x = 14, y = 8, fontface = "bold", size = 10, angle = -90, colour = "black") +
  draw_label("Receptor Agonists", x = 13.5, y = 8.1, fontface = "bold", size = 10, angle = -90, colour = "black") +
  draw_label("Other", x = 14, y = 2.9, fontface = "bold", size = 10, angle = -90, colour = "black") +
  draw_label("Normalized Animal Length (µm)", x = 0.4, y = 13.5, fontface = "bold", size = 15, angle = 90, colour = "black") +
  draw_label("Concentration (μM)", x = 6, y = 0.4, fontface = "bold", size = 15, angle = 0, colour = "black") +
  
  draw_label(expression(paste(bold("EC"["50"]), bold(": 4.74"))), x = 2.8, y = 21.3, size = 10, fontface = "bold", colour = "black") + #ABZ
  draw_label(expression(paste(bolditalic("h2"), bold(": 0.73"))), x = 2.7, y = 22.1, size = 10, fontface = "bold", colour = "black") + #ABZ
  draw_label(expression(paste(bolditalic("H2"),bold(": 0.83"))), x = 2.7, y = 21.7, size = 10, fontface = "bold", colour = "black") + #ABZ
  
  draw_label(expression(paste(bold("EC"["50"]), bold(": 6.09"))), x = 8.8, y = 21.3, size = 10, fontface = "bold", colour = "black") + #FBZ
  draw_label(expression(paste(bolditalic("h2"), bold(": 0.39"))), x = 8.7, y = 22.1, size = 10, fontface = "bold", colour = "black") + #FBZ
  draw_label(expression(paste(bolditalic("H2"), bold(": 0.71"))), x = 8.7, y = 21.7, size = 10, fontface = "bold", colour = "black") + #FBZ
  
  draw_label(expression(paste(bold("EC"["50"]), bold(": 6.8e-04"))), x = 3.0, y = 16.5, size = 10, fontface = "bold", colour = "black") + #IVM
  draw_label(expression(paste(bolditalic("h2"), bold(": 0.81"))), x = 2.7, y = 17.3, size = 10, fontface = "bold", colour = "black") + #IVM
  draw_label(expression(paste(bolditalic("H2"), bold(": 0.87"))), x = 2.7, y = 16.9, size = 10, fontface = "bold", colour = "black") + #IVM
  
  draw_label(expression(paste(bold("EC"["50"]), bold(": 0.51"))), x = 8.8, y = 16.5, size = 10, fontface = "bold", colour = "black") + #MILB
  draw_label(expression(paste(bolditalic("h2"), bold(": 0.69"))), x = 8.7, y = 17.3, size = 10, fontface = "bold", colour = "black") + #MILB
  draw_label(expression(paste(bolditalic("H2"), bold(": 0.82"))), x = 8.7, y = 16.9, size = 10, fontface = "bold", colour = "black") + #MILB
  
  draw_label(expression(paste(bold("EC"["50"]), bold(": 1.5e-03"))), x = 3.0, y = 11.7, size = 10, fontface = "bold", colour = "black") + #MOX
  draw_label(expression(paste(bolditalic("h2"), bold(": 0.26"))), x = 2.7, y = 12.5, size = 10, fontface = "bold", colour = "black") + #MOX
  draw_label(expression(paste(bolditalic("H2"), bold(": 0.37"))), x = 2.7, y = 12.1, size = 10, fontface = "bold", colour = "black") + #MOX
  
  draw_label(expression(paste(bold("EC"["50"]), bold(": 9.42"))), x = 2.8, y = 6.9, size = 10, fontface = "bold", colour = "black") + #LEV
  draw_label(expression(paste(bolditalic("h2"), bold(": 0.94"))), x = 2.7, y = 7.7, size = 10, fontface = "bold", colour = "black") + #LEV
  draw_label(expression(paste(bolditalic("H2"), bold(": 0.97"))), x = 2.7, y = 7.3, size = 10, fontface = "bold", colour = "black") + #LEV
  
  draw_label(expression(paste(bold("EC"["50"]), bold(": 201.3"))), x = 8.9, y = 6.9, size = 10, fontface = "bold", colour = "black") + #PYR
  draw_label(expression(paste(bolditalic("h2"), bold(": 0.59"))), x = 8.7, y = 7.7, size = 10, fontface = "bold", colour = "black") + #PYR
  draw_label(expression(paste(bolditalic("H2"), bold(": 0.75"))), x = 8.7, y = 7.3, size = 10, fontface = "bold", colour = "black") + #PYR
  
  draw_label(expression(paste(bold("EC"["50"]), bold(": 22.1"))), x = 2.8, y = 2.1, size = 10, fontface = "bold", colour = "black") + #CLO
  draw_label(expression(paste(bolditalic("h2"), bold(": 0.95"))), x = 2.7, y = 2.9, size = 10, fontface = "bold", colour = "black") + #CLO
  draw_label(expression(paste(bolditalic("H2"), bold(": 0.98"))), x = 2.7, y = 2.5, size = 10, fontface = "bold", colour = "black") + #CLO
  
  draw_label(expression(paste(bold("EC"["50"]), bold(": 0.48"))), x = 8.8, y = 2.1, size = 10, fontface = "bold", colour = "black") + #EMO
  draw_label(expression(paste(bolditalic("h2"), bold(": 0.99"))), x = 8.7, y = 2.9, size = 10, fontface = "bold", colour = "black") + #EMO
  draw_label(expression(paste(bolditalic("H2"), bold(": 0.99"))), x = 8.7, y = 2.5, size = 10, fontface = "bold", colour = "black")   #EMO
  
ggsave(filename= "paperV6", device = png, plot = final_plot, 
       path = "projects/Olivia/cbriggsae_drc/figures/final_plots", dpi = 300, width = 750, height = 1038, units = "px")

  
  
  
  
  
  
  


