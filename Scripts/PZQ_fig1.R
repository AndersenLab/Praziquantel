library(tidyverse)
library(cowplot)

###############################
# Download Praziquantel folder and as set working directory
setwd("~/Desktop/Praziquantel/")
source("Scripts/theme_PZQ.R")
###############################
assign("Dose_enant",get(load("Data/PZQ_Dose_enentiomers.Rda")))

# get average control phenotype for each strain
dosepruned <- Dose_enant %>%
  dplyr::group_by(strain, condition, dose, trait) %>%
  dplyr::filter(., dose == 0) %>%
  dplyr::mutate(mean_ctrl = mean(phenotype)) %>%
  dplyr::ungroup() %>%
  dplyr::select(strain, mean_ctrl,trait) %>% ##added trait
  dplyr::full_join(out_regressed, by = c("trait","strain")) %>% ## "trait" instead of "condition"
  dplyr::distinct() %>%
  dplyr::mutate(newpheno = phenotype - mean_ctrl) %>%
  dplyr::mutate(colours = paste(drug,strain, sep=""))

straincol <- c("N2" = "orange",
               "JU775" = "darkorchid4",
               "ECA601" = "grey",
               "ECA485" = "grey")

drugline <- c("PZQ" = "solid",
              "PZQ_p1" = "longdash",
              "PZQ_p2" = "dashed")  

tweaking_lab_curve =   theme(strip.background = element_blank(),
                             axis.text.x = ggplot2::element_text(size = 8, angle =0), #element_blank(), , hjust = -0.2
                             axis.text.y = ggplot2::element_text(size = 8, angle = 0), #, hjust = 0.5
                             axis.title.x.bottom = ggplot2::element_text(size = 12, face = "bold", color = "black"), #element_blank(), #, vjust = 0.4),
                             axis.title.y = ggplot2::element_text(size = 12, face = "bold", color = "black"), #ggplot2::element_blank(), #element_text(size = 12, face = "bold", color = "black"), 
                             strip.text.x = ggplot2::element_text(size = 8, face = "bold", color = "black"), 
                             strip.text.y = ggplot2::element_text(size = 8, face = "bold", color = "black"), 
                             panel.border = ggplot2::element_rect(fill = NA, color = "black", linewidth = 0.5),
                             legend.position = "none")  

averagePZQ <- dosepruned %>%
  dplyr::filter(strain %in% c("N2","JU775")) %>%
  dplyr::filter(trait == "median.norm.EXT") %>%
  dplyr::group_by(strain,drug,dose) %>%
  dplyr::summarise(avg = mean(newpheno), stdev = sd(newpheno))

avg_curves_dotted <- averagePZQ %>% 
  dplyr::mutate(dose = as.numeric(dose)) %>%
  dplyr::mutate(strain = factor(strain, levels = c("N2","JU775"))) %>%
  ggplot(aes(x=dose, y=avg, group = interaction(strain,drug), color = strain, fill = strain))+ #, fill = strain
  geom_errorbar(linewidth = 0.25,width=150, aes(ymin=avg-stdev, ymax=avg+stdev),position = position_dodge(width = 150)) +
  geom_line(aes(linetype = drug),linewidth = 0.5,position = position_dodge(width = 150))+
  geom_point(data = averagePZQ,fill = "white" ,alpha = 1, pch = 21, size  = 1.5, colour = "white",position = position_dodge(width = 150)) +
  geom_point(pch = 21, size  = 1.5, colour = "darkgrey",position = position_dodge(width = 150)) +
  scale_colour_manual(values = straincol)+
  scale_fill_manual(values = straincol)+ 
  scale_linetype_manual(values = drugline) +
  theme_PZQ +
  tweaking_lab_curve + 
  scale_x_continuous(breaks = c(0,250,500,1000,2000,3000), 
                     #labels = xlabel, 
                     labels = c("0","0.25","0.5","1","2","3"), 
                     limits = c(-75,3075),
                     NULL )+ #"Emodepside concentration (mM)") 
  scale_y_continuous("Phenotypic response", breaks = c(-0.4,0,0.4,0.8),  
                     labels = c("-0.4", "0","0.4","0.8"),expand = expansion(mult = c(0.05, 0.05))) 

final <- add_sub(avg_curves_dotted, "Drug concentration (mM)", y  = 0, vjust = 0, fontface = "bold", size =12)

ggsave2("OUT/PZQ_Figure1.png",final, height = 3.5, width = 5)



