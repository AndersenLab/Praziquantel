library(tidyverse)
library(data.table)
library(cowplot)

###############################
# Download Praziquantel folder and as set working directory
setwd("~/Desktop/Praziquantel/")
###############################

#theme
PZQ_fig_expr <- ggplot2::theme(panel.background = element_rect(fill = NA, colour = "black", linewidth = 0.5),
                               panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               legend.position="none", 
                               axis.text.y = element_text(size = 8, face = "bold"),
                               axis.text.x = element_text(size = 8, angle = 45, vjust = .5, face = "bold.italic"),
                               axis.title.x = element_text( size=12, face = "bold"),
                               axis.title.y = element_blank(),
                               title = element_text(size=12, face ="bold"),
                               strip.text = element_text(colour = 'black'),
                               axis.ticks = element_line(colour = "black", linewidth = 0.2), 
                               axis.ticks.length = unit(0.05, "cm"),
                               strip.background = element_rect(colour = "white", fill = "white"), 
                               strip.text.x = element_text(size=12, face = "bold.italic"))

### Figure A
PZQnew <- data.table::fread(file= "Data/RNAi.tsv", header = TRUE)

tukey_hsp16.2 <- as_tibble(PZQnew) %>%
  dplyr::filter(gene %in% c("hsp16.2")) %>%
  dplyr::mutate(condition = case_when(condition == "cct-8" ~ "cct8",
                                      condition == "act-1" ~ "act1",
                                      condition == "cct-8_act-1" ~ "cct8.act1",
                                      condition == "PZQ_act-1" ~ "PZQ.act1",
                                      condition != "cct-8" ~ condition ))%>%
  dplyr::ungroup() %>%
  dplyr::group_by(gene)  %>%
  rstatix::tukey_hsd(value ~ condition, ref.group = "control")

tukey_hsp70 <- as_tibble(PZQnew) %>%
  dplyr::filter(gene %in% c("hsp70")) %>%
  dplyr::mutate(condition = case_when(condition == "cct-8" ~ "cct8",
                                      condition == "act-1" ~ "act1",
                                      condition == "cct-8_act-1" ~ "cct8.act1",
                                      condition == "PZQ_act-1" ~ "PZQ.act1",
                                      condition != "cct-8" ~ condition ))%>%
  dplyr::ungroup() %>%
  dplyr::group_by(gene)  %>%
  rstatix::tukey_hsd(value ~ condition, ref.group = "control")

cond_col <- c("control" = "darkgrey",
              "act-1" = "darkgrey",
              "cct-8" = "darkgrey",
              "PZQ" = "darkgrey",
              "cct-8_act-1" = "darkgrey",
              "PZQ_act-1" = "darkgrey")  


RNAi_hsp70 <- as_tibble(PZQnew) %>%
  dplyr::mutate(condition = factor(condition, levels = c("control","act-1","cct-8","PZQ","cct-8_act-1","PZQ_act-1"))) %>%
  dplyr::mutate(gene = factor(gene, levels = c("hsp16.2","hsp16.11","hsp70"))) %>%
  dplyr::mutate(gene = case_when(gene == "hsp16.2" ~ "hsp-16.2",
                                 gene == "hsp70" ~ "hsp-70")) %>%
  dplyr::filter(gene %in% c("hsp-70")) %>%
  ggplot2::ggplot(aes(x = condition, y = value)) +
  geom_boxplot(aes(x = condition)) +
  facet_grid(~factor(gene, levels=c("hsp-16.2","hsp-70"))) +
  PZQ_fig_expr + 
  scale_colour_manual(values = cond_col) +
  theme(
    axis.text.x = element_text(size=10, face = "bold"),
    strip.background = element_rect(colour = "white", fill = "white")) +
  scale_y_continuous(name = NULL, limits = c(-2,110) ) +
  scale_x_discrete(name = NULL, labels = c("control","act-1","cct-8","PZQ","act-1 + cct-8","act-1 + PZQ")) 


RNAi_hsp16 <- as_tibble(PZQnew) %>%
  dplyr::mutate(condition = factor(condition, levels = c("control","act-1","cct-8","PZQ","cct-8_act-1","PZQ_act-1"))) %>%
  dplyr::mutate(gene = factor(gene, levels = c("hsp16.2","hsp16.11","hsp70"))) %>%
  dplyr::mutate(gene = case_when(gene == "hsp16.2" ~ "hsp-16.2",
                                 gene == "hsp70" ~ "hsp-70")) %>%
  dplyr::filter(gene %in% c("hsp-16.2")) %>%
  ggplot2::ggplot(aes(x = condition, y = value)) +
  geom_boxplot(aes(x = condition)) +
  facet_grid(~factor(gene, levels=c("hsp-16.2","hsp-70"))) +
  PZQ_fig_expr + 
  scale_colour_manual(values = cond_col) +
  theme(#strip.text.x = element_blank(),
    axis.text.x = element_text(size=10, face = "bold"),
    strip.background = element_rect(colour = "white", fill = "white")) +
  scale_y_continuous(name = NULL, limits = c(-2,220) )+#expand =c(0,50)) +
  scale_x_discrete(name = NULL, labels = c("control","act-1","cct-8","PZQ","act-1 + cct-8","act-1 + PZQ")) 


###FIGURE B
 ## MERGING ALL GENES (combining data from cct8, hsp70, hsp16 and rpl26_expression files)
 ## ALREADY MERGED AND NORMALIZED DATA BELOW #### (search DATAREADY) ####

assign("hsp16_full",get(load("Data/hsp16_expression.rda")))
assign("hsp70_full",get(load("Data/hsp70_expression.rda")))
assign("cct8_full",get(load("Data/cct8_expression.rda")))
assign("rpl26_full",get(load("Data/rpl26_expression.rda")))


ALL <- rbind(cct8_full,rpl26_full,hsp70_full,hsp16_full) %>%
  dplyr::filter(!sample %in% c("St1","St2","St3","St4","St5")) %>%
  dplyr::mutate(gene = case_when(gene == "hsp-70" ~ "hsp70",
                                 gene == "cct-8" ~ "cct8",
                                 !gene %in% c("hsp-70","cct-8") ~ gene )) %>%
  dplyr::mutate(sample = case_when(sample == "cal" & plate == "plate1" ~ "cal1",
                                   sample == "cal" & plate == "plate2" ~ "cal2",
                                   sample != "cal" ~ sample)) %>%
  dplyr::select(-new_mean,-log,-well)%>%
  tidyr::spread(key = gene, value = concentration)

cals <- ALL[ALL$sample %in% c("cal1","cal2"),] %>%
  dplyr::mutate(norm_cct8 = cct8/rpl26,
                norm_hsp70 = hsp70/rpl26,
                norm_hsp16 = hsp16_2/rpl26)

all_data <- ALL %>%
  dplyr::filter(! sample %in% c("cal1","cal2")) %>%
  dplyr::mutate(cal_cct8 = case_when(plate == "plate1" ~ as.character(cals[1,10]),
                                     plate == "plate2" ~ as.character(cals[2,10]))) %>%
  dplyr::mutate(cal_hsp70 = case_when(plate == "plate1" ~ as.character(cals[1,11]),
                                     plate == "plate2" ~ as.character(cals[2,11]))) %>%
  dplyr::mutate(cal_hsp16 = case_when(plate == "plate1" ~ as.character(cals[1,12]),
                                     plate == "plate2" ~ as.character(cals[2,12]))) %>%
  dplyr::mutate(cal_cct8 = as.numeric(cal_cct8),
                cal_hsp70 = as.numeric(cal_hsp70),
                cal_hsp16 = as.numeric(cal_hsp16)) %>%
  dplyr::mutate(norm_cct8 = cct8/rpl26/cal_cct8,
                norm_hsp70 = hsp70/rpl26/cal_hsp70,
                norm_hsp16 = hsp16_2/rpl26/cal_hsp16)


#### LOAD DATA ### DATAREADY
assign("all_data",get(load("Data/analysed_normalized_pzq.rda")))
    
  ######### STATS #############

### remove outliers
summary(all_data)
## cctt-8
cct8_no_outliers <- all_data[,-c(8,9,10,12,13,15,16)] %>%
  dplyr::group_by(strain,condition) %>%
  dplyr::summarize(mean = mean(norm_cct8), sd = sd(norm_cct8)) %>%
  dplyr::left_join(all_data[,-c(8,9,10,12,13,15,16)]) %>%
  dplyr::mutate(min = (mean - 2*sd),
                max = (mean + 2*sd)) %>%
  dplyr::filter(norm_cct8 > min & norm_cct8 < max) %>%
  dplyr::ungroup()

two.way.cct8_no_outliers <- aov(norm_cct8 ~ strain*condition + plate, data = cct8_no_outliers)
summary(two.way.cct8_no_outliers)

# analyse cct-8 PCR plates seperately because of plate effect
two.way.cct8.1.no <- aov(norm_cct8 ~ strain*condition, data = cct8_no_outliers[cct8_no_outliers$plate == "plate1",])
two.way.cct8.2.no <- aov(norm_cct8 ~ strain*condition, data = cct8_no_outliers[cct8_no_outliers$plate == "plate2",])

## hps-16
hsp16_no_outliers <- all_data[,-c(7,9,10,11,12,14,15)] %>%
  dplyr::group_by(strain,condition) %>%
  dplyr::summarize(mean = mean(norm_hsp16), sd = sd(norm_hsp16)) %>%
  dplyr::left_join(all_data[,-c(7,9,10,11,12,14,15)]) %>%
  dplyr::mutate(min = (mean - 2*sd),
                max = (mean + 2*sd)) %>%
  dplyr::filter(norm_hsp16 > min & norm_hsp16 < max) %>%
  dplyr::ungroup()

two.way.hsp16_no_outliers <- aov(norm_hsp16 ~ strain*condition + plate, data = hsp16_no_outliers)
summary(two.way.hsp16_no_outliers)

## hsp-70

hsp70_no_outliers <- all_data[,-c(7,8,10,11,13,14,16)] %>%
  dplyr::group_by(strain,condition) %>%
  dplyr::summarize(mean = mean(norm_hsp70), sd = sd(norm_hsp70)) %>%
  dplyr::left_join( all_data[,-c(7,8,10,11,13,14,16)]) %>%
  dplyr::mutate(min = (mean - 2*sd),
                max = (mean + 2*sd)) %>%
  dplyr::filter(norm_hsp70 > min & norm_hsp70 < max) %>%
  dplyr::ungroup()

two.way.hsp70_no_outliers <- aov(norm_hsp70 ~ strain*condition + plate, data = hsp70_no_outliers)
summary(two.way.hsp70_no_outliers)

## Plotting

plot_file <-  all_data %>%
  mutate(norm_cct8 = replace(norm_cct8, sample %in% c("9","10","13","30"), NA),
         norm_hsp16 = replace(norm_hsp16, sample %in% c("37"), NA),
         norm_hsp70 = replace(norm_hsp70, sample %in% c("37","43","10"), NA)) %>%
  dplyr::select(-c(7,8,9,10,11,12,13)) %>%
  tidyr::gather(key="gene", value = "norm_expr", norm_cct8, norm_hsp70 ,norm_hsp16)

#names for plot
gene_name <- c("norm_cct8" = "cct-8",
               "norm_hsp16" = "hsp-16.2",
               "norm_hsp70" = "hsp-70")  

expres_plot_cct <- plot_file %>%
  dplyr::mutate(gene = case_when(gene == "norm_cct8" ~ "cct-8",
                                 gene == "norm_hsp16" ~ "hsp-16.2",
                                 gene == "norm_hsp70" ~ "hsp-70")) %>%
  dplyr::mutate(group = paste(strain, condition, sep = "_")) %>%
  dplyr::mutate(strain =  factor(strain, levels = c("N2", "JU775"))) %>%
  dplyr::filter(gene ==  "cct-8") %>%
  ggplot(aes(x=strain,y=norm_expr, group = interaction(strain,condition), colour = condition)) +
  geom_boxplot(outlier.colour = NA) +
  facet_grid(~factor(gene, levels=c("hsp-16.2","hsp-70","cct-8"))) +
  geom_point(aes(x=strain,y=norm_expr, group = interaction(strain,condition), colour = condition, shape = plate), position =  position_jitterdodge(dodge.width=0.8))+
  PZQ_fig_expr + 
  scale_colour_manual(values = c("black","darkgrey")) + 
  theme(#strip.text.x = element_text(size=12, face = "bold.italic"),
    axis.text.x = element_text(size=10, face = "bold"),
    strip.background = element_rect(colour = "white", fill = "white")) +
  scale_x_discrete(name = NULL) + 
  scale_y_continuous(name = "Normalized expression", expand =c(0.02,1))

expres_plot_hsp16 <- plot_file %>%
  dplyr::mutate(gene = case_when(gene == "norm_cct8" ~ "cct-8",
                                 gene == "norm_hsp16" ~ "hsp-16.2",
                                 gene == "norm_hsp70" ~ "hsp-70")) %>%
  dplyr::mutate(group = paste(strain, condition, sep = "_")) %>%
  dplyr::mutate(strain =  factor(strain, levels = c("N2", "JU775"))) %>%
  dplyr::filter(gene ==  "hsp-16.2") %>%
  ggplot(aes(x=strain,y=norm_expr, group = interaction(strain,condition), colour = condition)) +
  geom_boxplot(outlier.colour = NA) +
  facet_grid(~factor(gene, levels=c("hsp-16.2","hsp-70","cct-8"))) +
  geom_point(aes(x=strain,y=norm_expr, group = interaction(strain,condition), colour = condition, shape = plate), position =  position_jitterdodge(dodge.width=0.8))+
  PZQ_fig_expr + 
  scale_colour_manual(values = c("black","darkgrey")) + 
  theme(#strip.text.x = element_text(size=12, face = "bold.italic"),
    axis.text.x = element_text(size=10, face = "bold"),
    strip.background = element_rect(colour = "white", fill = "white")) +
  scale_x_discrete(name = NULL) + 
  scale_y_continuous(name = "Normalized expression", expand =c(0.02,0.25))

expres_plot_hsp70 <- plot_file %>%
  dplyr::mutate(gene = case_when(gene == "norm_cct8" ~ "cct-8",
                                 gene == "norm_hsp16" ~ "hsp-16.2",
                                 gene == "norm_hsp70" ~ "hsp-70")) %>%
  dplyr::mutate(group = paste(strain, condition, sep = "_")) %>%
  dplyr::mutate(strain =  factor(strain, levels = c("N2", "JU775"))) %>%
  dplyr::filter(gene ==  "hsp-70") %>%
  ggplot(aes(x=strain,y=norm_expr, group = interaction(strain,condition), colour = condition)) +
  geom_boxplot(outlier.colour = NA) +
  facet_grid(~factor(gene, levels=c("hsp-16.2","hsp-70","cct-8"))) +
  geom_point(aes(x=strain,y=norm_expr, group = interaction(strain,condition), colour = condition, shape = plate), position =  position_jitterdodge(dodge.width=0.8))+
  PZQ_fig_expr + 
  scale_colour_manual(values = c("black","darkgrey")) + 
  theme(#strip.text.x = element_text(size=12, face = "bold.italic"),
    axis.text.x = element_text(size=10, face = "bold"),
    strip.background = element_rect(colour = "white", fill = "white")) +
  scale_x_discrete(name = NULL) + 
  scale_y_continuous(name = "Normalized expression", expand =c(0.02,0.25))




######
######
#MERGING A and B
######
######

#plots needed
expres_plot_hsp70
expres_plot_cct
expres_plot_hsp16
RNAi_hsp70
RNAi_hsp16

#Merge plots

align_r <- cowplot::align_plots(RNAi_hsp16,expres_plot_cct ,align = 'v')#, align = 'v', axis = 'l')

top_row <- plot_grid(align_r [[1]], NULL, RNAi_hsp70, ncol =3, rel_widths = c(1,0.05, 1))
bottom_row <- plot_grid(align_r [[2]],NULL, expres_plot_hsp16,NULL, expres_plot_hsp70, ncol =5, rel_widths = c(1,0.05,1,0.05, 1))

#full <- plot_grid(NULL, top_row,NULL, bottom_row, labels = c( "","A","","B"), label_size = 12, ncol = 1, rel_heights = c(0.01,1,0.01,1))

full <- plot_grid(top_row,NULL, bottom_row, ncol = 1, rel_heights = c(1,0.01,1))

label_space <- plot_grid(NULL, full, ncol = 2, rel_widths = c(.05,1))

#label_space

add_labels <- ggdraw(label_space) + 
  draw_label("A", color = "black", size = 12, fontface = "bold", angle = 0, x = 0.02, y = 0.97) +
  draw_label("B", color = "black", size = 12, fontface = "bold", angle = 0, x = 0.02, y = 0.47) +
  draw_label("Normalized espression", color = "black", size = 12, fontface = "bold", angle = 90, x = 0.03, y = 0.27) +
  draw_label("Normalized espression", color = "black", size = 12, fontface = "bold", angle = 90, x = 0.03, y = 0.78)

add_sign <- ggdraw(add_labels )+
  theme(plot.background = element_rect(fill="white", color = NA)) +
  
  ##cct-8
  
  draw_line(x = c(0.135,0.185),
            y = c(0.405), #0.87
            color = "black", size = .5
  ) +
  draw_line(x = c(0.25,0.30),
            y = c(0.405), #0.87
            color = "black", size = .5
  ) +
  draw_line(x = c(0.16), #c(0.125),
            y = c(0.405,.42), #0.87
            color = "black", size = .5
  ) +
  draw_line(x = c(0.275), #c(0.255),
            y = c(0.405,.42), #0.87
            color = "black", size = .5
  ) +
  draw_line(x = c(0.16,0.275),
            y = c(0.42), #0.87
            color = "black", size = .5) +
  
  ### hsp-70  
  
  draw_line(x = c(0.785,0.835),
            y = c(0.405), #0.87
            color = "black", size = .5
  ) +
  draw_line(x = c(0.895,0.945),
            y = c(0.405), #0.87
            color = "black", size = .5
  ) +
  draw_line(x = c(0.81), #c(0.785),
            y = c(0.405,.42), #0.87
            color = "black", size = .5
  ) +
  draw_line(x = c(0.92),# c(0.915),
            y = c(0.405,.42), #0.87
            color = "black", size = .5
  ) +
  draw_line(x = c(0.81,0.92),
            y = c(0.42), #0.87
            color = "black", size = .5) +
  
  ## hsp16.2 RNAi -control cct8
  draw_line(x = c(0.135), #c(0.785),
            y = c(0.915,0.925), #0.87
            color = "black", size = .5
  ) +
  draw_line(x = c(0.26),# c(0.915),
            y = c(0.915,0.925), #0.87
            color = "black", size = .5
  ) +
  draw_line(x = c(0.185), #c(0.785),
            y = c(0.895,0.905), #0.87
            color = "black", size = .5
  ) +
  draw_line(x = c(0.26),# c(0.915),
            y = c(0.895,0.905), #0.87
            color = "black", size = .5
  ) +
  draw_line(x = c(0.135,0.26),
            y = c(0.925), #0.87
            color = "black", size = .5) +
  draw_line(x = c(0.185,0.26),
            y = c(0.905), #0.87
            color = "black", size = .5) +
  draw_label("**", color = "black", size = 18, fontface = "bold", angle = 0, x = 0.2075, y = 0.424) +
  draw_label("***", color = "black", size = 18, fontface = "bold", angle = 0, x = 0.865, y = 0.424) +
  draw_label("*", color = "black", size = 18, fontface = "bold", angle = 0, x = 0.1975, y = 0.929) +
  draw_label("*", color = "black", size = 18, fontface = "bold", angle = 0, x = 0.2225, y = 0.909) 


#add_sign 
ggsave2("OUT/PZQ_Figure3.png", add_sign , height = 7.5, width = 7.5)

