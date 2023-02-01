library("tidyverse")
library("cowplot")
library("ggmsa")

###############################
# Download Praziquantel folder and as set working directory
setwd("~/Desktop/Praziquantel/")
###############################

##A -- Gene model
cct8_gene_info <- data.frame(feature = c("5'UTR", "Exon", "Intron", "Exon", "Intron", "Exon", "Intron", "Exon", "Intron", "Exon","Intron","Exon","Intron","Exon", "3'UTR"),
                             start = c(1093920,1093915,1093765,1093703,1093566,1093514,1093439,1092174,1091667,1090115,1089880,1089287,1088933,1088003,1087811),
                             stop = c(1093914,1093764,1093704,1093567,1093515,1093440,1092175,1091668,1090116,1089881,1089288,1088934,1088004,1087812,1087095),
                             color = c("blue","orange","gray","orange","gray","orange","gray","orange","gray","orange","gray","orange","gray","orange","blue"),
                             gene = "cct-8",
                             strand = "Watson")%>%
  dplyr::mutate(size = start-stop)

gene_model <- function(df, 
                       gene_color = "darkgrey",
                       intron_color = "black",
                       utr3_color = "lightgrey",
                       utr5_color = "gray60",
                       gene_alpha = 0.5){
  
  strand <- unique(df$strand)
  
  if(strand == "Watson"){
    
    exons <- df%>%
      dplyr::filter(feature == "Exon")
    
    introns <- df%>%
      dplyr::filter(feature == "Intron")%>%
      dplyr::mutate(halfpt = (start+stop )/2)
    
    startUTR <- df %>%
      dplyr::filter(feature == "5'UTR")
    
    endUTR <- df %>%
      dplyr::filter(feature == "3'UTR")%>%
      dplyr::select(start, stop)%>%
      tidyr::gather(loc, x)%>%
      dplyr::mutate(y = ifelse(loc == "start", 1, 0))
    
    endUTR_a <- dplyr::filter(endUTR, loc == "start")
    endUTR_a$y <- -1
    
    endUTR_pt <- bind_rows(endUTR, endUTR_a)
    
    ggplot(exons)+
      geom_rect( aes(xmin =  stop, xmax = start, ymin = -1 , ymax = 1), fill = gene_color, color = "black",alpha = gene_alpha)+
      geom_segment(aes(x = stop, y = 1, xend = halfpt, yend = 2), data = introns, color = intron_color)+
      geom_segment(aes(x = start, y = 1, xend = halfpt, yend = 2), data = introns, color = intron_color)+
      geom_rect(aes(xmin =  stop, xmax = start, ymin = -1 , ymax = 1), data = startUTR, fill = utr5_color, color = "black",alpha = gene_alpha)+
      geom_polygon(aes(x = x, y = y), fill = utr3_color,color = "black", data = endUTR_pt,alpha = gene_alpha)+
      theme_void()
  }}

cct8_genemodel <- gene_model(cct8_gene_info, intron_color = "black", gene_color = "darkgrey", utr5_color = "lightgrey", utr3_color = "lightgrey")+
  geom_segment(aes(x=1091858,xend=1091858, y=1,yend=1.32))+ 
  geom_label(aes(x=1091858,y=1.5,label = "G226V"))

##B -- Multiple sequence alignment
N2_JU_fasta <- "Data/N2_JU775.fasta"

#select fragment
N2JU <- ggmsa(N2_JU_fasta, 218, 229, char_width = 0.5, seq_name = T, color = "Clustal")+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12, face = "bold", colour= "black"))

##C -- response CRISP edited strains
assign("PZQ1000",get(load("Data/CRISPR_response_PZQ1000.Rda")))

PZQ1dose <- PZQ1000 %>%
  dplyr::mutate(strain = factor(strain, levels = c("N2","ECA485(G226V)","JU775","ECA601(V226G)")))%>%
  ggplot()+
  aes(x=strain, y=phenotype)+
  geom_boxplot(aes(fill=strain),outlier.shape = NA)+
  geom_jitter(width = 0.1)+
  theme_bw(12)+
  scale_fill_manual(values = c("N2"="orange","JU775"="purple","ECA485(G226V)"="grey","ECA601(V226G)"="grey"))+
  ylab("Normalized optical density")+
  xlab("Strain")+
  theme(legend.position = "None",
        panel.grid = element_blank(),
        axis.text.x = ggplot2::element_text(size = 10),
        axis.text.y = ggplot2::element_text(size = 10),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = ggplot2::element_text(size = 12, face ="bold", color = "black", vjust = 1)) +
  scale_y_continuous(breaks=c(-.5,-.25,0,.25,.5), limits = c(-0.5,0.70))


## Figure assembly
left <- cowplot::plot_grid(NULL, cct8_genemodel,NULL, N2JU,ncol =1, labels = c("A","","","B"), rel_heights = c(0.05,0.25,0.05,0.2))
total <- cowplot::plot_grid(left, NULL, PZQ1dose, ncol =3, labels = c("","","C") ,rel_widths = c(1,0.1,1.5)) +
  draw_label("ns", color = "black", size = 12, fontface = "italic", angle = 0, x = 0.635, y = 0.925) +
  draw_label("ns", color = "black", size = 12, fontface = "italic", angle = 0, x = 0.87, y = 0.925) +
  theme(plot.background = element_rect(fill="white", color =NA)) +
  draw_line(x = c(0.92),
            y = c(0.87,0.885), #0.87
            color = "black", size = .5
  ) +
  draw_line(x = c(0.82),
            y = c(0.87,0.885), #0.87
            color = "black", size = .5
  ) +
  draw_line(x = c(0.82,0.92),
            y = c(0.885), #0.87
            color = "black", size = .5
  ) +
  draw_line(x = c(0.685),
            y = c(0.87,0.885), #0.87
            color = "black", size = .5
  ) +
  draw_line(x = c(0.585),
            y = c(0.87,0.885), #0.87
            color = "black", size = .5
  ) +
  draw_line(x = c(0.585,0.685),
            y = c(0.885), #0.87
            color = "black", size = .5) 

ggsave(total , filename = "OUT/PZQ_figure4.png", device = "png", width = 7.5, height = 3.5, units = "in")

