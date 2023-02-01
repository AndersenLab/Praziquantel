library("tidyverse")
library("ggbeeswarm")
library("cowplot")

###############################
# Download Praziquantel folder and as set working directory
setwd("~/Desktop/Praziquantel/")
###############################


### FIGURE A
## QTL data
QTL_Chrom <- "IV"
QTL_Region_start <- 939925
QTL_Peak <- 1169239
QTL_Region_end <- 1334212

## Trait
trait_name <- "praziquantel.median.norm.EXT"   

# load independent tests result from NemaScan
total_independent_tests <- read.table("Data/total_independent_tests.txt", quote="\"", comment.char="", stringsAsFactors=FALSE)
independent_tests <- total_independent_tests[[1]]
independent_test_cutoff <- -log10(0.05/independent_tests)

###Manhattan Plot
processed_mapping <- read.delim("Data/processed_praziquantel_median_norm_EXT_AGGREGATE_mapping_inbred.tsv")

manhattan_plot <- processed_mapping %>%
  dplyr::filter(trait == unique(processed_mapping$trait)[1] & CHROM != "MtDNA") %>%
  dplyr::distinct(marker, .keep_all = T) %>%
  dplyr::mutate(EIGEN_CUTOFF = independent_test_cutoff) %>%
  dplyr::mutate(EIGEN_SIG = ifelse(log10p > BF, "1", 
                                   ifelse(log10p > EIGEN_CUTOFF, "2", "0")) )%>%
  ggplot2::ggplot(.) +
  ggplot2::aes(x = POS/1e6, y = log10p) +
  ggplot2::scale_color_manual(values = c("0" = "black", 
                                         "1" = "red",
                                         "2" = "black")) +
  ggplot2::geom_hline(ggplot2::aes(yintercept = BF),
                      color = "gray", 
                      alpha = .75,  
                      linewidth = .751) +
  ggplot2::geom_point( ggplot2::aes(color= factor(EIGEN_SIG)), size =1 ) +
  ggplot2::facet_grid( . ~ CHROM, scales = "free_x" , space = "free_x") +
  ggplot2::theme_bw(12) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(size = 10),
                 axis.text.y = ggplot2::element_text(size = 10),
                 axis.title.x = ggplot2::element_text(size = 12,color = "black", face = "bold", vjust = -0.3), 
                 axis.title.y = ggplot2::element_text(size = 12, face = "bold", color = "black", vjust = 1), 
                 strip.text.x = ggplot2::element_text(size = 12, face = "bold", color = "black"), 
                 strip.text.y = ggplot2::element_text(size = 12, face = "bold", color = "black"), 
                 panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(),
                 legend.position = "none",
                 strip.background = element_blank()) +
  scale_y_continuous(name = expression(bold(-log[10](italic(p)))), breaks=c(0,2,4,6), labels = c("0.0","0.2","0.4","0.6"))+ 
  ggplot2::labs(x = "Genomic position (Mb)",
                y = expression(bold(-log[10](italic(p)))))



### FIGURE B -- PxG Plot
pxg_df <- na.omit(processed_mapping) %>%
  dplyr::filter(CHROM==QTL_Chrom, peakPOS==QTL_Peak) %>%  
  dplyr::mutate(value=as.numeric(value)) %>% 
  group_by(allele) %>%
  dplyr::mutate(mean_pheno = mean(value, na.rm = T)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(straintype = case_when(strain == "N2" ~ "N2",
                                       strain == "JU775" ~ "JU775",
                                       TRUE  ~ "Other"))

pxg_plot <- pxg_df %>% 
  ggplot()+
  aes(x = factor(allele, levels = c(-1,1), labels = c("REF","ALT")))+
  geom_beeswarm(cex=3, priority='density',
                aes(y = value,fill = straintype),
                shape = 21, 
                size = 1.5)+
  theme_bw(12)+
  labs(y = "Response",
       x = "Genotype") +
  scale_fill_manual(values = c("N2"="orange","JU775"="purple","Other"="black"))+
  theme(axis.text.x = ggplot2::element_text(size = 10),
        axis.text.y = ggplot2::element_text(size = 10),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = ggplot2::element_text(size = 12, face ="bold", color = "black"), 
        axis.title.y = ggplot2::element_text(size = 12, face ="bold", color = "black", vjust = 1))+
  theme(legend.position = "none") 


### FIGURE C -- Fine Mapping Plot

genes_in_region <- data.table::fread("Data/praziquantel_median_norm_EXT_IV_939925-1334212_bcsq_genes_inbred.tsv")

gene_df <- genes_in_region %>%
  dplyr::filter(START_POS == unique(genes_in_region$START_POS)) %>%
  dplyr::distinct(WBGeneID, GENE_NAME, PEAK_MARKER, CHROM, STRAND, TRANSCRIPTION_START_POS, 
                  TRANSCRIPTION_END_POS, START_POS, END_POS, VARIANT_LOG10p) %>%
  dplyr::group_by(WBGeneID) %>%
  dplyr::mutate(VARIANT_LOG10p = max(VARIANT_LOG10p, na.rm = T)) %>%
  dplyr::distinct() %>%
  drop_na(STRAND)

peak_variant <- QTL_Peak

variant_df <- genes_in_region %>%
  dplyr::filter(START_POS == unique(genes_in_region$START_POS)) %>%
  dplyr::distinct(CHROM, POS, GENE_NAME, WBGeneID, VARIANT_LOG10p, VARIANT_IMPACT)

variant_df$VARIANT_IMPACT[is.na(variant_df$VARIANT_IMPACT)] <- "Intergenic"

xs <- unique(gene_df$START_POS)
xe <- unique(gene_df$END_POS)
cq <- unique(gene_df$CHROM)

max_logp <- unique(max(variant_df$VARIANT_LOG10p))/150

  ####Select data for Table S3
  table_s3 <- gene_df %>%
    dplyr::select(GENE_NAME, WBGeneID, VARIANT_LOG10p, TRANSCRIPTION_START_POS,TRANSCRIPTION_END_POS,STRAND) %>%
    dplyr::arrange(desc(VARIANT_LOG10p))%>%
    dplyr::mutate(Position = paste(TRANSCRIPTION_START_POS,"-",TRANSCRIPTION_END_POS)) %>%
    dplyr::rename(Gene = GENE_NAME,
                  "WormBase Gene ID" =  WBGeneID,
                  "-log10(p)" = VARIANT_LOG10p,
                  Strand = STRAND) %>%
    dplyr::select(Gene,"WormBase Gene ID","-log10(p)",Position,Strand)

  write_tsv(table_s3 , "~/Desktop/20230201_Gene_in_Region_PZQ.tsv")
  ####

gene_plot <- ggplot(gene_df) +
  geom_vline(aes(xintercept = peak_variant/1e6),
             linetype=3, color = "cyan")+
  geom_segment(aes(x = ifelse(STRAND == "+", TRANSCRIPTION_START_POS/1e6, TRANSCRIPTION_END_POS/1e6),
                   xend = ifelse(STRAND == "+", TRANSCRIPTION_END_POS/1e6, TRANSCRIPTION_START_POS/1e6),
                   y = VARIANT_LOG10p,
                   yend = VARIANT_LOG10p,
                   color = STRAND),
               size = 1.5) +
  geom_label(data=gene_df %>%filter(WBGeneID == "WBGene00021934"),
             aes(x = (TRANSCRIPTION_START_POS/1e6+TRANSCRIPTION_END_POS/1e6)/2,
                 y=VARIANT_LOG10p,label="cct-8"), nudge_x = -0.025,nudge_y =0.2,size=3.5,color="black",label.padding = unit(0.1,"lines"))+
  geom_segment(data = variant_df,
               aes(x = POS/1e6,
                   xend = POS/1e6,
                   y = VARIANT_LOG10p+max_logp,
                   yend = VARIANT_LOG10p-max_logp,
                   color = VARIANT_IMPACT)) +
  scale_color_manual(values = c("+" = "#3bb5c4",
                                "-" = "#341d8d",
                                "MODIFIER" = "gray50",
                                "LOW" = "gray30",
                                "MODERATE" = "orange",
                                "HIGH" = "red",
                                "Intergenic" = "gray80"),
                     name = "Variant Impact")+
  labs(x = "Genomic position (Mb)",
       y = expression(bold(-log[10](bold(italic(p)))))) +
  theme_bw(12)+
  ylim(2,NA)+
  theme(legend.position = "None",
        panel.grid = element_blank(),
        axis.text.x = ggplot2::element_text(size = 10),
        axis.text.y = ggplot2::element_text(size = 10),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = ggplot2::element_text(size = 12, face ="bold", color = "black", vjust = 1))


## Making figure

plots <- cowplot::align_plots(manhattan_plot, pxg_plot)#, align = 'v', axis = 'l')

# then build the bottom row
bottom_row <- plot_grid(plots[[2]], gene_plot, labels = c('B', 'C'), label_size = 12, rel_widths = c(.5, 1))

# then combine with the top row for final plot
full <- plot_grid(plots[[1]], bottom_row, labels = c('A', ''), label_size = 12, ncol = 1, rel_heights = c(3,2.5))

ggsave(full,filename = "OUT/PZQ_figure2.png", device = "png", width = 7.5, height = 5, units = "in")

