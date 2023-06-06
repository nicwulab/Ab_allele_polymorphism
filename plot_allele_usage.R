#R code
library(ggplot2)
library(scales)
library(RColorBrewer)
library(qualpalr)
library(readr)
library(tidyr)
library(reshape)
library(stringr)
library(dplyr)
library(data.table)
library(gridExtra)
library(stringr)
library(sinaplot)
library(ggforce)
require(cowplot)

plot_freq <- function(df,df_aa,pos,gene){
  print (c(gene, pos))
  df_aa <- select(df_aa,Id,pos) %>%
             filter(grepl(gene,Id)) %>%
             separate(Id, into=c("gene_id", "allele"), sep='\\*') 
  df <- filter(df, gene_id==gene) %>%
          inner_join(df_aa, id=c("gene_id", "allele")) %>%
          select(-X1)
  print (df)
  names(df) <- c('gene_id','allele','percentage','aa')
  df <- df %>%
          arrange(aa, allele) %>%
          mutate(allele=factor(allele,levels=allele)) %>%
          mutate(aa=paste(aa,pos,sep=''))
  ymax <- max(df$percentage)+10
  g_width  <- length(df$allele)/9.5+0.57
  #g_width  <- length(df$allele)/9.5+1
  palette  <- c(brewer.pal(3,"Accent"))
  textsize <- 8
  p <-  ggplot(df,aes(x=allele, y=percentage,fill=aa)) +
          geom_bar(stat='identity', width=0.6, position=position_dodge()) +
          theme_cowplot(12) +
          scale_fill_manual(values=palette,drop=FALSE) +
          theme(axis.text=element_text(size=textsize-1,face="bold",colour = 'black'),
                axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,colour = 'black'),
                axis.text.y=element_text(hjust=1,vjust=0.5,colour = 'black'),
                axis.title.x=element_text(size=textsize,face="bold"),
                axis.title.y=element_text(size=textsize,face="bold"),
                axis.line = element_line(colour = 'black', size = 0),
                panel.border = element_rect(colour = "black", fill=NA, size=1),
                legend.position = "top",
                legend.margin = margin(l = -8, t = 2, b = -10, unit = "pt"), 
                #legend.margin = margin(l = -8, unit = "pt"),
                legend.spacing.x = unit(3, 'pt'),
                legend.title    = element_blank(),
                legend.text=element_text(size=textsize,face="bold"),
                legend.justification='center',
                legend.key.size = unit(0.5,"line")) +
          ylab("Frequency (%)") +
          xlab(paste(gene,"allele",sep=' '))
          ylim(0,ymax)
  ggsave(paste('graph/allele_usage_',gene,'_resi',pos,'.png',sep=''),p,width=g_width, height=1.5, bg='white', dpi=1200)
}

df <- read_csv('results/baseline_variation/baseline_variation_table_1.csv')
df_aa_H <- read_csv('data/anarci_igv_output.csv_H.csv')
df_aa_KL <- read_csv('data/anarci_igv_output.csv_KL.csv')
plot_freq(df, df_aa_KL, '32', 'IGLV2-14')
plot_freq(df, df_aa_KL, '50', 'IGKV1-5')
plot_freq(df, df_aa_H, '50', 'IGHV1-2')
plot_freq(df, df_aa_H, '50', 'IGHV3-30')
plot_freq(df, df_aa_H, '52', 'IGHV3-30')
plot_freq(df, df_aa_H, '50', 'IGHV1-69')
plot_freq(df, df_aa_H, '52', 'IGHV3-23')
plot_freq(df, df_aa_H, '52', 'IGHV3-33')
plot_freq(df, df_aa_H, '50', 'IGHV4-4')
plot_freq(df, df_aa_H, '32', 'IGHV4-4')
plot_freq(df, df_aa_H, '31', 'IGHV3-33')
plot_freq(df, df_aa_KL, '50', 'IGLV3-21')
