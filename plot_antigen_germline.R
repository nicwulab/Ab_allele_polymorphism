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

plot_germline <- function(df, graphname){
  palette  <- c(brewer.pal(3,"Accent"))
  textsize <- 8
  p <- ggplot(data=df, aes(x=gene, y=DDG_binding, fill=gene_class, color=gene_class)) +
           geom_sina(pch=16, size=0.5, alpha=0.5, bin_limit=0.4, scale="width", maxwidth=0.8, color='grey') +
           geom_boxplot(width=0.3, color="black", outlier.shape=NA, alpha=0, lwd=0.3) +
           geom_abline(intercept=0 , slope=0, lty='11', col='black',lwd=0.3) +
           #scale_fill_manual(values=palette,drop=FALSE) +
           #scale_color_manual(values=palette,drop=FALSE) +
           theme_cowplot(12) +
           theme(axis.title=element_text(size=textsize,face="bold"),
                 axis.text=element_text(size=textsize,face="bold"),
                 axis.text.x=element_text(size=textsize,face="bold",angle=90,hjust=0.5,vjust=0.5),
                 legend.title=element_blank(),
                 legend.key.size=unit(0.1,'in'),
                 legend.spacing.x=unit(0.03, 'in'),
                 legend.text=element_text(size=textsize,face="bold"),
                 legend.position='none') +
           labs(x='', y=expression(bold('ΔΔG of binding (kcal/mol)')))
  ggsave(graphname,p,height=2,width=6,dpi=600, bg='white')
  }

plot_antigen <- function(df, graphname){
  palette  <- c(brewer.pal(7,"Set2"))
  textsize <- 8
  p <- ggplot(data=df, aes(x=antigen_abbrev, y=DDG_binding, color=antigen_class)) +
           geom_sina(pch=16, size=0.4, alpha=0.7, bin_limit=0.4, scale="width", maxwidth=0.8, color='grey') +
           geom_boxplot(width=0.3, color="black", outlier.shape=NA, alpha=0, lwd=0.3) +
           geom_abline(intercept=0 , slope=0, lty='11', col='black',lwd=0.3) +
           #scale_fill_manual(values=palette,drop=FALSE) +
           #scale_color_manual(values=palette,drop=FALSE) +
           theme_cowplot(12) +
           theme(axis.title=element_text(size=textsize,face="bold"),
                 axis.text=element_text(size=textsize,face="bold"),
                 axis.text.x=element_text(size=textsize,face="bold",angle=90,hjust=0.5,vjust=0.5),
                 legend.title=element_blank(),
                 legend.key.size=unit(0.1,'in'),
                 legend.spacing.x=unit(0.03, 'in'),
                 legend.text=element_text(size=textsize,face="bold"),
                 legend.position='none') +
           labs(x='', y=expression(bold('ΔΔG of binding (kcal/mol)')))
  ggsave(graphname,p,height=2.2,width=6,dpi=600, bg='white')
  }

factoring_antigen <- function(abbrev, antigen_levels){
  if (abbrev %in% antigen_levels){
    return (abbrev)
    }
  else {
    return ('Others')
    }
  }

plot_gene_specific_mut <- function(df, gene_name, graphname, antigen_count_table){
  palette  <- c(brewer.pal(7,"Set2"))
  textsize <- 8
  palette <- qualpal(n = 11, list(h = c(0, 360), s = c(0.4, 0.6), l = c(0.5, 0.85)))$hex
  #palette <- c(palette[1:3], palette[5:13])
  antigen_levels <- c('Sarbecovirus','HIV','Human','Influenza','Plasmodium','Yeast','HCV','Chicken','YFV','CMV','Others')
  df <- filter(df, gene==gene_name) %>%
          arrange(location, mut) %>%
          mutate(mut=factor(mut, levels=unique(mut))) %>%
          mutate(antigen_abbrev=mapply(factoring_antigen, as.character(antigen_abbrev),
                 MoreArgs=list(antigen_levels))) %>%
          mutate(antigen_abbrev=factor(antigen_abbrev, antigen_levels))
  if (gene == 'IGHV4-4'){
    df <- filter(df, mut!='S32W')
    }
  antigen_list <- as.character(df$antigen_abbrev)
  antigen_count_table <- c(antigen_count_table, antigen_list)
  g_width <- (length(unique(as.character(df$mut))))/6 + 0.57
  if (gene == 'IGHV3-23'){
    g_height <- 1.55
    }
  else{
    g_height <- 1.5
    }
  p <- ggplot(data=df, aes(x=mut, y=DDG_binding,color=antigen_abbrev,group=mut)) +
           geom_boxplot(width=0.3, color="black", outlier.shape=NA, alpha=0, lwd=0.2) +
           geom_sina(pch=16, size=0.5, alpha=0.8, bin_limit=0.4, scale="width", maxwidth=0.6) +
           geom_abline(intercept=0 , slope=0, lty='11', col='black',lwd=0.3) +
           #scale_fill_manual(values=palette,drop=FALSE) +
           scale_color_manual(values=palette,drop=FALSE) +
           theme_cowplot(12) +
           theme(plot.title=element_text(size=textsize,face="bold",hjust=0.5),
                 axis.title=element_text(size=textsize,face="bold"),
                 axis.text=element_text(size=textsize,face="bold"),
                 axis.text.x=element_text(size=textsize,face="bold",angle=90,hjust=0.5,vjust=0.5),
                 legend.title=element_blank(),
                 legend.key.size=unit(0.12,'in'),
                 legend.spacing.x=unit(0.03, 'in'),
                 legend.text=element_text(size=textsize,face="bold"),
                 legend.position='none') +
           labs(x='', y=expression(bold('ΔΔG of binding (kcal/mol)')), title=gene_name) +
           ylim(-5,10) +
           guides(color = guide_legend(override.aes = list(size = 1.2)))
  #ggsave(graphname,p,height=g_height, width=g_width,dpi=600, bg='white')
  ggsave(graphname,p,height=2, width=g_width,dpi=600, bg='white')
  return (antigen_count_table)
  }

set.seed(6)
antigen_class_levels <- c("Virus","Bacteria","Fungi","Parasite","Animal","Plant","Synthetic")
df <- read_csv('results/allele_var_info_table_final.csv') %>%
        mutate(gene_class=ifelse(grepl('IGHV',gene), 'IGHV', ifelse(grepl('IGKV',gene),'IGKV','IGLV'))) %>%
        mutate(mut=paste(amino_acid_original,location,variants,sep='')) %>%
        mutate(DDG_binding=if_else(DDG_binding > 10, 10, DDG_binding)) %>%
        mutate(DDG_binding=if_else(DDG_binding < -5, -5, DDG_binding)) %>%
        mutate(antigen_class=factor(antigen_class,levels=antigen_class_levels)) %>%
        arrange(antigen_class, antigen_abbrev) %>%
        mutate(antigen_abbrev=factor(antigen_abbrev,levels=unique(antigen_abbrev))) 

plot_germline(df, 'graph/DDG_vs_germline.png')
plot_antigen(df, 'graph/DDG_vs_antigen.png')

gene_list <- c('IGHV1-2', 'IGHV1-69', 'IGHV3-23', 'IGHV3-30', 'IGHV3-33','IGHV4-4',
               'IGHV4-34','IGHV4-39','IGHV4-59','IGLV2-14','IGKV1-5','IGLV3-21')
antigen_count_table <- c()
for (gene in gene_list){
  graphname <- paste('graph/DDG_vs_mut_',gene,'.png',sep='')
  antigen_count_table <- plot_gene_specific_mut(df, gene, graphname, antigen_count_table)
  }

#print (sort(table(antigen_count_table)))
gene_count_total <- length(unique(df$gene))
gene_count_DDG_0 <- length(unique(filter(df,DDG_binding>0)$gene))
gene_count_DDG_25 <- length(unique(filter(df,DDG_binding>2.5)$gene))
#print (c(gene_count_total, gene_count_DDG_0, gene_count_DDG_25))

df_virus <- filter(df, antigen_class == 'Virus')
df_human <- filter(df, antigen_abbrev == 'Human')
df_bacteria <- filter(df, antigen_class == 'Bacteria')
df_virus_DDG_0 <- length(filter(df_virus, DDG_binding>0)$DDG_binding)
df_human_DDG_0 <- length(filter(df_human, DDG_binding>0)$DDG_binding)
df_bacteria_DDG_0 <- length(filter(df_bacteria, DDG_binding>0)$DDG_binding)
df_virus_DDG_25 <- length(filter(df_virus, DDG_binding>2.5)$DDG_binding)
df_human_DDG_25 <- length(filter(df_human, DDG_binding>2.5)$DDG_binding)
df_bacteria_DDG_25 <- length(filter(df_bacteria, DDG_binding>2.5)$DDG_binding)
print (c(length(df_virus$DDG_binding), df_virus_DDG_0, df_virus_DDG_25))
print (c(length(df_human$DDG_binding), df_human_DDG_0, df_human_DDG_25))
print (c(length(df_bacteria$DDG_binding), df_bacteria_DDG_0, df_bacteria_DDG_25))
