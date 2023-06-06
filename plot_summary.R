#R code
library(ggplot2)
library(scales)
library(RColorBrewer)
library(qualpalr)
library(readr)
library(tidyr)
library(reshape)
library(stringr)
library(data.table)
library(dplyr)
library(gridExtra)
library(stringr)
require(cowplot)

plot_antigen_count <- function(df, graphname){
  textsize <- 7
  df <- df %>%
          group_by(antigen_abbrev) %>%
          summarise(count = n())
  p <-  ggplot(df,aes(x=antigen_abbrev, y=count)) +
          geom_bar(stat='identity', width=0.7, position=position_dodge()) +
          theme_cowplot(12) +
          theme(axis.text=element_text(size=textsize-1,face="bold",colour = 'black'),
                axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,colour = 'black'),
                axis.text.y=element_text(hjust=1,vjust=0.5,colour = 'black'),
                axis.title.x=element_text(size=textsize,face="bold"),
                axis.title.y=element_text(size=textsize,face="bold"),
                axis.line = element_line(colour = 'black', size = 0),
                panel.border = element_rect(colour = "black", fill=NA, size=1),
                legend.position = "top",
                legend.margin = margin(t = 2, b = -10, unit = "pt"),
                legend.spacing.x = unit(3, 'pt'),
                legend.title    = element_blank(),
                legend.text=element_text(size=textsize,face="bold"),
                legend.justification='center',
                legend.key.size = unit(0.5,"line")) +
          ylab("# of complex structures") +
          xlab('')
  ggsave(graphname,p,width=4.5, height=2, bg='white', dpi=600)
  }

plot_resolution <- function(df, graphname){
  #library(sinaplot)
  library(ggforce)
  print (max(df$resolution))
  print (min(df$resolution))
  print (median(df$resolution))
  textsize <- 8
  p <- ggplot() +
           geom_sina(data=df, aes(x=1, y=resolution), pch=16, size=0.5, alpha=0.5,
                     bin_limit=0.4, scale="width", maxwidth=0.8, color='black') +
           geom_boxplot(data=df, aes(x=1, y=resolution), width=0.3, color="black",
                        outlier.shape=NA, alpha=0, lwd=0.3) +
           theme_cowplot(12) +
           theme(axis.title=element_text(size=textsize,face="bold"),
                 axis.text=element_text(size=textsize,face="bold"),
                 axis.text.x=element_blank(),
                 axis.ticks.x = element_blank(),
                 legend.title=element_blank(),
                 legend.key.size=unit(0.1,'in'),
                 legend.spacing.x=unit(0.03, 'in'),
                 legend.text=element_text(size=textsize,face="bold"),
                 legend.position='none') +
           labs(x='', y=expression(bold('Resolution (Å)'))) +
           scale_y_continuous(breaks = seq(1,7), labels = seq(1,7))
  ggsave(graphname,p,height=2,width=2,dpi=600, bg='white')
  }

antigen_class_levels <- c("Virus","Bacteria","Fungi","Parasite","Animal","Plant","Synthetic")
df <- read_csv('results/allele_var_info_table_final.csv') %>%
        select(pdb_id, antigen_abbrev, antigen_class, resolution) %>%
        distinct() %>%
        mutate(antigen_class=factor(antigen_class,levels=antigen_class_levels)) %>%
        arrange(antigen_class, antigen_abbrev) %>%
        mutate(antigen_abbrev=factor(antigen_abbrev,levels=unique(antigen_abbrev))) 

plot_antigen_count(df, 'graph/summary_antigen.png')
plot_resolution(df, 'graph/summary_resolution.png')
print (df)
