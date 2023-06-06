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

plot_DDG_vs_DDG <- function(df, graphname){
    textsize <- 8
    p <- ggplot() +
           geom_point(data=df, aes(x=DDG, y=DDG_antibody,color=dssp_asa),
                      pch=16, size=1, alpha=0.7, stroke=0.2) +
           geom_abline(intercept=0 , slope=1, lty='11', col='black') +
           scale_color_gradientn(colours=c("#F7FCF0","#7BCCC4",'#084081'),
                values=rescale(c(0,0.5,1)),
                limits=c(0,1),
                breaks=c(0,0.5,1),
                labels=c('0','0.5','1'),
                guide="colorbar",
                na.value="grey") +
           theme_cowplot(12) +
           theme(axis.title=element_text(size=textsize,face="bold"),
                 axis.text=element_text(size=textsize,face="bold",hjust=0.5),
                 legend.title=element_blank(),
                 legend.key.size=unit(0.1,'in'),
                 legend.spacing.x=unit(0.03, 'in'),
                 legend.text=element_text(size=textsize,face="bold"),
                 legend.position='top') +
           labs(x=bquote(bold(Predicted~ΔΔG["complex"]~'(kcal/mol)')),y=bquote(bold(Predicted~ΔΔG["apo antibody"]~'(kcal/mol)')))
    ggsave(graphname,p,height=2.8,width=2.5,dpi=600, bg='white')
    }

plot_DDG_sina <- function(df, graphname){
  textsize <- 8
  p <- ggplot(data=df, aes(x='', y=DDG_binding, color=known, alpha=known, size=known,group=1)) +
           geom_sina(pch=16, method="counts", bin_limit=0.4, scale="width",
                     maxwidth=1) +
           scale_color_manual('known', values=c('grey30','red'),drop=FALSE) +
           scale_alpha_manual(name='known', values=c(0.5, 0.8)) +
           scale_size_manual(name='known', values=c(0.6, 1)) +
           theme_cowplot(12) +
           theme(axis.title=element_text(size=textsize,face="bold"),
                 axis.text=element_text(size=textsize,face="bold"),
                 axis.ticks.x = element_blank(),
                 legend.title=element_blank(),
                 legend.key.size=unit(0.1,'in'),
                 legend.spacing.x=unit(0.03, 'in'),
                 legend.text=element_text(size=textsize,face="bold"),
                 legend.position='none') +
           labs(x='', y=bquote(bold(Predicted~ΔΔG["binding"]~'(kcal/mol)')))
    ggsave(graphname,p,height=2.5,width=3,dpi=600, bg='white')
  }

plot_resolution_vs_DDG <- function(df, graphname){
  textsize <- 8
  p <- ggplot() +
	 geom_point(data=df, aes(x=resolution, y=DDG_binding),
		    pch=16, size=1, alpha=0.5, fill='black') +
	 theme_cowplot(12) +
	 theme(axis.title=element_text(size=textsize,face="bold"),
	       axis.text=element_text(size=textsize,face="bold",hjust=0.5),
	       legend.title=element_blank(),
	       legend.key.size=unit(0.1,'in'),
	       legend.spacing.x=unit(0.03, 'in'),
	       legend.text=element_text(size=textsize,face="bold"),
	       legend.position='none') +
	 labs(x=expression(bold('Resolution (Å)')),y=bquote(bold(Predicted~ΔΔG["binding"]~'(kcal/mol)')))
  ggsave(graphname,p,height=2,width=2,dpi=600, bg='white')
  }

set.seed(5)
df <- read_csv('results/allele_var_info_table_final.csv') %>%
        mutate(DDG=as.numeric(DDG)) %>%
        mutate(DDG_antibody=as.numeric(DDG_antibody)) %>%
        mutate(DDG=if_else(DDG > 10, 10, DDG)) %>%
        mutate(DDG=if_else(DDG < -5, -5, DDG)) %>%
        mutate(DDG_antibody=if_else(DDG_antibody > 10, 10, DDG_antibody)) %>%
        mutate(DDG_antibody=if_else(DDG_antibody < -5, -5, DDG_antibody))
plot_DDG_vs_DDG(df, 'graph/ddg_scatterplot.png')

df <- df %>%
        mutate(DDG_binding=if_else(DDG_binding > 10, 10, DDG_binding)) %>%
        mutate(DDG_binding=if_else(DDG_binding < -5, -5, DDG_binding)) %>%
        arrange(known)
plot_DDG_sina(df, 'graph/ddg_distribution.png')

print (cor(df$DDG_binding, df$resolution,method='spearman'))
print ("# of allelic variants:")
print (length(df$DDG))
print ("# of allelic variants with DDG > 2.5:")
print (length(filter(df, DDG_binding > 2.5)$DDG))
print ("# of allelic variants with 2.5 > DDG > 0:")
print (length(filter(df, DDG_binding > 0 & DDG_binding <= 2.5)$DDG))
print ("# of allelic variants with DDG > 0:")
print (length(filter(df, DDG_binding > 0)$DDG))
print ("% of allelic variants with DDG > 0:")
print (length(filter(df, DDG_binding > 2.5)$DDG)/length(df$DDG))
print ("% of allelic variants with 2.5 > DDG > 0:")
print (length(filter(df, DDG_binding > 0 & DDG_binding <= 2.5)$DDG)/length(df$DDG))
print ("Total # of PDBs with paratopic allelic variants")
print (length(unique(df$pdb_id)))
print ("Total # of PDBs with paratopic allelic variants with DDG > 2.5")
print (length(unique(filter(df, DDG_binding > 2.5)$pdb_id)))
print ("Total # of PDBs with paratopic allelic variants with 2.5 > DDG > 0")
print (length(unique(filter(df, DDG_binding > 0 & DDG_binding <= 2.5)$pdb_id)))

df_low_freq  <- filter(df, baseline_freq_original<50)
df_high_freq <- filter(df, baseline_freq_original>50)

print ('Total # of allelic variants with low freq:')
print (length(df_low_freq$DDG))
print ('Total # of allelic variants with high freq:')
print (length(df_high_freq$DDG))
print (length(filter(df_low_freq, DDG_binding > 2.5)$DDG)/length(df_low_freq$DDG))
print (length(filter(df_high_freq, DDG_binding > 2.5)$DDG)/length(df_high_freq$DDG))

plot_resolution_vs_DDG(df, 'graph/ddg_vs_resolution.png')
