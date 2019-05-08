#' ---
#' title: "Analysis of NP-PCR amplicon sequencing results"
#' author: "Chun-Hui Gao"
#' output: html_document
#' date: "`r Sys.Date()`"
#' ---

#' ## Import, Transform and Tidy
#' 
library(tidyverse)
library(ggpubr)
library(ggsci)
library(ggrepel)
library(cowplot)

# error rate and snp frequency were estimated using amplicon sequencing and mothur software
error_rate <- read_csv("data/error_rate.csv")
snp_freq <- read_csv("data/snp_freq.csv")

# read meta
meta <- read_tsv("data/meta.txt")

#+ label of plot
# NP Labels
np_labels <- c("CK"    = expression(paste(CK,"+")),
               "Fe2O3" = expression(paste(Fe[2],O[3])),
               "ZnO"   = expression(paste(ZnO)),
               "CeO2"  = expression(paste(Ce,O[2])),
               "Fe3O4" = expression(paste(Fe[3],O[4])),
               "Al2O3" = expression(paste(Al[2],O[3])),
               "CuO"   = expression(paste(CuO)),
               "TiO2"  = expression(paste(Ti,O[2])))
# color value
color_values <- pal_npg("nrc")(length(np_labels))
names(color_values) <- names(np_labels)

# DNA polymerase Labels
enzyme_labels <- c("extaq" = "Ex Taq","phusion" = "Phusion")

#+ function
# enzyme labeller
enzyme_labeller <- function(enzyme){
  return(as.character(enzyme_labels[enzyme]))
}
# Y-axis Labels
fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # remove +
  l <- gsub("\\+","",l)
  # return this as an expression
  parse(text=l)
}

meta$nps <- factor(meta$nps, levels = names(np_labels))

# join data
error_rate <- left_join(error_rate,meta)
snp_freq <- left_join(snp_freq,meta)

#' ## Data Visualization


#' Overall Error Rate - comparing with CK
#+ overall error rate, fig.width=6,fig.asp=0.618
plots <- lapply(c("extaq","phusion"), function(x){
  data <- filter(error_rate,enzyme==x)
  ggplot(data, aes(nps,error_rate,color=nps)) + 
    geom_boxplot(show.legend = F) +
    geom_jitter(show.legend=F) + 
    stat_compare_means(label = "p.signif",hide.ns = T,ref.group = "CK",show.legend = F,method = "t.test") +
    scale_x_discrete(labels=np_labels) +   
    scale_y_continuous(labels=fancy_scientific) +
    scale_color_manual(values = color_values,name="NPs") +
    labs(title=enzyme_labeller(x)) +
    xlab("Nanoparticles") + ylab("Overall Error Rate") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face = 'bold.italic'),
          axis.line.x = element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y = element_text(margin = margin(t=0,r=15,b=0,l=0)),
          legend.text.align = 0) # 0-left,1-right
})

plot_grid(plotlist = plots,labels = "AUTO")

ggsave("Figure 4.png",dpi = 300,path = "plots")

#' SNV Frequency - comparing with CK
#' For extaq
#+ extaq snv, fig.width=8,fig.asp=0.5
plots <- lapply(c("extaq","phusion"), function(x){
  data <- filter(snp_freq,enzyme==x)
  ggplot(data,aes(nps,freq,color=nps)) + 
    geom_boxplot(show.legend=F,outlier.alpha = 0) + 
    geom_jitter(show.legend=F) + 
    #stat_compare_means(label = "p.signif",hide.ns = T,ref.group = "CK",show.legend = F,method = "t.test") +
    facet_wrap(~snv,ncol=6)  + 
    scale_x_discrete(labels=np_labels) +
    scale_y_continuous(labels=fancy_scientific) +
    scale_color_manual(values = color_values,name="NPs") +
    labs(title=enzyme_labeller(x)) +
    xlab("Nanoparticles") + ylab("Frequency of SNVs") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face = 'bold.italic'),
          axis.line.x = element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y = element_text(margin = margin(t=0,r=15,b=0,l=0)),
          legend.text.align = 0) # 0-left,1-right
})

#' Put together
#+ cowplot, fig.width=8,fig.asp=1
plot_grid(plotlist = plots,ncol=1,align = "v",labels="AUTO")
ggsave("Figure 5.png",dpi = 300,path = "plots")