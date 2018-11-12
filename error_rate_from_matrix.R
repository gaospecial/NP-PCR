#' ---
#' title: "Analysis of NP-PCR amplicon sequencing results"
#' author: "Chun-Hui Gao"
#' output: html_document
#' date: "`r Sys.Date()`"
#' ---

#' ## Import, Transform and Tidy
#' 
require(tidyr)
require(dplyr)
require(tibble)
require(readr)

# error rate and snp frequency were estimated using amplicon sequencing and mothur software
error_rate <- read_csv("data/error_rate.csv")
snp_freq <- read_csv("data/snp_freq.csv")

# read meta
meta <- read_tsv("data/meta.txt")
meta$nps <- factor(meta$nps, levels = c("CK","ZnO","Fe2O3","CeO2","Fe3O4","Al2O3","CuO","TiO2"))
meta$enzyme <- factor(meta$enzyme,levels = c("extaq","phusion"))

# join data
error_rate <- left_join(error_rate,meta)
snp_freq <- left_join(snp_freq,meta)

#' ## Data Visualization

#+ load required gg environment, message=F
require(ggplot2)
require(ggpubr)
require(ggsci)

#+ label of plot
# NP Labels
np_labels <- c("Fe2O3" = expression(paste(Fe[2],O[3])),
               "ZnO"   = expression(paste(ZnO)),
               "CeO2"  = expression(paste(Ce,O[2])),
               "Fe3O4" = expression(paste(Fe[3],O[4])),
               "Al2O3" = expression(paste(Al[2],O[3])),
               "CuO"   = expression(paste(CuO)),
               "TiO2"  = expression(paste(Ti,O[2])),
               "CK"    = expression(paste(CK,"+")))
# DNA polymerase Labels
enzyme_labels <- c("extaq" = "TaKaRa Ex Taq","phusion" = "ThermoFisher Phusion")

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


#' Overall Error Rate - comparing with CK
#+ overall error rate, fig.width=6,fig.asp=0.618
p1 <- filter(error_rate,enzyme=="extaq") %>% ggplot(aes(nps,error_rate,color=nps)) + 
  geom_boxplot(show.legend=F,outlier.alpha = 0) + 
  geom_jitter(show.legend=F) + 
  stat_compare_means(label = "p.signif",hide.ns = T,ref.group = "CK",show.legend = F,method = "t.test") +
  scale_x_discrete(labels=np_labels) +   
  scale_y_continuous(labels=fancy_scientific) +
  scale_color_npg() +
  labs(title="Ex Taq") +
  xlab("Nanoparticles") + ylab("Overall Error Rate") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face = 'bold.italic'),
        axis.line.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y = element_text(margin = margin(t=0,r=15,b=0,l=0)),
        legend.text.align = 0) # 0-left,1-right
p2 <- filter(error_rate,enzyme=="phusion") %>% ggplot(aes(nps,error_rate,color=nps)) + 
  geom_boxplot(show.legend=F,outlier.alpha = 0) + 
  geom_jitter(show.legend=F) + 
  stat_compare_means(label = "p.signif",hide.ns = T,ref.group = "CK",show.legend = F,method = "t.test") +
  scale_x_discrete(labels=np_labels) +   
  scale_y_continuous(labels=fancy_scientific) +
  scale_color_npg() +
  labs(title="Phusion") +
  xlab("Nanoparticles") + ylab("Overall Error Rate") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face = 'bold.italic'),
        axis.line.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y = element_text(margin = margin(t=0,r=15,b=0,l=0)),
        legend.text.align = 0) # 0-left,1-right
cowplot::plot_grid(p1,p2,labels = "AUTO")
ggsave("Figure 4.png",dpi = 300,path = "plots")

#' SNV Frequency - comparing with CK
#' For extaq
#+ extaq snv, fig.width=8,fig.asp=0.5
p1 <- ggplot(filter(snp_freq,enzyme=="extaq"),aes(nps,freq,color=nps)) + 
  geom_boxplot(show.legend=F,outlier.alpha = 0) + 
  geom_jitter(show.legend=F) + 
  #stat_compare_means(label = "p.signif",hide.ns = T,ref.group = "CK",show.legend = F,method = "t.test") +
  facet_wrap(~snv,ncol=6)  + 
  scale_x_discrete(labels=np_labels) +
  scale_y_continuous(labels=fancy_scientific) +
  scale_color_npg() +
  labs(title="Ex Taq") +
  xlab("Nanoparticles") + ylab("Frequency of SNVs") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face = 'bold.italic'),
        axis.line.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y = element_text(margin = margin(t=0,r=15,b=0,l=0)),
        legend.text.align = 0) # 0-left,1-right

#' For Phusion
#+ phusion snv, fig.width=8,fig.asp=0.5
p2 <- ggplot(filter(snp_freq,enzyme=="phusion"),aes(nps,freq,color=nps)) + 
  geom_boxplot(show.legend=F,outlier.alpha = 0) + 
  geom_jitter(show.legend=F) + 
  #stat_compare_means(label = "p.signif",hide.ns = T,ref.group = "CK",show.legend = F,method = "t.test") +
  facet_wrap(~snv,ncol=6)  + 
  scale_x_discrete(labels=np_labels) +
  scale_y_continuous(labels=fancy_scientific) +
  scale_color_npg() +
  labs(title="Phusion") +
  xlab("Nanoparticles") + ylab("Frequency of SNVs") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face = 'bold.italic'),
        axis.line.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y = element_text(margin = margin(t=0,r=15,b=0,l=0)),
        legend.text.align = 0) # 0-left,1-right
#' Put together
#+ cowplot, fig.width=8,fig.asp=1
cowplot::plot_grid(p1,p2,ncol=1,align = "v",labels="AUTO")
ggsave("Figure 5.png",dpi = 300,path = "plots")