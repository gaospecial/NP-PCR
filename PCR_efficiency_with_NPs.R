#' ---
#' title: "NPs modulate PCR efficiency"
#' author: "Chun-Hui Gao"
#' date: "2018-10-12"
#' output: "html_document"
#' ---

#' # Read data
#'
#+ data, message=F
require(dplyr)
require(tibble)
require(readr)

results <- read_csv("data/RT-PCR_results.csv")
amplification <- read_csv("data/RT-PCR_amplification.csv")

#' # Visualization

#+ load required gg environment, message=F
require(ggplot2)
require(ggpubr)
require(ggsci)
library(ggrepel)

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
enzyme_labels <- c("extaq" = "Ex Taq","phusion" = "Phusion")

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

#' ## Amplication curve
#+ plot amplication curve, fig.width=10,fig.asp=0.618
amplification <- amplification %>%  mutate(label=as.character(np_labels[as.character(nps)]))


amplification_style <-  theme_bw() + theme(legend.position = "none")

p1 <- ggplot(filter(amplification,enzyme=="extaq"),aes(Cycle,mean,color=nps,group=nps)) + 
  geom_line(size=1,show.legend = F) + geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),show.legend = F) +
  labs(title = "Ex Taq") + geom_text_repel(mapping = aes(x=42,label=label,color=nps),
                                           data = filter(amplification,Cycle==40,enzyme=="extaq"),
                                           direction = "y",xlim=c(40,47),segment.alpha = 0,
                                           box.padding = 0.1,
                                           parse = T,show.legend = F,hjust=0) + 
  xlim(c(NA,47)) +  ylab("Rn") +  scale_y_continuous(labels=fancy_scientific) + 
  scale_color_npg(labels=np_labels,name="NPs") +
  amplification_style

p2 <- ggplot(filter(amplification,enzyme=="phusion"),aes(Cycle,mean,color=nps,group=nps)) + 
  geom_line(size=1,show.legend = F) + geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),show.legend = F) +
  labs(title = "Phusion") + geom_text_repel(mapping = aes(x=42,label=label),
                                            data = filter(amplification,Cycle==40,enzyme=="phusion"),
                                            direction = "y",xlim=c(40,47),segment.alpha = 0,
                                            box.padding = 0.1,
                                            parse = T,show.legend = F,hjust=0) + 
  xlim(c(NA,47)) +  ylab("Rn") +  scale_y_continuous(labels=fancy_scientific) + 
  scale_color_npg(labels=np_labels,name="NPs") +
  amplification_style

cowplot::plot_grid(p1,p2,labels = "AUTO")

ggsave("Figure 2.png",path = "plots",dpi=300)


#' ### Comparison of Ct values
#+ fig.width=6,fig.asp=0.618

# plot style
ct_plot_style <-  theme_bw() +
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,face = 'bold.italic'),
        axis.line.x = element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "none")

# left panel
p1 <- results %>% select(enzyme,nps,CT) %>% filter(enzyme=="extaq") %>% ggplot(aes(nps,CT,color=nps)) + geom_boxplot(show.legend = F,outlier.alpha = 0,alpha=I(1/2)) + 
  geom_jitter(show.legend = F,alpha=1/2,width = 0.25) +
  stat_compare_means(ref.group = "CK",label = "p.signif",label.y = 16.5,hide.ns = F,method = "t.test") + 
  #stat_compare_means(label.x = "CK",label.y = 20,method = "anova") +
  labs(title="Ex Taq") +  scale_color_npg() +
  scale_x_discrete(labels=np_labels) +   
  xlab("Nanoparticles") + ylab(expression(italic(C[t]))) + 
  ct_plot_style

# right panel
p2 <- results %>% select(enzyme,nps,CT) %>% filter(enzyme=="phusion") %>% ggplot(aes(nps,CT,color=nps)) + geom_boxplot(show.legend = F,outlier.alpha = 0,alpha=I(1/2)) + 
  geom_jitter(show.legend = F,alpha=1/2,width = 0.25) +
  stat_compare_means(ref.group = "CK",label = "p.signif",label.y = 19,hide.ns = F,method = "t.test") + 
  #stat_compare_means(label.x = "CK",label.y = 20,method = "anova") +  
  labs(title="Phusion") + scale_color_npg() +
  scale_x_discrete(labels=np_labels) +   
  xlab("Nanoparticles") + ylab(expression(italic(C[t]))) + 
  ct_plot_style


cowplot::plot_grid(p1,p2,labels = "AUTO")

ggsave("Figure 3.png",path = "plots",dpi=300)

