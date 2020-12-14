#TrAELseq data visualisation

#For line plots displaying read counts over specified regions

library(tidyverse)

setwd("path_to_Seqmonkfile")

#Seqmonk file contains read counts per window for each sample in separate columns
read_tsv("SeqmonkFile_withReadCounts") -> count.table

#Specify sample name for file output
#Set region positions to analyse e.g. chrom III nuc position 100-500kb
sample.name <- "All samples"
region.start <- 100
region.end <- 500
chrom <- "III"

#Line plot display for multiple samples
count.table %>% 
  mutate(Start = Start/1000, End = End/1000) %>% #make start position in kb coordinates
  mutate(Center = (Start+End)/2) %>% #center of quantification window
  filter(Chromosome == chrom, between(Center, region.start , region.end)) %>%
  ggplot(aes(x = Center)) +
  geom_line( aes(y = sample1, colour = "sample1"), size =0.4 , lineend = "round" ) +
  geom_line( aes(y = sample2, colour = "sample2"), size =0.4, lineend = "round" ) +
  geom_line( aes(y = sample3, colour = "sample3"), size =0.4, lineend = "round" ) +
  scale_colour_brewer( name = "sample", palette = "Set1")+
  labs(y= "Raw read counts") +
  scale_y_reverse(expand = c(0,0))+
  coord_cartesian(ylim = c(100,0)) + #change limits
  scale_x_continuous(position = "top") +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.margin = margin(0.3,0.3,0.3,0.3,"cm"),
    axis.line = element_line(size =0.5),
    axis.ticks.length = unit(1, "mm"),
    axis.text.y = element_text(size = 8, colour = "black"),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(vjust = 3, size = 9),
    aspect.ratio = 0.7
  )

#save to file according sample.name above
ggsave(paste(toString(sample.name), "ReadCountLinePlot.pdf"), 
       width = 6,
       height = 4,
       units = "cm",
       path = "path_to_file",
       device = "pdf")