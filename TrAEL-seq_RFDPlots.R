#TrAELseq data visualisation

#Creates RFD plots (replication fork directionality) using ggplot2 with Seqmonk files

library(tidyverse)

setwd("path_to_Seqmonkfiles")

#Select chromosome and region to analyse e.g. yeast chrom V 0-600kb
region.start <- 0
region.end <- 600
chrom <- "V"

#Read in Seqmonk file for forward counts per window
#Specify the center of each quantification window
read_tsv("SeqmonkFile_ForwardCounts") %>%
  mutate(Start = Start/1000, End = End/1000) %>% #make start position in kb sizes
  mutate(Center = (Start+End/2)) %>%
  filter(Chromosome == chrom, between(Center, region.start , region.end)) %>%
  rename(sample1.fwd = sample1, sample2.fwd = sample2, sample3.fwd = sample3) -> fwd.counts

#Read file for reverse counts and combine with forward table
#Calculate RFD by (rev-fwd)/(rev+ fwd)
read_tsv("SeqmonkFile_ReverseCounts") %>%
  mutate(Start = Start/1000, End = End/1000) %>% 
  mutate(Center = (Start+End/2)) %>%
  filter(Chromosome == chrom, between(Center, region.start , region.end)) %>%
  rename(sample1.rev = sample1, sample2.rev = sample2, sample3.rev = sample3) %>%
  left_join(fwd, by = c("Start", "End", "Chromosome")) %>%
  mutate(sample1.RFD = (sample1.rev - sample1.fwd)/
           (sample1.rev + sample1.fwd)) %>%
  mutate(sample2.RFD = (sample2.rev - sample2.fwd)/
           (sample2.rev + sample2.fwd)) %>%
  mutate(sample3.RFD = (sample3.rev - sample3.fwd)/
           (sample3.rev + sample3.fwd)) -> RFD.table

#For adding on ARS sites to plots in yeast
read_tsv("ARS_file.txt") %>%
  select(Chr, Start, End) %>%
  filter(Chr == chrom) %>%
  mutate(Site = ((Start+End)/2)/1000)-> ARS.sites

ARS.sites$Site -> sites

#Change sample file name and column to plot
sample <- RFD.table$sample1.RFD
sample.name <- "sample1"

#Adds column to dictate colour grouping
RFD.table %>%
  mutate(colour = ifelse(sample >0, "red", "blue"))-> RFD.table

#Displays a RFD dot plot with red/blue scale colouring
RFD.table %>%
  ggplot(aes(x = Center, colour = colour)) +
  scale_color_manual(values = c("dodgerblue2", "firebrick3")) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = sites, linetype = "dotted", size = 0.2, colour = "grey20") +
  geom_point(aes(y= sample), size =0.1) +
  labs(y= "RFD", x = paste("Chromosome", chrom, "position (kb)")) +
  scale_x_continuous(labels = comma, expand = c(0,0))+
  ylim(c(-1, 1)) +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.margin = margin(0.3,0.3,0.3,0.3,"cm"),
    axis.text = element_text(size = 8, colour = "black"),
    axis.title.x = element_text(vjust = -1, size = 9),
    axis.title.y = element_text(vjust = 1, size = 9),
    aspect.ratio = 0.3
  )

#Save to file according to input sample.name
ggsave(paste(toString(sample.name), "_RFDPlot.pdf"),
       width = 14,
       units = "cm",
       path = "path_to_file",
       device = "pdf")

#--------------------------------------------------------------

#Alternative Plot - plots RFD.table with multiple samples on one plot (as lines)

sample.name <- "All samples"

RFD.table %>%
  ggplot(aes(x = Center)) +
  geom_hline(yintercept = 0) +
  geom_line(aes(y= sample1.RFD), size =0.2, colour = "firebrick3") +
  geom_line(aes(y=sample2.RFD), size = 0.2, colour = "darkorange1") +
  geom_line(aes(y=sample3.RFD), size = 0.2, colour = "darkorchid4") +
  geom_vline(xintercept = sites, linetype = "dotted", size = 0.2, colour = "grey20") +
  labs(y= "RFD", x = paste("Chromosome", chrom, "position (kb)")) +
  scale_x_continuous(labels = comma, expand = c(0,0))+
  ylim(c(-1, 1)) +
  theme_classic() +
  theme(
    plot.margin = margin(0.3,0.3,0.3,0.3,"cm"),
    axis.line = element_line(size =0.3),
    axis.text = element_text(size = 8, colour = "black"),
    axis.title.x = element_text(vjust = -1, size = 9),
    axis.title.y = element_text(vjust = 1, size = 9),
    aspect.ratio = 0.3 
  )

#Save to file according to input sample.name
ggsave(paste(toString(sample.name), "_RFDLinePlot.pdf"),
       width = 14,
       units = "cm",
       path = "path_to_file",
       device = "pdf")
