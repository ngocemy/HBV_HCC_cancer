# Comparison of loop-related statistics between cancer and healthy libraries
# cmdoret, 20190708

library(tidyverse)

### LOAD DATA ###
args <- commandArgs(trailingOnly=T)
#loops <- args[1]
#samples <- args[2]
#res <- args[3]
loops <- read_tsv('data/output/loops/all_loops.tsv')
samples <- read_tsv('data/input/rnaseq/samples.tsv')
res <- 10000

### CLEAN DATA ###
# Add library metadata to loops table
loops <- loops %>%
    left_join(samples, by=c("library"="sample")) %>%
    mutate(condition = factor(condition))


### COMPUTE STATS ###
# Compute loop distance in Nb bins
loops <- loops %>%
    mutate(bin_dist = bin2_id - bin1_id)

### VISUALIZE ###
# Compare loop sizes
ggplot(data=loops %>% filter(library != 'PM25'), aes(x=condition, y=log10(bin_dist*res))) +
    geom_violin() +
    geom_boxplot(width=0.1)+
    theme_minimal() +
    ylab('Loop sizes [log10 bp]')

# Compare loop sets overlaps
ggplot
