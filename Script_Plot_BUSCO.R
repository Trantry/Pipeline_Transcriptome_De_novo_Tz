######################################
#
# BUSCO summary figure
# @version 4.0.0
# @since BUSCO 2.0.0
# 
# Copyright (c) 2016-2020, Evgeny Zdobnov (ez@ezlab.org)
# Licensed under the MIT license. See LICENSE.md file.
#
######################################

# Load the required libraries
library("ggplot2")
library("grid")

# !!! CONFIGURE YOUR PLOT HERE !!! 
# Output
my_output <- paste("busco_comparison_tetraripis_zetteli_corrected.png",sep="/") 
my_width <- 28
my_height <- 15
my_unit <- "cm"

# Colors
my_colors <- c("#56B4E9", "#3492C7", "#F0E442", "#F04442")
# Bar height ratio
my_bar_height <- 0.75

# Legend
my_title <- "BUSCO Assessment Results Comparison for Tetraripis zetteli"

# Font
my_family <- "sans"
my_size_ratio <- 1

# !!! SEE YOUR DATA HERE !!! 
# Data for old and new transcriptomes
my_species <- c( rep('Tetraripis zetteli (New)', 4),rep('Tetraripis zetteli (Old)', 4))
my_percentage <- c(37.1,38.6,7.1,17.2,34.9,30.5,6.1,28.5)
my_values <- c(930,970,177,433,876,766,152,716)
categories <- c('S', 'D', 'F', 'M', 'S', 'D', 'F', 'M')
my_species <- factor(my_species)
categories <- factor(categories, levels = c( 'S', 'D', 'F','M')) # set order of categories
######################################
######################################
######################################
# Code to produce the graph
labsize = 1
if (length(levels(my_species)) > 10){
  labsize = 0.66
}
print("Plotting the figure ...")
df = data.frame(my_species, my_percentage, my_values, categories)

figure <- ggplot() + 
  geom_bar(aes(y = my_percentage, x = my_species, fill = categories), position = position_stack(reverse = TRUE), data = df, stat="identity", width=my_bar_height) + 
  coord_flip() + 
  theme_gray(base_size = 8) + 
  scale_y_continuous(labels = c("0","20","40","60","80","100"), breaks = c(0,20,40,60,80,100)) + 
  scale_fill_manual(values = my_colors, labels =c("Complete (C) and single-copy (S)", 
                                                  "Complete (C) and duplicated (D)", 
                                                  "Fragmented (F)", 
                                                  "Missing (M)")) +   
  ggtitle(my_title) + 
  xlab("") + 
  ylab("\n%BUSCOs") + 
  theme(plot.title = element_text(family=my_family, hjust=0.5, colour = "black", size = rel(2.2)*my_size_ratio, face = "bold")) + 
  theme(legend.position="top",legend.title = element_blank()) + 
  theme(legend.text = element_text(family=my_family, size = rel(1.2)*my_size_ratio)) + 
  theme(panel.background = element_rect(color="#FFFFFF", fill="white")) + 
  theme(panel.grid.minor = element_blank()) + 
  theme(panel.grid.major = element_blank()) +
  theme(axis.text.y = element_text(family=my_family, colour = "black", size = rel(1.66)*my_size_ratio)) + 
  theme(axis.text.x = element_text(family=my_family, colour = "black", size = rel(1.66)*my_size_ratio)) + 
  theme(axis.line = element_line(size=1*my_size_ratio, colour = "black")) + 
  theme(axis.ticks.length = unit(.85, "cm")) + 
  theme(axis.ticks.y = element_line(colour="white", size = 0)) + 
  theme(axis.ticks.x = element_line(colour="#222222")) + 
  theme(axis.ticks.length = unit(0.4, "cm")) + 
  theme(axis.title.x = element_text(family=my_family, size=rel(1.2)*my_size_ratio)) + 
  guides(fill = guide_legend(override.aes = list(colour = NULL))) +
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

total_buscos_old <- sum(my_values[1:4])
total_buscos_new <- sum(my_values[5:8])

figure <- figure + 
  annotate("text", label=paste("C:", my_values[1] + my_values[2], " [S:", my_values[1], ", D:", my_values[2], "], F:", my_values[3], ", M:", my_values[4], ", n:", total_buscos_old, sep=""), 
           y=3, x = 1, size = labsize*4*my_size_ratio, colour = "black", hjust=0, family=my_family) +
  annotate("text", label=paste("C:", my_values[5] + my_values[6], " [S:", my_values[5], ", D:", my_values[6], "], F:", my_values[7], ", M:", my_values[8], ", n:", total_buscos_new, sep=""), 
           y=3, x = 2, size = labsize*4*my_size_ratio, colour = "black", hjust=0, family=my_family)

ggsave(figure, file=my_output, width = my_width, height = my_height, unit = my_unit)
print("Done")

figure


