# coeff plot for extinction risk under genome size for population trends

library(tidyverse)
library(dotwhisker)
library(patchwork)
library(fishualize)

output_wd <- "C:/Users/40274182/OneDrive - Queen's University Belfast/PhD/Projects/Genome_size_extinction_DPD/Output_dir"

setwd(output_wd)

coeffs <- read.csv("pop_trend_coeffs.csv")
options(scipen=999)

p1 <- ggplot(subset(coeffs, Model=="Climate"), aes(Mean, Predictor)) +
  geom_vline(aes(xintercept=0), size = .25, linetype = "dashed") +
  geom_errorbarh(aes(xmax = Upper, xmin = Lower), size = 1, height = .2, position = position_dodge(width = .5)) +
  geom_point(size = 6, position = position_dodge(width = .5)) +
  theme_classic() +
  scale_color_fish_d("Centropyge_loricula")  +
  scale_y_discrete(limits=c("Mean uvb ", "Precipitation seasonality ", "Annual temperature range", "Genome size"), labels=c("UVB","Pp seasonality", "Temp range", "Genome size")) +
  theme(axis.title.x = element_text(size = 22, color = "black"), axis.title.y = element_text(size = 22, color = "black"), text=element_text(size=16)) +
  theme(axis.text.y = element_text(size=18, color="black"), axis.text.x = element_text(size=18, color="black")) +
  labs(y= "", x="Model coefficient") +
  ggtitle("A") +
  scale_x_continuous(labels = function(x) format(x, scientific = TRUE))
p1


p2 <- ggplot(subset(coeffs, Model=="Life_history"), aes(Mean, Predictor)) +
  geom_vline(aes(xintercept=0), size = .25, linetype = "dashed") +
  geom_errorbarh(aes(xmax = Upper, xmin = Lower), size = 1, height = .2, position = position_dodge(width = .5)) +
  geom_point(size = 6, position = position_dodge(width = .5)) +
  theme_classic() +
  scale_color_fish_d("Centropyge_loricula")  +
  scale_y_discrete(limits=c("Clutch size", "Body mass", "Range size ","Genome Size"), labels=c("Clutch size", "Body size", "Range size", "Genome size")) +
  theme(axis.title.x = element_text(size = 22, color = "black"), axis.title.y = element_text(size = 22, color = "black"), text=element_text(size=16)) +
  theme(axis.text.y = element_text(size=18, color="black"), axis.text.x = element_text(size=18, color="black")) +
  labs(y= "", x="Model coefficient") +
  ggtitle("B") +
  scale_x_continuous(labels = function(x) format(x, scientific = TRUE))
p2


###############################################################################################################
########################################### export figure #####################################################

png(file=file.path(output_wd, 'Figure_3_pop_trend.png'),height=4000,width=9000,res=600)
(p1 + p2) + plot_annotation(title="IUCN Population rend",  theme = theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold")))
dev.off()

