# coeff plot for extinction risk under genome size

library(tidyverse)
library(dotwhisker)
library(patchwork)
library(fishualize)

output_wd <- "C:/Users/40274182/OneDrive - Queen's University Belfast/PhD/Projects/Genome_size_extinction_DPD/Output_dir"

setwd(output_wd)

bin_clim <- read.csv("Climate_binomial.csv")
bin_life <- read.csv("LH_binomial.csv")

# climate plot first

bin_clim <- bin_clim[-1,]

bin_clim$model <- "Binomial"

clim_coeffs <-bin_clim
colnames(clim_coeffs)[1] <- "Predictor"

# clim history plot


bin_life <- bin_life[-1,]

bin_life$model <- "Binomial"

life_coeffs <- bin_life
colnames(life_coeffs)[1] <- "Predictor"


# binomial plots
str(life_coeffs$Predictor)
life_coeffs$Predictor <- as.factor(life_coeffs$Predictor)

life_plot_bin <- ggplot(life_coeffs, aes(post.mean, Predictor)) +
  geom_vline(aes(xintercept=0), size = .25, linetype = "dashed") +
  geom_errorbarh(data = life_coeffs, aes(xmax = u.95..CI, xmin = l.95..CI), size = .6, height = .2, color = "grey12", position = position_nudge(y = 0)) +
  geom_point(data = life_coeffs, size = 6, position = position_nudge(y = 0)) +
  theme_classic() +
  scale_y_discrete(limits=c("ClutchSize", "Body_Mass", "RangeSize..km.2.","GenSizMed"), labels=c("Clutch size", "Body size", "Range size", "Genome size")) +
  theme(axis.title.x = element_text(size = 22, color = "black"), axis.title.y = element_text(size = 22, color = "black"), text=element_text(size=16)) +
  theme(axis.text.y = element_text(size=18, color="black"), axis.text.x = element_text(size=18, color="black")) +
  labs(y= "", x="") +
  ggtitle("B") 

life_plot_bin

# climate

clim_coeffs <- subset(clim_coeffs, model == "Binomial")

clim_plot_bin <- ggplot(clim_coeffs, aes(post.mean, Predictor)) +
  geom_vline(aes(xintercept=0), size = .25, linetype = "dashed") +
  geom_errorbarh(data = clim_coeffs, aes(xmax = u.95..CI, xmin = l.95..CI), size = .6, height = .2, color = "grey12", position = position_nudge(y = 0)) +
  geom_point(data = clim_coeffs, size = 6, position = position_nudge(y = 0)) +
  theme_classic() +
  scale_y_discrete(limits=c("mean.uvb", "PpSeasonality", "AnnTempRange", "GenSizMed"), labels=c("UVB","Pp seasonality", "Temp range", "Genome size")) +
  theme(axis.title.x = element_text(size = 22, color = "black"), axis.title.y = element_text(size = 22, color = "black"), text=element_text(size=16)) +
  theme(axis.text.y = element_text(size=18, color="black"), axis.text.x = element_text(size=18, color="black")) +
  labs(y= "", x="Model coefficient")+
  ggtitle("A") 

clim_plot_bin

##########################################################################################################
############################################ plot by order ##############################################


coeffs <- read.csv("Orders_coeffs_plot.csv")

coeffs$Predictor <- trimws(coeffs$Predictor)

names(coeffs)
table(coeffs$Predictor)

clim_p <- ggplot(subset(coeffs, Model == "Climate"), aes(Mean, Predictor, color = Order)) +
  geom_vline(aes(xintercept=0), size = .25, linetype = "dashed") +
  geom_errorbarh(aes(xmax = Upper, xmin = Lower), size = 1, height = .2, position = position_dodge(width = .5)) +
  geom_point(size = 6, position = position_dodge(width = .5)) +
  theme_classic() +
  scale_color_fish_d("Centropyge_loricula")  +
  scale_y_discrete(limits=c("Mean uvb", "Precipitation seasonality", "Annual temperature range", "Genome size"), labels=c("UVB","Pp seasonality", "Temp range", "Genome size")) +
  theme(axis.title.x = element_text(size = 22, color = "black"), axis.title.y = element_text(size = 22, color = "black"), text=element_text(size=16)) +
  theme(axis.text.y = element_text(size=18, color="black"), axis.text.x = element_text(size=18, color="black")) +
  labs(y= "", x="Model coefficient") +
  ggtitle("C") 
clim_p

life_p <- ggplot(subset(coeffs, Model == "Life_history"), aes(Mean, Predictor, color = Order)) +
  geom_vline(aes(xintercept=0), size = .25, linetype = "dashed") +
  geom_errorbarh(aes(xmax = Upper, xmin = Lower), size = 1, height = .2, position = position_dodge(width = .5)) +
  geom_point(size = 6, position = position_dodge(width = .5)) +
  theme_classic() +
  scale_color_fish_d("Centropyge_loricula")  +
  scale_y_discrete(limits=c("Clutch size", "Body mass", "Range size","Genome Size"), labels=c("Clutch size", "Body size", "Range size", "Genome size")) +
  theme(axis.title.x = element_text(size = 22, color = "black"), axis.title.y = element_text(size = 22, color = "black"), text=element_text(size=16)) +
  theme(axis.text.y = element_text(size=18, color="black"), axis.text.x = element_text(size=18, color="black")) +
  labs(y= "", x="Model coefficient") +
  ggtitle("D") 
life_p


###############################################################################################################
########################################### export figure #####################################################

png(file=file.path(output_wd, 'Figure_2.png'),height=6000,width=9000,res=600)
(clim_plot_bin + life_plot_bin) / (clim_p + life_p + plot_layout(guides="collect"))+ plot_annotation(title="IUCN Conservation Status",  theme = theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold")))
dev.off()
