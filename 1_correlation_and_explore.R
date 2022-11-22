###
# genome size and extinction risk in the worlda amphibians 
###

# libraries

library(dplyr)
library(ape)
library(MCMCglmm)
library(geiger)
library(phytools)
library(tibble)
library(tidyverse)

# working directories 

data_wd <- "C:/Users/40274182/OneDrive - Queen's University Belfast/PhD/Projects/Genome_size_extinction_DPD/Data_files"
output_wd <- "C:/Users/40274182/OneDrive - Queen's University Belfast/PhD/Projects/Genome_size_extinction_DPD/Output_dir"
phylo_wd <- "C:/Users/40274182/OneDrive - Queen's University Belfast/PhD/Phylogonies/Amphibian"

# load in data

setwd(data_wd)
df <- read.csv("DATASET_Funct Ecol_Jack_February.2022.csv")
uvb_df <- read.csv("UVB_data.csv")

uvb_df <- rename(uvb_df, Species = phylo.name)

df <- left_join(df, uvb_df, by="Species")

# standardisation function

myStd <- function(x) { (x-mean(x))/sd(x)}

# tiday data and produce corrlation matrix for life history and climate variables 

names(df)
summary(df$Body_Mass)
df$Body_Mass <- as.numeric(df$Body_Mass)
summary(df$Longevity)
summary(df$ClutchSize)
summary(df$AgeatMatur)
summary(df$EggHatching)
summary(df$Egg.Diameter)
summary(df$GenSizMed)
summary(df$RangeSize..km.2.)

lh_df <- select(df, c(Body_Mass, Longevity, ClutchSize, AgeatMatur, EggHatching, Egg.Diameter, GenSizMed, RangeSize..km.2.))


summary(df$Latitude)
summary(df$AnnMnTemp)
summary(df$AnnTempRange)
summary(df$AnnPrecipitation)
summary(df$PpSeasonality)
summary(df$NPP)

cl_df <- select(df, c(Latitude, AnnMnTemp, AnnTempRange, AnnPrecipitation, PpSeasonality, NPP, GenSizMed, mean.uvb))

library(corrmorant)

lh_corr <- ggcorrm(lh_df,
                   aes(col=.corr, fill=.corr),
                   bg_dia = "grey20",
                   rescale = "by_sd", 
                   corr_method = "spearman") +
  lotri_heatcircle(alpha = 1, col = 1) +
  utri_corrtext(nrow = 2, squeeze = 0.4) +
  dia_names(y_pos = 0.15, size = 3, color = "white") +
  dia_histogram(lower = 0.3, color = "grey80", fill = "grey60", size = .3) +
  scale_color_corr(aesthetics = c("fill", "color")) +
  theme(text = element_text(size = 14)) +
  theme(axis.text.y = element_text(size=0, color="white")) +
  theme(axis.text.x = element_text(size=0, color="white"))

lh_corr

cl_corr <- ggcorrm(cl_df,
                   aes(col=.corr, fill=.corr),
                   bg_dia = "grey20",
                   rescale = "by_sd", 
                   corr_method = "spearman") +
  lotri_heatcircle(alpha = 1, col = 1) +
  utri_corrtext(nrow = 2, squeeze = 0.4) +
  dia_names(y_pos = 0.15, size = 3, color = "white") +
  dia_histogram(lower = 0.3, color = "grey80", fill = "grey60", size = .3) +
  scale_color_corr(aesthetics = c("fill", "color")) +
  theme(text = element_text(size = 14)) +
  theme(axis.text.y = element_text(size=0, color="white")) +
  theme(axis.text.x = element_text(size=0, color="white"))
cl_corr

# save correlation matrix 

library(patchwork)

setwd(output_wd)
png(file = "correlation_matrix_genome_LH.png", height = 4000, width = 4000, res = 350)
lh_corr
dev.off()

png(file = "correlation_matrix_genome_CL.png", height = 4000, width = 4000, res = 350)
cl_corr
dev.off()

###
# density plots of genome size for each order
###

rm(df)
setwd(data_wd)
df <- read.csv("DATASET_Funct Ecol_Jack_February.2022.csv")

table(df$Order)
p_dens <- ggplot(df, aes(GenSizMed, group = Order, fill = Order)) +
  geom_density(adjust=1.5) +
  scale_fill_viridis_d(option = "D", alpha = 0.6) +
  theme_classic() +
  ylim(0, 0.25) +
  theme(axis.title.x = element_text(size = 22, color = "black"), axis.title.y = element_text(size = 22, color = "black"), text=element_text(size=16)) +
  theme(axis.text.y = element_text(size=18, color="black"), axis.text.x = element_text(size=18, color="black")) +
  labs(y= "Probabiltiy density", x="Genome size") 
p_dens

png(file = file.path(output_wd, "genome_size_density.png"), height = 2000, width = 3000, res = 350)
p_dens
dev.off()

###
# Plyogenetic Anovas for population trends 
###

table(df$IUCNTrend_2020)

#Decreasing Increasing   Stable 
# 173         11        216 

summary(df$IUCNTrend_2020)


p_trend <- ggplot(na.omit(df), aes(IUCNTrend_2020, GenSizMed, fill = IUCNTrend_2020)) +
  geom_violin(trim = F) +
  geom_boxplot(width = 0.1, fill = "white") +
  scale_fill_viridis_d(option = "D", alpha = 0.6) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 22, color = "black"), axis.title.y = element_text(size = 22, color = "black"), text=element_text(size=16)) +
  theme(axis.text.y = element_text(size=18, color="black"), axis.text.x = element_text(size=18, color="black")) +
  labs(y= "Genome size", x="Population trend")  +
  theme(legend.position="none")
p_trend


png(file = file.path(output_wd, "genome_size_pop_trend.png"), height = 2000, width = 3000, res = 350)
p_trend
dev.off()

