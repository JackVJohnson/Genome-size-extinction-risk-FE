
# libraries

library(tidyverse)
library(dplyr)
library(ape)
library(MCMCglmm)
library(geiger)
library(phytools)
library(tibble)

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

# tiday data 

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

table(df$IUCNTrend_2020)

df$pop_trend[df$IUCNTrend_2020 == "Decreasing"] <- 1
df$pop_trend[df$IUCNTrend_2020 == "Stable"] <- 2
df$pop_trend[df$IIUCNTrend_2020 == "Increasing"] <- 3


df3 <- df
names(df3)

# genome size with climate predictors for extinction risk 

df3 <- df3[!is.na(df3$GenSizMed),]
df3 <- df3[!is.na(df3$Latitude),]
df3 <- df3[!is.na(df3$AnnMnTemp),]
df3 <- df3[!is.na(df3$AnnTempRange),]
df3 <- df3[!is.na(df3$AnnPrecipitation),]
df3 <- df3[!is.na(df3$PpSeasonality),]
df3 <- df3[!is.na(df3$NPP),]
df3 <- df3[!is.na(df3$pop_trend),]

str(df3$Latitude)

df3$GenSizMed <- log1p(df3$GenSizMed)
df3$Latitude <- abs(df3$Latitude)
df3$Latitude <- myStd(df3$Latitude)
df3$AnnMnTemp <- myStd(df3$AnnMnTemp)
df3$AnnTempRange <- myStd(df3$AnnTempRange)
df3$AnnPrecipitation <- myStd(df3$AnnPrecipitation)
df3$PpSeasonality<- myStd(df3$PpSeasonality)
df3$NPP <- myStd(df3$NPP)
df3$mean.uvb <- myStd(df3$mean.uvb)

##########################################################################################################
########################################### climate model #################################################


df4 <- df

names(df4)

sum(is.na(df4$Longevity))
sum(is.na(df4$ClutchSize))
sum(is.na(df4$AgeatMatur))
sum(is.na(df4$EggHatching))
sum(is.na(df4$RangeSize..km.2.))

# for best resolution data and to avoid coliniarity use
# range size, clutch size, body mass, genome size

df4 <- df4[!is.na(df4$GenSizMed),]
df4 <- df4[!is.na(df4$ClutchSize),]
df4 <- df4[!is.na(df4$Body_Mass),]
df4 <- df4[!is.na(df4$RangeSize..km.2.),]

df4$GenSizMed <- log1p(df4$GenSizMed)
df4$ClutchSize <- log1p(df4$ClutchSize)
df4$Body_Mass <- log1p(df4$Body_Mass)
df4$RangeSize..km.2. <- log1p(df4$RangeSize..km.2.)

# 364 species

###################################################################
# climate regression 440 species 

library(parallel)

setwd(phylo_wd)

phylo <- read.nexus("Amphibia_Dated_NEXUS_7239Spp.txt")

has_rownames(df3)
df3 <- remove_rownames(df3)
df3 <- column_to_rownames(df3, "Species")

phylo <- as.phylo(phylo)
tmp <- name.check(phylo, df3)

newphy <- drop.tip(phylo, tip=tmp$tree_not_data)
newphy
name.check(newphy, df3)  


is.ultrametric(newphy)
newphy <- force.ultrametric(newphy)

inv.phylo <- inverseA(newphy)$Ainv

df3 <- rownames_to_column(df3, "PhyloName")


# prior is inverse-Gamma distribution with shape and scale parameters equal to 0.01, which is relatively canonical

prior<-list(G=list(G1=list(V=1,nu=0.02)),R=list(V=1,nu=0.02))

# run in parallel 

library(parallel)
parallel::detectCores()

setCores<-round(detectCores()*0.95) 
cl <- makeCluster(getOption("cl.cores",setCores))
cl.pkg <- clusterEvalQ(cl,library(MCMCglmm)) 
clusterExport(cl,"prior")
clusterExport(cl,"df3")
clusterExport(cl,"inv.phylo")


# exclude colinear variables. Use seasonality measures as they are superior predictors of extinction risk (Finn et al) 

model_clim<-parLapply(cl=cl,1,function(i){
  MCMCglmm(pop_trend~GenSizMed + PpSeasonality + AnnTempRange + mean.uvb, 
           random=~PhyloName, 
           family="ordinall",
           ginverse = list(PhyloName=inv.phylo),
           data=df3,
           prior = prior,
           nitt=500000, burnin=250000)
}
) 

summary.MCMCglmm(model_clim[[1]])

mod_1 <- model_clim[[1]]

out_1 <- summary.MCMCglmm(mod_1)
out_1 <- out_1$solutions

plot(mod_1$Sol)
plot(mod_1$VCV)

print(out_1)

# extract pagels lambda

lambda_1 <- mod_1$VCV[,'PhyloName']/
  (mod_1$VCV[,'PhyloName'] + mod_1$VCV[,'units'])
mean(lambda_1)

################################################################################################
##################################### life history #############################################

# life history model 
# data frame already tidy 

# work with phylo

rm(phylo)
setwd(phylo_wd)

phylo <- read.nexus("Amphibia_Dated_NEXUS_7239Spp.txt")

has_rownames(df4)
df4 <- remove_rownames(df4)
df4 <- column_to_rownames(df4, "Species")

phylo <- as.phylo(phylo)
tmp <- name.check(phylo, df4)

newphy <- drop.tip(phylo, tip=tmp$tree_not_data)
newphy
name.check(newphy, df4)   

is.ultrametric(newphy)
newphy <- force.ultrametric(newphy)

inv.phylo <- inverseA(newphy)$Ainv

df4 <- rownames_to_column(df4, "PhyloName")

# run the model

clusterExport(cl,"df4")
clusterExport(cl,"inv.phylo")

# binary model

names(df4)

model_LH<-parLapply(cl=cl,1,function(i){
  MCMCglmm(pop_trend~GenSizMed + Body_Mass + ClutchSize + RangeSize..km.2. , 
           random=~PhyloName, 
           family="ordinal",
           ginverse = list(PhyloName=inv.phylo),
           data=df4,
           prior = prior,
           nitt=500000, burnin=250000)
}
) 

summary.MCMCglmm(model_LH[[1]])

mod_LH <- model_LH[[1]]

out_LH <- summary.MCMCglmm(mod_LH)
out__LH <- out_LH$solutions

plot(mod_LH$Sol)
plot(mod_LH$VCV)
options(scipen=999)
print(out_LH)

# extract pagels lambda

lambda_LH <- mod_LH$VCV[,'PhyloName']/
  (modLH$VCV[,'PhyloName'] + modLH$VCV[,'units'])
mean(lambdaLH)
