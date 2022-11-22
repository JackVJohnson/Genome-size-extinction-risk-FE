
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


library(geiger)

df2 <- df
df2 <- df2[!is.na(df2$IUCNTrend_2020),]

setwd(phylo_wd)
phylo <- read.nexus("Amphibia_Dated_NEXUS_7239Spp.txt")

has_rownames(df2)
df2 <- remove_rownames(df2)
df2 <- column_to_rownames(df2, "Species")

tmp <- name.check(phylo, df2)
newphy <- drop.tip(phylo, tip=tmp$tree_not_data)
name.check(newphy, df2)

df2$IUCNTrend_2020 <- as.factor(df2$IUCNTrend_2020)

res <- df2$GenSizMed
grp <- df2$IUCNTrend_2020
names(res) = rownames(df2)
names(grp) = rownames(df2)

test_groups <- geiger::aov.phylo(res~grp, newphy, nsim = 1000)

phylosig(newphy, res, method = "lambda")

###
# bayesian phylogenetic regressions - binomial 
###

table(df$Order)
# Anura     Caudata Gymnophiona 
# 271         175          22 

df$IUCN_bin <- NA
df$IUCN_bin[df$IUCNStat_2020 == "CR"] <- 1
df$IUCN_bin[df$IUCNStat_2020 == "EN"] <- 1
df$IUCN_bin[df$IUCNStat_2020 == "VU"] <- 1
df$IUCN_bin[df$IUCNStat_2020 == "NT"] <- 0
df$IUCN_bin[df$IUCNStat_2020 == "LC"] <- 0
df$IUCN_bin[df$IUCNStat_2020 == "DD"] <- NA
table(df$IUCN_bin)

df$IUCN_rank <- NA
df$IUCN_rank[df$IUCNStat_2020 == "CR"] <- 5
df$IUCN_rank[df$IUCNStat_2020 == "EN"] <- 4
df$IUCN_rank[df$IUCNStat_2020 == "VU"] <- 3
df$IUCN_rank[df$IUCNStat_2020 == "NT"] <- 2
df$IUCN_rank[df$IUCNStat_2020 == "LC"] <- 1
df$IUCN_rank[df$IUCNStat_2020 == "DD"] <- NA

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
df3 <- df3[!is.na(df3$IUCN_bin),]

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


##########################################
# life history df
##########################################

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

str(df3$IUCN_bin)
df3$IUCN_bin <- as.factor(df3$IUCN_bin)
table(df3$IUCN_bin)
table(df3$IUCN_rank)
str(df3$IUCN_rank)

# exclude colinear variables. Use seasonality measures as they are superior predictors of extinction risk (Finn et al) 

model_bin<-parLapply(cl=cl,1,function(i){
  MCMCglmm(IUCN_bin~GenSizMed + PpSeasonality + AnnTempRange + mean.uvb, 
           random=~PhyloName, 
           family="categorical",
           ginverse = list(PhyloName=inv.phylo),
           data=df3,
           prior = prior,
           nitt=3000000, burnin=1000000)
}
) 

summary.MCMCglmm(model_bin[[1]])

mod_bin <- model_bin[[1]]

out_bin <- summary.MCMCglmm(mod_bin)
out_bin <- out_bin$solutions

plot(mod_bin$Sol)
plot(mod_bin$VCV)

print(out_bin)

# extract pagels lambda

lambda_bin <- mod_bin$VCV[,'PhyloName']/
  (mod_bin$VCV[,'PhyloName'] + mod_bin$VCV[,'units'])
mean(lambda_bin)

##########################################
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

model_bin_LH<-parLapply(cl=cl,1,function(i){
  MCMCglmm(IUCN_bin~GenSizMed + Body_Mass + ClutchSize + RangeSize..km.2. , 
           random=~PhyloName, 
           family="categorical",
           ginverse = list(PhyloName=inv.phylo),
           data=df4,
           prior = prior,
           nitt=3000000, burnin=1000000)
}
) 

summary.MCMCglmm(model_bin_LH[[1]])

mod_bin_LH <- model_bin_LH[[1]]

out_bin_LH <- summary.MCMCglmm(mod_bin_LH)
out_bin_LH <- out_bin_LH$solutions

plot(mod_bin_LH$Sol)
plot(mod_bin_LH$VCV)
options(scipen=999)
print(out_bin_LH)

# extract pagels lambda

lambda_bin_LH <- mod_bin_LH$VCV[,'PhyloName']/
  (mod_bin_LH$VCV[,'PhyloName'] + mod_bin_LH$VCV[,'units'])
mean(lambda_bin_LH)

# exprot models 


setwd(output_wd)
write.csv(out_bin, file = "Climate_binomial.csv")
write.csv(out_bin_LH, file = "LH_binomial.csv")


############################################
# models for each taxa
############################################

# anura climate 
names(df3)
dfanura_cl <- subset(df3, Order.x == "Anura")


rm(phylo)
setwd(phylo_wd)

phylo <- read.nexus("Amphibia_Dated_NEXUS_7239Spp.txt")

has_rownames(dfanura_cl)
dfanura_cl <- remove_rownames(dfanura_cl)
dfanura_cl <- column_to_rownames(dfanura_cl, "PhyloName")

phylo <- as.phylo(phylo)
tmp <- name.check(phylo, dfanura_cl)

newphy <- drop.tip(phylo, tip=tmp$tree_not_data)
newphy
name.check(newphy, dfanura_cl)   

is.ultrametric(newphy)
newphy <- force.ultrametric(newphy)

inv.phylo <- inverseA(newphy)$Ainv

dfanura_cl <- rownames_to_column(dfanura_cl, "PhyloName")

# run the model

clusterExport(cl,"dfanura_cl")
clusterExport(cl,"inv.phylo")

# binary model

names(dfanura_cl)

anura_climate<-parLapply(cl=cl,1,function(i){
  MCMCglmm(IUCN_bin~GenSizMed + PpSeasonality + AnnTempRange + mean.uvb , 
           random=~PhyloName, 
           family="categorical",
           ginverse = list(PhyloName=inv.phylo),
           data=dfanura_cl,
           prior = prior,
           nitt=3000000, burnin=1000000)
}
) 

summary.MCMCglmm(anura_climate[[1]])

mod_anura_cl <- anura_climate[[1]]

out_anura_cl <- summary.MCMCglmm(mod_anura_cl)
out_anura_cl <- out_anura_cl$solutions

plot(mod_anura_cl$Sol)
plot(mod_anura_cl$VCV)
options(scipen=999)
print(out_anura_cl)

# extract pagels lambda

lambda_anura_cl <- mod_anura_cl$VCV[,'PhyloName']/
  (mod_anura_cl$VCV[,'PhyloName'] + mod_anura_cl$VCV[,'units'])
mean(lambda_anura_cl)


################################
# anura life history 
################################

names(df4)
dfanura_lh <- subset(df4, Order.x == "Anura")


rm(phylo)
setwd(phylo_wd)

phylo <- read.nexus("Amphibia_Dated_NEXUS_7239Spp.txt")

has_rownames(dfanura_lh)
dfanura_lh <- remove_rownames(dfanura_lh)
dfanura_lh <- column_to_rownames(dfanura_lh, "Species")

phylo <- as.phylo(phylo)
tmp <- name.check(phylo, dfanura_lh)

newphy <- drop.tip(phylo, tip=tmp$tree_not_data)
newphy
name.check(newphy, dfanura_lh)   

is.ultrametric(newphy)
newphy <- force.ultrametric(newphy)

inv.phylo <- inverseA(newphy)$Ainv

dfanura_lh <- rownames_to_column(dfanura_lh, "PhyloName")

# run the model

clusterExport(cl,"dfanura_lh")
clusterExport(cl,"inv.phylo")

# binary model

names(dfanura_lh)

anura_lh<-parLapply(cl=cl,1,function(i){
  MCMCglmm(IUCN_bin~GenSizMed + Body_Mass + ClutchSize + RangeSize..km.2., 
           random=~PhyloName, 
           family="categorical",
           ginverse = list(PhyloName=inv.phylo),
           data=dfanura_lh,
           prior = prior,
           nitt=3000000, burnin=1000000)
}
) 

summary.MCMCglmm(anura_lh[[1]])

mod_anura_lh <- anura_lh[[1]]

out_anura_lh <- summary.MCMCglmm(mod_anura_lh)
out_anura_lh <- out_anura_lh$solutions

plot(mod_anura_lh$Sol)
plot(mod_anura_lh$VCV)
options(scipen=999)
print(out_anura_lh)

# extract pagels lambda

lambda_anura_lh <- mod_anura_lh$VCV[,'PhyloName']/
  (mod_anura_lh$VCV[,'PhyloName'] + mod_anura_lh$VCV[,'units'])
mean(lambda_anura_lh)


######################################
# caudata
######################################

# caudata climate 
names(df3)
dfcaudata_cl <- subset(df3, Order.x == "Caudata")


rm(phylo)
setwd(phylo_wd)

phylo <- read.nexus("Amphibia_Dated_NEXUS_7239Spp.txt")

has_rownames(dfcaudata_cl)
dfcaudata_cl <- remove_rownames(dfcaudata_cl)
dfcaudata_cl <- column_to_rownames(dfcaudata_cl, "PhyloName")

phylo <- as.phylo(phylo)
tmp <- name.check(phylo, dfcaudata_cl)

newphy <- drop.tip(phylo, tip=tmp$tree_not_data)
newphy
name.check(newphy, dfcaudata_cl)   

is.ultrametric(newphy)
newphy <- force.ultrametric(newphy)

inv.phylo <- inverseA(newphy)$Ainv

dfcaudata_cl <- rownames_to_column(dfcaudata_cl, "PhyloName")

# run the model

clusterExport(cl,"dfcaudata_cl")
clusterExport(cl,"inv.phylo")

# binary model

names(dfcaudata_cl)

caudata_climate<-parLapply(cl=cl,1,function(i){
  MCMCglmm(IUCN_bin~GenSizMed + PpSeasonality + AnnTempRange + mean.uvb , 
           random=~PhyloName, 
           family="categorical",
           ginverse = list(PhyloName=inv.phylo),
           data=dfcaudata_cl,
           prior = prior,
           nitt=3000000, burnin=1000000)
}
) 

summary.MCMCglmm(caudata_climate[[1]])

mod_caudata_cl <- caudata_climate[[1]]

out_caudata_cl <- summary.MCMCglmm(mod_caudata_cl)
out_caudata_cl <- out_caudata_cl$solutions

plot(mod_caudata_cl$Sol)
plot(mod_caudata_cl$VCV)
options(scipen=999)
print(out_caudata_cl)

# extract pagels lambda

lambda_caudata_cl <- mod_caudata_cl$VCV[,'PhyloName']/
  (mod_caudata_cl$VCV[,'PhyloName'] + mod_caudata_cl$VCV[,'units'])
mean(lambda_caudata_cl)



################################
# caudata life history 
################################

names(df4)
dfcaudata_lh <- subset(df4, Order.x == "Caudata")


rm(phylo)
setwd(phylo_wd)

phylo <- read.nexus("Amphibia_Dated_NEXUS_7239Spp.txt")

has_rownames(dfcaudata_lh)
dfcaudata_lh <- remove_rownames(dfcaudata_lh)
dfcaudata_lh <- column_to_rownames(dfcaudata_lh, "Species")

phylo <- as.phylo(phylo)
tmp <- name.check(phylo, dfcaudata_lh)

newphy <- drop.tip(phylo, tip=tmp$tree_not_data)
newphy
name.check(newphy, dfcaudata_lh)   

is.ultrametric(newphy)
newphy <- force.ultrametric(newphy)

inv.phylo <- inverseA(newphy)$Ainv

dfcaudata_lh <- rownames_to_column(dfcaudata_lh, "PhyloName")

# run the model

clusterExport(cl,"dfcaudata_lh")
clusterExport(cl,"inv.phylo")

# binary model

names(dfcaudata_lh)

caudata_lh<-parLapply(cl=cl,1,function(i){
  MCMCglmm(IUCN_bin~GenSizMed + Body_Mass + ClutchSize + RangeSize..km.2., 
           random=~PhyloName, 
           family="categorical",
           ginverse = list(PhyloName=inv.phylo),
           data=dfcaudata_lh,
           prior = prior,
           nitt=3000000, burnin=1000000)
}
) 

summary.MCMCglmm(caudata_lh[[1]])

mod_caudata_lh <- caudata_lh[[1]]

out_caudata_lh <- summary.MCMCglmm(mod_caudata_lh)
out_caudata_lh <- out_caudata_lh$solutions

plot(mod_caudata_lh$Sol)
plot(mod_caudata_lh$VCV)
options(scipen=999)
print(out_caudata_lh)

# extract pagels lambda

lambda_caudata_lh <- mod_caudata_lh$VCV[,'PhyloName']/
  (mod_caudata_lh$VCV[,'PhyloName'] + mod_caudata_lh$VCV[,'units'])
mean(lambda_caudata_lh)


########################################################
# Gymnophoina
########################################################


# gymnophoina climate 
names(df3)
table(df3$Order.x)
dfgymno_cl <- subset(df3, Order.x == "Gymnophiona")


rm(phylo)
setwd(phylo_wd)

phylo <- read.nexus("Amphibia_Dated_NEXUS_7239Spp.txt")

has_rownames(dfgymno_cl)
dfgymno_cl <- remove_rownames(dfgymno_cl)
dfgymno_cl <- column_to_rownames(dfgymno_cl, "PhyloName")

phylo <- as.phylo(phylo)
tmp <- name.check(phylo, dfgymno_cl)

newphy <- drop.tip(phylo, tip=tmp$tree_not_data)
newphy
name.check(newphy, dfgymno_cl)   

is.ultrametric(newphy)
newphy <- force.ultrametric(newphy)

inv.phylo <- inverseA(newphy)$Ainv

dfgymno_cl <- rownames_to_column(dfgymno_cl, "PhyloName")

# run the model

clusterExport(cl,"dfgymno_cl")
clusterExport(cl,"inv.phylo")

# binary model

names(dfgymno_cl)

gymno_climate<-parLapply(cl=cl,1,function(i){
  MCMCglmm(IUCN_bin~GenSizMed + PpSeasonality + AnnTempRange + mean.uvb , 
           random=~PhyloName, 
           family="categorical",
           ginverse = list(PhyloName=inv.phylo),
           data=dfgymno_cl,
           prior = prior,
           nitt=3000000, burnin=1000000)
}
) 

summary.MCMCglmm(gymno_climate[[1]])

mod_gymno_cl <- gymno_climate[[1]]

out_gymno_cl <- summary.MCMCglmm(mod_gymno_cl)
out_gymno_cl <- out_gymno_cl$solutions

plot(mod_gymno_cl$Sol)
plot(mod_gymno_cl$VCV)
options(scipen=999)
print(out_gymno_cl)

# extract pagels lambda

lambda_gymno_cl <- mod_gymno_cl$VCV[,'PhyloName']/
  (mod_gymno_cl$VCV[,'PhyloName'] + mod_gymno_cl$VCV[,'units'])
mean(lambda_gymno_cl)



################################
# gymnophoina life history 
################################

# this model does not converge so we cannot use in the results. 

names(df4)
dfgymno_lh <- subset(df4, Order.x == "Gymnophiona")


rm(phylo)
setwd(phylo_wd)

phylo <- read.nexus("Amphibia_Dated_NEXUS_7239Spp.txt")

has_rownames(dfgymno_lh)
dfgymno_lh <- remove_rownames(dfgymno_lh)
dfgymno_lh <- column_to_rownames(dfgymno_lh, "Species")

phylo <- as.phylo(phylo)
tmp <- name.check(phylo, dfgymno_lh)

newphy <- drop.tip(phylo, tip=tmp$tree_not_data)
newphy
name.check(newphy, dfgymno_lh)   

is.ultrametric(newphy)
newphy <- force.ultrametric(newphy)

inv.phylo <- inverseA(newphy)$Ainv

dfgymno_lh <- rownames_to_column(dfgymno_lh, "PhyloName")

# run the model

clusterExport(cl,"dfgymno_lh")
clusterExport(cl,"inv.phylo")

# binary model

names(dfgymno_lh)

gymno_lh<-parLapply(cl=cl,1,function(i){
  MCMCglmm(IUCN_bin~GenSizMed + Body_Mass + ClutchSize + RangeSize..km.2., 
           random=~PhyloName, 
           family="categorical",
           ginverse = list(PhyloName=inv.phylo),
           data=dfgymno_lh,
           prior = prior,
           nitt=16000000, burnin=8000000)
}
) 

summary.MCMCglmm(gymno_lh[[1]])

mod_gymno_lh <- gymno_lh[[1]]

out_gymno_lh <- summary.MCMCglmm(mod_gymno_lh)
out_gymno_lh <- out_gymno_lh$solutions

plot(mod_gymno_lh$Sol)
plot(mod_gymno_lh$VCV)
options(scipen=999)
print(out_gymno_lh)

# extract pagels lambda

lambda_gymno_lh <- mod_gymno_lh$VCV[,'PhyloName']/
  (mod_gymno_lh$VCV[,'PhyloName'] + mod_gymno_lh$VCV[,'units'])
mean(lambda_gymno_lh)
