# coeff plot for extinction risk under genome size

library(tidyverse)
library(dotwhisker)

output_wd <- "C:/Users/40274182/OneDrive - Queen's University Belfast/PhD/Projects/Genome_size_extinction_DPD/Output_dir"

setwd(output_wd)


genome_dir <- "C:/Users/40274182/OneDrive - Queen's University Belfast/PhD/Projects/Genome_size_extinction_DPD/Tree_figure"

setwd(genome_dir)


data <- read.csv("GenomeSize_JackNEWdataset.csv", row.names = 1, header = T)

data <- data[!is.na(data$IUCN_w),]
data <- data[!is.na(data$PopTrend_w),]
data <- data[!is.na(data$Bodysize),]
data <- data[!is.na(data$LNFecundity),]

hist(data$Genome)

names(data)

tree<-read.nexus(file.choose())

tmp <- name.check(tree, data)
tmp

str(tmp$tree_not_data)

newphy <- drop.tip(tree, tip=tmp$tree_not_data)
newphy
name.check(newphy, data)        

rm(tree)

data1 <- read.csv("GenomeSize_JackNEWdataset.csv")
data1 <- data1[!is.na(data1$IUCN_w),]
data1 <- data1[!is.na(data1$PopTrend_w),]
data1 <- data1[!is.na(data1$Bodysize),]
data1 <- data1[!is.na(data1$LNFecundity),]


p <- ggtree(newphy, layout = "fan", size = 0.1)               
p
p2 <- p +
  geom_fruit(
    data = data1,
    geom=geom_star,
    mapping=aes(y=X, fill=PopTrend_w),
    position = "identity",
    starshape=15,
    starstroke=0.2,
    size=1.5
  ) +
  scale_fill_manual(
    values=c("#CC0033", "#33FF00", "#3399FF")
  ) + 
  guides(fill=guide_legend(title="Population trend"))

p2

p3 <- p2 +
  new_scale_fill() +
  geom_fruit(
    data=data1,
    geom=geom_tile,
    mapping=aes(y=X, x=IUCN_n, fill=IUCN_w),
    offset=0.08,
    pwdith=0.25
  ) +
  scale_fill_discrete(breaks=c("CR", "EN", "VU", "NT", "LC")) +
  guides(fill=guide_legend(title="Red list category"))
p3

p4 <- p3 +
  geom_fruit(
    data=data1,
    geom=geom_bar,
    mapping=aes(y=X, x=Genome),
    pwdith=0.4,
    stat="identity",
    orientation="y",
    fill="gray80"
  ) 
p4

hist(data1$Bodysize)
hist(data1$LNBodysize)

p5 <- p4 +
  geom_fruit(
    data=data1,
    geom=geom_bar,
    mapping=aes(y=X, x=Bodysize),
    pwdith=0.4,
    stat="identity",
    orientation="y",
    fill="gray50"
  ) 

p5

hist(data1$Fecundity)
hist(data1$LNFecundity)

p6 <- p5 +
  geom_fruit(
    data=data1,
    geom=geom_bar,
    mapping=aes(y=X, x=LNFecundity),
    pwdith=0.4,
    stat="identity",
    orientation="y",
    fill="gray20"
  ) 

p6
