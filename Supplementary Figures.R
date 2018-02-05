#### Experimental design

# The following code was used in the analysis of data generated as part of a DOE funded project to explore the relationship between drought and microbial recruitment in Sorghum bicolor. 
# For detailed descriptions of experimental design, please see the associated publication. In brief, we planted two different sorghum cultivars (RTx430 and BTx642) within a randomized block design that accounted for treatments, genotypes and replication, with three replicate blocks in total. 
# From this field experiment, we collected a variety of plant phenotypes, soil measurements, and rhizosphere, root and soil samples for microbial community analysis and metatranscriptomics.
# All samples were collected weekly at the same time of day (between 10am and 1pm) and the same day of the week for seventeen weeks following seedling emergence (TP1 to TP17). 

##### Load all of the packages required for quality control, main figures, and supplementary figures.

library("phyloseq")
library("ggplot2")
library("scales")
library("grid")
theme_set(theme_bw())
library("DESeq2")
library("ape")
library("vegan")
library("data.table")
library("RColorBrewer")
library(colorRamps)
library("svglite")
library("plyr")
library("reshape2")
library(VennDiagram)

##### Fig.S4 Soil water content depletion measurements 

data<-read.table("figS4-soil-water-depletion.txt",header=T,sep="\t")
data <- data.frame(data)
data$Treatment=factor(data$Treatment,levels=c("Control","Pre-flowering","Post-flowering"))
data$Weeks<-factor(data$Weeks, levels=c("TP1","TP2","TP3","TP4","TP5","TP6","TP7","TP8","TP9","TP10","TP11","TP12","TP13","TP14","TP15","TP16","TP17"))
kearney<-subset(data,Field=="Kearney")
gill<-subset(data,Field=="Gill")

ggplot(kearney, aes(x=Weeks, y=mean, colour=Treatment, group=Treatment)) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.6) + xlab("Time point (in weeks)")+
  geom_point(shape=21,fill="white",size=3,stroke =1.5)+
  ylab("% Depletion Soil Moisture") + scale_color_manual(name="",values = c("#6AB187", "#DE7A22", "#F4CC70"))+
  scale_y_continuous(limits = c(0,100))+
   theme(axis.text.x=element_text(hjust=0.5,vjust=0.5,size=11,color="black",face="bold",angle=90), 
        axis.text.y=element_text(size=11,color="black",face="bold"), 
        axis.title=element_text(size=11,face="bold"),text=element_text(size=11))+
  theme(legend.text = element_text(colour="black", size = 11, face = "bold"))

ggplot(gill, aes(x=Weeks, y=mean, colour=Treatment, group=Treatment)) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.6) + xlab("")+
  geom_point(shape=21,fill="white",size=3,stroke =1.5)+
  ylab("% Depletion Soil Moisture") + scale_color_manual(name="",values = c("#6AB187", "#DE7A22"))+
  scale_y_continuous(limits = c(0,100))+
  theme(axis.text.x=element_text(hjust=0.5,vjust=0.5,size=11,color="black",face="bold",angle=90), 
        axis.text.y=element_text(size=11,color="black",face="bold"), 
        axis.title=element_text(size=11,face="bold"),text=element_text(size=11))+
  theme(legend.text = element_text(colour="black", size = 11, face = "bold"))

##### Fig.S5 Crop Water Stress Index analysis by treatment and timepoint

theme_set(theme_bw())
data<-read.table("figS5-water-stress-index.txt",header=T,sep="\t")
data <- data.frame(data)
data$Treatment=factor(data$Treatment,levels=c("Control","Pre-flowering","Post-flowering"))
data$Weeks<-factor(data$Weeks, levels=c("TP1","TP2","TP3","TP4","TP5","TP6","TP7","TP8","TP9","TP10","TP11","TP12","TP13","TP14","TP15","TP16","TP17"))

ggplot(data, aes(x=Weeks, y=CSWI, fill=Treatment, group=Treatment)) + 
  geom_bar(aes(fill=Treatment),stat="identity",position=position_dodge())+
  geom_errorbar(aes(ymin=CSWI-SE, ymax=CSWI+SE), width=0.2,position=position_dodge(0.9)) + xlab("Time point (in weeks)")+
  ylab("Crop Water Stress Index (CSWI)") + 
  scale_fill_manual(name="",values = c("#6AB187", "#DE7A22", "#F4CC70"))+
  scale_y_continuous(limits = c(-0.1,1.0),breaks = c(-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0))+
   theme(axis.text.x=element_text(hjust=0.5,vjust=0.5,size=11,color="black",face="bold",angle=90), 
        axis.text.y=element_text(size=11,color="black",face="bold"), 
        axis.title=element_text(size=11,face="bold"),text=element_text(size=11))+
  theme(legend.text = element_text(colour="black", size = 11, face = "bold"))

###### Fig.S6 Plant height

nba<-read.table("figS6-PlantHeight.txt",header = T,na.strings = "NA")
nba$Treatment=factor(nba$Treatment,levels=c("Control","Pre_flowering","Post_flowering"))
nba.melt<-melt(nba)

ggplot(nba.melt, aes(x = variable,color=Treatment,y=value))+
  geom_boxplot(outlier.colour=NA)+
  facet_grid(Hybrid~.)+
  theme(axis.text.x=element_text(hjust=1,vjust=0.5,size=11,color="black",angle=90), axis.text.y=element_text(size=11,color="black"), 
        axis.title=element_text(size=11,face="bold"),text=element_text(size=11)) + 
  scale_color_manual(values = c("#6AB187", "#DE7A22", "#F4CC70")) +
  ylab("Plant Height (cm)") + theme(panel.grid.major=element_line(colour=NA))+
  xlab("")+
  theme(legend.title = element_text(colour="black", size=11, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 11, face = "bold"))+
  theme(axis.text.x=element_text(hjust=1,vjust=0.5,size=11,color="black",angle=90,face="bold"), 
        axis.text.y=element_text(size=11,color="black",face="bold"), 
        axis.title=element_text(size=11,face="bold"),text=element_text(size=11,face="bold"))

##### Fig.S7 Grain harvest weight

rm(list = ls())
nba<-read.table("figS7-PlantGrainForage.txt",header = T,na.strings = "NA")
nba$Treatment=factor(nba$Treatment,levels=c("Control","Pre_flowering","Post_flowering"))
nba$Rep<-factor(nba$Rep)
nba<-nba[,-c(7,8)]
nba.melt<-melt(nba)

ggplot(nba, aes(x=Treatment,y=GrainHarvestWeight,colour=Hybrid))+
  geom_boxplot()+
  theme(axis.text.x=element_text(hjust=1,vjust=0.5,size=11,color="black",angle=90), 
        axis.text.y=element_text(size=11,color="black"), 
        axis.title=element_text(size=11,face="bold"),text=element_text(size=11)) + 
  scale_color_manual(values = c("#CB0000", "#50312F")) +
  ylab("Grain Harvest Weight") + theme(panel.grid.major=element_line(colour=NA))+
  xlab("")+
  geom_blank(aes(y = 0)) +
  geom_blank(aes(y = 25))

###### Fig.S8 PERMANOVA R-squared by experimental factor and timepoint

d <- read.table("figS8-permanova_figure.txt",header=T)
sub <- data.frame(d)
sub$Timepoint<-factor(sub$Timepoint, levels=c("TP1","TP2","TP3","TP4","TP5","TP6","TP7","TP8","TP9","TP10","TP11","TP12","TP13","TP14","TP15","TP16","TP17"))
sub$Conditions<-factor(sub$Conditions, levels=c("Control-Pre","Control-Post"))

plot <- ggplot(sub, aes(x=Timepoint, y=R_square,fill=Variation)) + 
  scale_fill_manual(values = c("#DE7A22","#6AB187","#F4CC70"))+
  geom_bar(stat = "identity", color = "NA")+facet_grid(~Conditions,scales="free_x", space = "free_x")+ 
  theme(axis.text.x=element_text(size=11,color="black",angle=90,hjust=1), axis.text.y=element_text(size=11,color="black"), 
        axis.title=element_text(size=14,face="bold"),text=element_text(size=14,face="bold"))+
  ylab("R Square")+xlab("") + theme(legend.title = element_blank())

##### Fig.S9a Bray Curtis distances for all control and pre-flowering drought rhizosphere samples

r = subset_samples(rar,SampleType == "Rhizosphere")
post <- subset_samples(r,Treatment!="Post_flowering")

plot_ordination(post, ordinate(post, "PCoA",distance="bray"), 
                     color = "Timepoint") + 
  scale_colour_manual(values=c("gray8","gray37",
                               brewer.pal(9,"YlGn")[c(9,8,7,6,4,3)],
                               brewer.pal(9,"YlGnBu")[c(3,4,6,7,8)],
                               brewer.pal(9,"RdPu")[c(9,8,7,6)]))+ 
  geom_point(size = 3)+ facet_wrap(~Treatment) +
  theme(axis.text.x=element_text(size=11,color="black",angle=90), axis.text.y=element_text(size=11,color="black"), 
        axis.title=element_text(size=11,face="bold"),text=element_text(size=11,face="bold")) 

##### Fig.S9b bray cutis distances for all control and pre-flowering drought soil samples

r = subset_samples(rar,SampleType == "Soil")
post <- subset_samples(r,Treatment!="Post_flowering")
plot_ordination(post, ordinate(post, "PCoA",distance="bray"), 
                color = "Timepoint") + 
  scale_colour_manual(values=c("gray8","gray37",
                               brewer.pal(9,"YlGn")[c(9,8,7,6,4,3)],
                               brewer.pal(9,"YlGnBu")[c(3,4,6,7,8)],
                               brewer.pal(9,"RdPu")[c(9,8,7,6)]))+ 
  geom_point(size = 3)+ facet_wrap(~Treatment) +
  theme(axis.text.x=element_text(size=11,color="black",angle=90), axis.text.y=element_text(size=11,color="black"), 
        axis.title=element_text(size=11,face="bold"),text=element_text(size=11,face="bold")) 

##### Fig.S10 Heat map of the mean Bray Curtis dissimilarity between pairs of adjacent time points for all replicates (y-axis) in soil, rhizosphere and root samples in the control treatment.

#Root Control compare timepoints
root <- subset_samples(rar,SampleType == "Root")
root_c <- subset_samples(root,Treatment=="Control")

x <- c("TP1","TP2","TP3","TP4","TP5","TP6","TP7","TP8","TP9","TP10","TP11","TP12","TP13","TP14","TP15","TP16")
rm <- {}
r <- subset_samples(root_c,SampleType=="Root")
for(i in x){
  n1 <- sapply(i, function(x) substr(x,3,4))
  n1 <- as.numeric(n1)
  n2 <- n1+1
  i1 <- paste("TP",n1,sep="")
  i2 <- paste("TP",n2,sep="")
  TP <- subset_samples(root_c,Timepoint==i1|Timepoint==i2)
  me <- merge_samples(TP,"Timepoint")
  bray <- phyloseq::distance(me, "bray")
  df <- as.matrix(bray)
  nba.m <- melt(df)
  sub <- subset(nba.m,value!=0)
  sub <- sub[1,]
  colnames(sub) <- c("Sample1","Sample2","Value")
  sub$Timepoint <- paste(i1,i2,sep="_")
  sub$SampleType <- "Root"
  rm <- rbind(rm,sub)
}

#Rhizosphere Control compare timepoints
root <- subset_samples(rar,SampleType == "Rhizosphere")
root_c <- subset_samples(root,Treatment=="Control")

x <- c("TP1","TP2","TP3","TP4","TP5","TP6","TP7","TP8","TP9","TP10","TP11","TP12","TP13","TP14","TP15","TP16")

r <- subset_samples(root_c,SampleType=="Rhizosphere")
for(i in x){
  n1 <- sapply(i, function(x) substr(x,3,4))
  n1 <- as.numeric(n1)
  n2 <- n1+1
  i1 <- paste("TP",n1,sep="")
  i2 <- paste("TP",n2,sep="")
  TP <- subset_samples(root_c,Timepoint==i1|Timepoint==i2)
  me <- merge_samples(TP,"Timepoint")
  bray <- phyloseq::distance(me, "bray")
  df <- as.matrix(bray)
  nba.m <- melt(df)
  sub <- subset(nba.m,value!=0)
  sub <- sub[1,]
  colnames(sub) <- c("Sample1","Sample2","Value")
  sub$Timepoint <- paste(i1,i2,sep="_")
  sub$SampleType <- "Rhizosphere"
  rm <- rbind(rm,sub)
}

#Soil Control compare timepoints
root <- subset_samples(rar,SampleType == "Soil")
root_c <- subset_samples(root,Treatment=="Control")

x <- c("TP1","TP2","TP3","TP4","TP5","TP6","TP7","TP8","TP9","TP10","TP11","TP12","TP13","TP14","TP15","TP16")

r <- subset_samples(root_c,SampleType=="Soil")
for(i in x){
  n1 <- sapply(i, function(x) substr(x,3,4))
  n1 <- as.numeric(n1)
  n2 <- n1+1
  i1 <- paste("TP",n1,sep="")
  i2 <- paste("TP",n2,sep="")
  TP <- subset_samples(root_c,Timepoint==i1|Timepoint==i2)
  me <- merge_samples(TP,"Timepoint")
  bray <- phyloseq::distance(me, "bray")
  df <- as.matrix(bray)
  nba.m <- melt(df)
  sub <- subset(nba.m,value!=0)
  sub <- sub[1,]
  colnames(sub) <- c("Sample1","Sample2","Value")
  sub$Timepoint <- paste(i1,i2,sep="_")
  sub$SampleType <- "Soil"
  rm <- rbind(rm,sub)
}

rm$Timepoint <- factor(rm$Timepoint, levels=c("TP1_TP2","TP2_TP3","TP3_TP4","TP4_TP5","TP5_TP6","TP6_TP7","TP7_TP8",
                                              "TP8_TP9","TP9_TP10","TP10_TP11","TP11_TP12","TP12_TP13","TP13_TP14",
                                              "TP14_TP15","TP15_TP16","TP16_TP17"))
rm$SampleType <- factor(rm$SampleType, levels=c("Soil","Rhizosphere","Root"))

ggplot(data = rm, aes(SampleType,Timepoint)) + 
  geom_tile(aes(fill = Value), colour = "white") +  theme(panel.border=element_blank())+
  facet_wrap(~SampleType,scales="free_x")+
  scale_fill_gradientn(name="Bray curtis distance",colours = terrain.colors(7))+
  theme(axis.text.x=element_text(size=10,color="black",angle=90,face="bold"), 
        axis.text.y=element_text(size=10,color="black",face="bold"), 
        axis.title=element_text(size=10,face="bold"),text=element_text(size=10))+xlab("")+ylab("")+
  theme(legend.title = element_text(colour="black", size=11, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 11, face = "bold"))+
  theme(strip.background = element_blank(),
        strip.text.x = element_blank())

##### Fig.S11 Phyla with significant differences in abundance across time points

data <- read.table("figS11-early-late-phylum.txt",header=T)
data$Phylum <- factor(data$Phylum,levels=c("Proteobacteria",
                                           "Bacteroidetes","Actinobacteria","Firmicutes",
                                           "Acidobacteria","Chloroflexi","Verrucomicrobia",
                                           "TM7","Gemmatimonadetes","Cyanobacteria","Planctomycetes","Armatimonadetes",
                                           "Nitrospirae","Chlorobi",
                                           "Spirochaetes",
                                           "Tenericutes","BRC1","Elusimicrobia","SPAM",
                                           "ZB2","CCM11b","Chlamydiae","GN02","SC3","BacteriaKI","SM2F11"))
data$Stage <- factor(data$Stage,levels=c("Early","Late"))

#Barplot
ggplot(data, aes(x=Phylum, y=log.baseMean, fill=Stage)) + 
  geom_bar(stat = "identity", color = "NA",position = "dodge")+ 
  scale_fill_manual(values = c("#EC96A4","#DFE166")) +
  theme(axis.text.x=element_text(size=11,color="black",angle=90,hjust=1), axis.text.y=element_text(size=11,color="black"), 
        axis.title=element_text(size=14,face="bold"),text=element_text(size=14,face="bold"))+
  ylab("log(baseMean)")+xlab("")

#Heatmap
ggplot(data = data, aes(Phylum,Stage)) + 
  geom_tile(aes(fill = log.padj), colour = "white") +  theme(panel.border=element_blank())+
  scale_fill_gradientn(name="padj",colours = terrain.colors(30))+
  theme(axis.text.x=element_text(size=10,color="black",angle=90,face="bold",hjust=1), 
        axis.text.y=element_text(size=10,color="black",face="bold"), 
        axis.title=element_text(size=10,face="bold"),text=element_text(size=10))+xlab("")+ylab("")+
  theme(legend.title = element_text(colour="black", size=11, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 11, face = "bold"))+
  theme(strip.background = element_blank(),
        strip.text.x = element_blank())+guides(fill=guide_legend(title="log(padj)"))

##### Fig.S12a-c Relative abundance patterns of the most abundant families within the proteobacteria

options(scipen=200)
actino <- subset_taxa(Bushman,Phylum == "Proteobacteria")
b1_Family <- tax_glom(actino, taxrank="Family")
c = subset_samples(b1_Family,Treatment=="Control")
x <- c("TP1","TP2","TP3","TP4","TP5","TP6","TP7","TP8","TP9","TP10","TP11","TP12","TP13","TP14","TP15","TP16","TP17")
rm <- {}
control <- {}
for(i in x){
  TP = subset_samples(c, Timepoint == i)
  merged_r = merge_samples(TP, "SampleType")
  merged_r_m = transform_sample_counts(merged_r , function(x) 100 * x/sum(x))
  physeqdf <- psmelt(merged_r_m)
  temp <- c("Soil","Rhizosphere","Root")
  physeqdf$Sample <- factor(physeqdf$Sample, levels = c(temp))
  n = nrow(physeqdf[order(physeqdf$Family),])
  x = sample(i,n,replace=TRUE)
  ST <- sample("Control",n,replace=TRUE)
  physeqdf[order(physeqdf$Family),]$Timepoint <- x 
  physeqdf[order(physeqdf$Family),]$Treatment <- ST
  physeqdf$Family <- as.character(physeqdf$Family)
  
  control <- rbind(control,physeqdf[order(physeqdf$Family),])
}

control <- data.frame(control)
pre1 <- subset(control,Timepoint=="TP1"|Timepoint=="TP2")
pre1$Treatment[pre1$Treatment == "Control"] <- "Pre_flowering"

pre2 <- {}
pre = subset_samples(b1_Family,Treatment=="Pre_flowering")
y <- c("TP3","TP4","TP5","TP6","TP7","TP8","TP9","TP10","TP11","TP12","TP13","TP14","TP15","TP16","TP17")
for(i in y){
  TP = subset_samples(pre, Timepoint == i)
  merged_r = merge_samples(TP, "SampleType")
  merged_r_m = transform_sample_counts(merged_r , function(x) 100 * x/sum(x))
  physeqdf <- psmelt(merged_r_m)
  temp <- c("Soil","Rhizosphere","Root")
  physeqdf$Sample <- factor(physeqdf$Sample, levels = c(temp))
  n = nrow(physeqdf[order(physeqdf$Family),])
  x = sample(i,n,replace=TRUE)
  ST <- sample("Pre_flowering",n,replace=TRUE)
  physeqdf[order(physeqdf$Family),]$Timepoint <- x 
  physeqdf[order(physeqdf$Family),]$Treatment <- ST
  physeqdf$Family <- as.character(physeqdf$Family)
  
  pre2 <- rbind(pre2,physeqdf[order(physeqdf$Family),])
}

add2 <- subset(control,Treatment=="Control")
add3 <- subset(add2,Timepoint=="TP1"|Timepoint=="TP2"|Timepoint=="TP3"|Timepoint=="TP4"|Timepoint=="TP5"|Timepoint=="TP6"|Timepoint=="TP7")
add3$Treatment[add3$Treatment == "Control"] <- "Post_flowering"

p <- {}
y <- c("TP8","TP9","TP10","TP11","TP12","TP13","TP14","TP15","TP16","TP17")
post = subset_samples(b1_Family,Treatment=="Post_flowering")
for(i in y){
  TP = subset_samples(post, Timepoint == i)
  merged_r = merge_samples(TP, "SampleType")
  merged_r_m = transform_sample_counts(merged_r , function(x) 100 * x/sum(x))
  physeqdf <- psmelt(merged_r_m)
  temp <- c("Soil","Rhizosphere","Root")
  physeqdf$Sample <- factor(physeqdf$Sample, levels = c(temp))
  n = nrow(physeqdf[order(physeqdf$Family),])
  x = sample(i,n,replace=TRUE)
  ST <- sample("Post_flowering",n,replace=TRUE)
  physeqdf[order(physeqdf$Family),]$Timepoint <- x 
  physeqdf[order(physeqdf$Family),]$Treatment <- ST
  physeqdf$Family <- as.character(physeqdf$Family)
  
  p <- rbind(p,physeqdf[order(physeqdf$Family),])
}
rm <- rbind(control,pre1,pre2,add3,p)
sum <- tapply(rm$Abundance,rm$Family, FUN=sum)
sum <- data.frame(sum)
sum <- sort(sum$sum,decreasing = T)
list <- rownames(sum[1:15])
a <- rm[rm$Family %in% list,]
b <- rm[!(rm$Family %in% list),]
b$Family <- "Other"
rm <- rbind(a,b)
rm$Treatment<-factor(rm$Treatment,levels=c("Control","Pre_flowering","Post_flowering"))
rm$Timepoint<-factor(rm$Timepoint, levels=c("TP1","TP2","TP3","TP4","TP5","TP6","TP7","TP8","TP9","TP10","TP11","TP12","TP13","TP14","TP15","TP16","TP17"))
col_SampleByTreatment<-c(colorRampPalette(brewer.pal(11,"Paired"))(8),colorRampPalette(brewer.pal(8,"Set1"))(8))
rm$Family <- factor(rm$Family,levels=c("Alteromonadaceae","Bradyrhizobiaceae","BurkholderialesOR","Comamonadaceae",
                                       "Haliangiaceae","Hyphomicrobiaceae","MyxococcalesOR","Oxalobacteraceae",
                                       "Rhizobiaceae","RhodocyclalesOR","Rhodospirillaceae","Sinobacteraceae",
                                       "Sphingomonadaceae","Syntrophobacteraceae","Xanthomonadaceae","Other" ))

sub <- subset(rm,Treatment=="Control")

ggplot(sub, aes(x=Timepoint, y=Abundance, fill=Family, )) + ylab("Relative Abundance") +
  xlab("")+
  geom_bar(stat = "identity", color = "NA")+facet_wrap(~Sample)+ 
  scale_fill_manual(values = col_SampleByTreatment) +
  theme(axis.text.x=element_text(size=14,color="black",angle=90), 
        axis.text.y=element_text(size=14,color="black"), 
        axis.title=element_text(size=14,face="bold"),
        text=element_text(size=14,face="bold")) 

##### Fig.S12d-e Relative abundance patterns of the most abundant families within the actinobacteria

options(scipen=200)
actino <- subset_taxa(Bushman,Phylum == "Actinobacteria")
b1_Family <- tax_glom(actino, taxrank="Family")
c = subset_samples(b1_Family,Treatment=="Control")
x <- c("TP1","TP2","TP3","TP4","TP5","TP6","TP7","TP8","TP9","TP10","TP11","TP12","TP13","TP14","TP15","TP16","TP17")
rm <- {}
control <- {}
for(i in x){
  TP = subset_samples(c, Timepoint == i)
  merged_r = merge_samples(TP, "SampleType")
  merged_r_m = transform_sample_counts(merged_r , function(x) 100 * x/sum(x))
  physeqdf <- psmelt(merged_r_m)
  temp <- c("Soil","Rhizosphere","Root")
  physeqdf$Sample <- factor(physeqdf$Sample, levels = c(temp))
  n = nrow(physeqdf[order(physeqdf$Family),])
  x = sample(i,n,replace=TRUE)
  ST <- sample("Control",n,replace=TRUE)
  physeqdf[order(physeqdf$Family),]$Timepoint <- x 
  physeqdf[order(physeqdf$Family),]$Treatment <- ST
  physeqdf$Family <- as.character(physeqdf$Family)
  
  control <- rbind(control,physeqdf[order(physeqdf$Family),])
}

control <- data.frame(control)
pre1 <- subset(control,Timepoint=="TP1"|Timepoint=="TP2")
pre1$Treatment[pre1$Treatment == "Control"] <- "Pre_flowering"

pre2 <- {}
pre = subset_samples(b1_Family,Treatment=="Pre_flowering")
y <- c("TP3","TP4","TP5","TP6","TP7","TP8","TP9","TP10","TP11","TP12","TP13","TP14","TP15","TP16","TP17")
for(i in y){
  TP = subset_samples(pre, Timepoint == i)
  merged_r = merge_samples(TP, "SampleType")
  merged_r_m = transform_sample_counts(merged_r , function(x) 100 * x/sum(x))
  physeqdf <- psmelt(merged_r_m)
  temp <- c("Soil","Rhizosphere","Root")
  physeqdf$Sample <- factor(physeqdf$Sample, levels = c(temp))
  n = nrow(physeqdf[order(physeqdf$Family),])
  x = sample(i,n,replace=TRUE)
  ST <- sample("Pre_flowering",n,replace=TRUE)
  physeqdf[order(physeqdf$Family),]$Timepoint <- x 
  physeqdf[order(physeqdf$Family),]$Treatment <- ST
  physeqdf$Family <- as.character(physeqdf$Family)
  
  pre2 <- rbind(pre2,physeqdf[order(physeqdf$Family),])
}

add2 <- subset(control,Treatment=="Control")
add3 <- subset(add2,Timepoint=="TP1"|Timepoint=="TP2"|Timepoint=="TP3"|Timepoint=="TP4"|Timepoint=="TP5"|Timepoint=="TP6"|Timepoint=="TP7")
add3$Treatment[add3$Treatment == "Control"] <- "Post_flowering"

p <- {}
y <- c("TP8","TP9","TP10","TP11","TP12","TP13","TP14","TP15","TP16","TP17")
post = subset_samples(b1_Family,Treatment=="Post_flowering")
for(i in y){
  TP = subset_samples(post, Timepoint == i)
  merged_r = merge_samples(TP, "SampleType")
  merged_r_m = transform_sample_counts(merged_r , function(x) 100 * x/sum(x))
  physeqdf <- psmelt(merged_r_m)
  temp <- c("Soil","Rhizosphere","Root")
  physeqdf$Sample <- factor(physeqdf$Sample, levels = c(temp))
  n = nrow(physeqdf[order(physeqdf$Family),])
  x = sample(i,n,replace=TRUE)
  ST <- sample("Post_flowering",n,replace=TRUE)
  physeqdf[order(physeqdf$Family),]$Timepoint <- x 
  physeqdf[order(physeqdf$Family),]$Treatment <- ST
  physeqdf$Family <- as.character(physeqdf$Family)
  
  p <- rbind(p,physeqdf[order(physeqdf$Family),])
}

rm <- rbind(control,pre1,pre2,add3,p)
sum <- tapply(rm$Abundance,rm$Family, FUN=sum)
sum <- data.frame(sum)
sum <- sort(sum$sum,decreasing = T)
list <- rownames(sum[1:10])
a <- rm[rm$Family %in% list,]
b <- rm[!(rm$Family %in% list),]
b$Family <- "Other"
rm <- rbind(a,b)
rm$Treatment<-factor(rm$Treatment,levels=c("Control","Pre_flowering","Post_flowering"))
rm$Timepoint<-factor(rm$Timepoint, levels=c("TP1","TP2","TP3","TP4","TP5","TP6","TP7","TP8","TP9","TP10","TP11","TP12","TP13","TP14","TP15","TP16","TP17"))
col_SampleByTreatment<-c(colorRampPalette(brewer.pal (11, "Set3" ))(11))
rm$Family <- factor(rm$Family,levels=c("AcidimicrobialesOR","Actinosynnemataceae","Glycomycetaceae","Micrococcaceae",
                                       "Micromonosporaceae","Nocardioidaceae","Promicromonosporaceae","Pseudonocardiaceae",
                                       "SolirubrobacteralesOR",
                                       "Streptosporangiaceae","Streptomycetaceae","Other"))

sub <- subset(rm,Treatment=="Control")

ggplot(sub, aes(x=Timepoint, y=Abundance, fill=Family, )) + 
  geom_bar(stat = "identity", color = "NA")+facet_wrap(~Sample)+ 
  scale_fill_manual(values = col_SampleByTreatment) +
  theme(axis.text.x=element_text(size=14,color="black",angle=90), 
        axis.text.y=element_text(size=14,color="black"), 
        axis.title=element_text(size=14,face="bold"),
        text=element_text(size=14,face="bold")) + ylab("Relative Abundance")+xlab("")

##### Fig.S13a Phylogenetic tree for all genera with significant differences in abundance between root and soil samples in the control treatment across TP3 to TP8

root <- subset_samples(Bushman,SampleType=="Root"|SampleType=="Soil")
tp8 = subset_samples(root,Timepoint=="TP3"|Timepoint=="TP4"|Timepoint=="TP5"|Timepoint=="TP6"|Timepoint=="TP7"|Timepoint=="TP8")
tp8_new = subset_samples(tp8,Treatment=="Control")
tp8_new= prune_taxa(taxa_sums(tp8_new)>=1, tp8_new)
# This converts to a DESeq2 format.
diagdds = phyloseq_to_deseq2(tp8_new, ~ SampleType)
# This helps to obtain the normalization factor for your dataset.
diagdds = estimateSizeFactors(diagdds)
# You also need to demonstrate that your dataset is effectively dispersed (?).
diagdds = estimateDispersions(diagdds)
# This will reduce sampling variance.
diagvst = getVarianceStabilizedData(diagdds)
dim(diagvst)

diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
res = results(diagdds, cooksCutoff = FALSE)
head(res,10)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(tp8_new)[rownames(sigtab), ], "matrix"))
#delete the OTUs which didn't assign a phylum name
new <- subset(sigtab,Phylum != "NA") 
Genus <- subset(new,Genus != "NA")
na <- setdiff(rownames(new),rownames(Genus))
na.complete <- new[na,]
#If the OTUs' class name is NA, then paste phylum name to class,order,family and genus
class <- na.complete[is.na(na.complete$Class),]
class$Class <- paste(class$Phylum,"PH",sep="")
class$Order <- paste(class$Phylum,"PH",sep="")
class$Family <- paste(class$Phylum,"PH",sep="")
class$Genus <- paste(class$Phylum,"PH",sep="")
order <- na.complete[is.na(na.complete$Order),]
order2 <- subset(order,Class!="NA")
#If the OTUs' order name is NA, then paste class name to order,family and genus
order2$Order <- paste(order2$Class,"OR",sep="")
order2$Family <- paste(order2$Class,"OR",sep="")
order2$Genus <- paste(order2$Class,"OR",sep="")

Family <- na.complete[is.na(na.complete$Family),]
Family2 <- subset(Family,Order!="NA")
Family2$Family <- paste(Family2$Order,"CL",sep="")
Family2$Genus <- paste(Family2$Order,"CL",sep="")

Ge <- na.complete[is.na(na.complete$Genus),]
Ge2 <- subset(Ge,Family!="NA")
Ge2$Genus <- paste(Ge2$Family,"FA",sep="")

all <- rbind(Ge2,Family2,order,class,Genus)
all$taxa <- paste(all$Phylum,all$Genus,sep="_")
mean <- tapply(all$log2FoldChange,all$taxa,mean)
as <- as.character(sort(as.numeric(row.names(all))))
all2 <- all[as,]
uniq <- all2[!duplicated(all2$taxa),]
uniq$name <- paste(rownames(uniq),uniq$log2FoldChange,uniq$taxa,sep="_")
write.csv(uniq,"fold_change.csv")

#input the data to iTOL for plotting the tree

##### 
#Fig.S13b
# Open the file containing the phyloseq object.
sorghum_data <- readRDS("Bushman.rds")

# Subset to all the control samples from time points 3-8 in roots and soil only.
rs38 = subset_samples(sorghum_data,Treatment=="Control")
rs38_new = subset_samples(rs38,SampleType=="Soil"|SampleType=="Root")
rs38_new = subset_samples(rs38_new, Timepoint=="TP3"|Timepoint=="TP4"|Timepoint=="TP5"|Timepoint=="TP6"|Timepoint=="TP7"|Timepoint=="rs38")
rs38_new= prune_taxa(taxa_sums(rs38_new)>=1, rs38_new)

###############################################################################
### Generating tables with readcounts, genera, etc. that we will use later. ###
###############################################################################

# Generate a table containing all of the summed readcounts for each phylum.
rs38_phylum <- tax_glom(rs38_new, taxrank = "Phylum")
phylum_tax <- (data.frame(tax_table(rs38_phylum)))
phylum_otu <- (data.frame(rowSums(otu_table(rs38_phylum))))
phylum_table <- cbind(phylum_tax[,2], phylum_otu)
colnames(phylum_table) <- c("Phylum", "Readcounts")
phylum_table$Phylum <- as.character(phylum_table$Phylum)

# Generate a table with all the genera with their readcounts.
rs38_genus <- tax_glom(rs38_new, taxrank = "Genus")
full_table <- cbind(tax_table(rs38_genus), data.frame(rowSums(otu_table(rs38_genus))))
full_table$Phylum <- as.character(full_table$Phylum)
full_table$Genus <- as.character(full_table$Genus)

# Generate a table with all the readcounts and taxonomy.
rs38_overall <- cbind(data.frame(tax_table(rs38_new)), data.frame(otu_table(rs38_new)))


################################
### Get the enriched genera. ###
################################

# This converts to a DESeq2 format.
diagdds = phyloseq_to_deseq2(rs38_new, ~ SampleType)
# This helps to obtain the normalization factor for your dataset.
diagdds = estimateSizeFactors(diagdds)
# You also need to demonstrate that your dataset is effectively dispersed (?).
diagdds = estimateDispersions(diagdds)
# This will reduce sampling variance.
diagvst = getVarianceStabilizedData(diagdds)

diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
res = results(diagdds, cooksCutoff = FALSE)

alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(rs38_new)[rownames(sigtab), ], "matrix"))
sigtab[sigtab==""] <- NA

sigtab$taxa <- paste(sigtab$Phylum, sigtab$Genus, sep = "_")
sigtab_mean <- tapply(sigtab$log2FoldChange,sigtab$taxa,mean)
as <- as.character(sort(as.numeric(row.names(sigtab))))
all2 <- sigtab[as,]
uniq <- all2[!duplicated(all2$taxa),]
uniq$name <- paste(rownames(uniq),uniq$log2FoldChange,uniq$taxa,sep="_")
uniq <- uniq[-which(is.na(uniq$Genus)),]
uniq$Phylum <- as.character(uniq$Phylum)
uniq$Genus <- as.character(uniq$Genus)
uniq_ordered <- uniq[order(uniq$log2FoldChange),]


# Get tables for only the enriched or depleted taxa.
root_enriched <- uniq[uniq$log2FoldChange < 0.00,]
soil_enriched <- uniq[uniq$log2FoldChange > 0.00,]


#########################################
### Start generating the final table. ###
#########################################

# 1. First step is to get all the phylum in the table. 

phylum_list <- (c(rep("Proteobacteria", 2), rep("Firmicutes", 2), rep("Bacteroidetes", 2), rep("Chloroflexi", 2),
                  rep("Verrucomicrobia", 2), rep("Actinobacteria", 2), rep("Cyanobacteria", 2), rep("Acidobacteria", 2)))
phylum_list <- data.frame(phylum_list)
colnames(phylum_list) <- c("Phylum")

# 2. Next, set up the root/soil enrichment status.

enrichment_list <- c(rep(c("Root", "Soil"),8))
enrichment_list <- data.frame(enrichment_list)
colnames(enrichment_list) <- c("Enrichment")

total_list <- cbind(phylum_list, enrichment_list)

# 3. Find out the number of genera that are enriched in soil or roots in all significantly enriched taxa. 
library(rlist)

enrichment <- data.frame(matrix(nrow=0,ncol=2))
colnames(enrichment) <- c("Genera_Sig", "Percent_Genera_Sig")

for(i in 1:nrow(total_list)){
  phylum <- total_list$Phylum[i]
  compartment <- total_list$Enrichment[i]
  if(compartment == "Root"){
    x <- sum(root_enriched$Phylum == phylum)
    y <- sum(root_enriched$Phylum == phylum) / sum(uniq$Phylum == phylum)
  }
  if(compartment == "Soil"){
    x <- sum(soil_enriched$Phylum == phylum)
    y <- sum(soil_enriched$Phylum == phylum) / sum(uniq$Phylum == phylum)
  }
  enrichment[i,1] <- x
  enrichment[i,2] <- y
}

total_list <- cbind(total_list, enrichment)

# 4. Find out the percent of genera present in the significantly enriched taxa compared with all taxa for a given phylum.

total_enrichment <- data.frame(matrix(nrow=0,ncol=1))
colnames(total_enrichment) <- c("Percent_Genera_Total")

for(i in 1:nrow(total_list)){
  phylum <- total_list$Phylum[i]
  x <- sum(uniq$Phylum == phylum) / sum(full_table$Phylum == phylum)
  total_enrichment[i,1] <- x
}

total_list <- cbind(total_list, total_enrichment)

# 5. Find out the readcounts represented by the significantly enriched lineages compared to the overall dataset.

# We originally generated the table with all the readcounts with taxonomies (rs38_overall). Use this.
# First subset it so that it only includes the rows that are included in uniq.

sig_readcounts <- as.data.frame(rs38_overall[rownames(rs38_overall) %in% rownames(uniq),])
sig_rowsums <- as.data.frame(rowSums(sig_readcounts[,7:ncol(sig_readcounts)]))
colnames(sig_rowsums) <- c("Sig_Rowsums")
sig_readcounts <- cbind(sig_readcounts, sig_rowsums)
sig_readcounts$Phylum <- as.character(sig_readcounts$Phylum)

total_readcounts <- data.frame(matrix(nrow=0,ncol=1))
colnames(total_readcounts) <- c("Percent_Readcounts")

for(i in 1:nrow(total_list)){
  phylum <- total_list$Phylum[i]
  sig_phycounts <- sum(sig_readcounts[which(sig_readcounts$Phylum == phylum),ncol(sig_readcounts)]) 
  all_phycounts <- sum(phylum_table[which(phylum_table$Phylum == phylum), 2])
  x <- sig_phycounts / all_phycounts
  total_readcounts[i,1] <- x
}

total_list <- cbind(total_list, total_readcounts)

write.table(total_list, "S13b_values.txt", sep = "\t")

# Now constructing the figure.

total_list$Phylum <- factor(total_list$Phylum,levels=c("Tenericutes","Bacteroidetes","TM7",
                                                       "Verrucomicrobia","Proteobacteria",
                                                       "Actinobacteria","Acidobacteria",
                                                       "Cyanobacteria","Firmicutes",
                                                       "Chloroflexi","Planctomycetes"))

plot <- ggplot(total_list, aes(x=Phylum, y=Percent_Genera_Sig, fill=Enrichment)) + 
  geom_bar(stat = "identity", color = "NA")+
  scale_fill_manual(values = c("#CD853F", "#000000")) +
  theme(axis.text.x=element_text(size=11,color="black"), 
        axis.text.y=element_text(size=11,color="black"), 
        axis.title=element_text(size=11,face="bold"),text=element_text(size=11,face="bold"))+
  ylab("Percentage")+xlab("")+ coord_flip()

ggsave("S13b.jpg", plot = plot, height = 6, width = 5)

######Fig.S14 qPCR for field samples

data <- read.table("figS14-actino1_field_qPCR.txt",header=T)

data$Timepoint <- factor(data$Timepoint,levels=c("TP4","TP8"))
data$Taxa <- factor(data$Taxa,levels=c("Bacteria","Actinobacteria","Firmicutes","Proteobacteria"))
#barplot
plot <- ggplot(data, aes(x=Taxa, y=Value, fill=Treatment)) + 
  geom_bar(stat = "identity", color = "NA",position = "dodge")+ 
  scale_fill_manual(values = c("#6AB187", "#DE7A22")) + facet_wrap(~Timepoint+Genotype+Normalization,ncol=6)+
  theme(axis.text.x=element_text(size=11,color="black",angle=90,hjust=1), axis.text.y=element_text(size=11,color="black"), 
        axis.title=element_text(size=14,face="bold"),text=element_text(size=14,face="bold"))+
  ylab("Relative Abundance")+xlab("") + geom_errorbar(aes(ymin=Value-SD, ymax=Value+SD),width=.2, 
                                                      position=position_dodge(.9),color="grey50")

ggsave(filename = "figS14-actino1_field_qPCR.jpg",plot = plot,width=15,height=10)

######Fig.S15

#Open the file containing the phyloseq object.
sorghum_data <- readRDS("Bushman.rds")

# Subset to the root samples from time point 8.
tp8 = subset_samples(sorghum_data,SampleType=="Root")
tp8_new = subset_samples(tp8,Timepoint=="TP8")
tp8_new = subset_samples(tp8_new, Treatment != "Post_flowering")
tp8_new= prune_taxa(taxa_sums(tp8_new)>=1, tp8_new)

###############################################################################
### Generating tables with readcounts, genera, etc. that we will use later. ###
###############################################################################

# Generate a table containing all of the summed readcounts for each phylum.
tp8_phylum <- tax_glom(tp8_new, taxrank = "Phylum")
phylum_tax <- (data.frame(tax_table(tp8_phylum)))
phylum_otu <- (data.frame(rowSums(otu_table(tp8_phylum))))
phylum_table <- cbind(phylum_tax[,2], phylum_otu)
colnames(phylum_table) <- c("Phylum", "Readcounts")
phylum_table$Phylum <- as.character(phylum_table$Phylum)


# Generate a table with all the genera with their readcounts.
tp8_genus <- tax_glom(tp8_new, taxrank = "Genus")
full_table <- cbind(tax_table(tp8_genus), data.frame(rowSums(otu_table(tp8_genus))))
full_table$Phylum <- as.character(full_table$Phylum)
full_table$Genus <- as.character(full_table$Genus)

# Generate a table with all the readcounts and taxonomy.
tp8_overall <- cbind(data.frame(tax_table(tp8_new)), data.frame(otu_table(tp8_new)))

################################
### Get the enriched genera. ###
################################

# This converts to a DESeq2 format.
diagdds = phyloseq_to_deseq2(tp8_new, ~ Treatment)
# This helps to obtain the normalization factor for your dataset.
diagdds = estimateSizeFactors(diagdds)
# You also need to demonstrate that your dataset is effectively dispersed (?).
diagdds = estimateDispersions(diagdds)
# This will reduce sampling variance.
diagvst = getVarianceStabilizedData(diagdds)

diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
res = results(diagdds, cooksCutoff = FALSE)

alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(tp8_new)[rownames(sigtab), ], "matrix"))
sigtab[sigtab==""] <- NA

sigtab$taxa <- paste(sigtab$Phylum, sigtab$Genus, sep = "_")
sigtab_mean <- tapply(sigtab$log2FoldChange,sigtab$taxa,mean)
as <- as.character(sort(as.numeric(row.names(sigtab))))
all2 <- sigtab[as,]
uniq <- all2[!duplicated(all2$taxa),]
uniq$name <- paste(rownames(uniq),uniq$log2FoldChange,uniq$taxa,sep="_")
uniq <- uniq[-which(is.na(uniq$Genus)),]
uniq$Phylum <- as.character(uniq$Phylum)
uniq$Genus <- as.character(uniq$Genus)
uniq_ordered <- uniq[order(uniq$log2FoldChange),]

# Get tables for only the enriched or depleted taxa.
control_enriched <- uniq[uniq$log2FoldChange < 0.00,]
drought_enriched <- uniq[uniq$log2FoldChange > 0.00,]

#########################################
### Start generating the final table. ###
#########################################

# 1. First step is to get all the phylum in the table. 

phylum_list <- (c(rep("Proteobacteria", 2), rep("Firmicutes", 2), rep("Bacteroidetes", 2), rep("Chloroflexi", 2),
                  rep("Verrucomicrobia", 2), rep("Actinobacteria", 2), rep("Cyanobacteria", 2), rep("Acidobacteria", 2)))
phylum_list <- data.frame(phylum_list)
colnames(phylum_list) <- c("Phylum")

# 2. Next, set up the root/soil enrichment status.

enrichment_list <- c(rep(c("Drought", "Control"),8))
enrichment_list <- data.frame(enrichment_list)
colnames(enrichment_list) <- c("Enrichment")

total_list <- cbind(phylum_list, enrichment_list)

# 3. Find out the number of genera that are enriched in soil or roots in all significantly enriched taxa. 
library(rlist)

enrichment <- data.frame(matrix(nrow=0,ncol=2))
colnames(enrichment) <- c("Genera_Sig", "Percent_Genera_Sig")

for(i in 1:nrow(total_list)){
  phylum <- total_list$Phylum[i]
  treatment <- total_list$Enrichment[i]
  if(treatment == "Drought"){
    x <- sum(drought_enriched$Phylum == phylum)
    y <- sum(drought_enriched$Phylum == phylum) / sum(uniq$Phylum == phylum)
  }
  if(treatment == "Control"){
    x <- sum(control_enriched$Phylum == phylum)
    y <- sum(control_enriched$Phylum == phylum) / sum(uniq$Phylum == phylum)
  }
  enrichment[i,1] <- x
  enrichment[i,2] <- y
}

total_list <- cbind(total_list, enrichment)

# 4. Find out the percent of genera present in the significantly enriched taxa compared with all taxa for a given phylum.

total_enrichment <- data.frame(matrix(nrow=0,ncol=1))
colnames(total_enrichment) <- c("Percent_Genera_Total")

for(i in 1:nrow(total_list)){
  phylum <- total_list$Phylum[i]
  x <- sum(uniq$Phylum == phylum) / sum(full_table$Phylum == phylum)
  total_enrichment[i,1] <- x
}

total_list <- cbind(total_list, total_enrichment)

# 5. Find out the readcounts represented by the significantly enriched lineages compared to the overall dataset.

# We originally generated the table with all the readcounts with taxonomies (tp8_overall). Use this.
# First subset it so that it only includes the rows that are included in uniq.
sig_readcounts <- as.data.frame(tp8_overall[rownames(tp8_overall) %in% rownames(uniq),])
sig_rowsums <- as.data.frame(rowSums(sig_readcounts[,7:ncol(sig_readcounts)]))
colnames(sig_rowsums) <- c("Sig_Rowsums")
sig_readcounts <- cbind(sig_readcounts, sig_rowsums)
sig_readcounts$Phylum <- as.character(sig_readcounts$Phylum)

total_readcounts <- data.frame(matrix(nrow=0,ncol=1))
colnames(total_readcounts) <- c("Percent_Readcounts")

for(i in 1:nrow(total_list)){
  phylum <- total_list$Phylum[i]
  sig_phycounts <- sum(sig_readcounts[which(sig_readcounts$Phylum == phylum),ncol(sig_readcounts)]) 
  all_phycounts <- sum(phylum_table[which(phylum_table$Phylum == phylum), 2])
  x <- sig_phycounts / all_phycounts
  total_readcounts[i,1] <- x
}

total_list <- cbind(total_list, total_readcounts)

write.table(total_list, "S14_values.txt", sep = "\t")

# Now constructing the figure.

total_list$Phylum <- factor(total_list$Phylum,levels=c("Bacteroidetes", "Verrucomicrobia",
                                                       "Proteobacteria", "Acidobacteria", "Chloroflexi",
                                                       "Cyanobacteria", "Actinobacteria", "Firmicutes"))

# We don't want Cyanobacteria or Acidobacteria.
total_list <- total_list[total_list$Phylum != "Cyanobacteria",]
total_list <- total_list[total_list$Phylum != "Acidobacteria",]

#######################

plot <- ggplot(total_list, aes(x=Phylum, y=Percent_Genera_Sig, fill=Enrichment)) + 
  geom_bar(stat = "identity", color = "NA")+
  scale_fill_manual(values = c("#6AB187","#DE7A22")) +
  theme(axis.text.x=element_text(size=11,color="black"), 
        axis.text.y=element_text(size=11,color="black"), 
        axis.title=element_text(size=11,face="bold"),text=element_text(size=11,face="bold"))+
  ylab("Percentage")+xlab("")+ coord_flip()

plot

ggsave("S15.jpg", plot = plot, height = 6, width = 5)

######Fig. S16

test<-read.table(file = "RootSoilRhizo9.txt",header = T)

test$SampleType<-factor(test$SampleType,levels=c("Root","Rhizosphere","Soil"))
test$Taxa<-factor(test$Taxa,levels=(test$Taxa)[rev(order(test$Order))])
test$PhylumRank<-as.factor(test$PhylumRank)

same<-subset(test,SignTest==0)
diff<-subset(test,SignTest==1)

#Fig 16A
ggplot(data=same,aes(SampleType,Taxa))+
  geom_tile(aes(fill=Score))+
  scale_fill_gradient2(low = "#6AB187", high = "#DE7A22")+
  theme(axis.text.x=element_text(hjust=1,vjust=0.5,size=10,color="black",angle=90), axis.text.y=element_text(size=2), axis.title=element_text(size=16,face="bold"),text=element_text(size=16))

#Fig 16B
ggplot(data=diff,aes(SampleType,Taxa))+
  geom_tile(aes(fill=Score))+
  scale_fill_gradient2(low = "#6AB187", high = "#DE7A22")+
  theme(axis.text.x=element_text(hjust=1,vjust=0.5,size=10,color="black",angle=90), axis.text.y=element_text(size=4), axis.title=element_text(size=16,face="bold"),text=element_text(size=16))


#Fig16A, Phyla
ggplot(data=same,aes(SampleType,Taxa))+
  geom_tile(aes(fill=PhylumRank),colour="white")+
  scale_fill_manual(values=c("lightskyblue", "lightblue1", "lightsalmon2", "cadetblue3", "cornsilk", "thistle","darkseagreen" ,"yellow", "slategray1", "lightpink", "plum"))+
  theme(axis.text.x=element_text(hjust=1,vjust=0.5,size=10,color="black",angle=90), axis.text.y=element_text(size=4), axis.title=element_text(size=16,face="bold"),text=element_text(size=16))

#Fig16A, Gram Status
ggplot(data=same,aes(SampleType,Taxa))+
  geom_tile(aes(fill=Status),colour="white")+
  scale_fill_manual(values=c("#669592","#CF9B54","#2C333F"))+
  theme(axis.text.x=element_text(hjust=1,vjust=0.5,size=10,color="black",angle=90), axis.text.y=element_text(size=4), axis.title=element_text(size=16,face="bold"),text=element_text(size=16))

#Fig 16B, Phylum
ggplot(data=diff,aes(SampleType,Taxa))+
  geom_tile(aes(fill=PhylumRank),colour="white")+
  scale_fill_manual(values=c("lightskyblue", "lightblue1", "lightsalmon2", "cadetblue3", "cornsilk","thistle", "darkseagreen","mistyrose2"))+
  theme(axis.text.x=element_text(hjust=1,vjust=0.5,size=10,color="black",angle=90), axis.text.y=element_text(size=4), axis.title=element_text(size=16,face="bold"),text=element_text(size=16))

#Fig 16B, Gram Status
ggplot(data=diff,aes(SampleType,Taxa))+
  geom_tile(aes(fill=Status),colour="white")+
  scale_fill_manual(values=c("#669592","#CF9B54","#2C333F"))+
  theme(axis.text.x=element_text(hjust=1,vjust=0.5,size=10,color="black",angle=90), axis.text.y=element_text(size=4), axis.title=element_text(size=16,face="bold"),text=element_text(size=16))

######Fig.S17b

#load R data and source itag_diversity scripts
sorghum_data <- readRDS("NewSorghumPilot_Phyloseq.rds")
sorghum_data <- rarefy_even_depth(sorghum_data)
source('/Users/dtnaylor/Desktop/Thesis/Extra_Analyses/Chapter_4_Gill/itag_diversity.R')

# Generate the phylum table
sorghum_phylum <- tax_glom(sorghum_data, taxrank="Phylum")
phylum_otus <- otu_table(sorghum_phylum)
rownames(phylum_otus) <- as.data.frame(tax_table(sorghum_phylum))$Phylum
phylum_otus <- t(phylum_otus)

condensed_OTU<-phylum_otus[,c("Actinobacteria","Chloroflexi","Cyanobacteria","Firmicutes","Bacteroidetes","Proteobacteria","Acidobacteria",
                              "Gemmatimonadetes","Nitrospirae","Planctomycetes","Tenericutes","TM7","Verrucomicrobia")]
others<-as.data.frame(rowSums(phylum_otus[,!(colnames(phylum_otus) %in% c("Actinobacteria","Chloroflexi","Cyanobacteria","Firmicutes","Bacteroidetes","Proteobacteria","Acidobacteria",
                                                                          "Gemmatimonadetes","Nitrospirae","Planctomycetes","Tenericutes","TM7","Verrucomicrobia"))]))
colnames(others)<-c("Other")

# Add in the metadata
env<-sample_data(sorghum_phylum)
env<-env[,1:4]
condensed_OTU<-cbind(env,condensed_OTU,others)
condensed_OTU$Timepoint <- as.character(condensed_OTU$Timepoint)
condensed_OTU$SampleType <- as.character(condensed_OTU$SampleType)

for(i in 1:nrow(condensed_OTU)){
  condensed_OTU$Timepoint[i] <- (paste("TP", condensed_OTU$Timepoint[i], sep = ""))
}

condensed_OTU$Timepoint <- factor(condensed_OTU$Timepoint, levels = c("TP3", "TP4", "TP5", "TP6", "TP7", "TP8", "TP9", "TP10", "TP11", "TP12", "TP13", "TP14"))
condensed_OTU$SampleType <- factor(condensed_OTU$SampleType, levels = c("Soil", "Rhizo", "Root"))

condensed_OTU_melt<-melt(data = condensed_OTU)
condensed_OTU_melt$variable<-factor(condensed_OTU_melt$variable,levels=levels(condensed_OTU_melt$variable))
#levels(condensed_OTU_melt$Timepoint) <- c("TP3", "TP4", "TP5", "TP6", "TP7", "TP8", "TP9", "TP10", "TP11", "TP12", "TP13", "TP14")
colnames(condensed_OTU_melt)[5] <- c("Phylum")

col_SampleByTreatment<-c("Acidobacteria"="thistle","Actinobacteria"="lightsalmon2","Armatimonadetes"="mediumpurple1","Bacteroidetes"="cornsilk","Chloroflexi"="lightskyblue","Cyanobacteria"="lightpink","Firmicutes"="lightblue1","Gemmatimonadetes"="darkseagreen","Nitrospirae"="mistyrose2",
                         "Planctomycetes"="gray69","Proteobacteria"="cadetblue3","SPAM"="darkorchid1","Tenericutes"="darkolivegreen2","TM7"="yellow","Verrucomicrobia"="slategray1","Other"="plum")



plot<-ggplot(data=condensed_OTU_melt,aes(x=Timepoint,y=value,fill=Phylum)) + 
  geom_bar(stat="identity",position = "fill") + 
  facet_wrap(SampleType~Treatment, ncol=2) +
  scale_fill_manual(values = col_SampleByTreatment) +
  ylab("Relative Abundance") +
  scale_y_continuous(labels = c("0", "25", "50", "75", "100")) +
  theme(axis.text.x=element_text(size=16,color="black",angle=90), axis.text.y=element_text(size=16,color="black"), axis.title=element_text(size=16,face="bold"),text=element_text(size=16)) 
ggsave("Revised_RelativeAbundance_EPICON.jpg",plot,width=7,height=8)

save.image(file = "SorghumPilot_RelativeAbundance.RData")

##### Fig.S19 Relative abundance of the top thirteen most abundant phyla for control (top panels) and post-flowering drought (bottom panels) treatments across different sample types

options(scipen=200)
b1_phylum <- tax_glom(Bushman, taxrank="Phylum")

c = subset_samples(b1_phylum,Treatment=="Control")
x <- c("TP1","TP2","TP3","TP4","TP5","TP6","TP7","TP8","TP9","TP10","TP11","TP12","TP13","TP14","TP15","TP16","TP17")
rm <- {}
control <- {}
for(i in x){
  TP = subset_samples(c, Timepoint == i)
  merged_r = merge_samples(TP, "SampleType")
  merged_r_m = transform_sample_counts(merged_r , function(x) 100 * x/sum(x))
  physeqdf <- psmelt(merged_r_m)
  temp <- c("Soil","Rhizosphere","Root")
  physeqdf$Sample <- factor(physeqdf$Sample, levels = c(temp))
  n = nrow(physeqdf[order(physeqdf$Phylum),])
  x = sample(i,n,replace=TRUE)
  ST <- sample("Control",n,replace=TRUE)
  physeqdf[order(physeqdf$Phylum),]$Timepoint <- x 
  physeqdf[order(physeqdf$Phylum),]$Treatment <- ST
  physeqdf$Phylum <- as.character(physeqdf$Phylum)
  control <- rbind(control,physeqdf[order(physeqdf$Phylum),])
}

control <- data.frame(control)
pre1 <- subset(control,Timepoint=="TP1"|Timepoint=="TP2")
pre1$Treatment[pre1$Treatment == "Control"] <- "Pre_flowering"

pre2 <- {}
pre = subset_samples(b1_phylum,Treatment=="Pre_flowering")
y <- c("TP3","TP4","TP5","TP6","TP7","TP8","TP9","TP10","TP11","TP12","TP13","TP14","TP15","TP16","TP17")
for(i in y){
  TP = subset_samples(pre, Timepoint == i)
  merged_r = merge_samples(TP, "SampleType")
  merged_r_m = transform_sample_counts(merged_r , function(x) 100 * x/sum(x))
  physeqdf <- psmelt(merged_r_m)
  temp <- c("Soil","Rhizosphere","Root")
  physeqdf$Sample <- factor(physeqdf$Sample, levels = c(temp))
  n = nrow(physeqdf[order(physeqdf$Phylum),])
  x = sample(i,n,replace=TRUE)
  ST <- sample("Pre_flowering",n,replace=TRUE)
  physeqdf[order(physeqdf$Phylum),]$Timepoint <- x 
  physeqdf[order(physeqdf$Phylum),]$Treatment <- ST
  physeqdf$Phylum <- as.character(physeqdf$Phylum)
  pre2 <- rbind(pre2,physeqdf[order(physeqdf$Phylum),])
}

add2 <- subset(control,Treatment=="Control")
add3 <- subset(add2,Timepoint=="TP1"|Timepoint=="TP2"|Timepoint=="TP3"|Timepoint=="TP4"|Timepoint=="TP5"|Timepoint=="TP6"|Timepoint=="TP7")
add3$Treatment[add3$Treatment == "Control"] <- "Post_flowering"

p <- {}
y <- c("TP8","TP9","TP10","TP11","TP12","TP13","TP14","TP15","TP16","TP17")
post = subset_samples(b1_phylum,Treatment=="Post_flowering")
for(i in y){
  TP = subset_samples(post, Timepoint == i)
  merged_r = merge_samples(TP, "SampleType")
  merged_r_m = transform_sample_counts(merged_r , function(x) 100 * x/sum(x))
  physeqdf <- psmelt(merged_r_m)
  temp <- c("Soil","Rhizosphere","Root")
  physeqdf$Sample <- factor(physeqdf$Sample, levels = c(temp))
  n = nrow(physeqdf[order(physeqdf$Phylum),])
  x = sample(i,n,replace=TRUE)
  ST <- sample("Post_flowering",n,replace=TRUE)
  physeqdf[order(physeqdf$Phylum),]$Timepoint <- x 
  physeqdf[order(physeqdf$Phylum),]$Treatment <- ST
  physeqdf$Phylum <- as.character(physeqdf$Phylum)
  p <- rbind(p,physeqdf[order(physeqdf$Phylum),])
}
rm <- rbind(control,pre1,pre2,add3,p)

sum <- tapply(rm$Abundance,rm$Phylum, FUN=sum)
sum <- data.frame(sum)
sum <- sort(sum$sum,decreasing = T)
list <- c("Acidobacteria","Actinobacteria","Bacteroidetes","Chloroflexi","Cyanobacteria",
          "Firmicutes","Gemmatimonadetes","Nitrospirae",
          "Planctomycetes","Proteobacteria","Tenericutes","TM7","Verrucomicrobia")
a <- rm[rm$Phylum %in% list,]
b <- rm[!(rm$Phylum %in% list),]
b$Phylum <- "Other"
rm <- rbind(a,b)
rm$Treatment<-factor(rm$Treatment,levels=c("Control","Pre_flowering","Post_flowering"))
rm$Timepoint<-factor(rm$Timepoint, levels=c("TP1","TP2","TP3","TP4","TP5","TP6","TP7","TP8","TP9","TP10","TP11","TP12","TP13","TP14","TP15","TP16","TP17"))

col_SampleByTreatment<-c("Acidobacteria"="thistle","Actinobacteria"="lightsalmon2","Armatimonadetes"="mediumpurple1","Bacteroidetes"="cornsilk","Chloroflexi"="lightskyblue","Cyanobacteria"="lightpink","Firmicutes"="lightblue1","Gemmatimonadetes"="darkseagreen","Nitrospirae"="mistyrose2",
                         "Planctomycetes"="gray69","Proteobacteria"="cadetblue3","SPAM"="darkorchid1","Tenericutes"="darkolivegreen2","TM7"="yellow","Verrucomicrobia"="slategray1","Other"="plum")


rm$Phylum <- factor(rm$Phylum,levels=c("Actinobacteria","Chloroflexi","Cyanobacteria","Firmicutes","Bacteroidetes","Proteobacteria",
                                       "Acidobacteria","Armatimonadetes","Gemmatimonadetes","Nitrospirae",
                                       "Planctomycetes","SPAM","Tenericutes","TM7","Verrucomicrobia","Other" ))

sub<-subset(rm,Treatment!="Pre_flowering")

ggplot(sub, aes(x=Timepoint, y=Abundance, fill=Phylum)) + 
  geom_bar(stat = "identity", color = "NA")+facet_wrap(~Sample+Treatment,ncol=2)+ scale_fill_manual(values = col_SampleByTreatment) +
  theme(axis.text.x=element_text(size=11,color="black",angle=90), axis.text.y=element_text(size=11,color="black"), 
        axis.title=element_text(size=14,face="bold"),text=element_text(size=14,face="bold"))+
  ylab("Relative Abundance")+xlab("")


##### Fig.S20 PCoA of Bray Curtis distance for root samples under control and pre-flowering drought

r = subset_samples(rar,SampleType == "Root")
r = subset_samples(r, Treatment != "Pre_flowering")
plot_ordination(r, ordinate(r, "MDS",distance="bray"), color = "Timepoint") + 
  scale_colour_manual(values=c("gray8","gray37",brewer.pal(9,"YlGn")[c(9,8,7,6,4,3)],brewer.pal(9,"YlGnBu")[c(3,4,6,7,8)],brewer.pal(9,"RdPu")[c(9,8,7,6)]))+ 
  geom_point(size = 3)+ facet_wrap(~Treatment) +
  theme(axis.text.x=element_text(size=11,color="black",angle=90), axis.text.y=element_text(size=11,color="black"), 
        axis.title=element_text(size=11,face="bold"),text=element_text(size=11,face="bold")) 

##### Fig.S21 Heat map of the Bray Curtis dissimilarity between treatments (x-axis) at different time points

#Control pre_flowering_drought
root <- subset_samples(rar,SampleType == "Root")
root_c_pre <- subset_samples(root,Treatment!="Post_flowering")
x <- c("TP3","TP4","TP5","TP6","TP7","TP8")
rm <- {}
for(i in x){
  TP <- subset_samples(root_c_pre,Timepoint==i)
  me <- merge_samples(TP,"Treatment")
  bray <- phyloseq::distance(me, "bray")
  df <- as.matrix(bray)
  nba.m <- melt(df)
  sub <- subset(nba.m,value!=0)
  sub <- sub[1,]
  colnames(sub) <- c("Sample1","Sample2","Value")
  sub$Timepoint <- i
  sub$Treatment <- "C-pre"
  rm <- rbind(rm,sub)
}

#Control post_flowering_drought
root <- subset_samples(rar,SampleType == "Root")
root_c_pre <- subset_samples(root,Treatment!="Pre_flowering")
x <- c("TP10","TP11","TP12","TP13","TP14","TP15","TP16","TP17")
for(i in x){
  TP <- subset_samples(root_c_pre,Timepoint==i)
  me <- merge_samples(TP,"Treatment")
  bray <- phyloseq::distance(me, "bray")
  df <- as.matrix(bray)
  nba.m <- melt(df)
  sub <- subset(nba.m,value!=0)
  sub <- sub[1,]
  colnames(sub) <- c("Sample1","Sample2","Value")
  sub$Timepoint <- i
  sub$Treatment <- "C-post"
  rm <- rbind(rm,sub)
}

rm$Timepoint <- factor(rm$Timepoint, levels=c("TP1","TP2","TP3","TP4","TP5","TP6","TP7","TP8","TP9","TP10","TP11","TP12","TP13","TP14","TP15","TP16","TP17"))

ggplot(data = rm, aes(Treatment,Timepoint)) + 
  geom_tile(aes(fill = Value), colour = "white") +  theme(panel.border=element_blank())+
  facet_wrap(~Treatment, scales="free")+
  scale_fill_gradientn(name="Unifrac distance",colours = terrain.colors(7))+
  theme(axis.text.x=element_text(size=10,color="black",angle=90,face="bold"), 
        axis.text.y=element_text(size=10,color="black",face="bold"), 
        axis.title=element_text(size=10,face="bold"),text=element_text(size=10))+xlab("")+ylab("")+
  theme(legend.title = element_text(colour="black", size=11, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 11, face = "bold"))+
  theme(strip.background = element_blank(),
        strip.text.x = element_blank())

######Fig.S22a Bar plot of phylum-level relative abundance for the differential expressed OTUs

d <- read.table("figS22a-pre_post_phylum_up_down.txt",header=T,sep="\t")
rm <- data.frame(d)
list <- c("Acidobacteria","Actinobacteria","Bacteroidetes","Chloroflexi","Cyanobacteria",
          "Firmicutes","Gemmatimonadetes","Nitrospirae",
          "Planctomycetes","Proteobacteria","Tenericutes","TM7","Verrucomicrobia")
a <- rm[rm$Phylum %in% list,]
b <- rm[!(rm$Phylum %in% list),]
b$Phylum <- "Other"
rm <- rbind(a,b)

col_SampleByTreatment<-c("Acidobacteria"="thistle","Actinobacteria"="lightsalmon2","Armatimonadetes"="mediumpurple1","Bacteroidetes"="cornsilk","Chloroflexi"="lightskyblue","Cyanobacteria"="lightpink","Firmicutes"="lightblue1","Gemmatimonadetes"="darkseagreen","Nitrospirae"="mistyrose2",
                         "Planctomycetes"="gray69","Proteobacteria"="cadetblue3","SPAM"="darkorchid1","Tenericutes"="darkolivegreen2","TM7"="yellow","Verrucomicrobia"="slategray1","Other"="plum")


rm$Phylum <- factor(rm$Phylum,levels=c("Actinobacteria","Chloroflexi","Cyanobacteria","Firmicutes","Bacteroidetes","Proteobacteria",
                                       "Acidobacteria","Armatimonadetes","Gemmatimonadetes","Nitrospirae",
                                       "Planctomycetes","SPAM","Tenericutes","TM7","Verrucomicrobia","Other" ))

ggplot(rm, aes(x=Phylum, y=otu_number, fill=Phylum, )) + 
  geom_bar(stat = "identity", color = "NA")+facet_wrap(~Treatment,ncol=2)+
  scale_fill_manual(values = col_SampleByTreatment) +
  theme(axis.text.x=element_text(size=11,color="black"), 
        axis.text.y=element_text(size=11,color="black"), 
        axis.title=element_text(size=14,face="bold"),text=element_text(size=14,face="bold"))+
  ylab("OTU number")+xlab("") + coord_flip()

# Fig.S22b Venn Diagram of enrichment overlap in Pre and Post flowering drought

library(VennDiagram)
grid.newpage()
draw.pairwise.venn(1059,1029,491, category = c("C-Post", "C-Pre"), 
                   lty = rep("blank",2), fill = c("#CCCC99","#666666"), 
                   alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.025, 2))

# Fig.S22c Log Fold Change in OTU Abundance in Pre and Post flowering drought colored by Phylum and Monoderm versus Diderm status

PrePostAbund<-read.table(file = "fig22c.txt",header = T)
library(ggplot2)

ggplot(data=PrePostAbund,aes(x=Pre,y=Post,color=Phylum_B))+
  geom_point(aes(size=OTUAbund))+facet_grid(~PreOrPostSign)+
  scale_y_continuous(limits=c(-10,10))+
  scale_x_continuous(limits=c(-10,10))+
  theme(axis.text.x=element_text(size=11,color="black",angle=90),
        axis.text.y=element_text(size=11,color="black"),
        axis.title=element_text(size=14,face="bold"),
        text=element_text(size=14,face="bold")) 

##### Fig.S23 Venn diagram indicating the number of differential expressed genes unique to and common between control and pre-flowering drought for rhizosphere and soil sample types at TP8 and TP9

grid.newpage()
draw.pairwise.venn(4685,13300,1457, category = c("Soil_tp8", "Rhi_tp8"), 
                   lty = rep("blank",2), fill = c("#CCCC99","#666666"), 
                   alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.025, 2))
grid.newpage()
draw.pairwise.venn(4619,5967,594, category = c("Soil_tp9", "Rhi_tp9"), 
                   lty = rep("blank",2), fill = c("#CCCC99","#666666"), 
                   alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.025, 2))


##### Fig.S24 Relative abundance for core gene set in metatrans data

biom_file <- "figS24_bacteria_count_core.biom"
map_file <- "figS24_metadata_metatrans.txt"

set.seed(1)
biomot=import_biom(biom_file,parseFunction = parse_taxonomy_greengenes)
bmsd = import_qiime_sample_data(map_file)
Bushman = merge_phyloseq(biomot,bmsd)
sample_data(Bushman)$SampleType<-factor(sample_data(Bushman)$SampleType, levels=c("Soil","Rhizosphere"))
sample_data(Bushman)$Treatment<-factor(sample_data(Bushman)$Treatment, levels=c("Control","Pre_flowering","Post_flowering"))
sample_data(Bushman)$Timepoint<-factor(sample_data(Bushman)$Timepoint, levels=c("TP8","TP9"))
Bushman = subset_taxa(Bushman, Phylum != "Streptophyta")

#Rhizosphere&soil_phylum_Bacteria
rhi = subset_samples(Bushman, SampleType=="Rhizosphere")
Phylum <- tax_glom(rhi, taxrank="Phylum")

x <- c("TP8","TP9")
rm <- {}
for(i in x){
  TP = subset_samples(Phylum, Timepoint == i)
  merged_r = merge_samples(TP, "Treatment")
  merged_r_m = transform_sample_counts(merged_r , function(x) 100 * x/sum(x))
  physeqdf <- psmelt(merged_r_m)
  temp <- c("Control","Pre_flowering","Post_flowering")
  physeqdf$Sample <- factor(physeqdf$Sample, levels = c(temp))
  n = nrow(physeqdf[order(physeqdf$Phylum),])
  x = sample(i,n,replace=TRUE)
  ST <- sample("Rhizosphere",n,replace=TRUE)
  physeqdf[order(physeqdf$Phylum),]$Timepoint <- x 
  physeqdf[order(physeqdf$Phylum),]$SampleType <- ST
  rm <- rbind(rm,physeqdf[order(physeqdf$Phylum),])
}

Soil = subset_samples(Bushman, SampleType=="Soil")
Phylum <- tax_glom(Soil, taxrank="Phylum")

x <- c("TP8","TP9")
rm2 <- {}
for(i in x){
  TP = subset_samples(Phylum, Timepoint == i)
  merged_r = merge_samples(TP, "Treatment")
  merged_r_m = transform_sample_counts(merged_r , function(x) 100 * x/sum(x))
  physeqdf <- psmelt(merged_r_m)
  temp <- c("Control","Pre_flowering","Post_flowering")
  physeqdf$Sample <- factor(physeqdf$Sample, levels = c(temp))
  n = nrow(physeqdf[order(physeqdf$Phylum),])
  x = sample(i,n,replace=TRUE)
  ST <- sample("Soil",n,replace=TRUE)
  physeqdf[order(physeqdf$Phylum),]$Timepoint <- x 
  physeqdf[order(physeqdf$Phylum),]$SampleType <- ST
  rm2 <- rbind(rm2,physeqdf[order(physeqdf$Phylum),])
}

rm3 <- rbind(rm,rm2)

sum <- tapply(rm3$Abundance,rm3$Phylum, FUN=sum)
sum <- data.frame(sum)
sum <- sort(sum$sum,decreasing = T)
list <- c("Acidobacteria","Actinobacteria","Bacteroidetes","Chloroflexi","Cyanobacteria",
          "Firmicutes","Gemmatimonadetes","Nitrospirae",
          "Planctomycetes","Proteobacteria","Tenericutes","TM7","Verrucomicrobia")

a <- rm3[rm3$Phylum %in% list,]
b <- rm3[!(rm3$Phylum %in% list),]
b$Phylum <- "Other"
rm <- rbind(a,b)
rm$Treatment<-factor(rm$Treatment,levels=c("Control","Pre_flowering","Post_flowering"))
rm$Timepoint<-factor(rm$Timepoint, levels=c("TP1","TP2","TP3","TP4","TP5","TP6","TP7","TP8","TP9","TP10","TP11","TP12","TP13","TP14","TP15","TP16","TP17"))

col_SampleByTreatment<-c("Acidobacteria"="thistle","Actinobacteria"="lightsalmon2","Armatimonadetes"="mediumpurple1","Bacteroidetes"="cornsilk","Chloroflexi"="lightskyblue","Cyanobacteria"="lightpink","Firmicutes"="lightblue1","Gemmatimonadetes"="darkseagreen","Nitrospirae"="mistyrose2",
                         "Planctomycetes"="gray69","Proteobacteria"="cadetblue3","SPAM"="darkorchid1","Tenericutes"="darkolivegreen2","TM7"="yellow","Verrucomicrobia"="slategray1","Other"="plum")

rm$Phylum <- factor(rm$Phylum,levels=c("Actinobacteria","Chloroflexi","Cyanobacteria","Firmicutes","Bacteroidetes","Proteobacteria",
                                       "Acidobacteria","Armatimonadetes","Gemmatimonadetes","Nitrospirae",
                                       "Planctomycetes","SPAM","Tenericutes","TM7","Verrucomicrobia","Other" ))

ggplot(rm, aes(x=Sample, y=Abundance, fill=Phylum, )) + 
  geom_bar(stat = "identity", color = "NA")+facet_wrap(~SampleType+Timepoint,ncol=4)+ 
  scale_fill_manual(values = col_SampleByTreatment) +
  theme(axis.text.x=element_text(size=11,color="black",angle=90), 
        axis.text.y=element_text(size=11,color="black"), 
        axis.title=element_text(size=11,face="bold"),
        text=element_text(size=11,face="bold"))+
  ylab("Relative Abundance")

##### Fig.S25 Gene Ontology (GO) enrichment analysis for all differentially expressed genes between drought and control treatments for both rhizosphere (left panel) and soils (right panels) at TP9

#TP8_input
data<-read.table("fig4b-1-new_heatmap_enrichment_TP8.txt",header=T,sep="\t")

tp8_rhi <- data[,c(1,3,4)]
ratio <- data.frame(cbind(tp8_rhi,tp8_rhi$Rhi/tp8_rhi$total))
dat <- ratio[,c(1,4)]
colnames(dat) <-c("Annotation","TP8_rhizosphere")
tp8_soil <- data[,c(1,2,4)]
ratio <- data.frame(tp8_soil$Soil/tp8_soil$total)
colnames(ratio) <- c("TP8_soil")
mer8 <- data.frame(cbind(dat,ratio))

#tp9 input
data<-read.table("fig4b-2-heatmap_enrich_tp9.txt",header=T,sep="\t")
tp8_rhi <- data[,c(1,3,4)]
ratio <- data.frame(cbind(tp8_rhi,tp8_rhi$rhi/tp8_rhi$total))
dat <- ratio[,c(1,4)]
colnames(dat) <- c("Annotation","TP9_rhizosphere")
tp8_soil <- data[,c(1,2,4)]
ratio <- data.frame(tp8_soil$soil/tp8_soil$total)
colnames(ratio) <- c("TP9_soil")
mer9 <- data.frame(cbind(dat,ratio))
rm <- merge(mer8,mer9,by="Annotation")
rm.m <- melt(rm)
rm.m$variable <- factor(rm.m$variable, levels=c("TP8_soil","TP9_soil",
                                                "TP8_rhizosphere","TP9_rhizosphere"))
dat <- rm.m[order(-rm.m[,3]),]
dat$Annotation <- factor(dat$Annotation, levels=c("Carbohydrate transport and metabolism",
                                                  "Amino acid transport and metabolism",
                                                  "Secondary metabolites biosynthesis, transport and catabolism",
                                                  "Cell cycle control, cell division, chromosome partitioning",
                                                  "Mobilome: prophages, transposons",
                                                  "Posttranslational modification, protein turnover, chaperones",
                                                  "Nucleotide transport and metabolism",
                                                  "Transcription",
                                                  "Replication, recombination and repair",
                                                  "Translation, ribosomal structure and biogenesis",
                                                  "Inorganic ion transport and metabolism",
                                                  "Energy production and conversion",
                                                  "Coenzyme transport and metabolism",
                                                  "Lipid transport and metabolism",
                                                  "Cell motility",
                                                  "Cell wall/membrane/envelope biogenesis",
                                                  "General function prediction only",
                                                  "Intracellular trafficking, secretion, and vesicular transport",
                                                  "Function unknown",
                                                  "Signal transduction mechanisms",
                                                  "Defense mechanisms",
                                                  "Cytoskeleton",
                                                  "Extracellular structures"))

samples_new <- sapply(as.character(dat$variable), function(x) strsplit(x, "[_]"))

Timepoint <- {}
for(i in 1:84){
  T <- samples_new[[i]][1]
  Timepoint <- append(Timepoint,T)
}

SampleType <- {}
for(i in 1:84){
  T <- samples_new[[i]][2]
  SampleType <- append(SampleType,T)
}

dat$SampleType = SampleType
dat$Timepoint = Timepoint
dat <- subset(dat,Timepoint=="TP9")

ggplot(dat, aes(x=value, y=Annotation,colour=Timepoint)) +
  scale_colour_manual(name="",values = c("#004445","#6FB98F","#9B4F0F","#C99E10"))+
  scale_shape_manual(name="",values=c(1,1,19,19))+facet_wrap(~SampleType,scales="free_x",ncol=2)+
  geom_point(size = 4,stroke = 1)+theme_bw()+
  geom_vline(aes(xintercept=1),colour="Black",size=1,linetype="dashed")+xlab("Enrichment ratio")+
  theme(axis.text.x=element_text(size=11,color="black",face="bold",angle=90), 
        axis.text.y=element_text(size=11,color="black",face="bold"), 
        axis.title=element_text(size=11,face="bold"),text=element_text(size=10))+
  theme(legend.text = element_text(colour="black", size = 11, face = "bold"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) +
  theme(strip.text = element_text(size=11,face="bold"))


##### Fig.S26 Barplot for G3P concentration for control and drought

data <- read.table("figS26-G3P_input.txt",header=T)
data$Genotype <- factor(data$Genotype,levels=c("RTx430","BTx642"))

ggplot(data, aes(x=Genotype, y=Mean, fill=Treatment)) + 
  geom_bar(stat = "identity", color = "NA",position = "dodge")+ 
  scale_fill_manual(values = c("#6AB187", "#DE7A22")) +
  theme(axis.text.x=element_text(size=11,color="black",angle=90,hjust=1), axis.text.y=element_text(size=11,color="black"), 
        axis.title=element_text(size=14,face="bold"),text=element_text(size=14,face="bold"))+
  ylab("umol/g")+xlab("") + geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD),width=.2, 
                                          position=position_dodge(.9),color="grey50")


##### Fig.S28 qPCR for taxa-specific primers (microbox)

data <- read.table("figS28_actino1_microbox_qPCR.txt",header=T)
data$Taxa <- factor(data$Taxa,levels=c("Sc1","Sc2","K"))

ggplot(data, aes(x=Taxa, y=Value, fill=Treatment)) + 
  geom_bar(stat = "identity", color = "NA",position = "dodge")+ 
  scale_fill_manual(values = c("#6AB187", "#DE7A22")) + facet_wrap(~Genotype+Normalization,ncol=6)+
  theme(axis.text.x=element_text(size=11,color="black",angle=90,hjust=1), axis.text.y=element_text(size=11,color="black"), 
        axis.title=element_text(size=14,face="bold"),text=element_text(size=14,face="bold"))+
  ylab("Relative Abundance")+xlab("") + geom_errorbar(aes(ymin=Value-SD, ymax=Value+SD),width=.2, 
                                                      position=position_dodge(.9),color="grey50")


##### Fig.S29 violin plots of the fluorescence intensity as measured with Sc2 using confocal fluorescence microscopy across control (green) and drought-treated (orange) root samples

rm(list = ls())
data<-read.table("fig5-Confocal.txt",header=T,sep="\t")
sai80 = subset(data,Strain=="SAI80")
sai80$Treatment <- as.factor(sai80$Treatment)

ggplot(sai80, aes(x=Treatment, y=Value, fill=Treatment)) +  
  geom_violin(trim = TRUE)+
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="pointrange", color="white")+
  scale_fill_manual(values=c("#6AB187", "#DE7A22"))

##### Fig.S30 violin plots of the fluorescence intensity as measured with Sc1 using confocal fluorescence microscopy across control (green) and PEG-treated (Orange) root samples.

rm(list = ls())
data<-read.table("fig5-Confocal.txt",header=T,sep="\t")
sai73 = subset(data,Strain=="SAI73")
sai73$Treatment <- as.factor(sai73$Treatment)
sai73_cd <- subset(sai73,Treatment != "Drought")

ggplot(sai73_cd, aes(x=Treatment, y=Value, fill=Treatment)) +  
  geom_violin(trim = TRUE)+
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="pointrange", color="white") +
  scale_fill_manual(values=c("#6AB187", "#DE7A22"))

##### Fig.S31 Violin plots of tissue fresh weight for roots (left panels) and shoots (right panels) under control (top) or drought (bottom) treatment. 

data<-read.table("figS31-fresh_weight.txt",header=T,sep="\t")
data$Treatment2 <- as.factor(data$Treatment2)


ggplot(data, aes(x=Strain, y=fresh_weight, fill=Strain)) +  
  geom_violin(trim = TRUE)+ facet_wrap(~Treatment2+SampleType,scales="free_y")+
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="pointrange", color="white") + xlab("")+ylab("Root fresh weight (mg)")+
  theme(axis.text.x=element_text(size=11,color="black",angle=90,vjust=0.5,hjust=0.5), 
        axis.text.y=element_text(size=11,color="black"), 
        axis.title=element_text(size=11,face="bold"),
        text=element_text(size=11,face="bold")) +
  scale_fill_manual(values=c("#DDDEDE", "#A5C05B", "#7BA4A8"))

##### Fig.S32
#Fig.S32a PCoA of Bray Curtis distance 
set.seed(1)
DNA <- readRDS("DNA.rds")
root = subset_taxa(DNA,Class != "Chloroplast")
#rar_s = rarefy_even_depth(root)
rar_s <- readRDS("rar_s.rds")

plot_ordination(rar_s, ordinate(rar_s, "MDS"), shape = "Treatment" ,color = "Replicate") + 
  geom_point(size = 4)+facet_wrap(~Timepoint) +
  scale_colour_manual(name="",values = c("#F77604","#B8D20B","#7F152E"))+
  scale_shape_manual(name="",values=c(18,19))+
  theme(axis.text.x=element_text(size=11,color="black",face="bold"), 
        axis.text.y=element_text(size=11,color="black",face="bold"), 
        axis.title=element_text(size=11,face="bold"),text=element_text(size=10))+
  theme(legend.text = element_text(colour="black", size = 11, face = "bold")) + 
  theme(strip.text = element_text(size=11,face="bold"))

#Fig.S32b boxplots of Shannon’s diversity

plot_richness(rar_s,x="Timepoint",color="Replicate",shape="Treatment",
              measures = "Shannon") + scale_color_manual(values=c("#F77604","#B8D20B","#7F152E"))+
  theme(legend.title = element_text(colour="black", size=11, face="bold"))+geom_point(size = 4)+
  theme(legend.text = element_text(colour="black", size = 11, face = "bold"))+
  theme(axis.text.x=element_text(hjust=1,vjust=0.5,size=11,
                                 color="black",angle=90,face="bold"), 
        axis.text.y=element_text(size=11,color="black",face="bold"), 
        axis.title=element_text(size=11,face="bold"),
        text=element_text(size=11,face="bold"))+
  scale_shape_manual(name="",values=c(18,19))+
  xlab("") + ylab("Shannon")+theme(plot.title = element_blank())

#Fig.S32c bar plots of phylum-level relative abundance 

merged_r_m = transform_sample_counts(root, function(x) 100 * x/sum(x))
physeqdf <- psmelt(merged_r_m)
n = nrow(physeqdf[order(physeqdf$Phylum),])
physeqdf$Phylum <- as.character(physeqdf$Phylum)
rm <- physeqdf
sum <- tapply(rm$Abundance,rm$Phylum, FUN=sum)
sum <- data.frame(sum)
sum <- sort(sum$sum,decreasing = T)
list <- c("Acidobacteria","Actinobacteria","Bacteroidetes","Chloroflexi","Cyanobacteria",
          "Firmicutes","Gemmatimonadetes","Nitrospirae",
          "Planctomycetes","Proteobacteria","Tenericutes","TM7","Verrucomicrobia")
a <- rm[rm$Phylum %in% list,]
b <- rm[!(rm$Phylum %in% list),]
b$Phylum <- "Other"
rm <- rbind(a,b)
rm$value <- rm$Abundance/2

col_SampleByTreatment<-c("Acidobacteria"="thistle","Actinobacteria"="lightsalmon2","Armatimonadetes"="mediumpurple1","Bacteroidetes"="cornsilk","Chloroflexi"="lightskyblue","Cyanobacteria"="lightpink","Firmicutes"="lightblue1","Gemmatimonadetes"="darkseagreen","Nitrospirae"="mistyrose2",
                         "Planctomycetes"="gray69","Proteobacteria"="cadetblue3","SPAM"="darkorchid1","Tenericutes"="darkolivegreen2","TM7"="yellow","Verrucomicrobia"="slategray1","Other"="plum")


rm$Phylum <- factor(rm$Phylum,levels=c("Actinobacteria","Chloroflexi","Cyanobacteria","Firmicutes","Bacteroidetes","Proteobacteria",
                                       "Acidobacteria","Armatimonadetes","Gemmatimonadetes","Nitrospirae",
                                       "Planctomycetes","SPAM","Tenericutes","TM7","Verrucomicrobia","Other" ))


ggplot(rm, aes(x=Replicate, y=value, fill=Phylum, )) + 
  geom_bar(stat = "identity", color = "NA")+  facet_grid(~Treatment+Timepoint)+
  theme(axis.text.x=element_text(size=11,color="black",angle=90), axis.text.y=element_text(size=11,color="black"), 
        axis.title=element_text(size=11,face="bold"),text=element_text(size=11,face="bold")) +
  scale_fill_manual(values = col_SampleByTreatment) + ylab("Relative Abundance") + xlab("")

##### End 