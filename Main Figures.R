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
library(VennDiagram)

##### Fig. 1a-c Shannon diversity for three treatments

richness_2 = estimate_richness(rar,measure="Shannon")
s <- data.frame(sample_data(rar))
alphadiv <- cbind(richness_2, s)
alphadiv <- data.frame(alphadiv)
alphadiv$Timepoint <- factor(alphadiv$Timepoint,levels=c("TP1","TP2","TP3","TP4","TP5","TP6","TP7","TP8","TP9","TP10","TP11","TP12","TP13","TP14","TP15","TP16","TP17"))
alphadiv$Treatment <- factor(alphadiv$Treatment,levels=c("Control","Pre_flowering","Post_flowering"))
alphadiv$SampleType <- factor(alphadiv$SampleType,levels=c("Soil","Rhizosphere","Root"))
tp1_2 <- subset(alphadiv,Timepoint=="TP1"|Timepoint=="TP2")
tp1_2$Treatment <- "Pre_flowering"
tp1_2$TreatmentByTimepoint <- paste(tp1_2$Treatment,tp1_2$Timepoint,sep="")
tp1_7 <- subset(alphadiv,Treatment == "Control")
tp1_7_2 <- subset(tp1_7,Timepoint=="TP1"|Timepoint=="TP2"|Timepoint=="TP3"|Timepoint=="TP4"|Timepoint=="TP5"|Timepoint=="TP6"|Timepoint=="TP7")
tp1_7_2$Treatment <- "Post_flowering"
tp1_7_2$TreatmentByTimepoint <- paste(tp1_7_2$Treatment,tp1_7_2$Timepoint,sep="")

data <- rbind(alphadiv,tp1_2,tp1_7_2)
input<-aggregate(Shannon ~ Treatment + Timepoint + SampleType, data, mean)
std <- aggregate(Shannon ~ Treatment + Timepoint + SampleType, data = data, FUN= sd)
input<-cbind(input,std$Shannon)
colnames(input)[5]<-"SD"
TreatmentSampleType<-with(input,interaction(Treatment,SampleType))
input<-cbind(input,TreatmentSampleType)
input$upper<-input$Shannon+input$SD
input$lower<-input$Shannon-input$SD
input<-data.table(input)
input[SampleType == "Soil",y_min := 7.5]
input[SampleType == "Soil",y_max := 5]
input[SampleType == "Rhizosphere",y_min := 4]
input[SampleType == "Rhizosphere",y_max := 6.5]
input[SampleType == "Root",y_min := 3]
input[SampleType == "Root",y_max := 5.5]

ggplot(input, aes(x = Timepoint, y = Shannon, color = Treatment, 
                       group = Treatment, shape = Treatment)) +
  geom_line(size = 0.8) + xlab("")+
  facet_wrap(~SampleType, scales = "free_y",ncol = 1) +
  scale_color_manual(values = c("##### 6AB187", "#DE7A22", "#F4CC70")) +
  scale_fill_manual(values = c("lightblue3","#FF8C00","#CDBE70"))+
  theme(legend.title = element_text(colour="black", size=11, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 11, face = "bold"))+
  theme(axis.text.x=element_text(hjust=1,vjust=0.5,size=10,color="black",angle=90,face="bold"), 
        axis.text.y=element_text(size=11,color="black",face="bold"), 
        axis.title=element_text(size=11,face="bold"),text=element_text(size=11,face="bold"))+
  geom_blank(aes(y = y_min)) +
  geom_blank(aes(y = y_max))+
  geom_ribbon(aes(ymin=lower,ymax=upper,fill=Treatment),alpha=0.4,colour=NA)+   
  scale_y_continuous(breaks=seq(3, 7.5, 0.5)) # make all y axis show 0.5 space


##### Fig. 1d PCoA of Bray Curtis distance for all samples

map_file <- "Meta_final_Figure1d.txt"
bmsd = import_qiime_sample_data(map_file)
rar_pre<-subset_samples(rar,Treatment!="Post_flowering")
sample_data(rar_pre) <- bmsd
sample_data(rar_pre)$Timepoint<-factor(sample_data(rar_pre)$Timepoint, levels=c("TP1", "TP2","TP3","TP4","TP5","TP6","TP7","TP8","TP9","TP10","TP11",
                                                                                "TP12","TP13","TP14","TP15","TP16","TP17"))
plot_ordination(rar_pre, ordinate(rar_pre, "MDS",distance="bray"),axes=1:2, color = "TreatmentByTimepoint5") + 
  scale_colour_manual(values=c(brewer.pal(9,"YlOrRd")[c(4,6,8)],brewer.pal(9,"YlGn")[c(4,6,8)],rep(brewer.pal(9,"YlGn")[c(8)],times=11),rep(brewer.pal(9,"YlOrRd")[c(8)],times=6),brewer.pal(9,"YlGn")[c(4,6,8)],rep(brewer.pal(9,"YlGn")[c(8)],times=6))) + #"green","green","green","green"brewer.pal(6,"YlOrBr")[c(4)],brewer.pal(11,"PiYG")[c(11)],brewer.pal(11,"PiYG")[c(8)],"black","grey")) + # scale_shape_manual(values=c(15,16,17)) +
  geom_point(size = 3,aes(shape=SampleType))+ #facet_wrap(~Treatment) +
  scale_shape_manual(values=c(15,16,17)) +
  theme(axis.text.x=element_text(size=11,color="black",angle=90), axis.text.y=element_text(size=11,color="black"), 
        axis.title=element_text(size=11,face="bold"),text=element_text(size=11,face="bold"))

##### Fig. 1e PCoA of Bray Curtis distance for root samples

r = subset_samples(rar_pre,SampleType == "Root")
plot_ordination(r, ordinate(r, "MDS",distance="bray"), color = "Timepoint") + 
  scale_colour_manual(values=c("gray8","gray37",brewer.pal(9,"YlGn")[c(9,8,7,6,4,3)],brewer.pal(9,"YlGnBu")[c(3,4,6,7,8)],brewer.pal(9,"RdPu")[c(9,8,7,6)]))+ 
  geom_point(size = 3)+ facet_wrap(~Treatment) +
  theme(axis.text.x=element_text(size=11,color="black",angle=90), axis.text.y=element_text(size=11,color="black"), 
        axis.title=element_text(size=11,face="bold"),text=element_text(size=11,face="bold")) 

##### Fig. 1f-h Heat map of the mean pairwise Bray Curtis dissimiliarity between all root sample replicates within the specified pair of treatments

data<-read.table("fig1f-h_heatmap_similarity_input.txt",header=T,sep="\t")
data$Var1<-factor(data$Var1, levels=c("ControlTP1",	"ControlTP2",	"ControlTP3",	"ControlTP4",	"ControlTP5",	"ControlTP6",	"ControlTP7",	"ControlTP8",	"ControlTP9",	"ControlTP10",	"ControlTP11",	"ControlTP12",	"ControlTP13",	"ControlTP14",	"ControlTP15",	"ControlTP16",	"ControlTP17",	"Pre_floweringTP3",	"Pre_floweringTP4",	"Pre_floweringTP5",	"Pre_floweringTP6",	"Pre_floweringTP7",	"Pre_floweringTP8",	"Pre_floweringTP9",	"Pre_floweringTP10",	"Pre_floweringTP11",	"Pre_floweringTP12",	"Pre_floweringTP13",	"Pre_floweringTP14",	"Pre_floweringTP15",	"Pre_floweringTP16",	"Pre_floweringTP17"))
data$Var2<-factor(data$Var2, levels=c("ControlTP1",	"ControlTP2",	"ControlTP3",	"ControlTP4",	"ControlTP5",	"ControlTP6",	"ControlTP7",	"ControlTP8",	"ControlTP9",	"ControlTP10",	"ControlTP11",	"ControlTP12",	"ControlTP13",	"ControlTP14",	"ControlTP15",	"ControlTP16",	"ControlTP17",	"Pre_floweringTP3",	"Pre_floweringTP4",	"Pre_floweringTP5",	"Pre_floweringTP6",	"Pre_floweringTP7",	"Pre_floweringTP8",	"Pre_floweringTP9",	"Pre_floweringTP10",	"Pre_floweringTP11",	"Pre_floweringTP12",	"Pre_floweringTP13",	"Pre_floweringTP14",	"Pre_floweringTP15",	"Pre_floweringTP16",	"Pre_floweringTP17"))
sub <- subset(data,Var1 != "ControlTP1")
sub <- subset(sub,Var1 != "ControlTP2")
sub <- subset(sub,Var2 != "ControlTP2")
sub <- subset(sub,Var2 != "ControlTP1")

ggplot(data = sub, aes(Var1,Var2)) + 
  geom_tile(aes(fill = value), colour = "white") + facet_wrap(~Variation, scales="free") + theme(panel.border=element_blank())+
  scale_fill_gradientn(name="Bray curtis distance",colours = terrain.colors(7))+
  theme(axis.text.x=element_text(size=10,color="black",angle=90,face="bold"), 
        axis.text.y=element_text(size=10,color="black",face="bold"), 
        axis.title=element_text(size=10,face="bold"),text=element_text(size=10))+xlab("")+ylab("")+
  theme(legend.title = element_text(colour="black", size=11, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 11, face = "bold"))+
  theme(strip.background = element_blank(),
        strip.text.x = element_blank())

##### Fig.2 Relative abundance for different phyla.

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

sub<-subset(rm,Treatment!="Post_flowering")

ggplot(sub, aes(x=Timepoint, y=Abundance, fill=Phylum)) + 
  geom_bar(stat = "identity", color = "NA")+facet_wrap(~Sample+Treatment,ncol=2)+ scale_fill_manual(values = col_SampleByTreatment) +
  theme(axis.text.x=element_text(size=11,color="black",angle=90), axis.text.y=element_text(size=11,color="black"), 
        axis.title=element_text(size=14,face="bold"),text=element_text(size=14,face="bold"))+
  ylab("Relative Abundance")+xlab("")

###### Fig. 3 Phylogenetic tree of all genera with significant differences in abundance between pre-flowering drought and control root samples at TP8.

TP <- subset_samples(Bushman, Timepoint =="TP8")
TP <- subset_samples(TP,Treatment != "Post_flowering")
TP <- subset_samples(TP,SampleType == "Root")

# This converts to a DESeq2 format.
diagdds = phyloseq_to_deseq2(TP, ~ Treatment)
# This helps to obtain the normalization factor for your dataset.
diagdds = estimateSizeFactors(diagdds)
# You also need to demonstrate that your dataset is effectively dispersed (?).
diagdds = estimateDispersions(diagdds)
# This will reduce sampling variance.
diagvst = getVarianceStabilizedData(diagdds)
dim(diagvst)
#norm=counts(diagdds, normalized=TRUE)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
resultsNames(diagdds)
res = results(diagdds)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(TP)[rownames(sigtab), ], "matrix"))


new <- subset(sigtab,Phylum != "NA")
Genus <- subset(new,Genus != "NA")
na <- setdiff(rownames(new),rownames(Genus))
na.complete <- new[na,]
class <- na.complete[is.na(na.complete$Class),]
class$Class <- paste(class$Phylum,"PH",sep="")
class$Order <- paste(class$Phylum,"PH",sep="")
class$Family <- paste(class$Phylum,"PH",sep="")
class$Genus <- paste(class$Phylum,"PH",sep="")

order <- na.complete[is.na(na.complete$Order),]
order2 <- subset(order,Class!="NA")
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

#hand over the data to iTOL for plotting
write.csv(uniq,"root_treatment_tp8.csv")


##### Fig.4a Relative abundance across the top thirteen phyla for all transcripts for which taxonomies could be assigned in the metatranscriptome data
 
biom_file <- "fig4a-1-bacteria_count.biom"
map_file <- "fig4a-2-metadata_metatrans.txt"

set.seed(1)
biomot=import_biom(biom_file,parseFunction = parse_taxonomy_greengenes)

saveRDS(biomot, "biomot.rds")
bmsd = import_qiime_sample_data(map_file)
Bushman = merge_phyloseq(biomot,bmsd)
saveRDS(Bushman, "Bushman.rds")
sample_data(Bushman)$SampleType<-factor(sample_data(Bushman)$SampleType, levels=c("Soil","Rhizosphere"))
sample_data(Bushman)$Treatment<-factor(sample_data(Bushman)$Treatment, levels=c("Control","Pre_flowering","Post_flowering"))
sample_data(Bushman)$Timepoint<-factor(sample_data(Bushman)$Timepoint, levels=c("TP8","TP9"))
Bushman = subset_taxa(Bushman, Phylum != "Streptophyta")

rhi = subset_samples(Bushman, SampleType=="Rhizosphere")
Phylum <- tax_glom(rhi, taxrank="Phylum")
saveRDS(Phylum, "rhi_soil_Phylum_Bacteria.rds")

x <- c("TP8","TP9")
rm <- {}
for(i in x){
  TP = subset_samples(Phylum, Timepoint == i)
  merged_r = merge_samples(TP, "Treatment")
  merged_r_m = transform_sample_counts(merged_r , function(x) 100 * x/sum(x))
  physeqdf <- psmelt(merged_r_m)
  #temp <- c("BT642Control","RT430Control","BT642Pre_flowering","RT430Pre_flowering","BT642Post_flowering","RT430Post_flowering")
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
saveRDS(Phylum, "soil_Phylum_Bacteria.rds")

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

##### Fig.4b Gene Ontology (GO) enrichment analysis for all differentially expressed genes

#TP8 input

data<-read.table("fig4b-1-new_heatmap_enrichment_TP8.txt",header=T,sep="\t")

tp8_rhi <- data[,c(1,3,4)]
ratio <- data.frame(cbind(tp8_rhi,tp8_rhi$Rhi/tp8_rhi$total))
dat <- ratio[,c(1,4)]
colnames(dat) <-c("Annotation","TP8_rhizosphere")

tp8_soil <- data[,c(1,2,4)]
ratio <- data.frame(tp8_soil$Soil/tp8_soil$total)
colnames(ratio) <- c("TP8_soil")

mer8 <- data.frame(cbind(dat,ratio))

#TP9 input
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
dat <- subset(dat,Timepoint=="TP8")

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

##### Fig.4c Relative abundance across the top thirteen phyla for all differentially expressed transcripts by treatment 

data<-read.table("fig4c-new_input_function_taxa.txt",header=T,sep="\t")
data <- data.frame(data)

sum <- tapply(data$Percentage,data$Phylum, FUN=sum)
sum <- data.frame(sum)
sum <- sort(sum$sum,decreasing = T)
list <- c("Acidobacteria","Actinobacteria","Bacteroidetes","Chloroflexi","Cyanobacteria",
          "Firmicutes","Gemmatimonadetes","Nitrospirae",
          "Planctomycetes","Proteobacteria","Tenericutes","TM7","Verrucomicrobia")

a <- data[data$Phylum %in% list,]
b <- data[!(data$Phylum %in% list),]
b$Phylum <- "Other"
rm <- rbind(a,b)

rm$Timepoint<-factor(rm$Timepoint, levels=c("TP8","TP9"))

col_SampleByTreatment<-c("Acidobacteria"="thistle","Actinobacteria"="lightsalmon2","Armatimonadetes"="mediumpurple1","Bacteroidetes"="cornsilk","Chloroflexi"="lightskyblue","Cyanobacteria"="lightpink","Firmicutes"="lightblue1","Gemmatimonadetes"="darkseagreen","Nitrospirae"="mistyrose2",
                         "Planctomycetes"="gray69","Proteobacteria"="cadetblue3","SPAM"="darkorchid1","Tenericutes"="darkolivegreen2","TM7"="yellow","Verrucomicrobia"="slategray1","Other"="plum")


rm$Phylum <- factor(rm$Phylum,levels=c("Actinobacteria","Chloroflexi","Cyanobacteria","Firmicutes","Bacteroidetes","Proteobacteria",
                                       "Acidobacteria","Armatimonadetes","Gemmatimonadetes","Nitrospirae",
                                       "Planctomycetes","SPAM","Tenericutes","TM7","Verrucomicrobia","Other" ))

rm$Category <- factor(rm$Category, levels=c("Carbohydrate transport and metabolism",
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

ggplot(rm, aes(x=Category, y=Percentage, fill=Phylum)) + 
  geom_bar(stat = "identity", color = "NA") + scale_fill_manual(values = col_SampleByTreatment) + 
  facet_wrap(~Timepoint+SampleType)+
  theme(axis.text.x=element_text(size=9,color="black",hjust=1,vjust=0,face="bold"), 
        axis.text.y=element_text(size=9,color="black",face="bold"), 
        axis.title=element_text(size=11,face="bold"),text=element_text(size=11,vjust=1,face="bold"))+
  theme(legend.title = element_text(colour="black", size=11, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 11, face = "bold"))+
  ylab("Relative Abundance") +
  coord_flip()

##### Fig. 5 violin plot

rm(list = ls())
data<-read.table("fig5-Confocal.txt",header=T,sep="\t")

sai73 = subset(data,Strain=="SAI73")
sai73$Treatment <- as.factor(sai73$Treatment)
sai73_cd <- subset(sai73,Treatment != "PEG")

ggplot(sai73_cd, aes(x=Treatment, y=Value, fill=Treatment)) +  
  geom_violin(trim = TRUE)+
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="pointrange", color="white") +
  scale_fill_manual(values=c("#6AB187", "#DE7A22"))