##### Experimental design

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

##### Metadata file: switch TP14_S22 with TP14_S20;switch TP14_R17 with TP14_R20;switch TP14_R18 with TP14_R19
set.seed(1)
Bushman = readRDS("Bushman.rds")

##### Throw away low quality samples
Bushman= subset_samples(Bushman,Sample_ID != "TP7_R20")
Bushman= subset_samples(Bushman,Sample_ID != "TP7_R14")
Bushman= subset_samples(Bushman,Sample_ID != "TP7_L20")
Bushman= subset_samples(Bushman,Sample_ID != "TP7_L14")
Bushman = subset_samples(Bushman,Sample_ID != "TP5_Z5")
Bushman = subset_samples(Bushman,Sample_ID != "TP8_S22")

##### Total read counts after quality filtering
a <- sample_sums(Bushman)
b = data.frame(a)
sum(b$a)

##### Removed OTUs without at least 5 reads in at least 3 samples
Bushman= prune_taxa(taxa_sums(Bushman)>=5, Bushman)
Bushman = filter_taxa(Bushman, function(x) sum(x >= 1) > 2, TRUE)

##### Remove samples that have less than 10000 reads count
Bushman= prune_samples(sample_sums(Bushman)>=10000, Bushman)
Bushman= prune_taxa(taxa_sums(Bushman)>=1, Bushman)
sample_data(Bushman)$SampleType<-factor(sample_data(Bushman)$SampleType, levels=c("Soil","Rhizosphere","Root"))
sample_data(Bushman)$Treatment<-factor(sample_data(Bushman)$Treatment, levels=c("Control","Pre_flowering","Post_flowering"))
sample_data(Bushman)$Timepoint<-factor(sample_data(Bushman)$Timepoint, levels=c("TP1","TP2","TP3","TP4","TP5","TP6","TP7","TP8","TP9","TP10","TP11","TP12","TP13","TP14","TP15","TP16","TP17"))
sample_data(Bushman)$TreatmentByTimepoint<-factor(sample_data(Bushman)$TreatmentByTimepoint, levels=c("ControlTP1",	"ControlTP2",	"ControlTP3",	"ControlTP4",	"ControlTP5",	"ControlTP6",	"ControlTP7",	"ControlTP8",	"ControlTP9",	"ControlTP10",	"ControlTP11",	"ControlTP12",	"ControlTP13",	"ControlTP14",	"ControlTP15",	"ControlTP16",	"ControlTP17",	"Pre_floweringTP3",	"Pre_floweringTP4",	"Pre_floweringTP5",	"Pre_floweringTP6",	"Pre_floweringTP7",	"Pre_floweringTP8",	"Pre_floweringTP9",	"Pre_floweringTP10",	"Pre_floweringTP11",	"Pre_floweringTP12",	"Pre_floweringTP13",	"Pre_floweringTP14",	"Pre_floweringTP15",	"Pre_floweringTP16",	"Pre_floweringTP17",	"Post_floweringTP10",	"Post_floweringTP11",	"Post_floweringTP12",	"Post_floweringTP13",	"Post_floweringTP14",	"Post_floweringTP15",	"Post_floweringTP16",	"Post_floweringTP17",	"Post_floweringTP8",	"Post_floweringTP9"	))

#rar = rarefy_even_depth(Bushman,sample.size=13000)
#rar = prune_taxa(taxa_sums(rar)>=1, rar)
rar<-readRDS("rar.rds")

##### Total reads count calculation
a <- sample_sums(rar)
b = data.frame(a)
sum(b$a)