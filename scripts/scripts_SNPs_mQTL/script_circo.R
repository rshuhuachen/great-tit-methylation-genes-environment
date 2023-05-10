# Script to filter output mQTL analysis and create circos from output mQTL analysis 
# Bernice Sepers, NIOO-KNAW, 2022

trans_raw<-read.table("/home/output2_trans.txt", sep = ",", header = TRUE)
nrow(trans_raw) #108,719,643
rm(trans_raw)
#filtered distance CpG-SNP
trans_raw<-read.table("/home/output2_trans_filtered.txt", sep = ",", header = TRUE)
nrow(trans_raw) #107,955,112

cis_raw<-read.table("/home/output2_cis.txt", sep = ",", header = TRUE)
nrow(cis_raw) #459,777

#check significant sites
head(cis_raw)
cis_raw$pvalue <- as.numeric(cis_raw$pvalue)
#cis p-value threshold: 0.05/459777 = 1.087484e-07
cis<-subset(cis_raw, cis_raw$pvalue < 1.087484e-07)
nrow(cis) #754
head(cis)
length(unique(cis$gene)) #680
length(unique(cis$snps)) #499
cis$X<-factor(cis$X)
write.table(cis, file = "/Users/bernicesepers/cis_sign.txt", append = FALSE, sep = ",", dec = ".",
            row.names = FALSE, col.names = TRUE)
rm(cis_raw)

head(trans_raw)
trans_raw$pvalue <- as.numeric(trans_raw$pvalue)
#trans p-value threshold: 0.05/107955112 = 4.631555e-10
trans<-subset(trans_raw, trans_raw$pvalue < 4.631555e-10)
nrow(trans) #4,202
length(unique(trans$gene)) #1,145
length(unique(trans$snps)) #1,836
head(trans)
trans$X<-factor(trans$X)
write.table(trans, file = "/Users/bernicesepers/trans_sign.txt", append = FALSE, sep = ",", dec = ".",
            row.names = FALSE, col.names = TRUE)
#remove X in first column in text files using nano in the terminal, 
#then use create_mQTL.py to create mqtl plot
rm(trans_raw)

library(dplyr)
trans_test <- trans %>%
  select(c("snps", "gene"))

cis_test <- cis %>%
  select(c("snps", "gene"))

df <- rbind(cis_test, trans_test)
nrow(df) #4956
nrow(cis)+nrow(trans) #4956
length(unique(df$gene)) #1600
length(unique(df$snps)) #2230
df_a<-as.data.frame(table(df$snps))
nrow(subset(df_a, df_a$Freq >= 2)) #871 snps that have an association with multiple CpGs

########## create circos
cis<-read.table("/Users/bernicesepers/cis_sign.txt", sep = ",", header = TRUE)
trans<-read.table("/Users/bernicesepers/trans_sign.txt", sep = ",", header = TRUE)

library(BioCircos)

myGenome = list("chr1" = 114059860,
                "chr2" = 	150265477,
                "chr3" = 	111636321,
                "chr4" = 	68030631,
                "chr4A" = 	19934923,
                "chr1A" = 	71365269,
                "chr5" = 	61875929,
                "chr6" = 	34342011,
                "chr7" = 	37653027,
                "chr8" = 	31324166,
                "chr9" = 	25106277,
                "chr10" = 	20202851,
                "chr11" = 	20315886,
                "chr12" = 	20466350,
                "chr13" = 	16480340,
                "chr14" = 	16193477,
                "chr15" = 	13820886,
                "chr17" = 	10486032,
                "chr18" = 	11572643,
                "chr19" = 	9871655,
                "chr20" = 	14661763,
                "chr21" = 	7693166,
                "chr22" = 	4276343,
                "chr23" = 	6655392,
                "chr24" = 	6808513,
                "chr25LG1" = 	1092960,
                "chr25LG2" = 	809223,
                "chr26" = 	6596997,
                "chr27" = 	4327975,
                "chr28" = 	5101010,
                "chrLGE22" = 	773534,
                "chrZ" = 	74514349,
                "scaffolds" = 2566323)

head(trans)
trans$chr<-as.factor(trans$snps)
library(tidyr)
trans <- trans %>% separate(chr, c("chr1","chr2", "snp_start"), sep = "([_])")
head(trans)
library(tidyverse)
trans <- add_column(trans, "snp_chr" = paste(trans$chr1, trans$chr2, sep = "_"))
trans$chr1 <- NULL
trans$chr2 <- NULL
trans$X<-NULL

trans$gene_chr<-as.factor(trans$gene)
trans <- trans %>% separate(gene_chr, c("chr1","chr2", "gene_start"), sep = "([_])")
trans <- add_column(trans, "gene_chr" = paste(trans$chr1, trans$chr2, sep = "_"))
trans$chr1 <- NULL
trans$chr2 <- NULL

chr_names<- read.table(file="/Users/bernicesepers/GCF_001522545.3_Parus_major1.1_assembly_report_chr_names.txt",sep=" ",na.strings = "na")
colnames(chr_names)<-c("Name","chr")
chr_names$chr<-as.factor(chr_names$chr)
chr_names$Name<-as.factor(chr_names$Name)
head(chr_names)
chr_names$chrName_scaffolds<-as.character(chr_names$Name)
Scaffolds<-grep(pattern = "Scaffold",x = chr_names$Name)
chr_names$chrName_scaffolds[Scaffolds]<-"scaffolds"
chr_names$chrName_scaffolds<-as.factor(chr_names$chrName_scaffolds)
head(chr_names)
tail(chr_names)
rm(Scaffolds)
chr_names$Name <-NULL
chr_names$Name <- chr_names$chrName_scaffolds
chr_names$chrName_scaffolds <- NULL

chr_names$gene_chr <- chr_names$chr
chr_names$chr <- NULL
head(trans)
head(chr_names)
trans = left_join(trans, chr_names,by = "gene_chr")
head(trans)
trans$gene_name <- trans$Name
trans$Name <- NULL

chr_names$snp_chr <- chr_names$gene_chr
chr_names$gene_chr <- NULL
trans = left_join(trans, chr_names,by = "snp_chr")
head(trans)
trans$snp_name <- trans$Name
trans$Name <- NULL

points_chromosomes<-as.character(trans$snp_name)
points_coordinates<-trans$snp_start
points_values<- -log10(trans$pvalue)

#tracklist with SNPs
tracklist = BioCircosSNPTrack('mySNPTrack', points_chromosomes, points_coordinates, 
                              points_values, colors = c("darkblue"), minRadius = 0.8, maxRadius = 0.9, size = 1)

#tracklist with SNPs and background
tracklist = tracklist + BioCircosBackgroundTrack("myBackgroundTrack", 
                                                 minRadius = 0.8, maxRadius = 0.9,
                                                 borderColors = "#AAAAAA", borderSize = 0.6, fillColors = "#B3E6FF") 
BioCircos(tracklist, genomeFillColor = "Blues", genome=myGenome,
          chrPad = 0.05, displayGenomeBorder = FALSE,
          genomeTicksDisplay = FALSE,  genomeLabelTextSize = 18, genomegenomeLabelDy = 30, genomeLabelOrientation = -90)

#add links
links_chromosomes_1 = as.character(trans$snp_name) # Chromosomes on which the links should start
links_pos_1 = as.numeric(trans$snp_start) #position where the links should start
links_chromosomes_2 = as.character(trans$gene_name) # Chromosomes on which the links should end
links_pos_2 = as.numeric(trans$gene_start) #position where the links should end

tracklist2 = tracklist + BioCircosLinkTrack('myLinkTrack', links_chromosomes_1, links_pos_1,
                                            links_pos_1, links_chromosomes_2, links_pos_2, links_pos_2,
                                            maxRadius = 0.80, width = "0.003em", color = "#99CCFF",  displayAxis = FALSE)

BioCircos(tracklist2, genomeFillColor = "Blues", genome=myGenome,
          chrPad = 0.05, displayGenomeBorder = FALSE,
          genomeTicksDisplay = FALSE,  genomeLabelTextSize = 25, genomeLabelDy = 55, genomeLabelOrientation = -90)


#DNMT3a
DNMT3a<-subset(trans, trans$snp_start == 117251)
links_chromosomes_1 = as.character(DNMT3a$snp_name) # Chromosomes on which the links should start
links_pos_1 = as.numeric(DNMT3a$snp_start) #position where the links should start
links_chromosomes_2 = as.character(DNMT3a$gene_name) # Chromosomes on which the links should end
links_pos_2 = as.numeric(DNMT3a$gene_start) #position where the links should end

tracklist3 = tracklist + BioCircosLinkTrack('myLinkTrack', links_chromosomes_1, links_pos_1,
                                             links_pos_1, links_chromosomes_2, links_pos_2, links_pos_2,
                                             maxRadius = 0.80, width = "0.1em", color = "#99CCFF", 
                                             displayAxis = FALSE) + BioCircosTextTrack(
                                              'CTNNA3', 'CTNNA3', size = "1.9em", 
                                             x = 0.03, y = 0.67, opacity = 0.65) + BioCircosTextTrack(
                                               'DNMT3a', 'DNMT3a', size = "1.9em", 
                                               x = -0.02, y = -0.67)

BioCircos(tracklist3, genomeFillColor = "Blues", genome=myGenome,
          chrPad = 0.05, displayGenomeBorder = FALSE,
          genomeTicksDisplay = FALSE,  genomeLabelTextSize = 25, genomeLabelDy = 55, genomeLabelOrientation = -90)

#check which SNPs have the most relationships
b<-as.data.frame(table(cis$snps)) #
b[which.max(b$Freq),] #NC_031771.1_1085160 has 13
subset(b, b$Freq >= 11) #NC_031770.1_2889149 has 11
#NC_031773.1_56904575 has 11
nrow(subset(b, b$Freq >= 2)) #108

a<-as.data.frame(table(trans$snps))
a[which.max(a$Freq),] #NC_031788.1_14006553 has 23
subset(a, a$Freq >= 18) #NC_031775.1_33814297 has 20
#NW_015379289.1_52389 has 18
subset(a, a$Freq >= 15)
nrow(subset(a, a$Freq >= 2)) #739

ab<-rbind(a,b)
ab[which.max(ab$Freq),] #still NC_031788.1_14006553

#NC_031788.1_14006553
NC_031788.1_14006553<-subset(trans, trans$snps == "NC_031788.1_14006553")
nrow(NC_031788.1_14006553)
links_chromosomes_1 = as.character(NC_031788.1_14006553$snp_name) # Chromosomes on which the links should start
links_pos_1 = as.numeric(NC_031788.1_14006553$snp_start) #position where the links should start
links_chromosomes_2 = as.character(NC_031788.1_14006553$gene_name) # Chromosomes on which the links should end
links_pos_2 = as.numeric(NC_031788.1_14006553$gene_start) #position where the links should end

tracklist3 = tracklist + BioCircosLinkTrack('myLinkTrack', links_chromosomes_1, links_pos_1,
                                             links_pos_1, links_chromosomes_2, links_pos_2, links_pos_2,
                                             maxRadius = 0.80, width = "0.05em", color = "darkblue", 
                                            displayAxis = FALSE) + BioCircosTextTrack('SLC9A8', 'SLC9A8', 
                                            size = "1.5em", x = -0.78, y = -0.08) + BioCircosTextTrack(
                                            'CPSF3L', 'CPSF3L',size = "1.5em", x = -0.68, y = -0.35, opacity = 0.65) + BioCircosTextTrack(
                                            'POLR1C', 'POLR1C', size = "1.5em", x = 0.4, y = -0.1, opacity = 0.65) + BioCircosTextTrack(
                                            'ST8SIA1', 'ST8SIA1', size = "1.5em", x = 0.05, y = 0.6, opacity = 0.65) + BioCircosTextTrack(
                                            'NID2', 'NID2', size = "1.5em", x = -0.03, y = 0.76, opacity = 0.65) + BioCircosTextTrack(
                                            'COL18A1', 'COL18A1', size = "1.5em", x = -0.3, y = 0.67, opacity = 0.65)

BioCircos(tracklist3, genomeFillColor = "Blues", genome=myGenome,
          chrPad = 0.05, displayGenomeBorder = FALSE,
          genomeTicksDisplay = FALSE,  genomeLabelTextSize = 25, genomeLabelDy = 55, genomeLabelOrientation = -90)

##############cis
head(cis)
cis$chr<-as.factor(cis$snps)
library(tidyr)
cis <- cis %>% separate(chr, c("chr1","chr2", "snp_start"), sep = "([_])")
head(cis)
library(tidyverse)
cis <- add_column(cis, "snp_chr" = paste(cis$chr1, cis$chr2, sep = "_"))
cis$chr1 <- NULL
cis$chr2 <- NULL
cis$X<-NULL

cis$gene_chr<-as.factor(cis$gene)
cis <- cis %>% separate(gene_chr, c("chr1","chr2", "gene_start"), sep = "([_])")
cis <- add_column(cis, "gene_chr" = paste(cis$chr1, cis$chr2, sep = "_"))
cis$chr1 <- NULL
cis$chr2 <- NULL

chr_names<- read.table(file="/Users/bernicesepers/GCF_001522545.3_Parus_major1.1_assembly_report_chr_names.txt",sep=" ",na.strings = "na")
colnames(chr_names)<-c("Name","chr")
chr_names$chr<-as.factor(chr_names$chr)
chr_names$Name<-as.factor(chr_names$Name)
head(chr_names)
chr_names$chrName_scaffolds<-as.character(chr_names$Name)
Scaffolds<-grep(pattern = "Scaffold",x = chr_names$Name)
chr_names$chrName_scaffolds[Scaffolds]<-"scaffolds"
chr_names$chrName_scaffolds<-as.factor(chr_names$chrName_scaffolds)
head(chr_names)
tail(chr_names)
rm(Scaffolds)
chr_names$Name <-NULL
chr_names$Name <- chr_names$chrName_scaffolds
chr_names$chrName_scaffolds <- NULL

chr_names$gene_chr <- chr_names$chr
chr_names$chr <- NULL
head(cis)
head(chr_names)
cis = left_join(cis, chr_names,by = "gene_chr")
head(cis)
cis$gene_name <- cis$Name
cis$Name <- NULL

chr_names$snp_chr <- chr_names$gene_chr
chr_names$gene_chr <- NULL
cis = left_join(cis, chr_names,by = "snp_chr")
head(cis)
cis$snp_name <- cis$Name
cis$Name <- NULL

points_chromosomes<-as.character(cis$snp_name)
points_coordinates<-cis$snp_start
points_values<- -log10(cis$pvalue)

#tracklist with SNPs
tracklist4 = BioCircosSNPTrack('mySNPTrack', points_chromosomes, points_coordinates, 
                               points_values, colors = c("darkblue"), minRadius = 0.5, maxRadius = 0.9, size = 1)

#tracklist with SNPs and background
tracklist5 = tracklist4 + BioCircosBackgroundTrack("myBackgroundTrack", 
                                                   minRadius = 0.5, maxRadius = 0.9,
                                                   borderColors = "#AAAAAA", borderSize = 0.6, fillColors = "#B3E6FF") 
BioCircos(tracklist5, genomeFillColor = "Blues", genome=myGenome,
          chrPad = 0.05, displayGenomeBorder = FALSE,
          genomeTicksDisplay = FALSE,  genomeLabelTextSize = 25, genomeLabelDy = 55, genomeLabelOrientation = -90)

####distance between SNPs and CpG sites in cis
identical(cis$snp_chr, cis$gene_chr) 
cis$distance <- as.numeric(cis$snp_start) - as.numeric(cis$gene_start)
head(cis)
tail(cis)
cis$distance
length(subset(cis$distance, cis$distance == 1)) #0
length(subset(cis$distance, cis$distance == -1)) #15
range(cis$distance) #-986975 - 990443
mean(cis$distance) #-8935.889
cis$distance2 <- as.numeric(gsub(cis$distance, pattern = "-", replacement = ""))
range(cis$distance2) #1 - 990443
mean(cis$distance2) #150007.7
median(cis$distance2) #278
nrow(cis) #754
length(subset(cis$distance2, cis$distance2 == 0)) #0
length(subset(cis$distance2, cis$distance2 == 1)) #15
length(subset(cis$distance2, cis$distance2 > 1)) #739
summary(cis$distance2)
c<-as.data.frame(table(cis$distance2))
c[which.max(c$Freq),]
subset(c, c$Freq == 23)

#distance SNP - CpG in tss, gene body and promoter region
CG<-read.table("/Users/bernicesepers/TSS_promo_genebody_CpGs_origin.txt", header = TRUE)
head(CG)
nrow(CG)
table(CG$region)
head(cis)
cis_CG = left_join(cis, CG,by = "gene")
head(cis_CG)
table(cis_CG$region)
cis_CG_promo<-subset(cis_CG, cis_CG$region == "promo" |
                       cis_CG$region == "tss")
head(cis_CG_promo)
table(cis_CG_promo$region)
cis_CG_promo$distance <- NULL
cis_CG_promo$distance2 <- NULL
cis_CG_promo$distance <- as.numeric(cis_CG_promo$snp_start) - as.numeric(cis_CG_promo$gene_start)
cis_CG_promo$distance2 <- as.numeric(gsub(cis_CG_promo$distance, pattern = "-", replacement = ""))
range(cis_CG_promo$distance2) #1-974622
mean(cis_CG_promo$distance2) #154211.7
median(cis_CG_promo$distance2) #159
nrow(cis_CG_promo) #97
length(unique(cis_CG_promo$snps)) #60
length(subset(cis_CG_promo$distance2, cis_CG_promo$distance2 == 0)) #0
length(subset(cis_CG_promo$distance2, cis_CG_promo$distance2 == 1)) #1
length(subset(cis_CG_promo$distance2, cis_CG_promo$distance2 > 1)) #96
summary(cis_CG_promo$distance2)

cis_CG_tss<-subset(cis_CG, cis_CG$region == "tss")
head(cis_CG_tss)
table(cis_CG_tss$region)
cis_CG_tss$distance <- NULL
cis_CG_tss$distance2 <- NULL
cis_CG_tss$distance <- as.numeric(cis_CG_tss$snp_start) - as.numeric(cis_CG_tss$gene_start)
cis_CG_tss$distance2 <- as.numeric(gsub(cis_CG_tss$distance, pattern = "-", replacement = ""))
range(cis_CG_tss$distance2) #3-509893
mean(cis_CG_tss$distance2) #96742.97
median(cis_CG_tss$distance2) #96
nrow(cis_CG_tss) #33
length(unique(cis_CG_tss$snps)) #15
length(subset(cis_CG_tss$distance2, cis_CG_tss$distance2 == 0)) #0
length(subset(cis_CG_tss$distance2, cis_CG_tss$distance2 == 1)) #0
length(subset(cis_CG_tss$distance2, cis_CG_tss$distance2 > 1)) #33
summary(cis_CG_tss$distance2)

cis_CG_body<-subset(cis_CG, cis_CG$region == "gene_body")
head(cis_CG_body)
table(cis_CG_body$region)
cis_CG_body$distance <- NULL
cis_CG_body$distance2 <- NULL
cis_CG_body$distance <- as.numeric(cis_CG_body$snp_start) - as.numeric(cis_CG_body$gene_start)
cis_CG_body$distance2 <- as.numeric(gsub(cis_CG_body$distance, pattern = "-", replacement = ""))
range(cis_CG_body$distance2) #1-990443
mean(cis_CG_body$distance2) #155755.4
median(cis_CG_body$distance2) #1065
nrow(cis_CG_body) #160
length(unique(cis_CG_body$snps)) #115
length(subset(cis_CG_body$distance2, cis_CG_body$distance2 == 0)) #0
length(subset(cis_CG_body$distance2, cis_CG_body$distance2 == 1)) #1
length(subset(cis_CG_body$distance2, cis_CG_body$distance2 > 1)) #159
summary(cis_CG_body$distance2)

head(trans)
trans_CG = left_join(trans, CG,by = "gene")
head(trans_CG)
table(trans_CG$region)