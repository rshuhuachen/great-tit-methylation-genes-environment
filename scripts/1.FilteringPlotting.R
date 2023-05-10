##### This script is for filtering data in methylKit, conducting PCA and conducting a DMA in DSS #####

#### Load all packages ####

library(data.table);library(stringr); library(dplyr);library(ggplot2)
library(devtools); library(reshape2); library(tibble); library(stringr)
library(ggdendro); library(tibble);library(RColorBrewer); library(plotly)
library(r2glmm); library(tidystats); library(plyr); library(robustlmm)
library(gaston); library(readxl); library(MuMIn); library(lmerTest)
library(sjmisc); library(renv); library(DSS);library(methylKit);require(bsseq)
library(GenomicFeatures); library(rtracklayer); library(genomation)
library(matrixStats); library("factoextra")


#### Set working directory ####
setwd("/home/nioo/rebeccash/2018_BS_manipulation/R_Scripts/complete_MS_pt1") #set this to the cloned github dir
# within the directory, you will find the scripts and the data in the data directory. 


## Other elements / vectors / dataframes to load #####

# Theme for ggplot
gtit <-  theme(panel.background = element_rect(fill="white", colour="white", size=0.5, linetype="solid", color="grey70")) +
  theme(panel.grid.major = element_line(colour="grey90", size=0.1, linetype="blank",lineend="butt",color="grey90"),
        panel.grid.minor = element_line(colour="grey90", size=0.1, linetype="blank",lineend="butt",color="grey90"),
        plot.title = element_text(face="bold", color="black", size=25, hjust = -0.1),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14))

# Reference genome
base::load(file = "data/reference_genome/GCF_001522545.3_Parus_major1.1_assembly_report_chr_names.RData")
#ref genome parus major v1.1 can be downloaded from https://www.ncbi.nlm.nih.gov/assembly/GCF_001522545.3

# Chicks parameters
chickdata <- read.csv("data/chickparameters.txt")

#### Load in raw unfiltered data from chicks ####
# this raw data can be downloaded from XXXXX , set paths accordingly
FileLocation.chicks <- "/scratch/VanOersgroup/2021_epiGBS_2018BSmanip_B1-B8/results/bismark_calling_new/" #adjust to your directory with raw data reads
FileNames.chicks <- list.files(path = FileLocation.chicks, pattern = "*CpG_report.txt")

list_files.chick <- chickdata$File
list_files.chick <- paste(FileLocation.chicks,list_files.chick)

list_files.chick <- str_replace_all(list_files.chick, fixed(" "), "")#to remove whitespaces
list_files.chick <- str_replace_all(list_files.chick, pattern = "/scratch/VanOersgroup/2021_epiGBS_2018BSmanip_B1-B8/results/bismark_calling_new", 
                                    replacement = "/scratch/VanOersgroup/2021_epiGBS_2018BSmanip_B1-B10/results/bismark_calling_new")
list_files.chick <- as.list(list_files.chick)
list_chicks <- as.list(chickdata$RNR)

# read in raw data to methylkit without filtering for coverage
setwd("/scratch/VanOersgroup/2021_epiGBS_2018BSmanip_B1-B10/results/bismark_calling_new")

myobj.chick.unfiltered <- methRead(list_files.chick, pipeline = "bismarkCytosineReport", 
                                   sample.id = list_chicks, assembly = "p.major1.1", 
                                   treatment = c(rep(1, each = length(list_files.chick))), context = "CpG", mincov = 1)

## Unite ##
meth.chick.unfiltered <- unite(myobj.chick.unfiltered, destrand = TRUE, min.per.group = 1L)
nrow(meth.chick.unfiltered) #2,768,598

#### Summary statistics ####
methData.chick.unfiltered <- getData(meth.chick.unfiltered)
sample_id.chick <- meth.chick.unfiltered@sample.ids
site.chick <- paste(methData.chick.unfiltered$chr, methData.chick.unfiltered$start, sep="_")
methDataNew.chick <- data.frame(site=site.chick) #dataframe with proportion methylated
methDataNew_Cov.chick <- data.frame(site=site.chick) #dataframe with coverage
methDataNew_Meth.chick <- data.frame(site=site.chick) #dataframe with methylated C's

ncol(methData.chick.unfiltered) #670 columns, 5th column is first individual, second one is 8th column and so on
#so for every individual (f), get coverage, methylated C's and proportion methylated.
for(f in seq(5,670,3)){
  cov <- methData.chick.unfiltered[,f]
  methCount <- methData.chick.unfiltered[,(f+1)]
  methDataNew.chick$methProp <- as.data.frame(methCount/cov)
  methDataNew_Cov.chick$cov <- as.data.frame(cov)
  methDataNew_Meth.chick$meth <- as.data.frame(methCount)
  no <- (f-2)/3 #number of samples, 286
  names(methDataNew.chick)[no+1] <- sample_id.chick[no] #change column names to sample names
  names(methDataNew_Cov.chick)[no+1] <- sample_id.chick[no] #change column names to sample names
  names(methDataNew_Meth.chick)[no+1] <- sample_id.chick[no] #change column names to sample names
}

row.names(methDataNew.chick) <- methDataNew.chick$site #sites will be row names in MethDataNew
methDataNew.chick$site <- NULL

row.names(methDataNew_Cov.chick) <- methDataNew_Cov.chick$site

methDataNew_Cov.chick$Mean <- rowMeans(methDataNew_Cov.chick[,-1], na.rm=T) #add column mean coverage
methDataNew_Cov.chick$Mean_meth <- rowMeans(methDataNew.chick[,-1], na.rm=T) #add column mean methylation
methDataNew_Cov.chick$SD_cov <- rowSds(as.matrix(methDataNew_Cov.chick[,-1]), na.rm=T) #add column SD coverage
methDataNew_Cov.chick$SD_meth <- rowSds(as.matrix(methDataNew.chick[,-1]), na.rm=T) #add column mean methylation
methDataNew_Cov.chick$NAs <- apply(methDataNew_Cov.chick, 1, function(x) sum(is.na(x))) #add column N
methDataNew_Cov.chick$N <- 222-methDataNew_Cov.chick$NAs

names(methDataNew_Cov.chick)[names(methDataNew_Cov.chick) == "Mean"] <- "Mean_coverage"
CovData.chick <- methDataNew_Cov.chick[,c(1,224, 225)] #table with mean coverage and methylation per site

row.names(methDataNew_Meth.chick) <- methDataNew_Meth.chick$site

ggplot(data=CovData.chick, aes(x=log(Mean_coverage), y=Mean_meth)) + 
  geom_point(size=0.7) + 
  theme_set(theme_grey()) #there's no obvious pattern between coverage and methylation level

#### PCA ####
#create dataset for PCA
methDataNew.chick.df <- data.frame(matrix(unlist(methDataNew.chick), nrow=nrow(methDataNew.chick)),stringsAsFactors=FALSE)
meth_PCA.chick <- methDataNew.chick[complete.cases(methDataNew.chick.df),] #because we only include complete.cases, we use unfiltered data as otherwise nothing will be left for PCA 
nrow(meth_PCA.chick) #31335

#do PCA for all CpG sites and added to that, individuals
#meth_PCA.chick.noout <- meth_PCA.chick[,c(1:64, 66:222)] #without outlier BD_83116 col 65
PCA.chick <- prcomp(t(meth_PCA.chick), center=F, scale.=F) # t() transposes the matrix meth_PCA to get one coordinate for each of the ring numbers 
biplot(PCA.chick,scale = T,center=F,cex=0.25)
str(PCA.chick)
plot(x=PCA.chick$x[,1],y=PCA.chick$x[,2],pch=1);text(x=PCA.chick$x[,1],y=PCA.chick$x[,2]+0.5,labels = rownames(PCA.chick$x),cex=0.5)

#Visualize eigenvalues (scree plot). Shows the percentage of variances explained by each principal component.
fviz_eig(PCA.chick)

eigs.chick <- PCA.chick$sdev^2
var.chick <- eigs.chick/sum(eigs.chick)
explained.chick<-100*eigs.chick/sum(eigs.chick)

head(explained.chick) #PC1 explains 94.8% of data, so apparent outlier BD_83116 is probably not problematic

# PCA where each point is an ID to identify clusters between broods, libraries or sexes
# add pc scores of pc1 and pc2 to dataframe chickdata, from PCA
chicks_PCs <- chickdata
chicks_PCs$PC1<-NA
chicks_PCs$PC2<-NA

for(i in rownames(PCA.chick$x)){
  chicks_PCs$PC1[chicks_PCs$RNR==i]<-PCA.chick$x[,1][grep(x=rownames(PCA.chick$x),pattern=i)]
  chicks_PCs$PC2[chicks_PCs$RNR==i]<-PCA.chick$x[,2][grep(x=rownames(PCA.chick$x),pattern=i)]
}

chicks_PCs$GeneticNK <- as.factor(chicks_PCs$GeneticNK)
chicks_PCs$Treatment <- as.factor(chicks_PCs$Treatment)
chicks_PCs$Mol_sex <- as.factor(chicks_PCs$Mol_sex)
chicks_PCs$Age <- as.factor(chicks_PCs$Age)
chicks_PCs$Library <- as.factor(chicks_PCs$Library)

plot(chicks_PCs$PC1~chicks_PCs$PC2,col=chicks_PCs$Mol_sex,pch=19)
legend(x = "topleft",        
       legend = levels(chicks_PCs$Mol_sex),
       col = 1:3,
       pch = 19)

plot(chicks_PCs$PC1~chicks_PCs$PC2,col=as.factor(chicks_PCs$Treatment),pch=19)
legend(x = "topleft",          # Position
       legend = levels(chicks_PCs$Treatment),  # Legend texts
       col = 1:3,       # Color of the squares
       pch = 19)

plot(chicks_PCs$PC1~chicks_PCs$PC2,col=as.factor(chicks_PCs$Library),pch=19)
legend(x = "topleft",          # Position
       legend = levels(chicks_PCs$Library),  # Legend texts
       col = 1:8,       # Color of the squares
       pch = 19)

plot(chicks_PCs$PC1~as.factor(chicks_PCs$Library))
plot(chicks_PCs$PC2~as.factor(chicks_PCs$Library))

# no obvious sex, treatment or library effects

#### Filter chick data using methylkit ####

#discard bases that have coverage below 10X and more than 99.9 percentile of coverage in each sample
myobj.filtered.chick <- filterByCoverage(myobj.chick.unfiltered,lo.count=10,lo.perc=NULL,
                                         hi.count=NULL,hi.perc=99.9)   

#normalise coverage
myobj.filtered.norm.chick <- normalizeCoverage(myobj.filtered.chick)

#unite samples
cov_80p.chick <- as.integer(round(0.8 * length(meth.chick.unfiltered@sample.ids)),0) 
meth.chick <- methylKit::unite(myobj.filtered.norm.chick, destrand=TRUE, min.per.group = cov_80p.chick)
nrow(meth.chick) #117525 CpG's that are shared between 80% of individuals/samples with a minimum coverage of 10 per ID 

#### Calculate methylation % per ID which will be used for model #### 

#make a loop to create methylation for each sample

meth.chick.df <- methylKit::getData(meth.chick)
for (i in 1:length(meth.chick@sample.ids)){
  meth.chick.df[[paste0("methylation_", i)]] <- meth.chick.df[[paste0("numCs", i)]] / meth.chick.df[[paste0("coverage", i)]]
}

#take out all CpG sites that have no variation in methylation - all 100% or all 0%
startmeth.chick <- 4+length(meth.chick@sample.ids)*3+1
endmeth.chick <- 4+length(meth.chick@sample.ids)*4
meth.chick.df <- meth.chick.df %>% mutate(summeth = rowSums(.[grep("methylation_", names(.))], na.rm = TRUE))
meth.chick.df$totalna <- apply(meth.chick.df[startmeth.chick:endmeth.chick], 1, function(x) sum(is.na(x))) #calculate na in only meth cols
meth.chick.df$if.100.or.0 <- case_when(meth.chick.df$summeth == length(meth.chick@sample.ids) - meth.chick.df$totalna ~ "100%",
                                       meth.chick.df$summeth == 0 ~"0%",
                                       meth.chick.df$summeth > 0 & meth.chick.df$summeth < length(meth.chick@sample.ids) - meth.chick.df$totalna ~ "Variation")
meth.chick.df$if.100.or.0 <- as.factor(meth.chick.df$if.100.or.0)
summary(meth.chick.df$if.100.or.0) #4 are 0%, 117521 are not 0%

meth.chick.df <- subset(meth.chick.df, meth.chick.df$if.100.or.0 == "Variation")
nrow(meth.chick.df) #117521 CpG sites

#### Adjusting chr names ####
#ensure to load reference genome with codes for corresponding chr numbers
meth.chick.adjchr <- meth.chick.df

meth.chick.adjchr$adjchr <- mapvalues(x=meth.chick.adjchr$chr, # x = the column with NW/NC/MT-codes that we need the corresponding chr names for ()
                                      from=levels(chr_names$Code),    # from = the matching NW/NC/MT codes in the chr_names dataframe
                                      to  =levels(chr_names$Name))    # to = the shortened chromosome names we want to add in the dateset


# we do not want to plot a category (equivalent to a whole chromosome) for each scaffold, so we will make one more column that contains the name "Scaffolds" for each scaffold
meth.chick.adjchr$chrName_scaffolds<-as.character(meth.chick.adjchr$adjchr)
Scaffolds<-grep(pattern = "Scaffold",x = meth.chick.adjchr$adjchr)
meth.chick.adjchr$chrName_scaffolds[Scaffolds]<-"scaffolds"
meth.chick.adjchr$chrName_scaffolds<-as.factor(meth.chick.adjchr$chrName_scaffolds)
meth.chick.adjchr$adjchr<-meth.chick.adjchr$chrName_scaffolds
summary(meth.chick.adjchr$adjchr) #5851 are scaffolds

meth.chick.adjchr$adjchr <- as.factor(meth.chick.adjchr$adjchr)
meth.chick.adjchr$strand <- as.factor(meth.chick.adjchr$strand)

#move around adjchr column to the front
meth.chick.adjchr <- meth.chick.adjchr[,c(1,896,2:892)]

#### How many CpG sites per ID for unfiltered and filtered data? ####

id_list.chick <- character()

for (i in myobj.chick.unfiltered) {
  out <- print(nrow(i))
  id_list.chick <- c(id_list.chick, out)
}

id_list.chick <- as.data.frame(id_list.chick)

id_list.chick <- data.frame(ID = getSampleID(myobj.chick.unfiltered), id_list.chick)
names(id_list.chick)[2]<- "n_CpG_unfiltered"

filtered.chick <- character()
for (i in myobj.filtered.chick) {
  out <- print(nrow(i))
  filtered.chick <- c(filtered.chick, out)
}
filtered.chick <- as.data.frame(filtered.chick)
#add together to id_list
id_list.chick <- cbind(id_list.chick, filtered.chick)
id_list.chick

#### Reorder dataframes ####
# Write out methylation data divided per ID 

#with adjusted chr names, per file nr C, nr T, cov and meth%
select_a = c(6:8) #adjust this to adjchr col
select_b = c(5+(length(meth.chick@sample.ids)*3)+1) #adjust this to adjchr col
select = c(select_a, select_b)

for (i in 1:length(meth.chick@sample.ids)){
  id_name <- id_list.chick[i,1]
  myfile <- file.path(getwd(), paste0(id_name, "_methylation_df",".csv"))
  assign(paste(id_name, "_methylation_df.csv",sep = ""),meth.chick.adjchr[,c(1:5,select)]) %>% write.table(file=myfile, sep = "\t", row.names = FALSE)
  select_a = select_a + 3
  select_b = select_b + 1
  select = c(select_a, select_b)
}

#only for DSS, write out per ID but with adjusted chr names, per file nr C, cov
select_a = c(6,7) 

for (i in 1:nrow(id_list.chick)){
  id_name <- id_list.chick[i,1]
  myfile <- file.path(getwd(), paste0(id_name, "_methylation_df_forDSS",".csv"))
  col_names <- c("chr", "pos", "N", "X")
  assign(paste(id_name, "_methylation_df_forDSS",sep = ""),meth.chick.adjchr[,c(1,3,select_a)]) %>% setNames(col_names)%>%write.table(file=myfile, sep = "\t", row.names = FALSE)
  select_a = select_a + 3
}

#### Load in data with DSS for DMA ####
setwd("/home/NIOO.INT/rebeccash/2018_BS_manipulation/R_Scripts/R_raw_files/library.chicks/")
#set this to the directory that you just outputted DSS files into

file.list.chick.dss <-paste("/home/nioo/rebeccash/2018_BS_manipulation/R_Scripts/R_raw_files/library.chicks/",
                            list.files("/home/nioo/rebeccash/2018_BS_manipulation/R_Scripts/R_raw_files/library.chicks/",
                                       pattern = "_methylation_df_forDSS.csv"),sep="")

individuals.dss.chick<-c()
for (i in 1:length(file.list.chick.dss)){
  individuals.dss.chick<-c(strsplit(strsplit(file.list.chick.dss[i],"library.chicks/")[[1]][2],"_methylation_df_forDSS.csv")[[1]][1],individuals.dss.chick)
}

file.list.chick.b = list.files(pattern="*df_forDSS.csv*")
for (i in 1:length(file.list.chick.b)) assign(file.list.chick.b[i], fread(file.list.chick.b[i], select = c(1:4), col.names = c("chr", "pos", "N", "X")))

DSSfiles<-grep("forDSS.csv",names(.GlobalEnv),value=TRUE)
DSSfiles_list<-do.call("list",mget(DSSfiles))

BSobj.chick <- makeBSseqData(dat = DSSfiles_list, 
                             sampleNames = c(rev(individuals.dss.chick)))

str(BSobj.chick)

### DMA
chicks.females <- subset(chickdata, Mol_sex == "Female" & Role == "Chick")
chicks.males <- subset(chickdata, Mol_sex == "Male" & Role == "Chick")

chicks.females <- as.data.frame(chicks.females)
chicks.males <- as.data.frame(chicks.males)

#group 1 females, group 2 males
dmlTest.chick.sex <- DMLtest(BSobj.chick, group1 = chicks.females[,2], 
                             group2 = chicks.males[,2], 
                             smoothing = FALSE)

head(dmlTest.chick.sex)

# make a volcano plot for visualisation
dmlTest.chick.sex.df <- as.data.frame(dmlTest.chick.sex)
fc.chick.sex=100*dmlTest.chick.sex.df$diff
sig.chick.sex=-log10(dmlTest.chick.sex.df$pval)
dmlTest.chick.sex.df<-cbind(dmlTest.chick.sex.df,fc.chick.sex,sig.chick.sex)
dmlTest.chick.sex.df$thre<-as.factor(abs(fc.chick.sex) >= 10 & dmlTest.chick.sex.df$fdr < 0.05)
ggplot(data=dmlTest.chick.sex.df, aes(x=fc.chick.sex, y=sig.chick.sex, color = thre)) +
  geom_point(alpha=1, size=0.75) +
  theme(legend.position = "none") +
  labs(subtitle =  "Volcano plot showing DMS between females and
       males in combined libraries in chicks", y=expression("-log"^"10"*"(p-value)"), x = expression(Delta*"(% methylation level)")) + 
  scale_x_continuous(breaks=c(-40,-30,-20, -10,0, 10, 20,30,40)) +gtit+
  scale_color_manual(values = c("grey60","red")) +
  geom_hline(yintercept = -log10(0.05/nrow(dmlTest.chick.sex.df)),color="black",linetype="dashed",alpha=0.7) +
  geom_vline(xintercept = -10,color="black",linetype="dashed",alpha=0.7) +
  geom_vline(xintercept = 10,color="black",linetype="dashed",alpha=0.7)

#p-value distribution
ggplot(dmlTest.chick.sex.df, aes(x=pval)) + geom_histogram(bins = 20) + 
  labs(subtitle = "p-value distribution Females vs Males 
in Chicks") + gtit + xlab("p-values") + ylab("Frequency")

#qqplot p-values
qqplot.pvalues(dmlTest.chick.sex.df$pval, col.abline = "red", col.CB = "gray80", CB=TRUE, CB.level = 0.95) 
#little bit underdispersed

## Call DML and DMR 
dmls.chick <- callDML(dmlTest.chick.sex, p.threshold = 0.05, delta = 0.1) #calling sites
nrow(dmls.chick) #0 significant sites
dmrs.chick <- callDMR(dmlTest.chick.sex, p.threshold = 0.05, minlen = 50) #calling regions
nrow(dmrs.chick) #15 significant regions

#### Reshape methylkit dataframe to a long table to make a model ####

# In order to run a lmer, we need one row per CpG site per ID (= huge df)

methpattern.chick.df <- meth.chick.df #make sure it's arranged well and lines up with colnrs

startmeth.chick <- 4+length(meth.chick@sample.ids)*3+1
endmeth.chick <- 4+length(meth.chick@sample.ids)*4

methpos.chick <- methpattern.chick.df[,c(1,2,startmeth.chick:endmeth.chick)]
methpattern.long.chick <- melt(methpos.chick, measure.vars = c(3:ncol(methpos.chick)), variable.name = "ID", value.name = "Methylation") #measure vars = meth cols
methpattern.long.chick$ID <- gsub(methpattern.long.chick$ID, pattern = "methylation_", replacement = "ID_")

cpos.chick <- seq(from = 6, to = (4+length(meth.chick@sample.ids)*3-1), by = 3) #check from and 4 is correct
methpattern.cpos.chick <- methpattern.chick.df[,c(1,2,cpos.chick)]
methpattern.long.cpos.chick <- melt(methpattern.cpos.chick, measure.vars = c(3:ncol(methpattern.cpos.chick)), variable.name = "ID", value.name = "nC")
methpattern.long.cpos.chick$ID <- gsub(methpattern.long.cpos.chick$ID, pattern = "numCs", replacement = "ID_")

tpos.chick <- seq(from = 7, to = (4+length(meth.chick@sample.ids)*3), by = 3)
methpattern.tpos.chick <- methpattern.chick.df[,c(1,2,tpos.chick)]
methpattern.long.tpos.chick <- melt(methpattern.tpos.chick, measure.vars = c(3:ncol(methpattern.tpos.chick)), variable.name = "ID", value.name = "nT")
methpattern.long.tpos.chick$ID <- gsub(methpattern.long.tpos.chick$ID, pattern = "numTs", replacement = "ID_")

methpattern.long.chick.join <- left_join(methpattern.long.chick, methpattern.long.cpos.chick, by = c("chr", "start", "ID"))
methpattern.long.chick.join <- left_join(methpattern.long.chick.join, methpattern.long.tpos.chick, by = c("chr", "start", "ID"))

#make df with id names and 'methylation nr'
ids.chick.ordered <- as.data.frame(names(methpos.chick[-c(1,2)])) #meth cols
ids.chick.ordered <- cbind(ids.chick.ordered,meth.chick@sample.ids)

ids.chick.ordered <- as.data.frame(ids.chick.ordered)
names(ids.chick.ordered) <- c("ID", "RNR")
ids.chick.ordered$ID <- gsub(ids.chick.ordered$ID, pattern = "methylation_", replacement = "ID_")

#join long table to include treatment and sex
ids.chick.ordered <- left_join(ids.chick.ordered,chickdata[,c(2,3,9,7,4,5,10,8,12)], by = "RNR")

methpattern.total.chick <- left_join(methpattern.long.chick.join, ids.chick.ordered, by = c("ID" = "ID"))

methpattern.total.chick <- add_column(methpattern.total.chick, "CHR_POS" = paste(methpattern.total.chick$chr, methpattern.total.chick$start, sep = "_"), .before = "chr")
methpattern.total.chick$ID <- as.factor(methpattern.total.chick$ID)
methpattern.total.chick$RNR <- as.factor(methpattern.total.chick$RNR)
names(methpattern.total.chick$Mol_sex) <- "Sex"

#### ANOVA to test for differences per library, treatment and sex #####
#in addition to visualisation from PCA
anova(lm(Methylation ~Library, data = methpattern.total.chick))
tukey.library <- TukeyHSD(aov(lm(Methylation ~Library, data = methpattern.total.chick)))
summary(tukey.library$Library[,2]) #despite significance, only very minimal with mean of -0.0016 difference in methylation %

anova(lm(Methylation ~Treatment, data = methpattern.total.chick))
tukey.treatment <- TukeyHSD(aov(lm(Methylation ~Treatment, data = methpattern.total.chick)))
summary(tukey.treatment$Treatment[,2])

anova(lm(Methylation ~Sex, data = methpattern.total.chick))
tukey.sex <- TukeyHSD(aov(lm(Methylation ~Sex, data = methpattern.total.chick)))
summary(tukey.sex$Sex[,2])


