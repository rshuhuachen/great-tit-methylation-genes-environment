#Annotation Parus major - Bernice Sepers (NIOO-KNAW)
#version 102: https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/9157/102/

#Below example of annotation CpG sites in TSS region
#TSS regions were defined as 300 bp upstream to 50 bp downstream of the annotated transcription 
#starting position of each gene (Laine et al., 2016 Nat. Commun.; Viitaniemi et al., 2019 Genome Biol. Evol.)

#resolved package dependencies
install.packages("remotes")
remotes::install_github("jasongraf1/JGmisc")
library(JGmisc)
detachAllPackages()

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("GenomicFeatures", "rtracklayer"))
library(GenomicFeatures)
library(rtracklayer)

BiocManager::install(c("BSgenome", "BiocGenerics"))
library(BiocGenerics)
library(BSgenome)

## --------------- load genome info --------------------------------------- ###

# load chr/sc length; !!! MT length not updated in new release - added it manually to chr.length.txt (NC_040875 16777) (https://www.ncbi.nlm.nih.gov/nuccore/NC_040875.1/) !!!
ParusM_1.1_genome <- read.table("chr.length.txt", header=FALSE, sep="\t") 
head(ParusM_1.1_genome)
#names(ParusM_1.1_genome) <- c("chrom","length")
colnames(ParusM_1.1_genome)<-c("chrom","length")
head(ParusM_1.1_genome)
# load PM .gff-file
gff <- makeTxDbFromGFF("GCF_001522545.3_Parus_major1.1_genomic.gff", format="gff3", organism="Parus major", chrominfo=ParusM_1.1_genome) 
## The 6 first orphan exons at least are immunological gene (C and V gene segments) so they might be a result of duplication or not

## --------------- get annotations ---------------------------------------- ###
TSS.laine <- promoters(gff, upstream=300, downstream=50, columns=c("tx_name", "gene_id")) # TSS as in Laine et al., 2016. Nature Communications
TSS.laine.t <- trim(TSS.laine, use.names=TRUE)

## --------------- export files ------------------------------------------- ###
export(TSS.laine.t, "TSS.laine.t.gff3")

####------Example annotation CpGs-----####
LMM_delta<-readRDS(file="Significant_CpG_sites.Rdata")

library(rtracklayer)
#BiocManager::install("GenomicFeatures")
library(GenomicFeatures)
#BiocManager::install("genomation")
library(genomation)
#BiocManager::install("GenomicRanges")
library(GenomicRanges)

#change object object into GRanges and change strand into both strands annotation
LMM_delta$strand="*"
LMM_delta$start <- LMM_delta$end
LMM_delta <- LMM_delta %>% relocate(CHR, .before = df)
LMM_delta <- LMM_delta %>% relocate(start, .before = df)
LMM_delta <- LMM_delta %>% relocate(end, .before = df)
LMM_delta <- LMM_delta %>% relocate(strand, .before = df)
LMM_delta_GR=as(LMM_delta, "GRanges")

##read the gene regions gff3 file

#TSS
TSS=gffToGRanges("TSS.laine.t.gff3")
TSS=unique(TSS)
#get table with "mapped" genes/transcripts(can be double)/locgenes
diff_TSS_gff=subsetByOverlaps(TSS, LMM_delta_GR)
#get table with "mapped" sites
diff_TSS_CG_gff=subsetByOverlaps(LMM_delta_GR, TSS)
write.table(x=diff_TSS_gff, file='diff_TSS_gff_delta',quote=FALSE)
write.table(x=diff_TSS_CG_gff, file='diff_TSS_CG_gff_delta',quote=FALSE)

#-----------Finally, associate sites in diff_TSS_CG_gff_delta to diff_TSS_gff_delta---------------

