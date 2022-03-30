#### Gene annotation ####

#### Load libraries #####
library(genomation); library(GenomicFeatures); library(rtracklayer);library(fuzzyjoin)
library(tibble); library(dplyr)

#### Load in dataframes #####

## Set working directory ##
setwd("/home/nioo/rebeccash/2018_BS_manipulation/R_Scripts/complete_MS_pt1") #set this to the cloned github dir

#Significant for geneticNK
base::load("data/SigCpGGenetic.RData") #sites.sig.geneticnk.M1.select.adjchr
#see script 2.Modelling how to get SigCpGGenetic.RData

#Significant for fosterNK
base::load("data/SigCpGFoster.RData") #sites.sig.fosternk.M1.select.adjchr
#see script 2.Modelling how to get SigCpGFoster.RData

#chrnames
base::load(file = "data/reference_genome/GCF_001522545.3_Parus_major1.1_assembly_report_chr_names.RData")
#ref genome parus major v1.1 can be downloaded from https://www.ncbi.nlm.nih.gov/assembly/GCF_001522545.3

## Annotation data ##
promoter=gffToGRanges("data/reference_genome/promoters_2k.t.gff3")
promoter=unique(promoter)
genes=gffToGRanges("data/reference_genome/genes.gff3")
genes=unique(genes)
TSS=gffToGRanges("data/reference_genome/TSS.laine.t.gff3")
TSS=unique(TSS)
exons_gene=gffToGRanges("data/reference_genome/exons.gene.gff3")
exons_gene=unique(exons_gene)
#exons_transcript=gffToGRanges("data/reference_genome/exons.transcript.gff3")
#exons_transcript=unique(exons_transcript)
introns=gffToGRanges("data/reference_genome/introns.transcripts.gff3")
introns=unique(introns)
downstream=gffToGRanges("data/reference_genome/downstream.laine.t.gff3")
downstream=unique(downstream)
upstream=gffToGRanges("data/reference_genome/upstream.laine.t.gff3")
upstream=unique(upstream)
threeUTR =gffToGRanges("data/reference_genome/threeUTRs.transcripts.unlist.gff3")
threeUTR=unique(threeUTR)
fiveUTR=gffToGRanges("data/reference_genome/fiveUTRs.transcripts.unlist.gff3")
fiveUTR=unique(fiveUTR)

#gff3 files can be downloaded at XXXXX
# From Laine et al. 2016

#### Genetic NK LMER ####
nrow(sites.sig.geneticnk.M1.select) #8315
names(sites.sig.geneticnk.M1.select)[3] <- "pvalue" #to format for GRanges
sites.sig.geneticnk.M1.select$end <- sites.sig.geneticnk.M1.select$start
sites.sig.geneticnk.M1.select$strand <- "+"
sites.sig.geneticnk.M1.select.GR <- as(sites.sig.geneticnk.M1.select, "GRanges")

#Promoters
annotated_M1.GeneticNK.promo <- subsetByOverlaps(promoter, sites.sig.geneticnk.M1.select.GR)
annotated_M1.GeneticNK.promo.df <- as.data.frame(annotated_M1.GeneticNK.promo)
annotated_M1.GeneticNK.promo.df <- fuzzy_inner_join(
  annotated_M1.GeneticNK.promo.df,sites.sig.geneticnk.M1.select,
  by = c("seqnames" = "chr","start" = "start","end" = "start"),
  match_fun = list(`==`, `<=`, `>=`))

annotated_M1.GeneticNK.promo.df <- annotated_M1.GeneticNK.promo.df %>% 
  arrange(annotated_M1.GeneticNK.promo.df$pvalue)

annotated_M1.GeneticNK.promo.GO <- annotated_M1.GeneticNK.promo.df[,"gene_id"] #for GOrilla

#Genes
annotated_M1.GeneticNK.genes <- subsetByOverlaps(genes, sites.sig.geneticnk.M1.select.GR)
annotated_M1.GeneticNK.genes.df <- as.data.frame(annotated_M1.GeneticNK.genes)
annotated_M1.GeneticNK.genes.df <- fuzzy_inner_join(
  annotated_M1.GeneticNK.genes.df,sites.sig.geneticnk.M1.select,
  by = c("seqnames" = "chr","start" = "start","end" = "start"),
  match_fun = list(`==`, `<=`, `>=`))

annotated_M1.GeneticNK.genes.df <- annotated_M1.GeneticNK.genes.df %>% 
  arrange(annotated_M1.GeneticNK.genes.df$pvalue)

annotated_M1.GeneticNK.genes.GO <- annotated_M1.GeneticNK.genes.df[,"gene_id"]

#TSS
annotated_M1.GeneticNK.TSS <- subsetByOverlaps(TSS, sites.sig.geneticnk.M1.select.GR)
annotated_M1.GeneticNK.TSS.df <- as.data.frame(annotated_M1.GeneticNK.TSS)
annotated_M1.GeneticNK.TSS.df <- fuzzy_inner_join(
  annotated_M1.GeneticNK.TSS.df,sites.sig.geneticnk.M1.select,
  by = c("seqnames" = "chr","start" = "start","end" = "start"),
  match_fun = list(`==`, `<=`, `>=`))

annotated_M1.GeneticNK.TSS.df <- annotated_M1.GeneticNK.TSS.df %>% 
  arrange(annotated_M1.GeneticNK.TSS.df$pvalue)

annotated_M1.GeneticNK.TSS.GO <- annotated_M1.GeneticNK.TSS.df[,"gene_id"]

#exon gene
annotated_M1.GeneticNK.exon.gene <- subsetByOverlaps(exons_gene, sites.sig.geneticnk.M1.select.GR)
annotated_M1.GeneticNK.exon.gene.df <- as.data.frame(annotated_M1.GeneticNK.exon.gene)
annotated_M1.GeneticNK.exon.gene.df <- fuzzy_inner_join(
  annotated_M1.GeneticNK.exon.gene.df,sites.sig.geneticnk.M1.select,
  by = c("seqnames" = "chr","start" = "start","end" = "start"),
  match_fun = list(`==`, `<=`, `>=`))

annotated_M1.GeneticNK.exon.gene.df <- annotated_M1.GeneticNK.exon.gene.df %>% 
  arrange(annotated_M1.GeneticNK.exon.gene.df$pvalue)

#annotated_M1.GeneticNK.exon.gene.GO <- annotated_M1.GeneticNK.exon.gene.df[,"gene_id"]

#exon transcript
annotated_M1.GeneticNK.exon.trans <- subsetByOverlaps(exons_transcript, sites.sig.geneticnk.M1.select.GR)
annotated_M1.GeneticNK.exon.trans.df <- as.data.frame(annotated_M1.GeneticNK.exon.trans)
annotated_M1.GeneticNK.exon.trans.df <- fuzzy_inner_join(
  annotated_M1.GeneticNK.exon.trans.df,sites.sig.geneticnk.M1.select,
  by = c("seqnames" = "chr","start" = "start","end" = "start"),
  match_fun = list(`==`, `<=`, `>=`))

annotated_M1.GeneticNK.exon.trans.df <- annotated_M1.GeneticNK.exon.trans.df %>% 
  arrange(annotated_M1.GeneticNK.exon.trans.df$pvalue)

#annotated_M1.GeneticNK.exon.trans.GO <- annotated_M1.GeneticNK.exon.trans.df[,"gene_id"]

#intron
annotated_M1.GeneticNK.intron <- subsetByOverlaps(introns, sites.sig.geneticnk.M1.select.GR)
annotated_M1.GeneticNK.intron.df <- as.data.frame(annotated_M1.GeneticNK.intron)
annotated_M1.GeneticNK.intron.df <- fuzzy_inner_join(
  annotated_M1.GeneticNK.intron.df,sites.sig.geneticnk.M1.select,
  by = c("seqnames" = "chr","start" = "start","end" = "start"),
  match_fun = list(`==`, `<=`, `>=`))

annotated_M1.GeneticNK.intron.df <- annotated_M1.GeneticNK.intron.df %>% 
  arrange(annotated_M1.GeneticNK.intron.df$pvalue)

#annotated_M1.GeneticNK.intron.GO <- annotated_M1.GeneticNK.intron.df[,"gene_id"]

#upstream
annotated_M1.GeneticNK.upstream <- subsetByOverlaps(upstream, sites.sig.geneticnk.M1.select.GR)
annotated_M1.GeneticNK.upstream.df <- as.data.frame(annotated_M1.GeneticNK.upstream)
annotated_M1.GeneticNK.upstream.df <- fuzzy_inner_join(
  annotated_M1.GeneticNK.upstream.df,sites.sig.geneticnk.M1.select,
  by = c("seqnames" = "chr","start" = "start","end" = "start"),
  match_fun = list(`==`, `<=`, `>=`))

annotated_M1.GeneticNK.upstream.df <- annotated_M1.GeneticNK.upstream.df %>% 
  arrange(annotated_M1.GeneticNK.upstream.df$pvalue)

annotated_M1.GeneticNK.upstream.GO <- annotated_M1.GeneticNK.upstream.df[,"gene_id"]

#downstream
annotated_M1.GeneticNK.downstream <- subsetByOverlaps(downstream, sites.sig.geneticnk.M1.select.GR)
annotated_M1.GeneticNK.downstream.df <- as.data.frame(annotated_M1.GeneticNK.downstream)
annotated_M1.GeneticNK.downstream.df <- fuzzy_inner_join(
  annotated_M1.GeneticNK.downstream.df,sites.sig.geneticnk.M1.select,
  by = c("seqnames" = "chr","start" = "start","end" = "start"),
  match_fun = list(`==`, `<=`, `>=`))

annotated_M1.GeneticNK.downstream.df <- annotated_M1.GeneticNK.downstream.df %>% 
  arrange(annotated_M1.GeneticNK.downstream.df$pvalue)

annotated_M1.GeneticNK.downstream.GO <- annotated_M1.GeneticNK.downstream.df[,"gene_id"]

# Add together for GOrilla
genes.for.GOrilla.M1.GeneticNK<-rbind(annotated_M1.GeneticNK.promo.df[,c("gene_id", "pvalue")], 
                                      annotated_M1.GeneticNK.genes.df[,c("gene_id", "pvalue")],
                                      annotated_M1.GeneticNK.TSS.df[,c("gene_id", "pvalue")],
                                      annotated_M1.GeneticNK.upstream.df[,c("gene_id", "pvalue")],
                                      annotated_M1.GeneticNK.downstream.df[,c("gene_id", "pvalue")])
genes.for.GOrilla.M1.GeneticNK <- genes.for.GOrilla.M1.GeneticNK %>% 
  arrange(genes.for.GOrilla.M1.GeneticNK$pvalue) #arrange by increasing p-value for GOrilla

genes.for.GOrilla.M1.GeneticNK.GO <- genes.for.GOrilla.M1.GeneticNK$gene_id 
write.table(genes.for.GOrilla.M1.GeneticNK.GO, row.names = FALSE, col.names=FALSE,
            quote = FALSE, "data/genes.for.GOrilla.GeneticNK.txt")

#### Count sites sig for GeneticNK per regulatory region ####
#should make loop for this
annotated_M1.GeneticNK.promo.CG <- subsetByOverlaps(sites.sig.geneticnk.M1.select.GR, promoter)
annotated_M1.GeneticNK.genes.CG <- subsetByOverlaps(sites.sig.geneticnk.M1.select.GR, genes)
annotated_M1.GeneticNK.TSS.CG <- subsetByOverlaps(sites.sig.geneticnk.M1.select.GR, TSS)
annotated_M1.GeneticNK.exons_gene.CG <- subsetByOverlaps(sites.sig.geneticnk.M1.select.GR, exons_gene)
annotated_M1.GeneticNK.introns.CG <- subsetByOverlaps(sites.sig.geneticnk.M1.select.GR, introns)
annotated_M1.GeneticNK.upstream.CG <- subsetByOverlaps(sites.sig.geneticnk.M1.select.GR, upstream)
annotated_M1.GeneticNK.downstream.CG <- subsetByOverlaps(sites.sig.geneticnk.M1.select.GR, downstream)
annotated_M1.GeneticNK.threeUTR.CG <- subsetByOverlaps(sites.sig.geneticnk.M1.select.GR, threeUTR)
annotated_M1.GeneticNK.fiveUTR.CG <- subsetByOverlaps(sites.sig.geneticnk.M1.select.GR, fiveUTR)

annotated_M1.GeneticNK.unplaced <- sites.sig.geneticnk.M1.select
annotated_M1.GeneticNK.unplaced <- add_column(annotated_M1.GeneticNK.unplaced, "CHR_POS" = paste(annotated_M1.GeneticNK.unplaced$chr, annotated_M1.GeneticNK.unplaced$start, sep = "_"))

annotated_M1.GeneticNK.promo.CG <- as.data.frame(annotated_M1.GeneticNK.promo.CG)
annotated_M1.GeneticNK.genes.CG <- as.data.frame(annotated_M1.GeneticNK.genes.CG)
annotated_M1.GeneticNK.TSS.CG <- as.data.frame(annotated_M1.GeneticNK.TSS.CG)
annotated_M1.GeneticNK.exons_gene.CG <- as.data.frame(annotated_M1.GeneticNK.exons_gene.CG)
annotated_M1.GeneticNK.introns.CG <- as.data.frame(annotated_M1.GeneticNK.introns.CG)
annotated_M1.GeneticNK.upstream.CG <- as.data.frame(annotated_M1.GeneticNK.upstream.CG)
annotated_M1.GeneticNK.downstream.CG <- as.data.frame(annotated_M1.GeneticNK.downstream.CG)
annotated_M1.GeneticNK.threeUTR.CG <- as.data.frame(annotated_M1.GeneticNK.threeUTR.CG)
annotated_M1.GeneticNK.fiveUTR.CG <- as.data.frame(annotated_M1.GeneticNK.fiveUTR.CG)

annotated_M1.GeneticNK.promo.CG <- add_column(annotated_M1.GeneticNK.promo.CG, "CHR_POS" = paste(annotated_M1.GeneticNK.promo.CG$seqnames, annotated_M1.GeneticNK.promo.CG$start, sep = "_"))
annotated_M1.GeneticNK.genes.CG <- add_column(annotated_M1.GeneticNK.genes.CG, "CHR_POS" = paste(annotated_M1.GeneticNK.genes.CG$seqnames, annotated_M1.GeneticNK.genes.CG$start, sep = "_"))
annotated_M1.GeneticNK.TSS.CG <- add_column(annotated_M1.GeneticNK.TSS.CG, "CHR_POS" = paste(annotated_M1.GeneticNK.TSS.CG$seqnames, annotated_M1.GeneticNK.TSS.CG$start, sep = "_"))
annotated_M1.GeneticNK.exons_gene.CG <- add_column(annotated_M1.GeneticNK.exons_gene.CG, "CHR_POS" = paste(annotated_M1.GeneticNK.exons_gene.CG$seqnames, annotated_M1.GeneticNK.exons_gene.CG$start, sep = "_"))
annotated_M1.GeneticNK.introns.CG <- add_column(annotated_M1.GeneticNK.introns.CG, "CHR_POS" = paste(annotated_M1.GeneticNK.introns.CG$seqnames, annotated_M1.GeneticNK.introns.CG$start, sep = "_"))
annotated_M1.GeneticNK.upstream.CG <- add_column(annotated_M1.GeneticNK.upstream.CG, "CHR_POS" = paste(annotated_M1.GeneticNK.upstream.CG$seqnames, annotated_M1.GeneticNK.upstream.CG$start, sep = "_"))
annotated_M1.GeneticNK.downstream.CG <- add_column(annotated_M1.GeneticNK.downstream.CG, "CHR_POS" = paste(annotated_M1.GeneticNK.downstream.CG$seqnames, annotated_M1.GeneticNK.downstream.CG$start, sep = "_"))
annotated_M1.GeneticNK.threeUTR.CG <- add_column(annotated_M1.GeneticNK.threeUTR.CG, "CHR_POS" = paste(annotated_M1.GeneticNK.threeUTR.CG$seqnames, annotated_M1.GeneticNK.threeUTR.CG$start, sep = "_"))
annotated_M1.GeneticNK.fiveUTR.CG <- add_column(annotated_M1.GeneticNK.fiveUTR.CG, "CHR_POS" = paste(annotated_M1.GeneticNK.fiveUTR.CG$seqnames, annotated_M1.GeneticNK.fiveUTR.CG$start, sep = "_"))

annotated_M1.GeneticNK.unplaced <- subset(annotated_M1.GeneticNK.unplaced, !CHR_POS %in% annotated_M1.GeneticNK.promo.CG$CHR_POS)
annotated_M1.GeneticNK.unplaced <- subset(annotated_M1.GeneticNK.unplaced, !CHR_POS %in% annotated_M1.GeneticNK.genes.CG$CHR_POS)
annotated_M1.GeneticNK.unplaced <- subset(annotated_M1.GeneticNK.unplaced, !CHR_POS %in% annotated_M1.GeneticNK.TSS.CG$CHR_POS)
annotated_M1.GeneticNK.unplaced <- subset(annotated_M1.GeneticNK.unplaced, !CHR_POS %in% annotated_M1.GeneticNK.exons_gene.CG$CHR_POS)
annotated_M1.GeneticNK.unplaced <- subset(annotated_M1.GeneticNK.unplaced, !CHR_POS %in% annotated_M1.GeneticNK.introns.CG$CHR_POS)
annotated_M1.GeneticNK.unplaced <- subset(annotated_M1.GeneticNK.unplaced, !CHR_POS %in% annotated_M1.GeneticNK.upstream.CG$CHR_POS)
annotated_M1.GeneticNK.unplaced <- subset(annotated_M1.GeneticNK.unplaced, !CHR_POS %in% annotated_M1.GeneticNK.downstream.CG$CHR_POS)
annotated_M1.GeneticNK.unplaced <- subset(annotated_M1.GeneticNK.unplaced, !CHR_POS %in% annotated_M1.GeneticNK.threeUTR.CG$CHR_POS)
annotated_M1.GeneticNK.unplaced <- subset(annotated_M1.GeneticNK.unplaced, !CHR_POS %in% annotated_M1.GeneticNK.fiveUTR.CG$CHR_POS)

distribution.M1.geneticNK <- data.frame("Functional.Region" = c("Promoter", "Genes", "TSS", "Exon", 
                                                                "Introns", "Upstream", "Downstream", 
                                                                "3-UTR", "5-UTR", "Unplaced", "All"), 
                                        "Nr.Sites" = c(nrow(annotated_M1.GeneticNK.promo.CG),
                                                       nrow(annotated_M1.GeneticNK.genes.CG),
                                                       nrow(annotated_M1.GeneticNK.TSS.CG),
                                                       nrow(annotated_M1.GeneticNK.exons_gene.CG),
                                                       nrow(annotated_M1.GeneticNK.introns.CG),
                                                       nrow(annotated_M1.GeneticNK.upstream.CG),
                                                       nrow(annotated_M1.GeneticNK.downstream.CG),
                                                       nrow(annotated_M1.GeneticNK.threeUTR.CG),
                                                       nrow(annotated_M1.GeneticNK.fiveUTR.CG),
                                                       length(unique(annotated_M1.GeneticNK.unplaced$CHR_POS)), 
                                                       nrow(sites.sig.geneticnk.M1.select)))

# make pie chart
distribution.M1Genetic.ggplot <- distribution.M1.geneticNK[-11,]
str(distribution.M1Genetic.ggplot)
distribution.M1Genetic.ggplot$Nr.Sites <- as.numeric(distribution.M1Genetic.ggplot$Nr.Sites)
distribution.M1Genetic.ggplot <- distribution.M1Genetic.ggplot %>% 
  arrange(Nr.Sites) %>%
  mutate(prop = Nr.Sites / sum(distribution.M1Genetic.ggplot$Nr.Sites) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )%>% 
  mutate(csum = rev(cumsum(rev(Nr.Sites))), 
         pos = Nr.Sites/2 + lead(csum, 1),
         pos = if_else(is.na(pos), Nr.Sites/2, pos))

library(ggplot2)
library(ggrepel)
library(tidyverse)

distribution.M1Genetic.ggplot[10,6] <- 2000

ggplot(distribution.M1Genetic.ggplot, aes(x = "" , y = Nr.Sites, fill = fct_inorder(Functional.Region))) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = c('#543005','#8c510a','#bf812d','#dfc27d','#f6e8c3','#c7eae5','#80cdc1','#35978f','#01665e','#003c30')) +
  geom_label_repel(data = distribution.M1Genetic.ggplot,
                   aes(y = pos, label = paste0(Nr.Sites)),
                   size = 8, nudge_x = 1, show.legend = FALSE, fill = "white", parse = TRUE) +
  guides(fill = guide_legend(title = "Functional region", reverse = T, label.position = c("left"))) +
  theme_void()+theme(text = element_text(size = 30))

#### Foster NK LMER ####
nrow(sites.sig.fosternk.M1.select) #101
names(sites.sig.fosternk.M1.select)[3] <- "pvalue"
names(sites.sig.fosternk.M1.select)[1] <- "chr"
names(sites.sig.fosternk.M1.select)[2] <- "start"
sites.sig.fosternk.M1.select$end <- sites.sig.fosternk.M1.select$start
sites.sig.fosternk.M1.select$strand <- "+"
sites.sig.fosternk.M1.select.GR <- as(sites.sig.fosternk.M1.select, "GRanges")

#Promoters
annotated_M1.FosterNK.promo <- subsetByOverlaps(promoter, sites.sig.fosternk.M1.select.GR)
annotated_M1.FosterNK.promo.df <- as.data.frame(annotated_M1.FosterNK.promo)
annotated_M1.FosterNK.promo.df <- fuzzy_inner_join(
  annotated_M1.FosterNK.promo.df,sites.sig.fosternk.M1.select,
  by = c("seqnames" = "chr","start" = "start","end" = "start"),
  match_fun = list(`==`, `<=`, `>=`))

annotated_M1.FosterNK.promo.df <- annotated_M1.FosterNK.promo.df %>% 
  arrange(annotated_M1.FosterNK.promo.df$pvalue)

annotated_M1.FosterNK.promo.GO <- annotated_M1.FosterNK.promo.df[,"gene_id"]

#Genes
annotated_M1.FosterNK.genes <- subsetByOverlaps(genes, sites.sig.fosternk.M1.select.GR)
annotated_M1.FosterNK.genes.df <- as.data.frame(annotated_M1.FosterNK.genes)
annotated_M1.FosterNK.genes.df <- fuzzy_inner_join(
  annotated_M1.FosterNK.genes.df,sites.sig.fosternk.M1.select,
  by = c("seqnames" = "chr","start" = "start","end" = "start"),
  match_fun = list(`==`, `<=`, `>=`))

annotated_M1.FosterNK.genes.df <- annotated_M1.FosterNK.genes.df %>% 
  arrange(annotated_M1.FosterNK.genes.df$pvalue)

annotated_M1.FosterNK.genes.GO <- annotated_M1.FosterNK.genes.df[,"gene_id"]

#TSS
annotated_M1.FosterNK.TSS <- subsetByOverlaps(TSS, sites.sig.fosternk.M1.select.GR)
annotated_M1.FosterNK.TSS.df <- as.data.frame(annotated_M1.FosterNK.TSS)
annotated_M1.FosterNK.TSS.df <- fuzzy_inner_join(
  annotated_M1.FosterNK.TSS.df,sites.sig.fosternk.M1.select,
  by = c("seqnames" = "chr","start" = "start","end" = "start"),
  match_fun = list(`==`, `<=`, `>=`))

annotated_M1.FosterNK.TSS.df <- annotated_M1.FosterNK.TSS.df %>% 
  arrange(annotated_M1.FosterNK.TSS.df$pvalue)

annotated_M1.FosterNK.TSS.GO <- annotated_M1.FosterNK.TSS.df[,"gene_id"]

#exon gene
annotated_M1.FosterNK.exon.gene <- subsetByOverlaps(exons_gene, sites.sig.fosternk.M1.select.GR)
annotated_M1.FosterNK.exon.gene.df <- as.data.frame(annotated_M1.FosterNK.exon.gene)
annotated_M1.FosterNK.exon.gene.df <- fuzzy_inner_join(
  annotated_M1.FosterNK.exon.gene.df,sites.sig.fosternk.M1.select,
  by = c("seqnames" = "chr","start" = "start","end" = "start"),
  match_fun = list(`==`, `<=`, `>=`))

annotated_M1.FosterNK.exon.gene.df <- annotated_M1.FosterNK.exon.gene.df %>% 
  arrange(annotated_M1.FosterNK.exon.gene.df$pvalue)

#annotated_M1.FosterNK.exon.gene.GO <- annotated_M1.FosterNK.exon.gene.df[,"gene_id"]

#exon transcript
annotated_M1.FosterNK.exon.trans <- subsetByOverlaps(exons_transcript, sites.sig.fosternk.M1.select.GR)
annotated_M1.FosterNK.exon.trans.df <- as.data.frame(annotated_M1.FosterNK.exon.trans)
annotated_M1.FosterNK.exon.trans.df <- fuzzy_inner_join(
  annotated_M1.FosterNK.exon.trans.df,sites.sig.fosternk.M1.select,
  by = c("seqnames" = "chr","start" = "start","end" = "start"),
  match_fun = list(`==`, `<=`, `>=`))

annotated_M1.FosterNK.exon.trans.df <- annotated_M1.FosterNK.exon.trans.df %>% 
  arrange(annotated_M1.FosterNK.exon.trans.df$pvalue)

#annotated_M1.FosterNK.exon.trans.GO <- annotated_M1.FosterNK.exon.trans.df[,"gene_id"]

#intron
annotated_M1.FosterNK.intron <- subsetByOverlaps(introns, sites.sig.fosternk.M1.select.GR)
annotated_M1.FosterNK.intron.df <- as.data.frame(annotated_M1.FosterNK.intron)
annotated_M1.FosterNK.intron.df <- fuzzy_inner_join(
  annotated_M1.FosterNK.intron.df,sites.sig.fosternk.M1.select,
  by = c("seqnames" = "chr","start" = "start","end" = "start"),
  match_fun = list(`==`, `<=`, `>=`))

annotated_M1.FosterNK.intron.df <- annotated_M1.FosterNK.intron.df %>% 
  arrange(annotated_M1.FosterNK.intron.df$pvalue)

#annotated_M1.FosterNK.intron.GO <- annotated_M1.FosterNK.intron.df[,"gene_id"]

#upstream
annotated_M1.FosterNK.upstream <- subsetByOverlaps(upstream, sites.sig.fosternk.M1.select.GR)
annotated_M1.FosterNK.upstream.df <- as.data.frame(annotated_M1.FosterNK.upstream)
annotated_M1.FosterNK.upstream.df <- fuzzy_inner_join(
  annotated_M1.FosterNK.upstream.df,sites.sig.fosternk.M1.select,
  by = c("seqnames" = "chr","start" = "start","end" = "start"),
  match_fun = list(`==`, `<=`, `>=`))

annotated_M1.FosterNK.upstream.df <- annotated_M1.FosterNK.upstream.df %>% 
  arrange(annotated_M1.FosterNK.upstream.df$pvalue)

annotated_M1.FosterNK.upstream.GO <- annotated_M1.FosterNK.upstream.df[,"gene_id"]

#downstream
annotated_M1.FosterNK.downstream <- subsetByOverlaps(downstream, sites.sig.fosternk.M1.select.GR)
annotated_M1.FosterNK.downstream.df <- as.data.frame(annotated_M1.FosterNK.downstream)
annotated_M1.FosterNK.downstream.df <- fuzzy_inner_join(
  annotated_M1.FosterNK.downstream.df,sites.sig.fosternk.M1.select,
  by = c("seqnames" = "chr","start" = "start","end" = "start"),
  match_fun = list(`==`, `<=`, `>=`))

annotated_M1.FosterNK.downstream.df <- annotated_M1.FosterNK.downstream.df %>% 
  arrange(annotated_M1.FosterNK.downstream.df$pvalue)

annotated_M1.FosterNK.downstream.GO <- annotated_M1.FosterNK.downstream.df[,"gene_id"]

# Add together for GOrilla
genes.for.GOrilla.M1.FosterNK<-rbind(annotated_M1.FosterNK.promo.df[,c("gene_id", "pvalue")], 
                                     annotated_M1.FosterNK.genes.df[,c("gene_id", "pvalue")],
                                     annotated_M1.FosterNK.TSS.df[,c("gene_id", "pvalue")],
                                     annotated_M1.FosterNK.upstream.df[,c("gene_id", "pvalue")],
                                     annotated_M1.FosterNK.downstream.df[,c("gene_id", "pvalue")])
genes.for.GOrilla.M1.FosterNK <- genes.for.GOrilla.M1.FosterNK %>% 
  arrange(genes.for.GOrilla.M1.FosterNK$pvalue)

genes.for.GOrilla.M1.FosterNK.GO <- genes.for.GOrilla.M1.FosterNK$gene_id 
write.table(genes.for.GOrilla.M1.FosterNK.GO, row.names = FALSE, col.names=FALSE,
            quote = FALSE, "data/genes.for.GOrilla.FosterNK.txt")

#### Count sites sig for FosterNK per regulatory region  ####
annotated_M1.fosternk.promo.CG <- subsetByOverlaps(sites.sig.fosternk.M1.select.GR, promoter)
annotated_M1.fosternk.genes.CG <- subsetByOverlaps(sites.sig.fosternk.M1.select.GR, genes)
annotated_M1.fosternk.TSS.CG <- subsetByOverlaps(sites.sig.fosternk.M1.select.GR, TSS)
annotated_M1.fosternk.exons_gene.CG <- subsetByOverlaps(sites.sig.fosternk.M1.select.GR, exons_gene)
annotated_M1.fosternk.introns.CG <- subsetByOverlaps(sites.sig.fosternk.M1.select.GR, introns)
annotated_M1.fosternk.upstream.CG <- subsetByOverlaps(sites.sig.fosternk.M1.select.GR, upstream)
annotated_M1.fosternk.downstream.CG <- subsetByOverlaps(sites.sig.fosternk.M1.select.GR, downstream)
annotated_M1.fosternk.threeUTR.CG <- subsetByOverlaps(sites.sig.fosternk.M1.select.GR, threeUTR)
annotated_M1.fosternk.fiveUTR.CG <- subsetByOverlaps(sites.sig.fosternk.M1.select.GR, fiveUTR)

annotated_M1.fosternk.unplaced <- sites.sig.fosternk.M1.select
annotated_M1.fosternk.unplaced <- add_column(annotated_M1.fosternk.unplaced, "CHR_POS" = paste(annotated_M1.fosternk.unplaced$chr, annotated_M1.fosternk.unplaced$start, sep = "_"))

annotated_M1.fosternk.promo.CG <- as.data.frame(annotated_M1.fosternk.promo.CG)
annotated_M1.fosternk.genes.CG <- as.data.frame(annotated_M1.fosternk.genes.CG)
annotated_M1.fosternk.TSS.CG <- as.data.frame(annotated_M1.fosternk.TSS.CG)
annotated_M1.fosternk.exons_gene.CG <- as.data.frame(annotated_M1.fosternk.exons_gene.CG)
annotated_M1.fosternk.introns.CG <- as.data.frame(annotated_M1.fosternk.introns.CG)
annotated_M1.fosternk.upstream.CG <- as.data.frame(annotated_M1.fosternk.upstream.CG)
annotated_M1.fosternk.downstream.CG <- as.data.frame(annotated_M1.fosternk.downstream.CG)
annotated_M1.fosternk.threeUTR.CG <- as.data.frame(annotated_M1.fosternk.threeUTR.CG)
annotated_M1.fosternk.fiveUTR.CG <- as.data.frame(annotated_M1.fosternk.fiveUTR.CG)

annotated_M1.fosternk.promo.CG <- add_column(annotated_M1.fosternk.promo.CG, "CHR_POS" = paste(annotated_M1.fosternk.promo.CG$seqnames, annotated_M1.fosternk.promo.CG$start, sep = "_"))
annotated_M1.fosternk.genes.CG <- add_column(annotated_M1.fosternk.genes.CG, "CHR_POS" = paste(annotated_M1.fosternk.genes.CG$seqnames, annotated_M1.fosternk.genes.CG$start, sep = "_"))
annotated_M1.fosternk.TSS.CG <- add_column(annotated_M1.fosternk.TSS.CG, "CHR_POS" = paste(annotated_M1.fosternk.TSS.CG$seqnames, annotated_M1.fosternk.TSS.CG$start, sep = "_"))
annotated_M1.fosternk.exons_gene.CG <- add_column(annotated_M1.fosternk.exons_gene.CG, "CHR_POS" = paste(annotated_M1.fosternk.exons_gene.CG$seqnames, annotated_M1.fosternk.exons_gene.CG$start, sep = "_"))
annotated_M1.fosternk.introns.CG <- add_column(annotated_M1.fosternk.introns.CG, "CHR_POS" = paste(annotated_M1.fosternk.introns.CG$seqnames, annotated_M1.fosternk.introns.CG$start, sep = "_"))
annotated_M1.fosternk.upstream.CG <- add_column(annotated_M1.fosternk.upstream.CG, "CHR_POS" = paste(annotated_M1.fosternk.upstream.CG$seqnames, annotated_M1.fosternk.upstream.CG$start, sep = "_"))
annotated_M1.fosternk.downstream.CG <- add_column(annotated_M1.fosternk.downstream.CG, "CHR_POS" = paste(annotated_M1.fosternk.downstream.CG$seqnames, annotated_M1.fosternk.downstream.CG$start, sep = "_"))
annotated_M1.fosternk.threeUTR.CG <- add_column(annotated_M1.fosternk.threeUTR.CG, "CHR_POS" = paste(annotated_M1.fosternk.threeUTR.CG$seqnames, annotated_M1.fosternk.threeUTR.CG$start, sep = "_"))
annotated_M1.fosternk.fiveUTR.CG <- add_column(annotated_M1.fosternk.fiveUTR.CG, "CHR_POS" = paste(annotated_M1.fosternk.fiveUTR.CG$seqnames, annotated_M1.fosternk.fiveUTR.CG$start, sep = "_"))

annotated_M1.fosternk.unplaced <- subset(annotated_M1.fosternk.unplaced, !CHR_POS %in% annotated_M1.fosternk.promo.CG$CHR_POS)
annotated_M1.fosternk.unplaced <- subset(annotated_M1.fosternk.unplaced, !CHR_POS %in% annotated_M1.fosternk.genes.CG$CHR_POS)
annotated_M1.fosternk.unplaced <- subset(annotated_M1.fosternk.unplaced, !CHR_POS %in% annotated_M1.fosternk.TSS.CG$CHR_POS)
annotated_M1.fosternk.unplaced <- subset(annotated_M1.fosternk.unplaced, !CHR_POS %in% annotated_M1.fosternk.exons_gene.CG$CHR_POS)
annotated_M1.fosternk.unplaced <- subset(annotated_M1.fosternk.unplaced, !CHR_POS %in% annotated_M1.fosternk.introns.CG$CHR_POS)
annotated_M1.fosternk.unplaced <- subset(annotated_M1.fosternk.unplaced, !CHR_POS %in% annotated_M1.fosternk.upstream.CG$CHR_POS)
annotated_M1.fosternk.unplaced <- subset(annotated_M1.fosternk.unplaced, !CHR_POS %in% annotated_M1.fosternk.downstream.CG$CHR_POS)
annotated_M1.fosternk.unplaced <- subset(annotated_M1.fosternk.unplaced, !CHR_POS %in% annotated_M1.fosternk.threeUTR.CG$CHR_POS)
annotated_M1.fosternk.unplaced <- subset(annotated_M1.fosternk.unplaced, !CHR_POS %in% annotated_M1.fosternk.fiveUTR.CG$CHR_POS)

distribution.M1.fosternk <- data.frame("Functional.Region" = c("Promoter", "Genes", "TSS", "Exon", 
                                                               "Introns", "Upstream", "Downstream", 
                                                               "3-UTR", "5-UTR", "Unplaced", "All"), 
                                       "Nr.Sites" = c(nrow(annotated_M1.fosternk.promo.CG),
                                                      nrow(annotated_M1.fosternk.genes.CG),
                                                      nrow(annotated_M1.fosternk.TSS.CG),
                                                      nrow(annotated_M1.fosternk.exons_gene.CG),
                                                      nrow(annotated_M1.fosternk.introns.CG),
                                                      nrow(annotated_M1.fosternk.upstream.CG),
                                                      nrow(annotated_M1.fosternk.downstream.CG),
                                                      nrow(annotated_M1.fosternk.threeUTR.CG),
                                                      nrow(annotated_M1.fosternk.fiveUTR.CG),
                                                      length(unique(annotated_M1.fosternk.unplaced$CHR_POS)), 
                                                      nrow(sites.sig.fosternk.M1.select)))

# make pie chart
distribution.M1Foster.ggplot <- distribution.M1.fosternk[-11,]
str(distribution.M1Foster.ggplot)
distribution.M1Foster.ggplot$Nr.Sites <- as.numeric(distribution.M1Foster.ggplot$Nr.Sites)
distribution.M1Foster.ggplot <- distribution.M1Foster.ggplot %>% 
  arrange(Nr.Sites) %>%
  mutate(prop = Nr.Sites / sum(distribution.M1Foster.ggplot$Nr.Sites) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )%>% 
  mutate(csum = rev(cumsum(rev(Nr.Sites))), 
         pos = Nr.Sites/2 + lead(csum, 1),
         pos = if_else(is.na(pos), Nr.Sites/2, pos))

distribution.M1Foster.ggplot[10,6] <- 25

ggplot(distribution.M1Foster.ggplot, aes(x = "" , y = Nr.Sites, fill = fct_inorder(Functional.Region))) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = c('#543005','#8c510a','#bf812d','#dfc27d','#f6e8c3','#c7eae5','#80cdc1','#35978f','#01665e','#003c30')) +
  geom_label_repel(data = distribution.M1Foster.ggplot,
                   aes(y = pos, label = paste0(Nr.Sites)),
                   size = 8, nudge_x = 1, show.legend = FALSE, fill = "white",col = "black", parse = TRUE) +
  guides(fill = guide_legend(title = "Functional region", reverse = T)) +
  theme_void()+theme(text = element_text(size = 30))

#### Genetic + Foster NK LMER ####
base::load("data/SigCpGGeneticFoster.RData") #sites.sig.fosternk.M1.select.adjchr

nrow(sites.sig.both.M1) #25
sites.sig.both.M1$end <- sites.sig.both.M1$start
sites.sig.both.M1$strand <- "+"
sites.sig.both.M1.GR <- as(sites.sig.both.M1, "GRanges")

#Promoters
annotated_M1.both.promo <- subsetByOverlaps(promoter, sites.sig.both.M1.GR)
annotated_M1.both.promo.df <- as.data.frame(annotated_M1.both.promo)
annotated_M1.both.promo.df <- fuzzy_inner_join(
  annotated_M1.both.promo.df,sites.sig.both.M1,
  by = c("seqnames" = "chr","start" = "start","end" = "start"),
  match_fun = list(`==`, `<=`, `>=`))

#Genes
annotated_M1.both.genes <- subsetByOverlaps(genes, sites.sig.both.M1.GR)
annotated_M1.both.genes.df <- as.data.frame(annotated_M1.both.genes)
annotated_M1.both.genes.df <- fuzzy_inner_join(
  annotated_M1.both.genes.df,sites.sig.both.M1,
  by = c("seqnames" = "chr","start" = "start","end" = "start"),
  match_fun = list(`==`, `<=`, `>=`))

#TSS
annotated_M1.both.TSS <- subsetByOverlaps(TSS, sites.sig.both.M1.GR)
annotated_M1.both.TSS.df <- as.data.frame(annotated_M1.both.TSS)
annotated_M1.both.TSS.df <- fuzzy_inner_join(
  annotated_M1.both.TSS.df,sites.sig.both.M1,
  by = c("seqnames" = "chr","start" = "start","end" = "start"),
  match_fun = list(`==`, `<=`, `>=`))

#exon gene
annotated_M1.both.exon.gene <- subsetByOverlaps(exons_gene, sites.sig.both.M1.GR)
annotated_M1.both.exon.gene.df <- as.data.frame(annotated_M1.both.exon.gene)
annotated_M1.both.exon.gene.df <- fuzzy_inner_join(
  annotated_M1.both.exon.gene.df,sites.sig.both.M1,
  by = c("seqnames" = "chr","start" = "start","end" = "start"),
  match_fun = list(`==`, `<=`, `>=`))

#exon transcript
annotated_M1.both.exon.trans <- subsetByOverlaps(exons_transcript, sites.sig.both.M1.GR)
annotated_M1.both.exon.trans.df <- as.data.frame(annotated_M1.both.exon.trans)
annotated_M1.both.exon.trans.df <- fuzzy_inner_join(
  annotated_M1.both.exon.trans.df,sites.sig.both.M1,
  by = c("seqnames" = "chr","start" = "start","end" = "start"),
  match_fun = list(`==`, `<=`, `>=`))

#intron
annotated_M1.both.intron <- subsetByOverlaps(introns, sites.sig.both.M1.GR)
annotated_M1.both.intron.df <- as.data.frame(annotated_M1.both.intron)
annotated_M1.both.intron.df <- fuzzy_inner_join(
  annotated_M1.both.intron.df,sites.sig.both.M1,
  by = c("seqnames" = "chr","start" = "start","end" = "start"),
  match_fun = list(`==`, `<=`, `>=`))

#upstream
annotated_M1.both.upstream <- subsetByOverlaps(upstream, sites.sig.both.M1.GR)
annotated_M1.both.upstream.df <- as.data.frame(annotated_M1.both.upstream)
annotated_M1.both.upstream.df <- fuzzy_inner_join(
  annotated_M1.both.upstream.df,sites.sig.both.M1,
  by = c("seqnames" = "chr","start" = "start","end" = "start"),
  match_fun = list(`==`, `<=`, `>=`))

#downstream
annotated_M1.both.downstream <- subsetByOverlaps(downstream, sites.sig.both.M1.GR)
annotated_M1.both.downstream.df <- as.data.frame(annotated_M1.both.downstream)
annotated_M1.both.downstream.df <- fuzzy_inner_join(
  annotated_M1.both.downstream.df,sites.sig.both.M1,
  by = c("seqnames" = "chr","start" = "start","end" = "start"),
  match_fun = list(`==`, `<=`, `>=`))

