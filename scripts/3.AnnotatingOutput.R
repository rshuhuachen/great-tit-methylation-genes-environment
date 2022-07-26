#### Gene annotation ####

#### Load libraries #####
library(genomation); library(GenomicFeatures); library(rtracklayer);library(fuzzyjoin)
library(tibble); library(dplyr); library(GenomicRanges)

#### Load in dataframes #####

## Set working directory ##
#setwd("/home/nioo/rebeccash/2018_BS_manipulation/R_Scripts/complete_MS_pt1/great-tit-methylation-genes-environment/") #set this to the cloned github dir

#All sites for background list
base::load("data/ModelOutput.RData") #complete model output

#Significant for geneticNK
base::load("data/SigCpGGenetic.RData") #sites.sig.geneticnk.M1.select.adjchr
base::load("data/SigCpGGenetic_overdisp.RData") #sites.sig.geneticnk.M1.select.overdisp.adjchr
#see script 2.Modelling how to get SigCpGGenetic.RData

#Significant for fosterNK
base::load("data/SigCpGFoster.RData") #sites.sig.fosternk.M1.select.adjchr
base::load("data/SigCpGFoster_overdisp.RData") #sites.sig.fosternk.M1.select.overdisp.adjchr
#see script 2.Modelling how to get SigCpGFoster.RData

#chrnames
base::load("~/2018_BS_manipulation/reference_genome/GCF_001522545.3_Parus_major1.1_assembly_report_chr_names.RData")
base::load(file = "/home/nioo/rebeccash/2018_BS_manipulation/R_Scripts/R_raw_files/GCF_001522545.3_Parus_major1.1_assembly_report_chr_names.RData")
#ref genome parus major v1.1 can be downloaded from https://www.ncbi.nlm.nih.gov/assembly/GCF_001522545.3

## Annotation data ##
promoter=gffToGRanges("~/2018_BS_manipulation/reference_genome/promoters_2k.t.gff3")
promoter=unique(promoter)
genes=gffToGRanges("~/2018_BS_manipulation/reference_genome/genes.gff3")
genes=unique(genes)
TSS=gffToGRanges("~/2018_BS_manipulation/reference_genome/TSS.laine.t.gff3")
TSS=unique(TSS)
exons_gene=gffToGRanges("~/2018_BS_manipulation/reference_genome/exons.gene.gff3")
exons_gene=unique(exons_gene)
exons_transcript=gffToGRanges("~/2018_BS_manipulation/reference_genome/exons.transcript.gff3")
exons_transcript=unique(exons_transcript)
introns=gffToGRanges("~/2018_BS_manipulation/reference_genome/introns.transcripts.gff3")
introns=unique(introns)
downstream=gffToGRanges("~/2018_BS_manipulation/reference_genome/downstream.laine.t.gff3")
downstream=unique(downstream)
upstream=gffToGRanges("~/2018_BS_manipulation/reference_genome/upstream.laine.t.gff3")
upstream=unique(upstream)
threeUTR =gffToGRanges("~/2018_BS_manipulation/reference_genome/threeUTRs.transcripts.unlist.gff3")
threeUTR=unique(threeUTR)
fiveUTR=gffToGRanges("~/2018_BS_manipulation/reference_genome/fiveUTRs.transcripts.unlist.gff3")
fiveUTR=unique(fiveUTR)

##
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

#### All sites for background list GOrilla ####
sites.all <- data_lmm_out[,c(3,4)] #select chr, start 

nrow(sites.all) #117515
sites.all$pvalue <- NA #to format for GRanges
sites.all$end <- sites.all$start
sites.all$strand <- "+"
sites.all.GR <- as(sites.all, "GRanges")

#Promoters
annotated_all.promo <- subsetByOverlaps(promoter, sites.all.GR)
annotated_all.promo <- as.data.frame(annotated_all.promo)
annotated_all.promo.GO <- annotated_all.promo[,"gene_id"] #for GOrilla

#Genes
annotated_all.genes <- subsetByOverlaps(genes, sites.all.GR)
annotated_all.genes <- as.data.frame(annotated_all.genes)
annotated_all.genes.GO <- annotated_all.genes[,"gene_id"]

#TSS
annotated_all.TSS <- subsetByOverlaps(TSS, sites.all.GR)
annotated_all.TSS <- as.data.frame(annotated_all.TSS)
annotated_all.TSS.GO <- annotated_all.TSS[,"gene_id"]

#exon gene
annotated_all.exon.gene <- subsetByOverlaps(exons_gene, sites.all.GR)
annotated_all.exon.gene <- as.data.frame(annotated_all.exon.gene)
annotated_all.exon.gene.GO <- annotated_all.exon.gene [,"ID"]

#exon transcript
annotated_all.exon.trans <- subsetByOverlaps(exons_transcript, sites.all.GR)
annotated_all.exon.trans <- as.data.frame(annotated_all.exon.trans)
annotated_all.exon.trans.GO <- annotated_all.exon.trans[,"ID"]

#intron
annotated_all.intron <- subsetByOverlaps(introns, sites.all.GR)
annotated_all.intron <- as.data.frame(annotated_all.intron)
annotated_all.intron.GO <- annotated_all.intron[,"ID"]

#upstream
annotated_all.upstream <- subsetByOverlaps(upstream, sites.all.GR)
annotated_all.upstream <- as.data.frame(annotated_all.upstream)
annotated_all.upstream.GO <- annotated_all.upstream[,"gene_id"]

#downstream
annotated_all.downstream <- subsetByOverlaps(downstream, sites.all.GR)
annotated_all.downstream <- as.data.frame(annotated_all.downstream)
annotated_all.downstream.GO <- annotated_all.downstream[,"gene_id"]

# Add together for GOrilla
genes.for.GOrilla.all<-c(annotated_all.promo.GO, annotated_all.genes.GO,
                             annotated_all.TSS.GO, annotated_all.exon.gene.GO,
                             annotated_all.exon.trans.GO, annotated_all.intron.GO,
                             annotated_all.upstream.GO, annotated_all.downstream.GO)

genes.for.GOrilla.all <- unique(genes.for.GOrilla.all)
write.table(genes.for.GOrilla.all, row.names = FALSE, col.names=FALSE,
            quote = FALSE, "data/all.genes.for.GOrilla.txt")

# make df one row per cpg site, then add a col to indicate promoter etc.
annotated_all.promo.CG <- as.data.frame(subsetByOverlaps(sites.all.GR, promoter))
annotated_all.promo.CG <- add_column(annotated_all.promo.CG, "CHR_POS" = paste(annotated_all.promo.CG$seqnames, annotated_all.promo.CG$start, sep = "_"))

annotated_all.genes.CG <- as.data.frame(subsetByOverlaps(sites.all.GR, genes))
annotated_all.genes.CG <- add_column(annotated_all.genes.CG, "CHR_POS" = paste(annotated_all.genes.CG$seqnames, annotated_all.genes.CG$start, sep = "_"))

annotated_all.TSS.CG <- as.data.frame(subsetByOverlaps(sites.all.GR, TSS))
annotated_all.TSS.CG <- add_column(annotated_all.TSS.CG, "CHR_POS" = paste(annotated_all.TSS.CG$seqnames, annotated_all.TSS.CG$start, sep = "_"))

annotated_all.exons_gene.CG <- as.data.frame(subsetByOverlaps(sites.all.GR, exons_gene))
annotated_all.exons_gene.CG <- add_column(annotated_all.exons_gene.CG, "CHR_POS" = paste(annotated_all.exons_gene.CG$seqnames, annotated_all.exons_gene.CG$start, sep = "_"))

annotated_all.introns.CG <- as.data.frame(subsetByOverlaps(sites.all.GR, introns))
annotated_all.introns.CG <- add_column(annotated_all.introns.CG, "CHR_POS" = paste(annotated_all.introns.CG$seqnames, annotated_all.introns.CG$start, sep = "_"))

annotated_all.upstream.CG <- as.data.frame(subsetByOverlaps(sites.all.GR, upstream))
annotated_all.upstream.CG <- add_column(annotated_all.upstream.CG, "CHR_POS" = paste(annotated_all.upstream.CG$seqnames, annotated_all.upstream.CG$start, sep = "_"))

annotated_all.downstream.CG <- as.data.frame(subsetByOverlaps(sites.all.GR, downstream))
annotated_all.downstream.CG <- add_column(annotated_all.downstream.CG, "CHR_POS" = paste(annotated_all.downstream.CG$seqnames, annotated_all.downstream.CG$start, sep = "_"))

annotated_all.threeUTR.CG <- as.data.frame(subsetByOverlaps(sites.all.GR, threeUTR))
annotated_all.threeUTR.CG <- add_column(annotated_all.threeUTR.CG, "CHR_POS" = paste(annotated_all.threeUTR.CG$seqnames, annotated_all.threeUTR.CG$start, sep = "_"))

annotated_all.fiveUTR.CG <- as.data.frame(subsetByOverlaps(sites.all.GR, fiveUTR))
annotated_all.fiveUTR.CG <- add_column(annotated_all.fiveUTR.CG, "CHR_POS" = paste(annotated_all.fiveUTR.CG$seqnames, annotated_all.fiveUTR.CG$start, sep = "_"))


list.annotated.all <- sites.all[,c(1,2)] 
list.annotated.all <- add_column(list.annotated.all, "CHR_POS" = paste(list.annotated.all$chr, list.annotated.all$start, sep = "_"))

## add column
list.annotated.all <- list.annotated.all%>% mutate(
  region = as.factor(case_when(
    list.annotated.all$CHR_POS %in% annotated_all.TSS.CG$CHR_POS ~ "TSS",
    list.annotated.all$CHR_POS %in% annotated_all.promo.CG$CHR_POS  ~ "Promoter",
    list.annotated.all$CHR_POS %in% annotated_all.genes.CG$CHR_POS  ~ "Gene body",
    list.annotated.all$CHR_POS %in% annotated_all.exons_gene.CG$CHR_POS  ~ "Gene body",
    list.annotated.all$CHR_POS %in% annotated_all.introns.CG$CHR_POS  ~ "Gene body",
    list.annotated.all$CHR_POS %in% annotated_all.threeUTR.CG$CHR_POS  ~ "Gene body",
    list.annotated.all$CHR_POS %in% annotated_all.fiveUTR.CG$CHR_POS  ~ "Gene body",
    list.annotated.all$CHR_POS %in% annotated_all.upstream.CG$CHR_POS  ~ "Up- or downstream",
    list.annotated.all$CHR_POS %in% annotated_all.downstream.CG$CHR_POS  ~ "Up- or downstream",
    is.na(list.annotated.all$region) ~ "Intergenic"
  )
))

summary(as.factor(list.annotated.all$region))

write.csv(list.annotated.all, "/home/nioo/rebeccash/2018_BS_manipulation/R_Scripts/complete_MS_pt1/great-tit-methylation-genes-environment/data/cpglist_all.csv")

# pie chart
list.annotated.all.ggplot <- as.data.frame(summary(list.annotated.all$region))
list.annotated.all.ggplot <- rownames_to_column(list.annotated.all.ggplot, "Functional.Region")
names(list.annotated.all.ggplot)[2] <- "Nr.Sites"
piepercent <- found(list.annotated.all.ggplot$prop,2)
list.annotated.all.ggplot <- list.annotated.all.ggplot %>% 
  mutate(prop = Nr.Sites / sum(list.annotated.all.ggplot$Nr.Sites) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )%>% 
  mutate(csum = rev(cumsum(rev(Nr.Sites))), 
         pos = Nr.Sites/2 + lead(csum, 1),
         pos = if_else(is.na(pos), Nr.Sites/2, pos))

#rearrange
list.annotated.all.ggplot <- list.annotated.all.ggplot[c(4,3,1,5,2),]

png(file = "/home/nioo/rebeccash/2018_BS_manipulation/R_Scripts/plots/piechart_allcpg.png")
pie(list.annotated.all.ggplot$Nr.Sites, labels = paste0(round(list.annotated.all.ggplot$prop,1),"%"), 
    col = c('#c7eae5','#80cdc1','#35978f','#01665e','#003c30'),clockwise=T)
legend(0.8,1.05, c(list.annotated.all.ggplot$Functional.Region), cex = 0.8, fill = c('#c7eae5','#80cdc1','#35978f','#01665e','#003c30'))
dev.off()

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
annotated_M1.GeneticNK.promo.CG <- as.data.frame(subsetByOverlaps(sites.sig.geneticnk.M1.select.GR, promoter))
annotated_M1.GeneticNK.genes.CG <- as.data.frame(subsetByOverlaps(sites.sig.geneticnk.M1.select.GR, genes))
annotated_M1.GeneticNK.TSS.CG <- as.data.frame(subsetByOverlaps(sites.sig.geneticnk.M1.select.GR, TSS))
annotated_M1.GeneticNK.exons_gene.CG <- as.data.frame(subsetByOverlaps(sites.sig.geneticnk.M1.select.GR, exons_gene))
annotated_M1.GeneticNK.introns.CG <- as.data.frame(subsetByOverlaps(sites.sig.geneticnk.M1.select.GR, introns))
annotated_M1.GeneticNK.upstream.CG <- as.data.frame(subsetByOverlaps(sites.sig.geneticnk.M1.select.GR, upstream))
annotated_M1.GeneticNK.downstream.CG <- as.data.frame(subsetByOverlaps(sites.sig.geneticnk.M1.select.GR, downstream))
annotated_M1.GeneticNK.threeUTR.CG <- as.data.frame(subsetByOverlaps(sites.sig.geneticnk.M1.select.GR, threeUTR))
annotated_M1.GeneticNK.fiveUTR.CG <- as.data.frame(subsetByOverlaps(sites.sig.geneticnk.M1.select.GR, fiveUTR))

#annotated_M1.GeneticNK.unplaced <- sites.sig.geneticnk.M1.select
#annotated_M1.GeneticNK.unplaced <- add_column(annotated_M1.GeneticNK.unplaced, "CHR_POS" = paste(annotated_M1.GeneticNK.unplaced$chr, annotated_M1.GeneticNK.unplaced$start, sep = "_"))

annotated_M1.GeneticNK.promo.CG <- add_column(annotated_M1.GeneticNK.promo.CG, "CHR_POS" = paste(annotated_M1.GeneticNK.promo.CG$seqnames, annotated_M1.GeneticNK.promo.CG$start, sep = "_"))
annotated_M1.GeneticNK.genes.CG <- add_column(annotated_M1.GeneticNK.genes.CG, "CHR_POS" = paste(annotated_M1.GeneticNK.genes.CG$seqnames, annotated_M1.GeneticNK.genes.CG$start, sep = "_"))
annotated_M1.GeneticNK.TSS.CG <- add_column(annotated_M1.GeneticNK.TSS.CG, "CHR_POS" = paste(annotated_M1.GeneticNK.TSS.CG$seqnames, annotated_M1.GeneticNK.TSS.CG$start, sep = "_"))
annotated_M1.GeneticNK.exons_gene.CG <- add_column(annotated_M1.GeneticNK.exons_gene.CG, "CHR_POS" = paste(annotated_M1.GeneticNK.exons_gene.CG$seqnames, annotated_M1.GeneticNK.exons_gene.CG$start, sep = "_"))
annotated_M1.GeneticNK.introns.CG <- add_column(annotated_M1.GeneticNK.introns.CG, "CHR_POS" = paste(annotated_M1.GeneticNK.introns.CG$seqnames, annotated_M1.GeneticNK.introns.CG$start, sep = "_"))
annotated_M1.GeneticNK.upstream.CG <- add_column(annotated_M1.GeneticNK.upstream.CG, "CHR_POS" = paste(annotated_M1.GeneticNK.upstream.CG$seqnames, annotated_M1.GeneticNK.upstream.CG$start, sep = "_"))
annotated_M1.GeneticNK.downstream.CG <- add_column(annotated_M1.GeneticNK.downstream.CG, "CHR_POS" = paste(annotated_M1.GeneticNK.downstream.CG$seqnames, annotated_M1.GeneticNK.downstream.CG$start, sep = "_"))
annotated_M1.GeneticNK.threeUTR.CG <- add_column(annotated_M1.GeneticNK.threeUTR.CG, "CHR_POS" = paste(annotated_M1.GeneticNK.threeUTR.CG$seqnames, annotated_M1.GeneticNK.threeUTR.CG$start, sep = "_"))
annotated_M1.GeneticNK.fiveUTR.CG <- add_column(annotated_M1.GeneticNK.fiveUTR.CG, "CHR_POS" = paste(annotated_M1.GeneticNK.fiveUTR.CG$seqnames, annotated_M1.GeneticNK.fiveUTR.CG$start, sep = "_"))

#make a list and add extra col

list.annotated.genetic <- sites.sig.geneticnk.M1.select[,c(1,2)] 
list.annotated.genetic <- add_column(list.annotated.genetic, "CHR_POS" = paste(list.annotated.genetic$chr, list.annotated.genetic$start, sep = "_"))

## add column
list.annotated.genetic <- list.annotated.genetic%>% mutate(
  region = as.factor(case_when(
    list.annotated.genetic$CHR_POS %in% annotated_M1.GeneticNK.TSS.CG$CHR_POS ~ "TSS",
    list.annotated.genetic$CHR_POS %in% annotated_M1.GeneticNK.promo.CG$CHR_POS  ~ "Promoter",
    list.annotated.genetic$CHR_POS %in% annotated_M1.GeneticNK.genes.CG$CHR_POS  ~ "Gene body",
    list.annotated.genetic$CHR_POS %in% annotated_M1.GeneticNK.exons_gene.CG$CHR_POS  ~ "Gene body",
    list.annotated.genetic$CHR_POS %in% annotated_M1.GeneticNK.introns.CG$CHR_POS  ~ "Gene body",
    list.annotated.genetic$CHR_POS %in% annotated_M1.GeneticNK.threeUTR.CG$CHR_POS  ~ "Gene body",
    list.annotated.genetic$CHR_POS %in% annotated_M1.GeneticNK.fiveUTR.CG$CHR_POS  ~ "Gene body",
    list.annotated.genetic$CHR_POS %in% annotated_M1.GeneticNK.upstream.CG$CHR_POS  ~ "Up- or downstream",
    list.annotated.genetic$CHR_POS %in% annotated_M1.GeneticNK.downstream.CG$CHR_POS  ~ "Up- or downstream",
    is.na(list.annotated.genetic$region) ~ "Intergenic")
  ))

list.annotated.genetic$region <- as.factor(list.annotated.genetic$region)

summary(list.annotated.genetic$region)
write.csv(list.annotated.genetic, "/home/nioo/rebeccash/2018_BS_manipulation/R_Scripts/complete_MS_pt1/great-tit-methylation-genes-environment/data/cpglist_geneticnk.csv")

# pie chart
list.annotated.genetic.ggplot <- as.data.frame(summary(list.annotated.genetic$region))
list.annotated.genetic.ggplot <- rownames_to_column(list.annotated.genetic.ggplot, "Functional.Region")
names(list.annotated.genetic.ggplot)[2] <- "Nr.Sites"
piepercent <- round(list.annotated.genetic.ggplot$prop,2)
list.annotated.genetic.ggplot <- list.annotated.genetic.ggplot %>% 
  mutate(prop = Nr.Sites / sum(list.annotated.genetic.ggplot$Nr.Sites) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )%>% 
  mutate(csum = rev(cumsum(rev(Nr.Sites))), 
         pos = Nr.Sites/2 + lead(csum, 1),
         pos = if_else(is.na(pos), Nr.Sites/2, pos))
#rearrange
list.annotated.genetic.ggplot<- list.annotated.genetic.ggplot[c(4,3,1,5,2),]

png(file = "/home/nioo/rebeccash/2018_BS_manipulation/R_Scripts/plots/piechart_geneticcpg.png")
pie(list.annotated.genetic.ggplot$Nr.Sites, labels = paste0(round(list.annotated.genetic.ggplot$prop,1),"%"), 
    col = c('#c7eae5','#80cdc1','#35978f','#01665e','#003c30'),clockwise=T)
legend(0.8,1.05, c(list.annotated.genetic.ggplot$Functional.Region), cex = 0.8, fill = c('#c7eae5','#80cdc1','#35978f','#01665e','#003c30'))
dev.off()

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
annotated_M1.fosternk.promo.CG <- as.data.frame(subsetByOverlaps(sites.sig.fosternk.M1.select.GR, promoter))
annotated_M1.fosternk.genes.CG <- as.data.frame(subsetByOverlaps(sites.sig.fosternk.M1.select.GR, genes))
annotated_M1.fosternk.TSS.CG <- as.data.frame(subsetByOverlaps(sites.sig.fosternk.M1.select.GR, TSS))
annotated_M1.fosternk.exons_gene.CG <- as.data.frame(subsetByOverlaps(sites.sig.fosternk.M1.select.GR, exons_gene))
annotated_M1.fosternk.introns.CG <- as.data.frame(subsetByOverlaps(sites.sig.fosternk.M1.select.GR, introns))
annotated_M1.fosternk.upstream.CG <- as.data.frame(subsetByOverlaps(sites.sig.fosternk.M1.select.GR, upstream))
annotated_M1.fosternk.downstream.CG <- as.data.frame(subsetByOverlaps(sites.sig.fosternk.M1.select.GR, downstream))
annotated_M1.fosternk.threeUTR.CG <- as.data.frame(subsetByOverlaps(sites.sig.fosternk.M1.select.GR, threeUTR))
annotated_M1.fosternk.fiveUTR.CG <- as.data.frame(subsetByOverlaps(sites.sig.fosternk.M1.select.GR, fiveUTR))

#annotated_M1.fosternk.unplaced <- sites.sig.fosternk.M1.select
#annotated_M1.fosternk.unplaced <- add_column(annotated_M1.fosternk.unplaced, "CHR_POS" = paste(annotated_M1.fosternk.unplaced$chr, annotated_M1.fosternk.unplaced$start, sep = "_"))

annotated_M1.fosternk.promo.CG <- add_column(annotated_M1.fosternk.promo.CG, "CHR_POS" = paste(annotated_M1.fosternk.promo.CG$seqnames, annotated_M1.fosternk.promo.CG$start, sep = "_"))
annotated_M1.fosternk.genes.CG <- add_column(annotated_M1.fosternk.genes.CG, "CHR_POS" = paste(annotated_M1.fosternk.genes.CG$seqnames, annotated_M1.fosternk.genes.CG$start, sep = "_"))
annotated_M1.fosternk.TSS.CG <- add_column(annotated_M1.fosternk.TSS.CG, "CHR_POS" = paste(annotated_M1.fosternk.TSS.CG$seqnames, annotated_M1.fosternk.TSS.CG$start, sep = "_"))
annotated_M1.fosternk.exons_gene.CG <- add_column(annotated_M1.fosternk.exons_gene.CG, "CHR_POS" = paste(annotated_M1.fosternk.exons_gene.CG$seqnames, annotated_M1.fosternk.exons_gene.CG$start, sep = "_"))
annotated_M1.fosternk.introns.CG <- add_column(annotated_M1.fosternk.introns.CG, "CHR_POS" = paste(annotated_M1.fosternk.introns.CG$seqnames, annotated_M1.fosternk.introns.CG$start, sep = "_"))
annotated_M1.fosternk.upstream.CG <- add_column(annotated_M1.fosternk.upstream.CG, "CHR_POS" = paste(annotated_M1.fosternk.upstream.CG$seqnames, annotated_M1.fosternk.upstream.CG$start, sep = "_"))
annotated_M1.fosternk.downstream.CG <- add_column(annotated_M1.fosternk.downstream.CG, "CHR_POS" = paste(annotated_M1.fosternk.downstream.CG$seqnames, annotated_M1.fosternk.downstream.CG$start, sep = "_"))
annotated_M1.fosternk.threeUTR.CG <- add_column(annotated_M1.fosternk.threeUTR.CG, "CHR_POS" = paste(annotated_M1.fosternk.threeUTR.CG$seqnames, annotated_M1.fosternk.threeUTR.CG$start, sep = "_"))
annotated_M1.fosternk.fiveUTR.CG <- add_column(annotated_M1.fosternk.fiveUTR.CG, "CHR_POS" = paste(annotated_M1.fosternk.fiveUTR.CG$seqnames, annotated_M1.fosternk.fiveUTR.CG$start, sep = "_"))


#make a list and add extra col

list.annotated.foster <- sites.sig.fosternk.M1.select[,c(1,2)] 
list.annotated.foster <- add_column(list.annotated.foster, "CHR_POS" = paste(list.annotated.foster$chr, list.annotated.foster$start, sep = "_"))

## add column
list.annotated.foster <- list.annotated.foster%>% mutate(
  region = as.factor(case_when(
    list.annotated.foster$CHR_POS %in% annotated_M1.fosternk.TSS.CG$CHR_POS ~ "TSS",
    list.annotated.foster$CHR_POS %in% annotated_M1.fosternk.promo.CG$CHR_POS  ~ "Promoter",
    list.annotated.foster$CHR_POS %in% annotated_M1.fosternk.genes.CG$CHR_POS  ~ "Gene body",
    list.annotated.foster$CHR_POS %in% annotated_M1.fosternk.exons_gene.CG$CHR_POS  ~ "Gene body",
    list.annotated.foster$CHR_POS %in% annotated_M1.fosternk.introns.CG$CHR_POS  ~ "Gene body",
    list.annotated.foster$CHR_POS %in% annotated_M1.fosternk.threeUTR.CG$CHR_POS  ~ "Gene body",
    list.annotated.foster$CHR_POS %in% annotated_M1.fosternk.fiveUTR.CG$CHR_POS  ~ "Gene body",
    list.annotated.foster$CHR_POS %in% annotated_M1.fosternk.upstream.CG$CHR_POS  ~ "Up- or downstream",
    list.annotated.foster$CHR_POS %in% annotated_M1.fosternk.downstream.CG$CHR_POS  ~ "Up- or downstream",
    is.na(list.annotated.foster$region) ~ "Intergenic")
  ))

list.annotated.foster$region <- as.factor(list.annotated.foster$region)

summary(list.annotated.foster$region)
write.csv(list.annotated.foster, "/home/nioo/rebeccash/2018_BS_manipulation/R_Scripts/complete_MS_pt1/great-tit-methylation-genes-environment/data/cpglist_fosternk.csv")

# pie chart
list.annotated.foster.ggplot <- as.data.frame(summary(list.annotated.foster$region))
list.annotated.foster.ggplot <- rownames_to_column(list.annotated.foster.ggplot, "Functional.Region")
names(list.annotated.foster.ggplot)[2] <- "Nr.Sites"
piepercent <- found(list.annotated.foster.ggplot$prop,2)
list.annotated.foster.ggplot <- list.annotated.foster.ggplot %>% 
  mutate(prop = Nr.Sites / sum(list.annotated.foster.ggplot$Nr.Sites) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )%>% 
  mutate(csum = rev(cumsum(rev(Nr.Sites))), 
         pos = Nr.Sites/2 + lead(csum, 1),
         pos = if_else(is.na(pos), Nr.Sites/2, pos))

#rearrange
list.annotated.foster.ggplot<-list.annotated.foster.ggplot[c(4,3,1,5,2),]
png(file = "/home/nioo/rebeccash/2018_BS_manipulation/R_Scripts/plots/piechart_fostercpg.png")
pie(list.annotated.foster.ggplot$Nr.Sites, labels = paste0(round(list.annotated.foster.ggplot$prop,1),"%"), 
    col = c('#c7eae5','#80cdc1','#35978f','#01665e','#003c30'),clockwise=T)
legend(0.8,1.05, c(list.annotated.foster.ggplot$Functional.Region), cex = 0.8, fill = c('#c7eae5','#80cdc1','#35978f','#01665e','#003c30'))
dev.off()

## venn diagram
library(VennDiagram)
draw.pairwise.venn(area1=nrow(sites.sig.geneticnk.M1.select), area2=nrow(sites.sig.fosternk.M1.select),
                 cross.area=nrow(sites.sig.both.M1), fill = c("#01665e", "#c7eae5"), col = c("#01665e", "#c7eae5"), 
                 cex=1, cat.cex = 0, category=c("Brood of Origin", "Brood of Rearing"),
                 head = "(a)                                                                  ") #to edit

### combine graphs in one


png(file = "/home/nioo/rebeccash/2018_BS_manipulation/R_Scripts/plots/piechart_combined.png")
par(mfrow=c(2,2))
draw.pairwise.venn(area1=nrow(sites.sig.geneticnk.M1.select), area2=nrow(sites.sig.fosternk.M1.select),
                   cross.area=nrow(sites.sig.both.M1), fill = c("#01665e", "#c7eae5"), col = c("#01665e", "#c7eae5"), 
                   cex=1, cat.cex = 0) #to edit

pie(list.annotated.all.ggplot$Nr.Sites, labels = paste0(round(list.annotated.all.ggplot$prop,1),"%"), 
    col = c('#c7eae5','#80cdc1','#35978f','#01665e','#003c30'),clockwise=T)

pie(list.annotated.genetic.ggplot$Nr.Sites, labels = paste0(round(list.annotated.genetic.ggplot$prop,1),"%"), 
    col = c('#c7eae5','#80cdc1','#35978f','#01665e','#003c30'),clockwise=T)

pie(list.annotated.foster.ggplot$Nr.Sites, labels = paste0(round(list.annotated.foster.ggplot$prop,1),"%"), 
    col = c('#c7eae5','#80cdc1','#35978f','#01665e','#003c30'),clockwise=T)
#legend(0.8,1.05, c(list.annotated.foster.ggplot$Functional.Region), cex = 0.8, fill = c('#c7eae5','#80cdc1','#35978f','#01665e','#003c30'))
dev.off()

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

########### GO TERMS FOR SITES FILTERED BY OVERDISPERSION #############

#### Genetic NK LMER ####
nrow(sites.sig.geneticnk.M1.select.overdisp.adjchr) #6054
names(sites.sig.geneticnk.M1.select.overdisp.adjchr)[3] <- "pvalue" #to format for GRanges
sites.sig.geneticnk.M1.select.overdisp.adjchr$end <- sites.sig.geneticnk.M1.select.overdisp.adjchr$start
sites.sig.geneticnk.M1.select.overdisp.adjchr$strand <- "+"
sites.sig.geneticnk.M1.select.overdisp.adjchr.GR <- as(sites.sig.geneticnk.M1.select.overdisp.adjchr, "GRanges")

#Promoters
annotated_geneticnk_overdisp.promo <- subsetByOverlaps(promoter, sites.sig.geneticnk.M1.select.overdisp.adjchr.GR)
annotated_geneticnk_overdisp.promo <- as.data.frame(annotated_geneticnk_overdisp.promo)
annotated_geneticnk_overdisp.promo.GO <- annotated_geneticnk_overdisp.promo[,"gene_id"] #for GOrilla

#Genes
annotated_geneticnk_overdisp.genes <- subsetByOverlaps(genes, sites.sig.geneticnk.M1.select.overdisp.adjchr.GR)
annotated_geneticnk_overdisp.genes <- as.data.frame(annotated_geneticnk_overdisp.genes)
annotated_geneticnk_overdisp.genes.GO <- annotated_geneticnk_overdisp.genes[,"gene_id"]

#TSS
annotated_geneticnk_overdisp.TSS <- subsetByOverlaps(TSS, sites.sig.geneticnk.M1.select.overdisp.adjchr.GR)
annotated_geneticnk_overdisp.TSS <- as.data.frame(annotated_geneticnk_overdisp.TSS)
annotated_geneticnk_overdisp.TSS.GO <- annotated_geneticnk_overdisp.TSS[,"gene_id"]

#exon gene
annotated_geneticnk_overdisp.exon.gene <- subsetByOverlaps(exons_gene, sites.sig.geneticnk.M1.select.overdisp.adjchr.GR)
annotated_geneticnk_overdisp.exon.gene <- as.data.frame(annotated_geneticnk_overdisp.exon.gene)
annotated_geneticnk_overdisp.exon.gene.GO <- annotated_geneticnk_overdisp.exon.gene [,"ID"]

#exon transcript
annotated_geneticnk_overdisp.exon.trans <- subsetByOverlaps(exons_transcript, sites.sig.geneticnk.M1.select.overdisp.adjchr.GR)
annotated_geneticnk_overdisp.exon.trans <- as.data.frame(annotated_geneticnk_overdisp.exon.trans)
annotated_geneticnk_overdisp.exon.trans.GO <- annotated_geneticnk_overdisp.exon.trans[,"ID"]

#intron
annotated_geneticnk_overdisp.intron <- subsetByOverlaps(introns, sites.sig.geneticnk.M1.select.overdisp.adjchr.GR)
annotated_geneticnk_overdisp.intron <- as.data.frame(annotated_geneticnk_overdisp.intron)
annotated_geneticnk_overdisp.intron.GO <- annotated_geneticnk_overdisp.intron[,"ID"]

#upstream
annotated_geneticnk_overdisp.upstream <- subsetByOverlaps(upstream, sites.sig.geneticnk.M1.select.overdisp.adjchr.GR)
annotated_geneticnk_overdisp.upstream <- as.data.frame(annotated_geneticnk_overdisp.upstream)
annotated_geneticnk_overdisp.upstream.GO <- annotated_geneticnk_overdisp.upstream[,"gene_id"]

#downstream
annotated_geneticnk_overdisp.downstream <- subsetByOverlaps(downstream, sites.sig.geneticnk.M1.select.overdisp.adjchr.GR)
annotated_geneticnk_overdisp.downstream <- as.data.frame(annotated_geneticnk_overdisp.downstream)
annotated_geneticnk_overdisp.downstream.GO <- annotated_geneticnk_overdisp.downstream[,"gene_id"]

# Add together for GOrilla
genes.for.GOrilla_overdisp.downstream<-c(annotated_geneticnk_overdisp.promo.GO, annotated_geneticnk_overdisp.genes.GO,
                         annotated_geneticnk_overdisp.TSS.GO, annotated_geneticnk_overdisp.exon.gene.GO,
                         annotated_geneticnk_overdisp.exon.trans.GO, annotated_geneticnk_overdisp.intron.GO,
                         annotated_geneticnk_overdisp.upstream.GO, annotated_geneticnk_overdisp.downstream.GO)

genes.for.GOrilla_overdisp.downstream <- unique(genes.for.GOrilla_overdisp.downstream)
write.table(genes.for.GOrilla_overdisp.downstream, row.names = FALSE, col.names=FALSE,
            quote = FALSE, "data/genes.for.GOrilla.GeneticNK.woverdisp.txt")


#### foster NK LMER ####
nrow(sites.sig.fosternk.M1.select.overdisp.adjchr) #101
names(sites.sig.fosternk.M1.select.overdisp.adjchr)[3] <- "pvalue" #to format for GRanges
sites.sig.fosternk.M1.select.overdisp.adjchr$end <- sites.sig.fosternk.M1.select.overdisp.adjchr$start
sites.sig.fosternk.M1.select.overdisp.adjchr$strand <- "+"
sites.sig.fosternk.M1.select.overdisp.adjchr.GR <- as(sites.sig.fosternk.M1.select.overdisp.adjchr, "GRanges")

#Promoters
annotated_fosternk_overdisp.promo <- subsetByOverlaps(promoter, sites.sig.fosternk.M1.select.overdisp.adjchr.GR)
annotated_fosternk_overdisp.promo <- as.data.frame(annotated_fosternk_overdisp.promo)
annotated_fosternk_overdisp.promo.GO <- annotated_fosternk_overdisp.promo[,"gene_id"] #for GOrilla

#Genes
annotated_fosternk_overdisp.genes <- subsetByOverlaps(genes, sites.sig.fosternk.M1.select.overdisp.adjchr.GR)
annotated_fosternk_overdisp.genes <- as.data.frame(annotated_fosternk_overdisp.genes)
annotated_fosternk_overdisp.genes.GO <- annotated_fosternk_overdisp.genes[,"gene_id"]

#TSS
annotated_fosternk_overdisp.TSS <- subsetByOverlaps(TSS, sites.sig.fosternk.M1.select.overdisp.adjchr.GR)
annotated_fosternk_overdisp.TSS <- as.data.frame(annotated_fosternk_overdisp.TSS)
annotated_fosternk_overdisp.TSS.GO <- annotated_fosternk_overdisp.TSS[,"gene_id"]

#exon gene
annotated_fosternk_overdisp.exon.gene <- subsetByOverlaps(exons_gene, sites.sig.fosternk.M1.select.overdisp.adjchr.GR)
annotated_fosternk_overdisp.exon.gene <- as.data.frame(annotated_fosternk_overdisp.exon.gene)
annotated_fosternk_overdisp.exon.gene.GO <- annotated_fosternk_overdisp.exon.gene [,"ID"]

#exon transcript
annotated_fosternk_overdisp.exon.trans <- subsetByOverlaps(exons_transcript, sites.sig.fosternk.M1.select.overdisp.adjchr.GR)
annotated_fosternk_overdisp.exon.trans <- as.data.frame(annotated_fosternk_overdisp.exon.trans)
annotated_fosternk_overdisp.exon.trans.GO <- annotated_fosternk_overdisp.exon.trans[,"ID"]

#intron
annotated_fosternk_overdisp.intron <- subsetByOverlaps(introns, sites.sig.fosternk.M1.select.overdisp.adjchr.GR)
annotated_fosternk_overdisp.intron <- as.data.frame(annotated_fosternk_overdisp.intron)
annotated_fosternk_overdisp.intron.GO <- annotated_fosternk_overdisp.intron[,"ID"]

#upstream
annotated_fosternk_overdisp.upstream <- subsetByOverlaps(upstream, sites.sig.fosternk.M1.select.overdisp.adjchr.GR)
annotated_fosternk_overdisp.upstream <- as.data.frame(annotated_fosternk_overdisp.upstream)
annotated_fosternk_overdisp.upstream.GO <- annotated_fosternk_overdisp.upstream[,"gene_id"]

#downstream
annotated_fosternk_overdisp.downstream <- subsetByOverlaps(downstream, sites.sig.fosternk.M1.select.overdisp.adjchr.GR)
annotated_fosternk_overdisp.downstream <- as.data.frame(annotated_fosternk_overdisp.downstream)
annotated_fosternk_overdisp.downstream.GO <- annotated_fosternk_overdisp.downstream[,"gene_id"]

# Add together for GOrilla
genes.for.GOrilla_fosternk_overdisp.downstream<-c(annotated_fosternk_overdisp.promo.GO, annotated_fosternk_overdisp.genes.GO,
                                         annotated_fosternk_overdisp.TSS.GO, annotated_fosternk_overdisp.exon.gene.GO,
                                         annotated_fosternk_overdisp.exon.trans.GO, annotated_fosternk_overdisp.intron.GO,
                                         annotated_fosternk_overdisp.upstream.GO, annotated_fosternk_overdisp.downstream.GO)

genes.for.GOrilla_fosternk_overdisp.downstream <- unique(genes.for.GOrilla_fosternk_overdisp.downstream)
write.table(genes.for.GOrilla_fosternk_overdisp.downstream, row.names = FALSE, col.names=FALSE,
            quote = FALSE, "data/genes.for.GOrilla.FosterNK.woverdisp.txt")
