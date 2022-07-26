### Here we will construct the model per CpG site ###

## Load libraries ##
library(MuMIn); library(lmerTest)
library(sjmisc); library(data.table);library(stringr); 
library(dplyr); library(insight)
library(pbnm); library(pbkrtest); library(r2glmm); library(tidystats)
library(performance); library(parallel)

## Set working directory ##
#setwd("/home/nioo/rebeccash/2018_BS_manipulation/R_Scripts/complete_MS_pt1") #set this to the cloned github dir

#dataframe chicks
#load methpattern.total.chick dataframe 
base::load(file = "data/ChickData.RData") #=methpattern.total.chick
#see script 1.FilteringPlotting how to get ChickData.RData

base::load(file = "data/reference_genome/GCF_001522545.3_Parus_major1.1_assembly_report_chr_names.RData")
#ref genome parus major v1.1 can be downloaded from https://www.ncbi.nlm.nih.gov/assembly/GCF_001522545.3

#define overdispersion function

overdisp.lmer_fun <- function(model) {
  vpars <- function(m) {
    nrow(m)*(nrow(m)+1)/2
  }
  model.df <- sum(sapply(VarCorr(model),vpars))+length(fixef(model))
  rdf <- nrow(model.frame(model))-model.df
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}
#function taken from http://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#overdispersion accessed latest 29/03/2022

#run model with treatment and sex as fixed factor

# make list to run glmms in parallel 
data <- methpattern.total.chick %>% group_split(CHR_POS)
# write function for fitting glmms:
function_lmm <- function(df) {
  
  # fit the model and null-model, do lrt, and get p-value
  df <- as.data.frame(df)
  CHR_POS <- as.character(df[1,2])
  null <- lmerTest::lmer(Methylation ~ 1 + (1|GeneticNK)+ (1|FosterNK), df)
  model <- lmerTest::lmer(Methylation ~ Sex + (1|GeneticNK)+ (1|FosterNK), df)
  
  lrt <- anova(null, model)
  pvalue=lrt$"Pr(>Chisq)"[2]
  deviance <- lrt$deviance[2]
  Df <- lrt$Df[2]
  Chisq <- lrt$Chisq[2]
  GeneticNK.Intercept <- tidy_stats(model)$effects$random_effects$groups[[1]]$terms[[1]]$statistics$var
  GeneticNK.SD <- tidy_stats(model)$effects$random_effects$groups[[1]]$terms[[1]]$statistics$SD
  GeneticNK.LRT <- ranova(model)[4][[1]][2]
  GeneticNK.pval <- ranova(model)[6][[1]][2]
  FosterNK.Intercept <- tidy_stats(model)$effects$random_effects$groups[[2]]$terms[[1]]$statistics$var
  FosterNK.SD <- tidy_stats(model)$effects$random_effects$groups[[2]]$terms[[1]]$statistics$SD
  FosterNK.LRT <- ranova(model)[4][[1]][3]
  FosterNK.pval <- ranova(model)[6][[1]][3]
  Residual.Var <- tidy_stats(model)$effects$random_effects$groups[[3]]$statistics$var
  Residual.SD <- tidy_stats(model)$effects$random_effects$groups[[3]]$statistics$SD
  RsqM <- r.squaredGLMM(model)[1]
  RsqC <- r.squaredGLMM(model)[2]
  Dispersion.chisq <- as.data.frame(overdisp.lmer_fun(model)[1])
  Dispersion.ratio <- as.data.frame(overdisp.lmer_fun(model)[2])
  Dispersion.rdf <- as.data.frame(overdisp.lmer_fun(model)[3])
  Dispersion.pval <- as.data.frame(overdisp.lmer_fun(model)[4])
  isSingular <- isSingular(model)
  ICC_GeneticNK <-icc(model, by_group = TRUE, tolerance = 0)[1,2]
  ICC_FosterNK <- icc(model, by_group = TRUE, tolerance = 0)[2,2]
  return(data.frame(CHR_POS=CHR_POS, 
                    df=Df, 
                    Chisq=Chisq, 
                    pval_lrt=pvalue, 
                    GeneticNK.Intercept = GeneticNK.Intercept,
                    GeneticNK.SD = GeneticNK.SD,
                    GeneticNK.LRT = GeneticNK.LRT,
                    GeneticNK.pval = GeneticNK.pval,
                    FosterNK.Intercept = FosterNK.Intercept,
                    FosterNK.SD = FosterNK.SD,
                    FosterNK.LRT = FosterNK.LRT,
                    FosterNK.pval = FosterNK.pval,
                    Residual.Var = Residual.Var,
                    Residual.SD = Residual.SD,
                    RsqM = RsqM,
                    RsqC = RsqC,
                    Dispersion.chisq = Dispersion.chisq,
                    Dispersion.ratio = Dispersion.ratio,
                    Dispersion.rdf = Dispersion.rdf,
                    Dispersion.pval = Dispersion.pval,
                    isSingular = isSingular,
                    ICC_GeneticNK = ICC_GeneticNK,
                    ICC_FosterNK = ICC_FosterNK))
}

## run model using mclapply to spread across multiple cores
lmm_out <- mclapply(data, function_lmm, mc.cores=40)

#then save as dataframe rather than list
data_lmm_out <- do.call(rbind.data.frame, lmm_out)

#correct the chr_pos names
chr_pos <- sapply(data,'[',1,1)
chr_pos <- as.data.frame(unlist(chr_pos))
names(chr_pos)[1] <- "CHR_POS"

names(data_lmm_out)[1] <- "CHR"
data_lmm_out <- cbind(data_lmm_out, chr_pos)
data_lmm_out <- data_lmm_out %>% relocate("CHR_POS", .after = "CHR") 

#multiple testing correction for p-values using FDR
data_lmm_out$GeneticNK.qval <- p.adjust(data_lmm_out$GeneticNK.pval, method="fdr", n = nrow(data_lmm_out))
data_lmm_out$FosterNK.qval <- p.adjust(data_lmm_out$FosterNK.pval, method="fdr", n = nrow(data_lmm_out))

## count how many sig
nrow(subset(data_lmm_out, GeneticNK.qval < 0.05 & ICC_GeneticNK > 0.15))
nrow(subset(data_lmm_out, FosterNK.qval < 0.05 & ICC_FosterNK > 0.15))
nrow(subset(data_lmm_out, GeneticNK.qval<0.05 & ICC_GeneticNK > 0.15 & FosterNK.qval<0.05 & ICC_FosterNK > 0.15))

data_lmm_out_notreat <- data_lmm_out
base::save(data_lmm_out_notreat, file = "data/ModelOutput_noTreatment.RData")
load("data/ModelOutput.RData") #this file can be found in github directory
load("data/ModelOutput_noTreatment.RData") #this file can be found in github directory

data_lmm_out$Var.Perc.Gen
data_glmm_out_sign_genetic <- subset(data_lmm_out, GeneticNK.qval<0.05 & Var.Perc.Gen > 15)
data_glmm_out_sign_foster <- subset(data_lmm_out, FosterNK.qval<0.05 & Var.Perc.Fos > 15)
data_glmm_out_sign_overlap <- subset(data_lmm_out, GeneticNK.qval<0.05 & Var.Perc.Gen > 0.15 & FosterNK.qval<0.05 & Var.Perc.Fos > 0.15)

data_glmm_out_sign_genetic_icc <- subset(data_lmm_out, GeneticNK.qval<0.05 & ICC_GeneticNK > 0.15)
data_glmm_out_sign_foster_icc <- subset(data_lmm_out, FosterNK.qval<0.05 & ICC_FosterNK > 0.15)
data_glmm_out_sign_overlap_icc <- subset(data_lmm_out, GeneticNK.qval<0.05 & ICC_GeneticNK > 0.15 & ICC_FosterNK<0.05 & Var.Perc.Fos > 0.15)

data_glmm_out_sign_genetic2 <- subset(data_lmm_out2, GeneticNK.qval<0.05 & ICC_GeneticNK > 0.15)
data_glmm_out_sign_foster2 <- subset(data_lmm_out2, FosterNK.qval<0.05 & ICC_FosterNK > 0.15)
data_glmm_out_sign_overlap2 <- subset(data_lmm_out2, GeneticNK.qval<0.05 & ICC_GeneticNK > 0.15 & FosterNK.qval<0.05 & ICC_FosterNK > 0.15)

data_glmm_out_sign_genetic_notreat <- subset(data_lmm_out_notreat, GeneticNK.qval<0.05 & ICC_GeneticNK > 0.15)
data_glmm_out_sign_foster_notreat <- subset(data_lmm_out_notreat, FosterNK.qval<0.05 & ICC_FosterNK > 0.15)
data_glmm_out_sign_overlap_notreat <- subset(data_lmm_out_notreat, GeneticNK.qval<0.05 & ICC_GeneticNK > 0.15 & FosterNK.qval<0.05 & ICC_FosterNK > 0.15)

summary(data_glmm_out_sign_genetic_notreat$CHR_POS %in% data_glmm_out_sign_genetic$CHR_POS)
summary(data_glmm_out_sign_foster_notreat$CHR_POS %in% data_glmm_out_sign_foster$CHR_POS)
summary(data_glmm_out_sign_overlap_notreat$CHR_POS %in% data_glmm_out_sign_overlap$CHR_POS)
