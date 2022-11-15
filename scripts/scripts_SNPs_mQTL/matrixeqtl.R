# Autor: Matrix eQTL by Andrey A. Shabalin
# URL link: http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/
# Date: 26-9-2022 by Bernice Sepers & Michelle Memelink
# Before running this code, change the input and output pathways to the 
# appropriate directories. 

library(MatrixEQTL)
library(png)

## Settings
# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR; 

# Genotype file name
SNP_file_name = paste("/home/SNP_filt.txt", sep="");
snps_location_file_name = paste("/home/snpsloc_filt.txt", sep="");

# methylation file name
expression_file_name = paste("/home/GE_noSNPs3.txt", sep="");
gene_location_file_name = paste("/home/geneloc_noSNPs3.txt", sep="");

# Covariates file name
covariates_file_name = paste("/home/Covariates.txt", sep="");

# Output file name
output_file_name_cis = tempfile();
output_file_name_tra = tempfile();

# Only associations significant at this level will be saved
pvOutputThreshold_cis = 1;
pvOutputThreshold_tra = 1;

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();

# max distance for local CpG-SNP pairs
cisDist = 1e6;

## Load genotype data
snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

## Load methylation data
gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

## Load covariates
cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}

## Run the analysis
snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);

me = Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name     = output_file_name_tra,
  pvOutputThreshold     = pvOutputThreshold_tra,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos,
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = 100,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);

#6685 CpGs
#16332 SNPs
#100.00% done, 459,777 cis-eQTLs, 108,719,643 trans-eQTLs

me2 = Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name     = output_file_name_tra,
  pvOutputThreshold     = pvOutputThreshold_tra,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos,
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);

unlink(output_file_name_tra);
unlink(output_file_name_cis);

## Results:
cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
png(filename="/home/histogram.png")
plot(me)
dev.off()

cat('Analysis done in: ', me2$time.in.sec, ' seconds', '\n');
png(filename="/home/q_q.png")
plot(me2)
dev.off()

write.csv(me$trans$eqtls, "/home/output2_trans.txt")
write.csv(me$cis$eqtls, "/home/output2_cis.txt")
