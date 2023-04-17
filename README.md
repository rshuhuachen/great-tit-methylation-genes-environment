# great-tit-methylation-genes-environment
This is the full workflow for Sepers et al. (2023) "Variation in DNA methylation in avian nestlings is largely determined by genetic effects" https://doi.org/10.1093/molbev/msad086

The workflow is divided up in 4 scripts:

1. Filtering and plotting data using methylKit, conducting a PCA and DMA using DSS
2. Fitting a model per CpG site for variance partitioning of genetic brood (brood of origin) and foster brood (brood of rearing)
3. Annotating the significant CpG sites and creating pie charts
4. SNPs


Several data files can be found within /data
1. chickparameters.txt: Information per nestling ID including file names, brood names, molecular sex, etc.
2. ModelOutput.RData: This is the RData file containing the complete model output as fitted in script number 2, with each row containing model output per CpG site
3. genes.for.GOrilla.GeneticNK.txt and genes.for.GOrilla.FosterNK.txt: These files contain lists of genes significant for brood of origin (GeneticNK) and brood of rearing (FosterNK), respectively. These datafiles were used as input for GOrilla, and are sorted by descending significance (i.e. lowest p-value at the top).

In addition, raw data can be found in the following public repositories:
1. Raw reads can be found on NCBI under BioProject ID PRJNA208335 with accession numbers X
2. Reference genome parus major v1.1 can be downloaded from https://www.ncbi.nlm.nih.gov/assembly/GCF_001522545.3

