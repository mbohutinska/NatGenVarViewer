# NatGenVarViewer
## The natural genetic variability viewer

This pipeline allows to screen any set of genes to explore natural genetic variability in these genes in a given species. It is primarily designed to screen for natural variability in Arabidopsis arenosa, which is ecologically diverse species of Arabidopsis. The pipeline consists of three steps:

### STEP 1. visualization of allele frequencies 
of single nucleotide polymorphisms within a given set of genes across 74 A. arenosa populations

Example output for one gene, population sorted by altitude:
![image](https://user-images.githubusercontent.com/40301863/218489804-db5cb2c6-98d6-434c-8e0d-d094d848f902.png)



### STEP 2. simple allele frequency difference-based selection scan 
across given set of genes

Example result of the output list of outlier genes which are candidates fr positive selection:
![image](https://user-images.githubusercontent.com/40301863/218490352-2783be8c-352c-4f10-b62a-f4dd5b2f690e.png)


### STEP 3. visualization of outlier genes and SNPs 
identified in STEP 2

Example visualization of candidate SNPs, differentiated between high and low altitude populations: 
![image](https://user-images.githubusercontent.com/40301863/218491704-d6843f0b-19da-4a3a-af98-7bce0476ccd5.png)




All the steps are well described in the script thePipeline.R.

Note that for a mere exploration you can also only run STEP 1.


## Instructions: 

- Download all files from the repository, open the R script thePipeline.R and run it. All steps are clearly described within the script. It will require to load certain tables (popInfo.txt, geneSet.txt, etc.), examples of which are all given in this repository. I recommend to place them into a separate /data directory. 

- Results will be output into a separate /results directory, which the user should create within the working directory.

- User needs to supply information about a total genetic variability of the studied species. For A. arenosa, this table (ALL613.table.recode.txt) consists of > 3.1M rows and 623 columns. Thus I only provide the example 10k rows here (ALL613.table.recode.example.txt.gz). You can download the total table from my account on Metacentrum /storage/pruhonice1-ibot/home/holcovam/ScanTools/VCF_arenosa613.merged.annotated_DP4.M0.3/ALL613.table.recode.txt (ecolgen team & associates) or contact me through e-mail. This table was filtered for minimum individual read depth = 4, maximum fraction of filtered genotype = 0.3 and minimum minor allele frequency across all populations > 0.05.






