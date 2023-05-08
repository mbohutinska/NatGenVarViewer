############# The natural genetic variability viewer #############
############# Extension for genome-wide scans #############

setwd("/home/aa/BetyProject/claraMS/cuticle/thePipeline/")

# Libraries
library(data.table)
library(dplyr)
library(gplots)
library(RColorBrewer)

### PART 0. DEFINE PARAMETERS ### 
### this part is common for all other parts
# User-supplied parameters 
adaptation<-"Altitude" # select environment of interest. It has to correspond to one of columns in 'popInfo.txt'.
project<-"GuillaumePlasticity" # project name will appear in the output tables name
tresholdSNP<-0.99 # outlier SNPs above this threshold, only for Part 2 & 3
tresholdGene<-0.25 # outlier SNPs density per gene cutoff, takes top tresholdGene genes. (0.25 takes top 25%). 
popInfo<-fread("data/popInfo.txt",h=T,sep="\t") # names and environmental variables for populations
typeOfSNPs<-c("upstream_gene_variant","5_prime_UTR_variant","5_prime_UTR_premature_start_codon_gain_variant") # "all" or anything else

derivedEnv<-c("VEL","LAC","TKO","SCH")#populations used for selections scan. In this case alpine populaions. only for Part 2 & 3
ancestralEnv<-c("SUB","ING","HRA","DRA") #populations used for selections scan. In this case foothill populaions. only for Part 2 & 3
#Alternatively define based on environmental values from popInfo.txt:
#derivedEnv<-subset(popInfo,popInfo$Altitude>1600)$pop
#ancestralEnv<-subset(popInfo,popInfo$Altitude<900 & popInfo$Altitude>400)$pop


### PART 4: GENOME-WIDE SELECTION SCAN ###
totAF<-fread("/home/aa/BetyProject/envAssoc/Data/WGDataAF.txt",h=T)

if (typeOfSNPs == "all") {} else
{totAF<-subset(totAF,totAF$ann %in% typeOfSNPs)}

afANC <- as.data.frame(totAF %>% dplyr:: select(starts_with(ancestralEnv)))
afDER <- as.data.frame(totAF %>% dplyr:: select(starts_with(derivedEnv)))

# Compute allele frequency difference (AFD)
totAF$ancAF<-apply(afANC,1,mean,na.rm=T)
totAF$derAF<-apply(afDER,1,mean,na.rm=T) #runs a bit (1 minute?)
totAF$AFD<-abs(totAF$ancAF-totAF$derAF)

hist(totAF$AFD,breaks = 100)
summary(totAF$AFD)
# Perform the outlier scan
outSNPs<-subset(totAF,totAF$AFD >= quantile(totAF$AFD,probs = tresholdSNP,na.rm = T))
print(paste0("selecting SNPs with AFD > ", quantile(totAF$AFD,probs = tresholdSNP,na.rm = T)))
write.table(x = outSNPs, file = paste("results/",adaptation,project,".outSNPs.",tresholdSNP,".txt",sep=""),quote = F,row.names = F,sep="\t")

# Distribution of N outlier SNPs/gene, given the size of region covered by data.
outGenes<-as.data.frame(table(outSNPs$ID))
outGenes$start<-""
outGenes$end<-""
for (id in outGenes$Var1) { # id = "AL4G18690"
  idInfo<-subset(totAF, totAF$ID %in% id)
  outGenes$start[which(outGenes$Var1 %in% id)]<-unlist(c(idInfo[1,"start"]))[[1]]
  outGenes$end[which(outGenes$Var1 %in% id)] <-idInfo[nrow(idInfo),"start"]
}
outGenes$length<-as.numeric(outGenes$end)-as.numeric(outGenes$start)
outGenes$density<-as.numeric(outGenes$Freq)/as.numeric(outGenes$length)
hist(outGenes$density,breaks = 100)
# Outlier SNP density outliers
outGenes<-outGenes[ order(outGenes[,"density"],decreasing = T), ]
outGenes<-outGenes[1:(nrow(outGenes)*tresholdGene),]
outGenes<-subset(outGenes,Freq>2) 
print(nrow(outGenes))
outGenes <- apply(outGenes,2,as.character)
write.table(x = as.data.frame(outGenes), file = paste("results/",adaptation,project,".outGenes.",tresholdSNP,".",tresholdGene,".txt",sep=""),quote = F,sep="\t", row.names = F)


### PART 5: VISUALIZATION OF CANDIDATE GENES SCAN ###
# follows up directly after Part 4. 
ann<-fread("data/LyV2_TAIR11orth_des_20171231.txt") # gene annotation and description
# Plot exploratory heatmap
outGenes<-fread(paste("results/",adaptation,project,".outGenes.",tresholdSNP,".",tresholdGene,".txt",sep=""),h=T)
outSNPs<-fread(paste("results/",adaptation,project,".outSNPs.",tresholdSNP,".txt",sep=""))

popInfo<-subset(popInfo,popInfo$pop %in% colnames(totAF)[5:c(length(colnames(totAF)) - 3)])
popInfo<-popInfo[ order(popInfo %>% dplyr:: select(all_of(adaptation)),decreasing = T)]

# Final tables
plottingPol<-dplyr::select(totAF,"ID","aas","ann","start",popInfo$pop)

plottingPolOutl<-subset(plottingPol,paste(plottingPol$ID,plottingPol$start) %in% paste(outSNPs$ID, outSNPs$start))

# Repolarise 
anc<-dplyr::select(plottingPolOutl,all_of(ancestralEnv))
anc$af<-apply(anc,1,mean,na.rm=T)
der<-dplyr::select(plottingPolOutl,all_of(derivedEnv))
der$af<-apply(X = der,1,mean,na.rm=T)

for (i in  1:nrow(plottingPolOutl)){ # i=1
  if (anc$af[i]>der$af[i]) 
  {plottingPolOutl[i,5:(length(plottingPolOutl[i,]))]<-(1-plottingPolOutl[i,5:(length(plottingPolOutl[i,]))])
  plottingPolOutl[i,2]<-paste("R.",plottingPolOutl[i,2],sep = "")
  print(i)
  } else {}}

my_palette <- colorRampPalette(c("gold", "forestgreen", "blue2"))(n = 100) #
my_paletteAltit <- colorRampPalette(c("grey", "black"))(n = ncol(plottingPolOutl)-4)

# Plot gene by gene
pdf(paste("results/",adaptation,project,".heatmapCandidSNPs",tresholdSNP,".",tresholdGene,".pdf",sep=""),height = 11,width = 13)
for (id in unique(outGenes$Var1)) { # id = outGenes$Var1[1]
  annId<-subset(ann,ann$AL %in% id)
  plottingCanddf<-subset(plottingPolOutl,plottingPolOutl$ID %in% id & paste(plottingPolOutl$ID,plottingPolOutl$start) %in% paste(outSNPs$ID,outSNPs$start))
  plottingCanddf<-as.matrix(plottingCanddf[,5:ncol(plottingCanddf)],rownames = plottingCanddf$aas)
  
  heatmap.2(x = plottingCanddf,dendrogram = "none",Colv="NA", Rowv="NA",key = F,col=my_palette,colsep= c(0,11,ncol(plottingCanddf)),sepcolor= c("black"),sepwidth = c(0.02),trace="none",ColSideColors = my_paletteAltit,labCol=popInfo$pop,colCol= "black",offsetRow = c(0.05),offsetCol= c(0.05),cexRow = 1,cexCol = 1.6,margins = c(7,7),lhei = c(0.2,9),lwid = c(0.03,0.8),xlab = paste(id,substr(annId$AT,1,9),annId$annShort,sep=" - "))
  
}
dev.off()





