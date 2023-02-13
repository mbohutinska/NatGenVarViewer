############# The natural genetic variability viewer #############

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
project<-"cuticle26Jan" # project name will appear in the output tables name
tresholdSNP<-0.99 # outlier SNPs above this threshold, only for Part 2 & 3
tresholdGene<-0.5 # outlier SNPs density per gene cutoff, only for Part 2 & 3

# User-supplied datasets
name<-fread("data/PopKey613.csv",h=F) # get individual names
name<-c(name[2:nrow(name),1])[[1]] 
geneSet<-read.table("geneSet.txt",h=F) # user supplied set of genes for the analysis. My exampe involves set of cutitule forming genes.
thaliana<-T # genes defined using A. thaliana identifiers (for A. lyrata set to FALSE)
popInfo<-fread("data/popInfo.txt",h=T,sep="\t") # names and environmental variables for populations

derivedEnv<-c("BAL","LAC","WIL","ZEP","SUB","SCH","TKO","TRT") #populations used for selections scan. In this case alpine populaions. only for Part 2 & 3
ancestralEnv<-c("DRA","TIS","KAS","SUB","BAB","HRA","SPI","ING") #populations used for selections scan. In this case foothill populaions. only for Part 2 & 3
#Alternatively define based on environmental values from popInfo.txt:
derivedEnv<-subset(popInfo,popInfo$Altitude>1600)$pop
ancestralEnv<-subset(popInfo,popInfo$Altitude<900 & popInfo$Altitude>400)$pop

ann<-fread("data/LyV2_TAIR11orth_des_20171231.txt") # gene annotation and description
geneSetDataGenrated<-T # if you already run this pipeline previously and generated data for the subset of your genes, set this to T to save time
repolarise<-T # False saves time but makes exploratory heatmaps more difficult to read. If the gene set is reasonably big, I recomment to set it to T
x="data/ALL613.table.recode.txt" #All variants of 613 individuals of A. arenosa. Filtered for MAF = 0.05, min depth per genotype > 4 and mafimum fraction of filtered genotypes = 0.4



if (geneSetDataGenrated!=T) {
  tot<-fread(x,h=F,na.strings = "-9",nThread = 3)
  # Add header 
  colnames(tot) <- c("pop","ploidy","scaff","start","AN","DP","ID","ann","aas",name,"nan")
  # Convert A. thaliana to A. lyrata
  if (thaliana==T) {
    dict<-fread("data/ALATdict.txt")
    geneSet<-subset(dict, substr(dict$AT,1,9) %in% geneSet$V1)
    geneSet<-geneSet[ order(geneSet$AL)]
    geneSet<-geneSet[!duplicated(geneSet$AL)]
    geneSet<-geneSet$AL
  } else {
    geneSet<-geneSet$V1[ order(geneSet$V1)]
    geneSet<-geneSet[!duplicated(geneSet)]
  }
  # Subset data from genes of interest
  geneSetData<-subset(tot,tot$ID %in% geneSet)
  write.table(geneSetData,paste("results/",adaptation,"_",project,"rawGT.txt",sep = ""),row.names = F,quote = F,sep = "\t")
  rm(tot)
} else {geneSetData<-fread(paste("results/",adaptation,"_",project,"rawGT.txt",sep = ""),h=T)} #read geneSetData


### PART 1: EXPLORATION OF A GIVEN SET OF GENES ###
# Meaningful for any set of genes, gives visualization of allele frequencies.
# Calculate AF for each lineage


popnames<-unique(substr(name,1,3))
for (pop in popnames) { # pop="BAB"
  #define ploidy
  if (substr(subset(name,name %like% pop)[1],7,7)=="d") {
    ploidy=2
  } else {ploidy=4}
  afh1 <- geneSetData %>% dplyr:: select(starts_with(pop))
  geneSetData$ACh1<-rowSums(afh1,na.rm = T)
  geneSetData$NAh1<-apply(is.na(afh1), 1, sum)
  geneSetData$ANh1<-(ncol(afh1)-geneSetData$NAh1)*ploidy   
  geneSetData$POP<-geneSetData$ACh1/geneSetData$ANh1
  colnames(geneSetData)[which(colnames(geneSetData) %in% "POP")]<-pop}
# Selecting order for plotting
popInfo<-subset(popInfo,popInfo$pop %in% popnames)
popInfo<-popInfo[ order(popInfo %>% dplyr:: select(adaptation),decreasing = T)]

# Final tables
plottingPol<-dplyr::select(geneSetData,"ID","aas","ann","start",popInfo$pop)

# Repolarise 
if (repolarise == T) {
  for (i in  1:nrow(plottingPol)){ # i=1
    if (mean(as.numeric(c(plottingPol[i,5:ncol(plottingPol)])),na.rm = T) > 0.5) 
    {plottingPol[i,5:(length(plottingPol[i,]))]<-(1-plottingPol[i,5:(length(plottingPol[i,]))])
    plottingPol[i,2]<-paste("R.",plottingPol[i,2],sep = "")
    print(i)} else {}}
} else{print("no repolarisation of exploratory heatmap")}

write.table(plottingPol,paste("results/",adaptation,"_",project,"AF.txt",sep = ""),row.names = F,quote = F,sep = "\t")
rm(geneSetData)

# Plot exploratory heatmap
my_palette <- colorRampPalette(c("gold", "forestgreen", "blue2"))(n = 100)
my_paletteAltit <- colorRampPalette(c("grey", "black"))(n = ncol(plottingPol)-4)

# Plot
pdf(paste("results/",adaptation,project,".heatmapExplor.pdf",sep=""),height = 11,width = 13)
for (id in unique(plottingPol$ID)) { # id = "AL1G10770"
  annId<-subset(ann,ann$AL %in% id)
  plottingPoldf<-subset(plottingPol, plottingPol$ID %in% id)
  if (nrow(plottingPoldf)>1) {
    plottingPoldfAAS<-subset(plottingPoldf, plottingPoldf$ann %in% "missense_variant")
    plottingPoldf<-as.matrix(plottingPoldf[,5:ncol(plottingPoldf)],rownames = plottingPoldf$aas)
    if (nrow(plottingPoldfAAS)>1){
      plottingPoldfAAS<-as.matrix(plottingPoldfAAS[,5:ncol(plottingPoldfAAS)],rownames = plottingPoldfAAS$aas)}
    # Plot amino acid substitutions  
    if (nrow(plottingPoldfAAS)>1) {
      heatmap.2(x = plottingPoldfAAS,dendrogram = "none",Colv="NA", Rowv="NA",key = F,col=my_palette,colsep= c(0,11,ncol(plottingPoldf)),sepcolor= c("black"),sepwidth = c(0.02),trace="none",ColSideColors = my_paletteAltit,labCol=popInfo$pop,offsetRow = c(0.05),offsetCol= c(0.05),cexRow = 1.1,cexCol = 1.6,margins = c(7,7),lhei = c(0.2,9),lwid = c(0.03,0.8),xlab = paste(id,substr(annId$AT,1,9),annId$annShort,sep=" - "))} 
    else {print(paste("skipping",id))}
    # Plot all  
    heatmap.2(x = plottingPoldf,dendrogram = "none",Colv="NA", Rowv="NA",key = F,col=my_palette,colsep= c(0,11,ncol(plottingPoldf)),sepcolor= c("black"),sepwidth = c(0.02),trace="none",ColSideColors = my_paletteAltit,labCol=popInfo$pop,colCol= "black",offsetRow = c(0.05),offsetCol= c(0.05),cexRow = 1,cexCol = 1.6,margins = c(7,7),lhei = c(0.2,9),lwid = c(0.03,0.8),xlab = paste(id,substr(annId$AT,1,9),annId$annShort,sep=" - ")) 
  } else {print(paste("skipping",id))}
} 
dev.off()
#

#}




### PART 2: SELECTION SCAN ###
# this is only meaningful to run given sufficiently big set of genes. I recommend at least 40, optimally more than 100. 
geneSetData<-fread(paste("results/",adaptation,"_",project,"rawGT.txt",sep = ""),h=T)

# Calculate AF for each lineage
geneSetData$ancAC<-0 #ancestral
geneSetData$ancAN<-0
for (pop in ancestralEnv) { # pop=ancestralEnv[1]
  #define ploidy
  if (substr(subset(name,name %like% pop)[1],7,7)=="d") {
    ploidy=2
  } else {ploidy=4}
  afh1 <- geneSetData[,1:length(name)] %>% dplyr:: select(starts_with(pop))
  geneSetData$ancAC<-as.numeric(geneSetData$ancAC)+rowSums(afh1,na.rm = T)
  geneSetData$NAh1<-apply(is.na(afh1), 1, sum)
  geneSetData$ancAN<-as.numeric(geneSetData$ancAN)+((ncol(afh1)-geneSetData$NAh1)*ploidy)}

geneSetData$derAC<-0 #derived
geneSetData$derAN<-0
for (pop in derivedEnv) {
  #define ploidy
  if (substr(subset(name,name %like% pop)[1],7,7)=="d") {
    ploidy=2
  } else {ploidy=4}
  afh1 <- geneSetData[,1:length(name)] %>% dplyr:: select(starts_with(pop))
  geneSetData$derAC<-as.numeric(geneSetData$derAC)+rowSums(afh1,na.rm = T)
  geneSetData$NAh1<-apply(is.na(afh1), 1, sum)
  geneSetData$derAN<-as.numeric(geneSetData$derAN)+((ncol(afh1)-geneSetData$NAh1)*ploidy)}
# Compute allele frequency difference (AFD)
geneSetData$AFD<-abs((geneSetData$derAC/geneSetData$derAN) - (geneSetData$ancAC/geneSetData$ancAN))
hist(geneSetData$AFD,breaks = 100)
summary(geneSetData$AFD)
# Perform the outlier scan
outSNPs<-subset(geneSetData,geneSetData$AFD >= quantile(geneSetData$AFD,probs = tresholdSNP,na.rm = T))
write.table(x = outSNPs, file = paste("results/",adaptation,project,".outSNPs.",tresholdSNP,".txt",sep=""),quote = F,row.names = F,sep="\t")

# Distribution of N outlier SNPs/gene
outGenes<-as.data.frame(table(outSNPs$ID))
outGenes$start<-c()
outGenes$end<-c()
for (id in outGenes$Var1) { # id = droplevels(outGenes[1,1])
  idInfo<-subset(geneSetData, geneSetData$ID %in% id)
  outGenes$start[which(outGenes$Var1 %in% id)]<-idInfo[1,"start"]
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


### PART 3: VISUALIZATION OF CANDIDATE GENES ###
# follows up directly after Part 2. 
# Plot exploratory heatmap
outGenes<-fread(paste("results/",adaptation,project,".outGenes.",tresholdSNP,".",tresholdGene,".txt",sep=""),h=T)
outSNPs<-fread(paste("results/",adaptation,project,".outSNPs.",tresholdSNP,".txt",sep=""))
plottingPol<-fread(paste("results/",adaptation,"_",project,"AF.txt",sep = ""),h=T)
plottingPolOutl<-subset(plottingPol,paste(plottingPol$ID,plottingPol$start) %in% paste(outSNPs$ID, outSNPs$start))

# Repolarise 
anc<-dplyr::select(plottingPolOutl,ancestralEnv)
anc$af<-apply(anc,1,mean,na.rm=T)
der<-dplyr::select(plottingPolOutl,derivedEnv)
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


# Plot all candidates together: excluding synonymous
plottingPolOutl2<-subset(plottingPolOutl,!plottingPolOutl$ann %in% "synonymous_variant")

if (nrow(plottingPolOutl)<20) {
  pdf(paste("results/",adaptation,project,".heatmapCandidSNPsAll",tresholdSNP,".",tresholdGene,".pdf",sep=""),height = 11,width = 14)
} else{pdf(paste("results/",adaptation,project,".heatmapCandidSNPsall",tresholdSNP,".",tresholdGene,".pdf",sep=""),height = nrow(plottingPolOutl)/4,width = 14)}

  plottingCanddf<-as.matrix(plottingPolOutl2[,5:ncol(plottingPolOutl2)],rownames = paste(substr(plottingPolOutl2$ann,1,8),plottingPolOutl2$ID,sep = "."))
  heatmap.2(x = plottingCanddf,dendrogram = "none",Colv="NA", Rowv="NA",key = F,col=my_palette,colsep= c(0,11,ncol(plottingCanddf)),sepcolor= c("black"),sepwidth = c(0.02),trace="none",ColSideColors = my_paletteAltit,labCol=popInfo$pop,colCol= "black",offsetRow = c(0.05),offsetCol= c(0.05),cexRow = 1,cexCol = 1.6,margins = c(7,12),lhei = c(0.2,9),lwid = c(0.03,0.8),xlab = paste("All candidate SNPs, excluding synonymous"))
dev.off()


# Plot all candidates together: excluding synonymous
plottingPolOutl2<-subset(plottingPolOutl,plottingPolOutl$ann %in% "missense_variant")

if (nrow(plottingPolOutl2)<20) {
  pdf(paste("results/",adaptation,project,".heatmapCandidSNPsMissense",tresholdSNP,".",tresholdGene,".pdf",sep=""),height = 11,width = 14)
} else{pdf(paste("results/",adaptation,project,".heatmapCandidSNPsMissense",tresholdSNP,".",tresholdGene,".pdf",sep=""),height = nrow(plottingPolOutl2)/4,width = 14)}
par(mar=c(5,4,4,20) + 0.1)
plottingCanddf<-as.matrix(plottingPolOutl2[,5:ncol(plottingPolOutl2)],rownames = paste(plottingPolOutl2$aas,plottingPolOutl2$ID,sep = "."))
heatmap.2(x = plottingCanddf,dendrogram = "none",Colv="NA", Rowv="NA",key = F,col=my_palette,colsep= c(0,11,ncol(plottingCanddf)),sepcolor= c("black"),sepwidth = c(0.02),trace="none",ColSideColors = my_paletteAltit,labCol=popInfo$pop,colCol= "black",offsetRow = c(0.05),offsetCol= c(0.05),cexRow = 1,cexCol = 1.6,margins = c(7,12),lhei = c(0.2,9),lwid = c(0.03,0.8),xlab = paste('All candidate nonsynonymous SNPs'))
dev.off()





aaa<-aaa[ order(aaa$AL)]
aaa<-aaa[!duplicated(aaa$AL)]
write.table(aaa, "geneTableAnnCutecleFinal.txt", quote =F, sep = "\t", row.names = F)






