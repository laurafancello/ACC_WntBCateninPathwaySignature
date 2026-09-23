library(GEOquery)
library(ggplot2)
library(ggrepel)
library(SummarizedExperiment)
library(msigdbr)
library(fgsea)
options(digits=15)
source("./Scripts/Miscellaneous_Functions.R")
inPath <- "./ACC_datasets_formatted/MicroarrayAnnotations/"
out <- "./DEanalyses_ACC/"
if (!(file.exists(out))){
  dir.create(file.path(out))
}

ACC_signature = c("AFF3","AXIN2","FSCN1","LEF1","MYC","PDE2A","SP5","TBX3","TNFRSF19")
dataset_name <- "Assie"
cutoffType <- "median"

### Set whether to use LEF1 or signature expression
#marker <- "LEF1"
marker <- "signature"

if(marker == "LEF1"){
  outPath <- paste0(out, "LEF1_median/")
  if (!(file.exists(outPath))){
    dir.create(file.path(outPath))
  }
}
if(marker == "signature"){
  outPath <- paste0(out, "signature_median/")
  if (!(file.exists(outPath))){
    dir.create(file.path(outPath))
  }
}


#  IMPORT DATA FROM GEO ---------------------------------------
# Load series and platform data from GEO
gset <- getGEO("GSE49278", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

### Log2 transformation, if necessary 
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { 
  ex[which(ex <= 0)] <- NaN # data were already log2 transformed: no need to log transform
}

# Log2 norm data with probe and gene ID annotation 
microarray_annotation <- read.table(file=paste0(inPath, "/GPL16686.an.txt"), sep="\t", header=T, fill=T, quote="") # probeID to geneSymbol file from https://gemma.msl.ubc.ca/arrays/showArrayDesign.html?id=787
colnames(microarray_annotation)[c(1,2)] <- c("FROM","TO")
dataset <- cbind(rownames(ex),ex)
colnames(dataset)[1] <- "FROM"
dataset <- merge(dataset, microarray_annotation[,c('FROM','TO')], by="FROM", all=F)
dataset <- dataset[,c(1,ncol(dataset),2:(ncol(dataset)-1))]
colnames(dataset)[1:2] <- c("Probe","Gene")
saveRDS(file="./ACC_datasets_formatted/Assie/Log2NormData_all.RDS", dataset)

# SPLIT INTO GROUPS BASED ON LEF1/SIGNATURE MEDIAN CUTOFF ----------------------------------
dataset <- readRDS("./ACC_datasets_formatted/Assie/Log2NormData_all.RDS")
if(marker == "LEF1"){
  select <- dataset[dataset$Gene %in% c(marker),]
  select <- as.matrix(select[,-which(colnames(select) %in% c("Gene","Probe"))]) # Remove character columns Gene and Probe
  select <- apply(select, 2, as.numeric)
  vec <- as.data.frame(cbind(names(select), select))
  colnames(vec) <- c("TumorID",marker)
  vec[,marker] <- as.numeric(as.vector(vec[,marker]))
  cutoff <- round(median(vec[,marker]), digits=2)
  print(paste0("LEF1 median expression in ACC Assie (after norm and log2 transformation): ", cutoff))
  vec$Level <- NA
  vec[vec$LEF1 > cutoff, ]$Level <- "high"
  vec[vec$LEF1 <= cutoff, ]$Level <- "low"
}

if(marker=="signature"){
  select <- dataset[dataset$Gene %in% ACC_signature,]
  select <- as.matrix(select[,-which(colnames(select) %in% c("Gene","Probe"))])
  select <- apply(select, 2, as.numeric)
  
  vec <- apply(select, 2, mean)
  vec <- as.data.frame(cbind(names(vec), vec))
  colnames(vec) <- c("TumorID",marker)
  vec[,marker] <- as.numeric(as.vector(vec[,marker]))
  cutoff <- round(median(vec[,marker]), digits=2)
  print(paste0("Signature median expression in ACC Assie (after norm and log2 transformation): ", cutoff))
  vec$Level <- NA
  vec[vec$signature > cutoff, ]$Level <- "high"
  vec[vec$signature <= cutoff, ]$Level <- "low"
}


#  PREPARE FOR DE ANALYSIS ------------------------------------------------
dataset <- readRDS(file="./ACC_datasets_formatted/Assie/Log2NormData_all.RDS")
probes <- dataset$Probe
dataset <- as.matrix(dataset[,-which(colnames(dataset) %in% c("Gene","Probe"))])
dataset <- apply(dataset, 2, as.numeric)
rownames(dataset) <- probes
vec <- vec[match(colnames(dataset), vec$TumorID), ] # Order metadata as columns order in affyExpr object
vec$condition <- ifelse(vec$Level == "high", "altered", "wt")
vec$condition <- factor(vec$condition, levels=c("wt","altered"))

if(identical(colnames(dataset), vec$TumorID)){ # Check that counts and metadata have same sample order
  
  ### DE analysis
  DE_result <- DE_microarray(expr=dataset, metadata=vec, condition=vec$condition, microarray_annotation=microarray_annotation)
  #saveRDS(file=paste0(outPath, dataset_name, "_DE_LEF1", cutoffType, "_fullDEoutputObject.RDS"), DE_result)
  res <- as.data.frame(DE_result$affyDEResults)
  write.csv(file=paste0(outPath, dataset_name, "_DE_", marker, cutoffType, ".csv"), res)
  ### Volcano plots
  outFile <- paste0(dataset_name, "_DE_", marker, cutoffType, "_volcanoPlot_L2FC05_pFDR01")
  dataset <- dataset_name
  analysis <- cutoffType
  res <- read.csv(file=paste0(outPath, dataset_name, "_DE_", marker, cutoffType, ".csv"))
  colnames(res) <- c("X","PROBEID","gene","p.value","statistic","logFC","avgExpr","logFCSE","sampleSize","pFDR")
  plots <- volcanoPlots_BcatTargets(outPath, res, dataset, analysis)
  png(paste0(outPath, "Volcano_",outFile, "_BcatSignificant.png"), width=900, height = 900)
  print(plots$BcatSignificant)
  dev.off()
  png(paste0(outPath, "Volcano_",outFile, "_9selectedBcatTargets.png"), width=900, height = 900)
  print(plots$NineSignatureTargets)
  dev.off()
  pdf(paste0(outPath, "Volcano_",outFile, "_BcatSignificant.pdf"), width=10, height = 10, useDingbats = FALSE)
  print(plots$BcatSignificant)
  dev.off()
  #pdf(paste0(outPath, "Volcano_",outFile, "_9selectedBcatTargets.pdf"), width=10, height = 10, useDingbats = FALSE)
  #print(plots$NineSignatureTargets)
  #dev.off()
  
  ### Output for iRegulon
  if (!(file.exists(paste0(outPath, "/SignificantUpAndDownDEgenes_iRegulon/")))){
    dir.create(file.path(paste0(outPath, "/SignificantUpAndDownDEgenes_iRegulon/")))
  }
  sign <- res[res$pFDR<0.1,]
  sign <- sign[order(sign$pFDR),]
  
  sum(sign$logFC > 0, na.rm = TRUE)
  sum(sign$logFC < 0, na.rm = TRUE)
  if (!(file.exists(paste0(outPath, "/SignificantUpAndDownDEgenes_iRegulon/")))){
    dir.create(file.path(paste0(outPath, "/SignificantUpAndDownDEgenes_iRegulon/")))
  }
  write.table(file=paste0(outPath, "/SignificantUpAndDownDEgenes_iRegulon/SignificantDEgenesUp_", dataset_name, "_", marker, cutoffType, ".txt"), sign[sign$logFC>0,'gene'], col.names = F, row.names = F, quote = F, sep="\n")
  write.table(file=paste0(outPath, "/SignificantUpAndDownDEgenes_iRegulon/SignificantDEgenesDown_", dataset_name, "_", marker, cutoffType, ".txt"), sign[sign$logFC<0,'gene'], col.names = F, row.names = F, quote = F, sep="\n")
  
}else{
  print("Counts and metadata sample order is not the same. Correct prior to DE analysis")
}


# GSEAs ----------------------------
if (!(file.exists(paste0(outPath, "/GSEAs/")))){
  dir.create(file.path(paste0(outPath, "/GSEAs/")))
}
analysis <- cutoffType
colNameGenes <- "ID"
RankedListType <- "FoldChangeRanked"

### Generate ranked list of genes (by fold change)
res <- read.csv(file=paste0(outPath, dataset_name, "_DE_", marker, cutoffType,".csv")) 
res$log2FoldChange <- res$logFC 
RankedList <- sign(res$log2FoldChange)*(-(log10(res$p.value)))
names(RankedList) <- res[,colNameGenes]
if(length( which(is.na(RankedList)))>0){
  RankedList <- RankedList[-(which(is.na(RankedList)))]
}
RankedList <- sort(RankedList)
RankedList <- as.data.frame(cbind(names(RankedList),RankedList))
colnames(RankedList) <- c("GeneName","RankedList")
write.table(file=paste0(outPath, "/GSEAs/", marker, analysis, "_",dataset_name, "_", RankedListType,".rnk"), RankedList, col.names = T, row.names = T, sep="\t")

### Run GSEAs
values <- as.numeric(as.vector(RankedList$RankedList))
names(values) <- RankedList$GeneName
RankedList <- values
# Transfrom -Inf into zeros (when calculating log2FC no pseudocount was used)
RankedList[which(RankedList== -Inf)] <- 0

print("Running GSEA on Hallmarks")
outPathHallmarks <- GSEA_MSigDBHallmark_v2(bulkType=dataset_name, analysis=paste0(marker,cutoffType), RankedList=RankedList, RankedListType=RankedListType, outPath=paste0(outPath, "GSEAs/"), geneIdType="gene_symbol")
#write.table(file=paste0(outPath, "GSEAs/", dataset_name, RankedListType,"_", marker, cutoffType, "_Hallmarks.txt"), outPathHallmarks[order(outPathHallmarks$padj),c("pathway","pval","padj","log2err","ES","NES","size")], col.names = T, row.names = F, quote=F, sep="\t")
saveRDS(file=paste0(outPath, "GSEAs/", dataset_name, RankedListType, "_", marker, cutoffType, "_Hallmarks.RDS"), outPathHallmarks[order(outPathHallmarks$padj),])

### Enrichment plots of pathways of interest
h_gene_sets = msigdbr(species = "Homo sapiens", category = "H")
msigdbr_list_H = split(x = h_gene_sets$gene_symbol, f = h_gene_sets$gs_name)
msigdbr_list_H_nr <- lapply(msigdbr_list_H, unique)

selectPathway <- "HALLMARK_WNT_BETA_CATENIN_SIGNALING"
data <- readRDS(file=paste0(outPath, "GSEAs/", dataset_name, RankedListType, "_", marker, cutoffType, "_Hallmarks.RDS"))
sub <- data[data$pathway == selectPathway,]
pval <- signif(sub$pval, digits=3)
padj <- signif(sub$padj, digits=3)
ES <- round(sub$ES, digits=2)
NES <- round(sub$NES, digits=2)
Ranks <- read.table(file=paste0(outPath, "/GSEAs/", marker, analysis, "_",dataset_name, "_", RankedListType,".rnk"), header = T, sep="\t")
Ranks_ok <- as.numeric(as.vector(Ranks$RankedList))
names(Ranks_ok) <- as.character(as.vector(Ranks$GeneName))
png(paste0(outPath, "/GSEAs/", dataset_name, RankedListType, "_", marker, cutoffType, "_HALLMARK_WNT_BETA_CATENIN_SIGNALING.png"))
print(plotEnrichment(msigdbr_list_H_nr[[selectPathway]],
                     Ranks_ok) + labs(title=ggtitle(paste0(selectPathway, ", ", dataset_name, "\n p=", pval, ", padj=", padj, ", ES=", ES, ", NES=", NES))))
dev.off()
