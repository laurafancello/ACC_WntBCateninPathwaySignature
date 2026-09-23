library(ggplot2)
library(ggrepel)
library(DESeq2)
library(fgsea)
library(writexl)
library(msigdbr)
library(gridExtra)
library(stringr)
options(digits=15)

source("./Scripts/Miscellaneous_Functions.R")

#  Settings ------------------------------------------------------------
inPath <- "./ACC_datasets_formatted/"
out <- "./DEanalyses_ACC/"
if (!(file.exists(out))){
  dir.create(file.path(out))
}
ACC_signature = c("AFF3","AXIN2","FSCN1","LEF1","MYC","PDE2A","SP5","TBX3","TNFRSF19")
dataset_name <- "tcga"
cutoffType <- "Median" # Set cutoff type

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

# Generate normalized counts and separate by LEF1 or signature level ----------------------------------------------
rna <- read.table(file=paste0(inPath, "/tcga/RNAseq_STAR_RawCounts_fromTCGABiolinks_GeneSymbols.txt"), sep="\t", header=T)
colnames(rna) <- str_replace_all(colnames(rna), '\\.', "-")

### Get normalized counts
sample_info <- as.data.frame(colnames(rna))
colnames(sample_info) <- "patientId" 
sample_info$patientId <- factor(sample_info$patientId)
dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = rna,
  colData = sample_info,
  design = ~ patientId)
dds <- estimateSizeFactors(dds)
deseqNorm <- counts(dds, normalized=TRUE)
norm_log2 <- as.data.frame(apply(deseqNorm, 2, function(x) log2(x+1)))
GeneSymbol <- rownames(deseqNorm)
norm_log2 <- as.data.frame(cbind(GeneSymbol,norm_log2))

### Split samples into groups based on cutoff
if(marker == "LEF1"){
  select <- norm_log2[norm_log2$GeneSymbol %in% c(marker),]
  select <- as.matrix(select[,-1])
  select <- apply(select, 2, as.numeric)
  vec <- as.data.frame(cbind(names(select), select))
  colnames(vec) <- c("TumorID",marker)
  vec[,marker] <- as.numeric(as.vector(vec[,marker]))
  cutoff <- round(median(vec[,marker]), digits=2)
  print(paste0("LEF1 median expression in ACC TCGA (after DESeq norm and log2 transformation): ", cutoff))
  vec$Level <- NA
  vec[vec$LEF1 > cutoff, ]$Level <- "high"
  vec[vec$LEF1 <= cutoff, ]$Level <- "low"
} 
if(marker == "signature"){
  select <- norm_log2[norm_log2$GeneSymbol %in% ACC_signature,]
  select <- as.matrix(select[,-which(colnames(select) %in% c("GeneSymbol"))])
  select <- apply(select, 2, as.numeric)
  
  vec <- apply(select, 2, mean)
  vec <- as.data.frame(cbind(names(vec), vec))
  colnames(vec) <- c("TumorID",marker)
  vec[,marker] <- as.numeric(as.vector(vec[,marker]))
  cutoff <- round(median(vec[,marker]), digits=2)
  print(paste0("Signature median expression in ACC tcga (after norm and log2 transformation): ", cutoff))
  vec$Level <- NA
  vec[vec$signature > cutoff, ]$Level <- "high"
  vec[vec$signature <= cutoff, ]$Level <- "low"
}

# Differential expression analysis ----------------------------
rna <- read.table(file=paste0(inPath, "/tcga/RNAseq_STAR_RawCounts_fromTCGABiolinks_GeneSymbols.txt"), sep="\t", header=T)
colnames(rna) <- stringr::str_replace_all(colnames(rna), '\\.', "-")
rna <- as.matrix(rna)   # Format counts matrix for DE analysis
  
sample_info <- as.data.frame(colnames(rna))
colnames(sample_info) <- "TumorID"
sample_info <- dplyr::left_join(sample_info, vec, by="TumorID")
sample_info$Level <- factor(sample_info$Level)
sample_info <- sample_info[match(colnames(rna), sample_info$TumorID), ] # Order metadata as col order in counts object
  
if(identical(colnames(rna), sample_info$TumorID)){ # Check that counts and sample info objects have same patient order
    
  ### DE analysis
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = rna,
    colData = sample_info,
    design = ~ Level)
  dds$Level <- relevel(dds$Level, ref = "low")
  dds <- DESeq(dds)
  res <- results(dds)
  res <- res[!(is.na(res$"padj")),]
  write.csv(file=paste0(outPath, dataset_name, "_DE_", marker, cutoffType, ".csv"), res)
  
  ### Volcano plots
  outFile <- paste0(dataset_name, "_DE_", marker, cutoffType, "_volcanoPlot_L2FC05_pFDR01")
  dataset <- dataset_name
  analysis <- cutoffType
  res <- read.csv(file=paste0(outPath, dataset_name, "_DE_", marker, cutoffType, ".csv"))
  colnames(res) <- c("gene","baseMean","logFC","lfcSE","stat","p.value","pFDR")
  plots <- volcanoPlots_BcatTargets(outPath, res, dataset, analysis)
  png(paste0(outPath, "Volcano_",outFile, "_SignificantBCateninTargets.png"), width=900, height = 900)
  print(plots$BcatSignificant)
  dev.off()
  # png(paste0(outPath, "Volcano_",outFile, "_9SignatureBcatTargets.png"), width=900, height = 900)
  # print(plots$NineSignatureTargets)
  # dev.off()
  pdf(paste0(outPath, "Volcano_",outFile, "_SignificantBCateninTargets.pdf"), width=10, height = 10, useDingbats = FALSE)
  print(plots$BcatSignificant)
  dev.off()
  #pdf(paste0(outPath, "Volcano_",outFile, "_9SignatureBcatTargets.pdf"), width=10, height = 10, useDingbats = FALSE)
  #print(plots$NineSignatureTargets)
  #dev.off()
  
  ### Output sign up and down genes for iRegulon
  if (!(file.exists(paste0(outPath, "/SignificantUpAndDownDEgenes_iRegulon/")))){
    dir.create(file.path(paste0(outPath, "/SignificantUpAndDownDEgenes_iRegulon/")))
  }
  sign <- res[res$pFDR<0.1,]
  sign <- sign[order(sign$pFDR),]
  if (!(file.exists(paste0(outPath, "/SignificantUpAndDownDEgenes_iRegulon/")))){
    dir.create(file.path(paste0(outPath, "/SignificantUpAndDownDEgenes_iRegulon/")))
  }
  write.table(file=paste0(outPath, "/SignificantUpAndDownDEgenes_iRegulon/SignificantDEgenesUp_", dataset_name, "_", marker, cutoffType, ".txt"), sign[sign$logFC>0,'gene'], col.names = F, row.names = F, quote = F, sep="\n")
  write.table(file=paste0(outPath, "/SignificantUpAndDownDEgenes_iRegulon/SignificantDEgenesDown_", dataset_name, "_", marker, cutoffType, ".txt"), sign[sign$logFC<0,'gene'], col.names = F, row.names = F, quote = F, sep="\n")
  
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
colnames(res) <- c("ID","baseMean","logFC","lfcSE","stat","p.value","pFDR")

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
