library(ggplot2)
library(ggpubr)
library(gridExtra)
library(ggExtra)
library(ggrepel)

out <- "./SignatureAndLEF1_inACC/output/"

# Read input --------------------------------------------------------------

# Read signature genes
ACC_BcatSignature <- c("AFF3","AXIN2","FSCN1","LEF1","MYC","PDE2A","SP5","TBX3","TNFRSF19")

### Read input --------------------------------------------------------------
### Read Z-score transformed log2 norm counts 
all_data <- readRDS("C:/Users/LF260934/Documents/WORKING_FOLDER/Results/ACC_datasets_integration/Counts_Zscore_allDatasets.RDS")

### Read metadata
Assie_metadata <- read.csv("./ACC_datasets_formatted/Assie/Metadata/Metadata_Assie.csv")
Assie_metadata$BcatStatus <- Assie_metadata$CTNNB1_ZNRF3_Alteration
Assie_metadata <- Assie_metadata[!(is.na(Assie_metadata$BcatStatus)),]
### Separate ZNRF3 deletions or inact mutations for Assie --------------------
Assie_metadata$BcatStatus <- as.character(as.vector(Assie_metadata$BcatStatus))
Assie_metadata$BcatStatus <- stringr::str_replace_all(Assie_metadata$BcatStatus, "ZNRF3", "ZNRF3_mut")
Assie_metadata$BcatStatus <- stringr::str_replace_all(Assie_metadata$BcatStatus, "CTNNB1", "CTNNB1_mut")
Assie_metadata[Assie_metadata$ZNRF3 == "hdel",]$BcatStatus <- "ZNRF3del"
Assie_metadata$BcatStatus <- factor(Assie_metadata$BcatStatus, levels=c("CTNNB1_mut","ZNRF3_mut","ZNRF3del","wt"))

Heaton_metadata <- read.csv("./ACC_datasets_formatted/Heaton/Metadata/Metadata_Heaton.csv")
Heaton_metadata <- Heaton_metadata[Heaton_metadata$Histotype %in% c("ACC"),]
Heaton_metadata <- Heaton_metadata[!(is.na(Heaton_metadata$BetaCateninStaining)),]
Heaton_metadata$BcatStatus <- factor(Heaton_metadata$BetaCateninStaining, levels=c("Nuclear","Membrane"))
Heaton_metadata <- Heaton_metadata %>% mutate(OS5years = ifelse(OS_Years <= 5 & OS_Status == "dead", 1, 0))
Heaton_metadata[Heaton_metadata$OS_Years %in% c("unknown"),]$OS5years <- NA

tcga_metadata <- read.csv("./ACC_datasets_formatted/tcga/Metadata/Metadata_tcga.csv")
tcga_Bcat <- readRDS(file="./ACC_datasets_formatted/tcga/Metadata/tcga_sampleIDs_CTNNB1mut_ZNRF3mut_wt.RDS")
colnames(tcga_Bcat)[1] <- "TumorID"
tcga_Bcat_ZNRF3details <- read.table(file="./ACC_datasets_formatted/tcga/Metadata/CuratedStatus_MUTandCN.txt", header=T, sep="\t")
tcga_Bcat_ZNRF3details[((tcga_Bcat_ZNRF3details$Manual.curation.boxplot=="ZNRF3_mut")&(is.na(tcga_Bcat_ZNRF3details$gene_proteinChange))),]$Manual.curation.boxplot <- "ZNRF3del"
tcga_Bcat_ZNRF3details <- tcga_Bcat_ZNRF3details[,c("sampleId","Manual.curation.boxplot")]
colnames(tcga_Bcat_ZNRF3details) <- c("TumorID","Alteration")
tcga_Bcat_ZNRF3details_ok <- rbind(tcga_Bcat[tcga_Bcat$Alteration=="wt",], tcga_Bcat_ZNRF3details)
tcga_metadata <- merge(tcga_metadata, tcga_Bcat_ZNRF3details_ok, all=F, by="TumorID")
tcga_metadata$BcatStatus <- factor(tcga_metadata$Alteration, levels=c("CTNNB1_mut","ZNRF3_mut","ZNRF3del","wt"))
removeCols <- which(colnames(tcga_metadata) %in% c("proteinChange","gene_proteinChange")) # remove columns proteinChange and gene_proteinChange which only refer to TP53 mutations and cause duplcation of "TCGA-OR-A5JA" entry (because harboring 2 TP53 mutations)
tcga_metadata <- unique(tcga_metadata[,-removeCols])

# Function to plot signature & LEF1 expression with known pathway activation status and aggressiveness --------------------
plot_BcatSignature_LEF1_5yearsOS <- function(dataset, dataset_name, signature, signatureType, corrType, 
                                             LEF1cutoffLow, LEF1cutoffHigh, LEF1_1cutoffEM, Signature_1cutoffEM,
                                             metadata, label=NULL){
  
  print(dataset_name)
  
  select <- dataset[dataset$GeneSymbol %in% signature,]
  select <- as.matrix(select[,-1])
  select <- apply(select, 2, as.numeric)
  
  mean <- apply(select, 2, mean)
  mean <- as.data.frame(cbind(names(mean), mean))
  colnames(mean) <- c("TumorID","BcatSignatureMean")
  
  LEF1 <- dataset[dataset$GeneSymbol %in% c("LEF1"),-1]
  LEF1 <- t(LEF1)
  LEF1 <- as.data.frame(cbind(rownames(LEF1), LEF1))
  colnames(LEF1) <- c("TumorID","LEF1expr")
  LEF1$LEF1expr <- as.numeric(as.vector(LEF1$LEF1expr))
  
  all <- merge(mean, LEF1, by="TumorID")
  all$LEF1expr <- as.numeric(as.vector(all$LEF1expr))
  all$BcatSignatureMean <- as.numeric(as.vector(all$BcatSignatureMean))
  
  all <- merge(all, metadata, by="TumorID", all=F)
  if(length(intersect(unique(all$BcatStatus), c("CTNNB1","ZNRF3","wt")))==3){
   all$BcatStatus <- factor(all$BcatStatus, levels=c("CTNNB1","ZNRF3","wt"))
   colors_manual <- c('red','orange', 'green4')
  }
  if((length(unique(all$BcatStatus))==3)&(length(intersect(unique(all$BcatStatus), c("CTNNB1_mut","ZNRF3_mut","wt")))==3)){
    all$BcatStatus <- factor(all$BcatStatus, levels=c("CTNNB1_mut","ZNRF3_mut","wt"))
    colors_manual <- c('red','orange', 'green4')
  }
  if((length(unique(all$BcatStatus))==4)&(length(intersect(unique(all$BcatStatus), c("CTNNB1_mut","ZNRF3_mut","ZNRF3del","wt")))==4)){
    all$BcatStatus <- factor(all$BcatStatus, levels=c("CTNNB1_mut","ZNRF3_mut","ZNRF3del","wt"))
    colors_manual <- c('red','orange', 'yellow1', 'green4')
  }
  if(length(intersect(unique(all$BcatStatus), c("Nuclear","Membrane")))==2){
    all$BcatStatus <- factor(all$BcatStatus, levels=c("Nuclear","Membrane"))
    colors_manual <- c('red', 'green4')
  }
  
  if(nrow(all[is.na(all$OS5years),])>0){
    all[is.na(all$OS5years),]$OS5years <- "unknown"
    all$OS5years <- factor(all$OS5years, levels=c("0","1","unknown"))
  }else{
    all$OS5years <- factor(all$OS5years, levels=c("0","1"))
  }
  
  all$label <- ""
  if(!(is.null(label))){
    all[all$TumorID %in% c(label),]$label <- all[all$TumorID %in% c(label),]$TumorID
  }
  
  cor <- cor.test(all$LEF1expr, all$BcatSignatureMean, method=corrType)
  if((is.null(LEF1cutoffLow))&(is.null(LEF1cutoffHigh))&(is.null(LEF1_1cutoffEM))&(is.null(Signature_1cutoffEM))){
    p <- ggplot(data=all, aes(x=LEF1expr,y=BcatSignatureMean)) +
      geom_point(aes(color=BcatStatus, shape=OS5years),size=3) +
      geom_vline(xintercept=LEF1cutoffLow, linetype=2, colour="gray40") + geom_vline(xintercept=LEF1cutoffHigh, linetype=2, colour="gray40")  +
      geom_vline(xintercept=LEF1_1cutoffEM, linetype=2,colour="gray40") +
      geom_hline(yintercept=Signature_1cutoffEM, linetype=2,colour="gray40") +
      scale_color_manual(values=colors_manual) +
      geom_smooth(data=all, aes(x=LEF1expr,y=BcatSignatureMean),method="loess") +
      geom_text_repel(aes(label=label), min.segment.length = 0.2, max.overlaps = 20) +
      theme_bw() + theme(legend.position="bottom") +
      ggtitle(paste0(dataset_name,",",signatureType,"\n",corrType,"=",round(cor$estimate,digits=2),",p=", signif(cor$p.value, digits=2)))
  }else{
    p <- ggplot(data=all, aes(x=LEF1expr,y=BcatSignatureMean)) +
      geom_point(aes(color=BcatStatus, shape=OS5years),size=3) +
      scale_color_manual(values=colors_manual) +
      geom_smooth(data=all, aes(x=LEF1expr,y=BcatSignatureMean),method="loess") +
      geom_text_repel(aes(label=label), min.segment.length = 0.2, max.overlaps = 20) +
      theme_bw() + theme(legend.position="bottom") +
      ggtitle(paste0(dataset_name,",",signatureType,"\n",corrType,"=",round(cor$estimate,digits=2),",p=", signif(cor$p.value, digits=2)))
  }
  ggMarginal(p, type="densigram", fill = "slateblue")
}


# Plot signature-LEF1 expression correlation with known pathway activation status and aggressiveness ----------------
### Plot ----------------
### Do not plot any cutoff
LEF1cutoffLow <- NULL
LEF1cutoffHigh <- NULL
LEF1_1cutoffEM <- NULL
Signature_1cutoffEM <- NULL

signature_select <- ACC_BcatSignature
signatureType <- "ACC_signature"

for(dataset_name in c("Assie","Heaton","tcga")){
  print(dataset_name)
  
  genes <- all_data$GeneSymbol
  dataset <- all_data[,grep(dataset_name,colnames(all_data))]
  colnames(dataset) <- stringr::str_replace_all(colnames(dataset), paste0(dataset_name, "_"), "")
  dataset$GeneSymbol <- genes
  
  if(dataset_name == "Assie"){
    metadata_obj <- Assie_metadata
  }
  if(dataset_name == "Heaton"){
    metadata_obj <- Heaton_metadata
  }
  if(dataset_name == "tcga"){
    metadata_obj <- tcga_metadata
  }
  
  png(paste0(out,"/BcateninSignatures_LEF1expr_",signatureType,"_",dataset_name,"_NoCutoffShown_5yearsOS_LOESS_AvgProbes.png"), width=500, height=500)
  print(plot_BcatSignature_LEF1_5yearsOS(dataset=dataset, dataset_name=dataset_name, signatureType=signatureType, corrType="pearson", signature=signature_select,
                                         LEF1cutoffLow=LEF1cutoffLow, LEF1cutoffHigh=LEF1cutoffHigh, LEF1_1cutoffEM=LEF1_1cutoffEM, Signature_1cutoffEM=Signature_1cutoffEM, metadata=metadata_obj) )
  dev.off()
}

