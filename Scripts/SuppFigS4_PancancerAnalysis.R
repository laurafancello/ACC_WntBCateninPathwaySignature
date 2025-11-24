library(ggplot2)
library(ggpubr)
library(gridExtra)
library(ggExtra)
library(cutpointr)
library(ggrepel)
library(survival)
library(survminer)
library(ggsurvfit)
library(GEOquery)
library(DESeq2)

source("./Scripts/Miscellaneous_Functions.R")

inpath <- "./Pancancer/input/"
out <- "./Pancancer/output/"

# SETTINGS ----------------------------------------------------------------
analysis <- "onlyCTNNB1MutWt" # compare tumors with CTNNNB1 mutations vs tumors wt for CTNNB1. Remove ZNRF3-altered tumors.
metric <- "sens_constrain"

# DE ANALYSIS & HEATMAP FC ALL B-CATENIN TARGETS ---------------------
if (!(file.exists(paste0(out, "/DEanalysis/")))){
  dir.create(file.path(paste0(out, "/DEanalysis/")))
}
# DE analysis CTNNB1 mut vs wt --------------------------------------------
for(cancer in c("LIHC","UCEC","PAAD","SKCM","STAD","COAD")){
  print(cancer)
  
  ### Prepare input for DE analysis ----------------------------------
  rna_tumor  <- readRDS(file=paste0(inpath, "TCGA-",cancer,"_RNA_STARcountsMatrix_GeneSymbol_onlyTumors.RDS")) # raw counts from TCGA
  altered_vs_wt <- read.table(file=paste0(inpath, cancer, "_sample_matrix.txt"), header=T, sep="\t")
  altered_vs_wt$sum <- rowSums(altered_vs_wt[,-1])  
  altered_vs_wt$BcatMut <- ifelse(altered_vs_wt$sum>0, "mut", "wt")
  altered_vs_wt$patientId <- matrix(unlist(stringr::str_split(altered_vs_wt$studyID.sampleId, ":")), byrow=T, ncol=2)[,2]
  CTNNB1mut <- altered_vs_wt[altered_vs_wt$CTNNB1>0,]
  wt <- altered_vs_wt[altered_vs_wt$BcatMut == "wt",]
  CTNNB1mut_allWt <- rbind(CTNNB1mut,wt)
  altered_vs_wt <- CTNNB1mut_allWt
  colnames(rna_tumor) <- apply(matrix(unlist(stringr::str_split(colnames(rna_tumor), "-")), ncol=7, byrow=T)[,1:4], 1, paste, collapse="-" )
  colnames(rna_tumor) <- substr(colnames(rna_tumor),1,nchar(colnames(rna_tumor))-1)
  
  mut_cna_rna_samples <- intersect(colnames(rna_tumor), altered_vs_wt$patientId)
  
  rna_tumor <- rna_tumor[,which(colnames(rna_tumor) %in% mut_cna_rna_samples)]
  altered_vs_wt <- altered_vs_wt[which(altered_vs_wt$patientId %in% mut_cna_rna_samples),]
  
  sample_info <- as.data.frame(altered_vs_wt[,c("patientId","BcatMut")])
  
  ### Check that counts and sample info have same sample order ----------------
  matrix <- rna_tumor[,match(as.character(sample_info$patientId),colnames(rna_tumor))]
  identical(sample_info$patientId, colnames(matrix))
  
  ### DE analysis -----------------------------------------------------
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = matrix,
    colData = sample_info,
    design = ~ BcatMut )
  dds$BcatMut <- relevel(dds$BcatMut, ref = "wt")
  dds <- estimateSizeFactors(dds)
  deseqNormCounts <- counts(dds, normalized=TRUE)
  log2counts <- deseqNormCounts + 1
  log2counts <- as.data.frame(apply(log2counts, 2, log2))
  log2counts <- as.data.frame(cbind(rownames(log2counts), log2counts))
  colnames(log2counts)[1] <- "GeneSymbol"
  if (!(file.exists(paste0(out, "/DEanalysis/")))){
    dir.create(file.path(paste0(out, "/DEanalysis/")))
  }
  saveRDS(file=paste0(out, "DEanalysis/TCGA-",cancer,"_DESeqLog2counts_", analysis,".RDS"), log2counts)
  dds <- DESeq(dds) 
  res <- results(dds)
  res <- res[!(is.na(res$padj)),]
  saveRDS(file=paste0(out, "DEanalysis/TCGA-",cancer,"_DE_",analysis,".RDS"), res)
}

knownTargets <- c("ABCB1", "AFF3", "AXIN2", "BCL2L2", "BIRC5", "CCND1", "CDC25A", "CDKN2A", "CDX1", "CLDN1", "CTLA4", "DKK1", "EDN1", "ENAH", "ENC1", "FGF18", 
                  "FGF4", "FGFBP1", "FOSL1", "FSCN1", "FST", "FZD7", "GBX2", "HES1", "HNF1A", "ID2", "JAG1", "JUN", "KRT5", "L1CAM","LAMC2", "LEF1", "LGR5", "MMP14", 
                  "MMP7", "MYC", "MYCBP", "NEDD9", "NEUROD1", "NEUROG1", "NOS2", "NOTCH2", "NRCAM", "PDE2A", "PLAU", 
                  "PPARD", "PTGS2", "S100A4", "SGK1", "SMC3", "SP5", "SUZ12", "TBX3", "TBXT", "TCF4", "TERT", "TNC", "TNFRSF19", "VCAN", "VEGFA", 
                  "YY1AP1") # known B-catenin targets (those reported on Nusse lab website + Herbst et al 2014 paper + PDE2A + AFF3)
all <- matrix(ncol=5)
colnames(all) <- c("GeneSymbol","log2FoldChange","padj", "significant","cancer")
analysis <- "onlyCTNNB1MutWt"
for(cancer in c("LIHC","UCEC","PAAD","SKCM","STAD","COAD")){
  
  res <- readRDS(file=paste0(out, "/DEanalysis/TCGA-",cancer,"_DE_",analysis,".RDS"))
  select <- res[which(rownames(res) %in% knownTargets),]
  select <- as.data.frame(cbind(rownames(select), select$log2FoldChange, select$padj))
  colnames(select) <- c("GeneSymbol","log2FoldChange","padj")
  select$log2FoldChange <- as.numeric(as.vector(select$log2FoldChange))
  select$padj <- as.numeric(as.vector(select$padj))
  select$significant <- ifelse(select$padj < 0.1, "*", "")
  select$significant <- factor(select$significant, levels=c("*","") )
  select$cancer <- cancer
  all <- rbind(all, select)
  
}
all <- all[-1,]
# Add DE analysis for ACC, already performed and available in input. NOTE THAT FOR ACC TUMORS CTNNB1 MUTATIONS
# WERE MANUALLY INSPECTED: two patients present a CTNNB1 truncating mutation and were thus not classified
# among patients with CTNNB1 activating mutations.
ACC_DE <- read.table(file="./ACC_datasets_formatted/tcga/CTNNB1mutOnly_vs_wt_allResults.txt", header=T, sep="\t") 
ACC_DE$GeneSymbol <- rownames(ACC_DE)
ACC_DE <- ACC_DE[ACC_DE$GeneSymbol %in% knownTargets,]
ACC_DE$significant <- ifelse(ACC_DE$padj<0.1, "*", "")
ACC_DE$cancer <- "ACC"
ACC_DE <- ACC_DE[,c("GeneSymbol","log2FoldChange","padj","significant","cancer")]
all <- rbind(all, ACC_DE)
png(paste0(out, "Heatmap_DE_FC_padj_allTargets_",analysis,".png"), height=900)
ggplot(data=all, aes(y=GeneSymbol, x=cancer, fill=log2FoldChange)) + geom_tile() +
  theme_bw() + ggtitle(analysis) + scale_fill_gradient2(low = "green",
                                                        mid = "white",
                                                        high = "red") +
  geom_text(aes(label=significant), color="black")
dev.off()


# READ SIGNATURES ---------------------------------------------
pancancer3genes <- c("AXIN2","SP5","TNFRSF19")
LEF1 <- c("LEF1")
ACC_BcatSignature <- c("AFF3","AXIN2","FSCN1","LEF1","MYC","PDE2A","SP5","TBX3","TNFRSF19")
list_signatures <- list(ACC_BcatSignature=ACC_BcatSignature,
                        pancancer3genes=pancancer3genes,
                        LEF1=LEF1)

# ACC CANCER  ----------------------------------------------------------
cancer <- "ACC"

# Read sample classification in CTNNB1 altered, ZNRF3 altered, wild type for both CTNNB1 and ZNRF3 (Assie and tcga) or nuclear vs membrane B-Catenin accumulation (Heaton)
Assie_Bcat <- readRDS(file="./ACC_datasets_formatted/Assie/Metadata/Assie_sampleIDs_CTNNB1mut_ZNRF3mut_wt.RDS")
Heaton_Bcat <- readRDS(file="./ACC_datasets_formatted/Heaton/Metadata/Heaton_sampleIDs_BcatNuclearMembrane_onlyACCs.RDS")
tcga_Bcat <- readRDS(file="./ACC_datasets_formatted/tcga/Metadata/tcga_sampleIDs_CTNNB1mut_ZNRF3mut_wt.RDS")
colnames(tcga_Bcat)[1] <- "TumorID"

# Read metadata
Assie_metadata <- read.csv("./ACC_datasets_formatted/Assie/Metadata/Metadata_Assie.csv")
Assie_metadata <- Assie_metadata[!(is.na(Assie_metadata$BcatStatus)),]
Assie_metadata$BcatStatus <- factor(Assie_metadata$BcatStatus, levels=c("CTNNB1","ZNRF3","wt"))

Heaton_metadata <- read.csv("./ACC_datasets_formatted/Heaton/Metadata/Metadata_Heaton.csv")
Heaton_metadata <- Heaton_metadata[Heaton_metadata$Histotype %in% c("ACC"),]
Heaton_metadata <- Heaton_metadata[!(is.na(Heaton_metadata$BetaCateninStaining)),]
Heaton_metadata$BcatStatus <- factor(Heaton_metadata$BetaCateninStaining, levels=c("Nuclear","Membrane"))

tcga_metadata <- read.csv("./ACC_datasets_formatted/tcga/Metadata/Metadata_tcga.csv")
tcga_metadata <- merge(tcga_metadata, tcga_Bcat, all=F, by="TumorID")
tcga_metadata$BcatStatus <- factor(tcga_metadata$Alteration, levels=c("CTNNB1_mut","ZNRF3_mut","wt"))

# Read Z-score transformed log2 normalized expression (average if multiple probes) and metadata -----------
ACC_Zscored <- readRDS(file="./ACC_datasets_integration/Counts_Zscore_allDatasets.RDS")
genes <- ACC_Zscored$GeneSymbol

Assie <- ACC_Zscored[,grep("Assie",colnames(ACC_Zscored))]
colnames(Assie) <- stringr::str_replace_all(colnames(Assie), paste0("Assie", "_"), "")
Assie$GeneSymbol <- genes
Heaton <- ACC_Zscored[,grep("Heaton",colnames(ACC_Zscored))]
colnames(Heaton) <- stringr::str_replace_all(colnames(Heaton), paste0("Heaton", "_"), "")
Heaton$GeneSymbol <- genes
tcga <- ACC_Zscored[,grep("tcga",colnames(ACC_Zscored))]
colnames(tcga) <- stringr::str_replace_all(colnames(tcga), paste0("tcga", "_"), "")
tcga$GeneSymbol <- genes

list_datasets <- list(tcga=tcga,Heaton=Heaton,Assie=Assie)

# Plots with AUC & co values in title ---------------------------------
heatmapAUC <- as.data.frame(matrix(ncol=5))
colnames(heatmapAUC) <- c("dataset", "cancer", "signature","analysis","AUC")
for(d in 1:length(list_datasets)){
  dataset_name <- names(list_datasets)[d]
  dataset <- list_datasets[[d]]
  print(dataset_name)
  
  if(dataset_name == "Assie"){
    Assie_Bcat <- as.data.frame(Assie_Bcat)
    Assie_Bcat <- Assie_Bcat[Assie_Bcat$Alteration %in% c("CTNNB1", "wt"),]
    Bcat_obj <- Assie_Bcat
    alterations <- c("CTNNB1","wt")}
  if(dataset_name == "Heaton"){
    Bcat_obj <- Heaton_Bcat
    alterations <- c("Nuclear","Membrane")}
  if(dataset_name == "tcga"){
    tcga_Bcat <- as.data.frame(tcga_Bcat)
    tcga_Bcat <- tcga_Bcat[tcga_Bcat$Alteration %in% c("CTNNB1_mut", "wt"),]
    Bcat_obj <- tcga_Bcat
    alterations <- c("CTNNB1_mut","wt")}
  
  for(s in 1:length(list_signatures)){
    signature_name <- names(list_signatures)[s]
    signature_genes <- list_signatures[[s]]
    print(signature_name)
    
    ### AUC, sensitivity, specificity values
    if(length(signature_genes)>1){
      index <- which(colnames(dataset) == "GeneSymbol")
      signature <- dataset[which(dataset$GeneSymbol %in% signature_genes),-index]
      mean <- apply(as.matrix(signature), 2, mean)
      mean <- as.data.frame(cbind(names(mean),mean))
      colnames(mean) <- c("patientId","Signature")
      colnames(Bcat_obj)[1] <- "patientId"
      input <- merge(mean, Bcat_obj, by="patientId")
      input <- input[input$Alteration %in% alterations,]
      input$Signature <- as.numeric(as.vector(input$Signature))
    }else{
      index <- which(colnames(dataset) == "GeneSymbol")
      singleGene <- dataset[which(dataset$GeneSymbol==signature_genes),-index]
      singleGene <- t(rbind(names(singleGene), singleGene))
      colnames(singleGene) <- c("patientId","Signature")
      colnames(Bcat_obj)[1] <- "patientId"
      input <- merge(singleGene, Bcat_obj, by="patientId")
      input <- input[input$Alteration %in% alterations,]
      input$Signature <- as.numeric(as.vector(input$Signature))
    }
    
    if(metric == "sum_sens_spec"){
      cp <- cutpointr(input, Signature, Alteration,
                      method = maximize_metric, metric = sum_sens_spec)
    }
    if(metric == "sens_constrain"){
      cp <- cutpointr(input, Signature, Alteration,
                      method = maximize_metric, metric = sens_constrain)
    }
    saveRDS(file=paste0(out,"Cutpointr_", metric, "_", dataset_name, "_", cancer, "_", signature_name, "_", analysis, ".RDS"), summary(cp))
    
    ### Plot
    plot_out <- plot_Density_byBcatStatus(dataset=dataset, dataset_name=dataset_name, signature=signature_genes, signatureName=signature_name, Bcat=input[,c('patientId','Alteration')])
    png(paste0(out,"Density_perBcatStatus_",signature_name, "_", analysis, "_", cancer, "_", dataset_name, "_AvgProbes_",metric,".png"))
    print(plot_out +
            ggtitle(paste0(cancer, ", ", dataset_name,  ", AUC=",round(cp$AUC, digits=2), "\nsensitivity=", round(cp$sensitivity,digits=2), ",specificity=", round(cp$specificity,digits=2))))
    dev.off()
    
    vec <- c(dataset_name, cancer, signature_name, analysis, round(cp$AUC, digits=3)) 
    names(vec) <- c("dataset", "cancer", "signature","analysis","AUC")
    heatmapAUC <- rbind(heatmapAUC, vec)
  }
}
heatmapAUC <- heatmapAUC[-1,]
heatmapAUC$AUC <- as.numeric(as.vector(heatmapAUC$AUC))
heatmapAUC <- heatmapAUC[heatmapAUC$signature %in% c("LEF1", "ACC_BcatSignature", "pancancer3genes"),]
heatmapAUC$signature <- factor(heatmapAUC$signature, levels=c("LEF1", "ACC_BcatSignature", "pancancer3genes"))
saveRDS(file=paste0(out, "/AUCs_per_signature_AssieHeatonTCGA_", cancer, "_", analysis, "_", metric, ".RDS"), heatmapAUC)


# OTHER CANCER TYPES, USING TCGA DATA ----------------------------------------------------------
heatmapAUC <- as.data.frame(matrix(ncol=5))
colnames(heatmapAUC) <- c("dataset", "cancer", "signature","analysis","AUC")
dataset_name <- "TCGA"
for(cancer in c("LIHC","UCEC","SKCM","STAD","PAAD","COAD")){
  
  print(cancer)
  
  for(s in 1:length(list_signatures)){
    signature_name <- names(list_signatures)[s]
    signature_genes <- list_signatures[[s]]
    print(signature_name)
    
    dataset <- readRDS(file=paste0(out, "DEanalysis/TCGA-",cancer,"_DESeqLog2counts_", analysis,".RDS"))
    index <- which(colnames(dataset) == "GeneSymbol")
    expr <- as.matrix(dataset[,-index])
    # Apply Z-score transformation
    # Need to remove genes with st dev = 0 (for which Z-scores is NAs) 
    check_sd <- apply(expr, 1, sd)
    toRemove <- which(check_sd == 0)
    if(length(toRemove)>0){ # it will be always 0 if only most variable genes selected. Otherwise it depends on the dataset
      expr <- expr[-(toRemove),]
      geneNames <- rownames(expr)
    }else{
      geneNames <- rownames(expr)
    }
    # Z-score normalization on a per gene basis
    expr_Z <- as.data.frame(t(scale(t(expr))))
    expr_Z$GeneSymbol <- geneNames
    dataset <- expr_Z
    
    altered_vs_wt <- read.table(file=paste0(inpath, "/" ,cancer, "_sample_matrix.txt"), header=T, sep="\t")
    altered_vs_wt$sum <- rowSums(altered_vs_wt[,-1])
    altered_vs_wt$BcatMut <- ifelse(altered_vs_wt$sum>0, "mut", "wt")
    altered_vs_wt$patientId <- matrix(unlist(stringr::str_split(altered_vs_wt$studyID.sampleId, ":")), byrow=T, ncol=2)[,2]
    if(analysis == "onlyCTNNB1MutWt"){
      CTNNB1mut <- altered_vs_wt[altered_vs_wt$CTNNB1>0,]
      wt <- altered_vs_wt[altered_vs_wt$BcatMut == "wt",]
      CTNNB1mut_allWt <- rbind(CTNNB1mut,wt)
      dataset <- dataset[,which(colnames(dataset) %in% CTNNB1mut_allWt$patientId)]
      altered_vs_wt <- CTNNB1mut_allWt
      dataset <- cbind(rownames(dataset), dataset)
      colnames(dataset)[1] <- "GeneSymbol"
    }
    alterations <- c("mut","wt")
    
    if(length(signature_genes)>1){
      index <- which(colnames(dataset) == "GeneSymbol")
      signature <- dataset[which(dataset$GeneSymbol %in% signature_genes),-index]
      mean <- apply(as.matrix(signature), 2, mean)
      mean <- as.data.frame(cbind(names(mean),mean))
      colnames(mean) <- c("patientId","Signature")
      input <- merge(mean, altered_vs_wt, by="patientId")
      input$Signature <- as.numeric(as.vector(input$Signature))
      input$Alteration <- factor(input$BcatMut)
    }else{
      index <- which(colnames(dataset) == "GeneSymbol")
      singleGene <- dataset[which(dataset$GeneSymbol==signature_genes),-index]
      singleGene <- t(rbind(names(singleGene), singleGene))
      colnames(singleGene) <- c("patientId","Signature")
      input <- merge(singleGene, altered_vs_wt, by="patientId")
      input$Alteration <- input$BcatMut
      input <- input[input$Alteration %in% alterations,]
      input$Signature <- as.numeric(as.vector(input$Signature))
    }
    
    ### Estimate sensitivity, specificity, AUC
    if(metric == "sum_sens_spec"){
      cp <- cutpointr(input, Signature, Alteration,
                      method = maximize_metric, metric = sum_sens_spec)
    }
    if(metric == "sens_constrain"){
      cp <- cutpointr(input, Signature, Alteration,
                      method = maximize_metric, metric = sens_constrain)
    }
    saveRDS(file=paste0(out,"Cutpointr_", metric, "_", dataset_name, "_", cancer, "_", signature_name, "_", analysis, ".RDS"), summary(cp))
    
    ### Plot
    plot_out <- plot_Density_byBcatStatus(dataset=dataset, dataset_name=dataset_name, signature=signature_genes, signatureName=signature_name, Bcat=input[,c('patientId','Alteration')])
    png(paste0(out,"Density_perBcatStatus_",signature_name, "_", analysis, "_", cancer, "_", dataset_name, "_AvgProbes_",metric,".png"))
    print(plot_out +
            ggtitle(paste0(cancer, ", ", dataset_name,  ", AUC=",round(cp$AUC, digits=2), "\nsensitivity=", round(cp$sensitivity,digits=2), ",specificity=", round(cp$specificity,digits=2))))
    dev.off()
    
    vec <- c(dataset_name, cancer, signature_name, analysis, round(cp$AUC, digits=3)) 
    names(vec) <- c("dataset", "cancer", "signature","analysis","AUC")
    heatmapAUC <- rbind(heatmapAUC, vec)
  }
}
heatmapAUC <- heatmapAUC[-1,]
heatmapAUC$AUC <- as.numeric(as.vector(heatmapAUC$AUC))
heatmapAUC <- heatmapAUC[heatmapAUC$signature %in% c("LEF1", "ACC_BcatSignature", "pancancer3genes"),]
heatmapAUC$signature <- factor(heatmapAUC$signature, levels=c("LEF1", "ACC_BcatSignature", "pancancer3genes"))
saveRDS(file=paste0(out, "/AUCs_per_signature_", dataset_name, "_OtherCancers_", analysis, "_", metric, ".RDS"), heatmapAUC)


# UCEC SUPP DATA (NON-TCGA)  -----------------------------------------
# Prepare input -----------------------------------------------------------
ids <- read.table(paste0(inpath, "/CPTAC3_UCEC/mRNAFileID_cases.submitter_id_fromTCGABiolinks.txt"), header=T, sep="\t") ### read correspondance fileID - sampleID (to merge rna and mut data)
ids$file_nameBIS <- stringr::str_replace(ids$file_name, ".rna_seq.augmented_star_gene_counts.tsv", "")

rna <- read.table(file=paste0(inpath, "/CPTAC3_UCEC/all_genes_mRNAexpr_unstranded_matrix.txt"), header=T, sep="\t") ### read input mRNA
colnames(rna)[-1] <- stringr::str_replace_all(stringr::str_replace_all(colnames(rna)[-1], '\\.', "-"),'^X',"")
genes <- rna$GeneSymbol
rna <- as.matrix(rna[,-which(colnames(rna) == "GeneSymbol")])
rna <- apply(rna, 2, as.numeric)
rownames(rna) <- genes

# Perform DESeq normalization
sample_info <- as.data.frame(colnames(rna))
colnames(sample_info) <- "patientId" 
sample_info$patientId <- factor(sample_info$patientId)
dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = rna,
  colData = sample_info,
  design = ~ patientId)
dds <- estimateSizeFactors(dds)
deseqNorm_UCEC <- counts(dds, normalized=TRUE)

# Apply Z-score transformation
# Need to remove genes with st dev = 0 (for which Z-scores is NAs) 
check_sd <- apply(deseqNorm_UCEC, 1, sd)
toRemove <- which(check_sd == 0)
if(length(toRemove)>0){ # it will be always 0 if only most variable genes selected. Otherwise it depends on the dataset
  deseqNorm_UCEC <- deseqNorm_UCEC[-(toRemove),]
  geneNames <- rownames(deseqNorm_UCEC)
}else{
  geneNames <- rownames(deseqNorm_UCEC)
}
# Z-score normalization on a per gene basis
expr_Z <- as.data.frame(t(scale(t(deseqNorm_UCEC))))
rownames(expr_Z) <- geneNames

# Metadata (CTNNB1 alterations)
mut <- read.table(file=paste0(inpath, "/CPTAC3_UCEC/UCEC_CTNNB1mut_sample_matrix_from_cBioportal.txt"), header=T, sep="\t") ### read input mutations and CNAs manually downloaded from cBioportal
mut$studyID.sampleId <- stringr::str_replace_all(mut$studyID.sampleId, "uec_cptac_gdc:", "")
colnames(mut)[1] <- "cases.submitter_id"
m <- matrix(unlist(stringr::str_split(mut$cases.submitter_id, "-")), ncol=3, byrow = T)[,1:2]
mut$cases.submitter_id <- apply(m, 1, function(x) paste(x, collapse="-"))

mut_id <- merge(mut, ids, by="cases.submitter_id", all=F)
mut_id <- mut_id[,c("CTNNB1", "cases.submitter_id", "file_nameBIS")]

# Plots with AUC & co values in title ---------------------------------
dataset_name <- "CPTAC3GDC"
cancer <- "UCEC"
metadata <- mut_id # read CTNNB1 status
metadata$CTNNB1.mutation <- metadata$CTNNB1
metadata$Sample <- metadata$file_nameBIS

heatmapAUC <- as.data.frame(matrix(ncol=5))
colnames(heatmapAUC) <- c("dataset","cancer","signature","analysis","AUC")
for(s in 1:length(list_signatures)){
  
  print(s)
  signature_name <- names(list_signatures)[s]
  signature_genes <- list_signatures[[s]]
  print(signature_name)
  
  CTNNB1mut <- metadata[!(metadata$CTNNB1.mutation == 0),]
  CTNNB1mut$BcatStatus <- "mut"
  wt <- metadata[metadata$CTNNB1.mutation == 0,]
  wt$BcatStatus <-"wt"
  CTNNB1mut_allWt <- rbind(CTNNB1mut,wt)
  dataset <- expr_Z
  genes <- rownames(expr_Z)
  dataset <- dataset[,which(colnames(dataset) %in% CTNNB1mut_allWt$Sample)]
  rownames(dataset) <- genes
  altered_vs_wt <- CTNNB1mut_allWt
  alterations <- c("mut","wt")
  
  if(length(signature_genes)>1){
    signature <- as.matrix(dataset[which(rownames(dataset) %in% signature_genes),])
    signature <- apply(signature, 2, as.numeric)
    mean <- apply(signature, 2, mean)
    mean <- as.data.frame(cbind(names(mean),mean))
    colnames(mean) <- c("Sample","Signature")
    input <- merge(mean, altered_vs_wt, by="Sample")
    input$Signature <- as.numeric(as.vector(input$Signature))
    input$BcatStatus <- factor(input$BcatStatus, levels=c("mut","wt"))
    input$Alteration <- input$BcatStatus
  }else{
    singleGene <- as.matrix(dataset[which(rownames(dataset) %in% signature_genes),])
    singleGene <- as.data.frame(cbind(rownames(singleGene), singleGene))
    colnames(singleGene) <- c("Sample","Signature")
    input <- merge(singleGene, altered_vs_wt, by="Sample")
    input$Alteration <- input$BcatStatus
    input <- input[input$Alteration %in% alterations,]
    input$Signature <- as.numeric(as.vector(input$Signature))
  }
  
  ### Estimate sensitivity, specificity, AUC
  if(metric == "sum_sens_spec"){
    cp <- cutpointr(input, Signature, Alteration,
                    method = maximize_metric, metric = sum_sens_spec)
  }
  if(metric == "sens_constrain"){
    cp <- cutpointr(input, Signature, Alteration,
                    method = maximize_metric, metric = sens_constrain)
  }
  saveRDS(file=paste0(out,"Cutpointr_", metric, "_", dataset_name, "_", cancer, "_", signature_name, "_", analysis, ".RDS"), summary(cp))
  
  ### Plot
  dataset <- expr_Z
  GeneSymbol <- rownames(dataset)
  dataset <- apply(dataset, 2, as.numeric)
  dataset <- apply(dataset, 2, function(x) log2(x+1))
  dataset <- as.data.frame(cbind(GeneSymbol, dataset))
  colnames(dataset)[1] <- "GeneSymbol"
  plot_out <- plot_Density_byBcatStatus(dataset=dataset, dataset_name=dataset_name, signature=signature_genes, signatureName=signature_name, Bcat=input[,c('Sample','BcatStatus')])
  png(paste0(out,"Density_perBcatStatus_",signature_name, "_", analysis, "_", cancer, "_", dataset_name, "_", metric,".png"))
  print(plot_out + 
          ggtitle(paste0(dataset_name, ", AUC=",round(cp$AUC, digits=2), "\nsensitivity=", round(cp$sensitivity,digits=2), ",specificity=", round(cp$specificity,digits=2))))
  dev.off()
  
  vec <- c(dataset_name, cancer, signature_name, analysis, round(cp$AUC, digits=3)) 
  names(vec) <- c("dataset","cancer","signature","analysis","AUC")
  heatmapAUC <- rbind(heatmapAUC, vec)
}
heatmapAUC <- heatmapAUC[-1,]
heatmapAUC$AUC <- as.numeric(as.vector(heatmapAUC$AUC))
heatmapAUC <- heatmapAUC[heatmapAUC$signature %in% c("LEF1", "ACC_BcatSignature", "pancancer3genes"),]
heatmapAUC$signature <- factor(heatmapAUC$signature, levels=c("LEF1", "ACC_BcatSignature", "pancancer3genes"))
saveRDS(file=paste0(out, "/AUCs_per_signature_", dataset_name, "_", cancer, "_", analysis, "_", metric, ".RDS"), heatmapAUC)


# LIHC SUPP DATA (NON-TCGA)  -----------------------------------------
cancer <- "LIHC"
dataset_name <- "Chiang"
### Read inputs
metadata <- read.csv(file=paste0(inpath, "LIHC_Chiang/SupplementaryTable5.csv")) # read Bcat status
dataset <- readRDS(file=paste0(inpath, "LIHC_Chiang/Chiang_norm_counts_AvgMultipleProbes.RDS")) # read counts
genes <- dataset$GeneSymbol
dataset <- as.matrix(dataset[,-which(colnames(dataset) == "GeneSymbol")])
dataset <- apply(dataset, 2, as.numeric)
rownames(dataset) <- genes

# Apply Z-score transformation
# Need to remove genes with st dev = 0 (for which Z-scores is NAs) 
check_sd <- apply(dataset, 1, sd)
toRemove <- which(check_sd == 0)
if(length(toRemove)>0){ # it will be always 0 if only most variable genes selected. Otherwise it depends on the dataset
  dataset <- dataset[-(toRemove),]
  geneNames <- rownames(dataset)
}else{
  geneNames <- rownames(dataset)
}
# Z-score normalization on a per gene basis
expr_Z <- as.data.frame(t(scale(t(dataset))))
rownames(expr_Z) <- geneNames
saveRDS(file=paste0(inpath, "/LIHC_Chiang/Chiang_norm_counts_AvgMultipleProbes_Zscore.RDS"), expr_Z)

heatmapAUC <- as.data.frame(matrix(ncol=5))
colnames(heatmapAUC) <- c("dataset","cancer","signature","analysis","AUC")
for(index in 1:length(list_signatures)){
  
  print(index)
  signature_name <- names(list_signatures)[index]
  signature_genes <- list_signatures[[index]]
  print(signature_name)
  
  metadata <- metadata[!(metadata$CTNNB1.mutation == ""),] # remove NA for CTNNB1.mutation
  CTNNB1mut <- metadata[(!(metadata$CTNNB1.mutation == 0) & (!(metadata$CTNNB1.mutation == ""))),]
  CTNNB1mut$BcatStatus <- "mut"
  wt <- metadata[metadata$CTNNB1.mutation == 0,]
  wt$BcatStatus <-"wt"
  CTNNB1mut_allWt <- rbind(CTNNB1mut,wt)
  genes <- rownames(expr_Z)
  dataset <- expr_Z[,which(colnames(expr_Z) %in% CTNNB1mut_allWt$Sample)]
  rownames(dataset) <- genes
  altered_vs_wt <- CTNNB1mut_allWt
  alterations <- c("mut","wt")
  
  if(length(signature_genes)>1){
    signature <- as.matrix(dataset[which(rownames(dataset) %in% signature_genes),])
    signature <- apply(signature, 2, as.numeric)
    mean <- apply(signature, 2, mean)
    mean <- as.data.frame(cbind(names(mean),mean))
    colnames(mean) <- c("Sample","Signature")
    input <- merge(mean, altered_vs_wt, by="Sample")
    input$Signature <- as.numeric(as.vector(input$Signature))
    input$BcatStatus <- factor(input$BcatStatus, levels=c("mut","wt"))
    input$Alteration <- input$BcatStatus
    input <- input[input$Alteration %in% alterations,]
  }else{
    singleGene <- as.matrix(dataset[which(rownames(dataset) %in% signature_genes),])
    singleGene <- t(rbind(colnames(singleGene), singleGene))
    colnames(singleGene) <- c("Sample","Signature")
    input <- merge(singleGene, altered_vs_wt, by="Sample")
    input$Alteration <- input$BcatStatus
    input <- input[input$Alteration %in% alterations,]
    input$Signature <- as.numeric(as.vector(input$Signature))
  }
  
  ### Estimate sensitivity, specificity, AUC
  if(metric == "sum_sens_spec"){
    cp <- cutpointr(input, Signature, Alteration,
                    method = maximize_metric, metric = sum_sens_spec)
  }
  if(metric == "sens_constrain"){
    cp <- cutpointr(input, Signature, Alteration,
                    method = maximize_metric, metric = sens_constrain)
  }
  saveRDS(file=paste0(out,"Cutpointr_", metric, "_", dataset_name, "_", cancer, "_", signature_name, "_", analysis, ".RDS"), summary(cp))
  
  ### Plot
  dataset <- readRDS(file=paste0(inpath, "/LIHC_Chiang/Chiang_norm_counts_AvgMultipleProbes_Zscore.RDS"))
  GeneSymbol <- rownames(dataset)
  dataset <- apply(dataset, 2, as.numeric)
  dataset <- as.data.frame(cbind(GeneSymbol, dataset))
  colnames(dataset)[1] <- "GeneSymbol"
  plot_out <- plot_Density_byBcatStatus(dataset=dataset, dataset_name=dataset_name, signature=signature_genes, signatureName=signature_name, Bcat=input[,c('Sample','BcatStatus')])
  png(paste0(out,"Density_perBcatStatus_",signature_name, "_", analysis, "_", cancer, "_", dataset_name, "_", metric,".png"))
  print(plot_out +
          ggtitle(paste0(dataset_name, ", AUC=",round(cp$AUC, digits=2), "\nsensitivity=", round(cp$sensitivity,digits=2), ",specificity=", round(cp$specificity,digits=2))))
  dev.off()
  
  vec <- c(dataset_name, cancer, signature_name, analysis, round(cp$AUC, digits=3)) 
  names(vec) <- c("dataset", "cancer","signature","analysis","AUC")
  heatmapAUC <- rbind(heatmapAUC, vec)
}
heatmapAUC <- heatmapAUC[-1,]
heatmapAUC$AUC <- as.numeric(as.vector(heatmapAUC$AUC))
heatmapAUC <- heatmapAUC[heatmapAUC$signature %in% c("LEF1", "ACC_BcatSignature", "pancancer3genes"),]
heatmapAUC$signature <- factor(heatmapAUC$signature, levels=c("LEF1", "ACC_BcatSignature", "pancancer3genes"))
saveRDS(file=paste0(out, "/AUCs_per_signature_", dataset_name, "_", cancer, "_", analysis, "_", metric, ".RDS"), heatmapAUC)


# PLOT HEAMAPS AUCs ------------------------------------------------------
ACCs <- readRDS(file=paste0(out, "/AUCs_per_signature_AssieHeatonTCGA_ACC_", analysis, "_", metric, ".RDS"))
OtherCancers <- readRDS(file=paste0(out, "/AUCs_per_signature_TCGA_OtherCancers_", analysis, "_", metric, ".RDS"))
LIHC_supp <- readRDS(file=paste0(out, "/AUCs_per_signature_Chiang_LIHC_", analysis, "_", metric, ".RDS"))
UCEC_supp <- readRDS(file=paste0(out, "/AUCs_per_signature_CPTAC3GDC_UCEC_", analysis, "_", metric, ".RDS"))

all <- rbind(ACCs, OtherCancers, LIHC_supp, UCEC_supp)
all$dataset <- stringr::str_replace_all(all$dataset, "tcga", "TCGA")
all$ID <- paste(all$dataset, all$cancer, sep="-")

heatmapAUC <- all[all$ID %in% c("TCGA-ACC","TCGA-LIHC","TCGA-UCEC","TCGA-SKCM","TCGA-STAD","TCGA-PAAD","TCGA-COAD", "Assie-ACC", "Heaton-ACC", "Chiang-LIHC","CPTAC3GDC-UCEC"),]
heatmapAUC <- heatmapAUC[heatmapAUC$signature %in% c("pancancer3genes","ACC_BcatSignature"),]
heatmapAUC$AUC <- round(heatmapAUC$AUC, digits=2)

png(file=paste0(out,"/HeatmapAUCs_onlyCTNNB1MutWt_MultiCancer3gene_ACCsignature_",metric,"_TrainAndValidation.png"), width=500, height = 200)
ggplot(data=heatmapAUC, aes(x=factor(ID, levels=c("TCGA-ACC","TCGA-LIHC","TCGA-UCEC","TCGA-SKCM","TCGA-STAD","TCGA-PAAD","TCGA-COAD", "Assie-ACC", "Heaton-ACC", "Chiang-LIHC","CPTAC3GDC-UCEC")),
                            y=factor(signature, levels=c("pancancer3genes","ACC_BcatSignature")), fill=AUC)) + geom_tile() +
  geom_text(aes(label=AUC)) + theme_bw() + xlab("") + ylab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_gradient2(low="navy", mid="white", high="red", midpoint = 0.6)
dev.off()


