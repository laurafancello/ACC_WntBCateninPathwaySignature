library(ggplot2)
library(ggpubr)
library(ggExtra)
library(gridExtra)
library(cutpointr)
library(ggrepel)
library(survival)
library(survminer)
library(ggsurvfit)
library(patchwork)

out <- "./SignatureAndLEF1_inACC/output/"

source("./Scripts/Miscellaneous_Functions.R")

#  Read signatures ---------------------------------------------
signature <- c("AFF3","AXIN2","FSCN1","LEF1","MYC","PDE2A","SP5","TBX3","TNFRSF19") # home-made signature of Wnt/B-Catenin pathway activation in ACC
Nusse_PDE2A <- c("ABCB1", "AFF3", "AXIN2", "BCL2L2", "BIRC5", "CCND1", "CDC25A", "CDKN2A", "CDX1", "CLDN1", "CTLA4", "DKK1", "EDN1", "ENAH", "ENC1", "FGF18", 
                 "FGF4", "FGFBP1", "FOSL1", "FSCN1", "FST", "FZD7", "GBX2", "HES1", "HNF1A", "ID2", "JAG1", "JUN", "KRT5", "L1CAM","LAMC2", "LEF1", "LGR5", "MMP14", 
                 "MMP7", "MYC", "MYCBP", "NEDD9", "NEUROD1", "NEUROG1", "NOS2", "NOTCH2", "NRCAM", "PDE2A", "PLAU", 
                 "PPARD", "PTGS2", "S100A4", "SGK1", "SMC3", "SP5", "SUZ12", "TBX3", "TBXT", "TCF4", "TERT", "TNC", "TNFRSF19", "VCAN", "VEGFA", 
                 "YY1AP1") # known B-catenin targets (those reported on Nusse lab website + Herbst et al 2014 paper + PDE2A + AFF3)
KeggWntPathwayGenes <- read.table(file="./SignatureAndLEF1_inACC/input/MSigDB_C2CPCurated_KEGG_WNT_SIGNALING_PATHWAY.v2024.1.Hs.txt", header=T, sep="\t") # KEGG genes for Wnt signaling pathway
KeggWntPathwayGenes <- KeggWntPathwayGenes$Gene_symbol
pancancer3gene <- c("AXIN2","SP5","TNFRSF19") # home-made 3-gene signature of Wnt/B-Catenin pathway activation for multiple cancer types
                 
genesOfInterest <- c(signature, "CTNNB1", "ZNRF3", "APC", "LGR5", "MYC", "TCF7", "TCF7L1", "TCF7L2")
list_signatures <- list(CTNNB1=c("CTNNB1"), ZNRF3=c("ZNRF3"),APC=c("APC"),
                        AXIN2=c("AXIN2"), LGR5=c("LGR5"), SP5=c("SP5"),# other canonical targets of Wnt/B-catenin pathway
                        MYC=c("MYC"),AFF3=c("AFF3"),PDE2A=c("PDE2A"), # other non canonical targets of Wnt/B-catenin pathway
                        LEF1=c("LEF1"), TCF7=c("TCF7"), TCF7L1=c("TCF7L1"),TCF7L2=c("TCF7L2"), # LEF1 and other members of TCF/LEF complex
                        signature=signature, 
                        Nusse_PDE2A=Nusse_PDE2A,
                        KeggWntPathwayGenes=KeggWntPathwayGenes, 
                        pancancer3gene=c("AXIN2","SP5","TNFRSF19") 
) 

# Read inputs ----------------------------

### Read Z-score transformed log2 norm counts (average if multiple probes)
all_data <- readRDS("./ACC_datasets_integration/Counts_Zscore_allDatasets.RDS")

### Read metadata
Assie_metadata <- read.csv("./ACC_datasets_formatted/Assie/Metadata/Metadata_Assie.csv")
Assie_metadata <- Assie_metadata[!(is.na(Assie_metadata$BcatStatus)),]
Assie_metadata$BcatStatus <- factor(Assie_metadata$BcatStatus, levels=c("CTNNB1","ZNRF3","wt"))

Heaton_metadata <- read.csv("./ACC_datasets_formatted/Heaton/Metadata/Metadata_Heaton.csv")
Heaton_metadata <- Heaton_metadata[Heaton_metadata$Histotype %in% c("ACC"),]
Heaton_metadata <- Heaton_metadata[!(is.na(Heaton_metadata$BetaCateninStaining)),]
Heaton_metadata$BcatStatus <- factor(Heaton_metadata$BetaCateninStaining, levels=c("Nuclear","Membrane"))

tcga_Bcat <- readRDS(file="./ACC_datasets_formatted/tcga/Metadata/tcga_sampleIDs_CTNNB1mut_ZNRF3mut_wt.RDS") # manual annotation of TCGA samples with CTNNB1 activating mutation, ZNRF3 inactivating alteration or wt fro CTNNB1 and ZNRF3
colnames(tcga_Bcat)[1] <- "TumorID"
tcga_metadata <- read.csv("./ACC_datasets_formatted/tcga/Metadata/Metadata_tcga.csv")
tcga_metadata <- merge(tcga_metadata, tcga_Bcat, all=F, by="TumorID")
tcga_metadata$BcatStatus <- factor(tcga_metadata$Alteration, levels=c("CTNNB1_mut","ZNRF3_mut","wt"))

Assie_Bcat <- readRDS(file="./ACC_datasets_formatted/Assie/Metadata/Assie_sampleIDs_CTNNB1mut_ZNRF3mut_wt.RDS")
Heaton_Bcat <- readRDS(file="./ACC_datasets_formatted/Heaton/Metadata/Heaton_sampleIDs_BcatNuclearMembrane_onlyACCs.RDS")
tcga_Bcat <- readRDS(file="./ACC_datasets_formatted/TCGA/Metadata/tcga_sampleIDs_CTNNB1mut_ZNRF3mut_wt.RDS")


metric <- "sens_constrain"
#metric <- "sum_sens_spec"

# Single PNG plots per gene/signature per dataset ---------------------------------
for(dataset_name in c("Assie","Heaton","tcga")){
  print(dataset_name)
  
  genes <- all_data$GeneSymbol
  dataset <- all_data[,grep(dataset_name,colnames(all_data))]
  colnames(dataset) <- stringr::str_replace_all(colnames(dataset), paste0(dataset_name, "_"), "")
  dataset$GeneSymbol <- genes
  
  if(dataset_name == "Assie"){
    Bcat_obj <- Assie_Bcat
    alterations <- c("CTNNB1","wt")}
  if(dataset_name == "Heaton"){
    Bcat_obj <- Heaton_Bcat
    alterations <- c("Nuclear","Membrane")}
  if(dataset_name == "tcga"){
    Bcat_obj <- tcga_Bcat
    alterations <- c("CTNNB1_mut","wt")}
  if(dataset_name == "Lefevre"){
    Bcat_obj <- Lefevre_Bcat
    alterations <- c("expressed","repressed")}
  
  for(s in 1:length(list_signatures)){
    signatureName <- names(list_signatures)[s]
    signature_genes <- list_signatures[[s]]
    print(signatureName)
    
    ### AUC, sensitivity, specificity values
    if(length(signature_genes)>1){
      index <- which(colnames(dataset) == "GeneSymbol")
      signature <- dataset[which(dataset$GeneSymbol %in% signature_genes),-index]
      mean <- apply(as.matrix(signature), 2, mean)
      mean <- as.data.frame(cbind(names(mean),mean))
      colnames(mean) <- c("Sample","Signature")
      input <- merge(mean, Bcat_obj, by="Sample")
      input <- input[input$Alteration %in% alterations,]
      input$Signature <- as.numeric(as.vector(input$Signature))
    }else{
      index <- which(colnames(dataset) == "GeneSymbol")
      singleGene <- dataset[which(dataset$GeneSymbol==signature_genes),-index]
      singleGene <- t(rbind(names(singleGene), singleGene))
      colnames(singleGene) <- c("Sample","Signature")
      input <- merge(singleGene, Bcat_obj, by="Sample")
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
    saveRDS(file=paste0(out,"/Cutpointr_", metric, "_", dataset_name, "_", signatureName, ".RDS"), cp)
    
    ### Plot
    plot_out <- plot_Density_byBcatStatus(dataset=dataset, dataset_name=dataset_name, signature=signature_genes, signatureName=signatureName, Bcat=Bcat_obj)
    png(paste0(out,"/GeneSignaturesDensityDistribution_perCTNNB1Status_",dataset_name,"_",signatureName,"_AvgProbes_", metric,".png"))
    print(plot_out +
            ggtitle(paste0(dataset_name, ", AUC=",round(cp$AUC, digits=2), "\nsensitivity=", round(cp$sensitivity,digits=2), ",specificity=", round(cp$specificity,digits=2))))
    dev.off()
  }
}

# Unique plot per dataset for supplementary  ---------------------------------
ordered_signatures <- c(
  "Nusse_PDE2A", "KeggWntPathwayGenes", "LEF1", "TCF7L2", "TCF7L1", 
  "TCF7", "SP5", "LGR5", "AXIN2", "APC", "ZNRF3", "CTNNB1"
)

for(dataset_name in c("Assie","Heaton","tcga")){
  print(dataset_name)
  
  genes <- all_data$GeneSymbol
  dataset <- all_data[,grep(dataset_name,colnames(all_data))]
  colnames(dataset) <- stringr::str_replace_all(colnames(dataset), paste0(dataset_name, "_"), "")
  dataset$GeneSymbol <- genes
  
  if(dataset_name == "Assie"){
    Bcat_obj <- Assie_Bcat
    alterations <- c("CTNNB1","wt")}
  if(dataset_name == "Heaton"){
    Bcat_obj <- Heaton_Bcat
    alterations <- c("Nuclear","Membrane")}
  if(dataset_name == "tcga"){
    Bcat_obj <- tcga_Bcat
    alterations <- c("CTNNB1_mut","wt")}
  if(dataset_name == "Lefevre"){
    Bcat_obj <- Lefevre_Bcat
    alterations <- c("expressed","repressed")}
  
  # only keep available signatures
  sig_list <- list_signatures[ordered_signatures]
  
  all_plots <- list()   # store plots here
  plot_index <- 1
  
  for(signatureName in ordered_signatures){
    
    signature_genes <- sig_list[[signatureName]]
    if(is.null(signature_genes)) next
    
    print(signatureName)
    
    ### Compute signature scores
    if(length(signature_genes) > 1){
      index <- which(colnames(dataset) == "GeneSymbol")
      signature <- dataset[which(dataset$GeneSymbol %in% signature_genes), -index]
      mean <- apply(as.matrix(signature), 2, mean)
      mean <- data.frame(Sample = names(mean), Signature = as.numeric(mean))
      input <- merge(mean, Bcat_obj, by="Sample")
      input <- input[input$Alteration %in% alterations,]
    } else {
      index <- which(colnames(dataset) == "GeneSymbol")
      singleGene <- dataset[which(dataset$GeneSymbol==signature_genes), -index]
      singleGene <- data.frame(Sample = names(singleGene),
                               Signature = as.numeric(singleGene))
      input <- merge(singleGene, Bcat_obj, by="Sample")
      input <- input[input$Alteration %in% alterations,]
    }
    
    ### Run cutpointr
    if(metric == "sum_sens_spec"){
      cp <- cutpointr(input, Signature, Alteration,
                      method = maximize_metric, metric = sum_sens_spec)
    }
    if(metric == "sens_constrain"){
      cp <- cutpointr(input, Signature, Alteration,
                      method = maximize_metric, metric = sens_constrain)
    }
    
    saveRDS(file=paste0(out,"/Cutpointr_", metric, "_", dataset_name, "_", signatureName, ".RDS"), cp)
    
    ### Create the plot
    p <- plot_Density_byBcatStatus(
      dataset=dataset,
      dataset_name=dataset_name,
      signature=signature_genes,
      signatureName=signatureName,
      Bcat=Bcat_obj
    ) +
      ggtitle(paste0(dataset_name,
                     ", AUC=", round(cp$AUC,2),
                     "\nsensitivity=", round(cp$sensitivity,2),
                     ", specificity=", round(cp$specificity,2)))
    
    all_plots[[plot_index]] <- p
    plot_index <- plot_index + 1
  }
  
  ### -------------------------------------------------------------
  ### COMBINE ALL SIGNATURE PLOTS INTO ONE PAGE (3 COLUMNS)
  ### -------------------------------------------------------------
  combined_plot <- wrap_plots(all_plots, ncol = 3)
  
  ### SAVE TO A SINGLE PAGE
  #pdf(file = paste0(out, "/GeneSignaturesDensityDistribution_", dataset_name, "_AvgProbes_", metric, ".pdf"), width = 14, height = 12, useDingbats=FALSE)
  png(file=paste0(out,"/GeneSignaturesDensityDistribution_perCTNNB1Status_",dataset_name,"_allGenesSignatures_AvgProbes_", metric,".png"), width=1200, height = 1200)
  print(combined_plot)
  
  dev.off()
}


# HEATMAP AUC ROC values of genes/signatures of interest -------------------------------------------------------
heatmapAUC <- as.data.frame(matrix(ncol=3))
colnames(heatmapAUC) <- c("dataset","signature","AUC")
for(dataset_name in c("Assie","Heaton","tcga")){
  
  for(signatureName in c("CTNNB1","ZNRF3","APC","AXIN2","LGR5","SP5","TCF7","TCF7L1","TCF7L2","KeggWntPathwayGenes", "Nusse_PDE2A", "LEF1", "signature","pancancer3gene")){
    cp <- readRDS(file=paste0(out,"/Cutpointr_", metric, "_", dataset_name, "_", signatureName, ".RDS"))
    vec <- c(dataset_name, signatureName, round(cp$AUC, digits=3)) 
    names(vec) <- c("dataset","signature","AUC")
    heatmapAUC <- rbind(heatmapAUC, vec)
  }
  
}
heatmapAUC <- heatmapAUC[-1,]
heatmapAUC$AUC <- as.numeric(as.vector(heatmapAUC$AUC))
saveRDS(file=paste0(out, "HeatmapAUCs_dataset_signature_",metric,".RDS"), heatmapAUC)

heatmapAUC <- readRDS(file=paste0(out, "HeatmapAUCs_dataset_signature_",metric,".RDS"))
heatmapAUC <- heatmapAUC[heatmapAUC$signature %in% c("CTNNB1","ZNRF3","APC","AXIN2","LGR5","SP5","TCF7","TCF7L1","TCF7L2","KeggWntPathwayGenes", "Nusse_PDE2A", "LEF1", "signature"),]
heatmapAUC$signature <- factor(heatmapAUC$signature, levels=c("CTNNB1","ZNRF3","APC","AXIN2","LGR5","SP5","TCF7","TCF7L1","TCF7L2","LEF1","KeggWntPathwayGenes", "Nusse_PDE2A", "signature"))
png(paste0(out, "HeatmapAUCs_dataset_GenesSignature_ForSuppFig_",metric, ".png"))
ggplot(data=heatmapAUC, aes(x=dataset,y=signature,fill=AUC)) + geom_tile() +
  geom_text(aes(label=AUC)) + theme_bw() +
  scale_fill_gradient2(low="navy", mid="white", high="red", midpoint = 0.6)
dev.off()
