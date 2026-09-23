library(tidyr)
library(ggplot2)
library(stats)
library(stringr)
library(gridExtra)
library(ggpubr)
library(smplot2)
library(magrittr)
library(dplyr)

outPath <- "./BoxplotsGeneExpression/"

# 1. GENE EXPRESSION CTNNB1 MUT vs WT TUMORS --------------------------

### Settings ------------------------
# Plot samples with ZNRF3 alterations or not
useZNRF3 <- TRUE
#useZNRF3 <- FALSE
# Genes of which to plot expression

### Read log2 norm counts  --------
in_Assie <- readRDS(file="./ACC_datasets_formatted/Assie/Log2NormData_all_SampleID_ok_AvgMultipleProbesAnyGene.RDS")
in_Heaton <- readRDS(file="./ACC_datasets_formatted/Heaton/Log2NormData_all_SampleID_ok_AvgMultipleProbesAnyGene.RDS")
in_tcga <- readRDS(file="./ACC_datasets_formatted/tcga/Log2Pseudocount1_DESeqNormCounts_tcga.RDS")
in_Lefevre <- readRDS(file="./ACC_datasets_formatted/Lefevre/Log2NormData_all_SampleID_ok_AvgMultipleProbesAnyGene.RDS")

# Keep only CTNNB1 mut or both CTNNB1 mut and ZNRF3 mut
if(useZNRF3 == FALSE){
  # Identify all cases NOT ZNRF3mut and NOT NA (only keep CTNNB1 mut/nuclear and wt/membrane) ----------------------------------------------
  samples_BcatMutWt_TCGA <- readRDS(file="./ACC_datasets_formatted/tcga/Metadata/tcga_sampleIDs_CTNNB1mut_ZNRF3mut_wt.RDS")
  samples_BcatMutWt_TCGA <- samples_BcatMutWt_TCGA[samples_BcatMutWt_TCGA$Alteration %in% c("CTNNB1_mut","wt"),]
  samples_BcatMutWt_Assie <- as.data.frame(readRDS(file="./ACC_datasets_formatted/Assie/Metadata/Assie_sampleIDs_CTNNB1mut_ZNRF3mut_wt.RDS"))
  samples_BcatMutWt_Assie <- samples_BcatMutWt_Assie[samples_BcatMutWt_Assie$Alteration %in% c("CTNNB1","wt"),]
  samples_BcatMutWt_Heaton <- readRDS(file="./ACC_datasets_formatted/Heaton/Metadata/Heaton_sampleIDs_BcatNuclearMembrane_onlyACCs.RDS")
  samples_BcatMutWt_Heaton <- samples_BcatMutWt_Heaton[samples_BcatMutWt_Heaton$Alteration %in% c("Membrane","Nuclear"),]
  
  in_tcga_select <- in_tcga[,c(1,which(colnames(in_tcga) %in% samples_BcatMutWt_TCGA$Sample))]
  in_Assie_select <- in_Assie[,c(1,which(colnames(in_Assie) %in% samples_BcatMutWt_Assie$Sample))]
  in_Heaton_select <- in_Heaton[,c(1,which(colnames(in_Heaton) %in% samples_BcatMutWt_Heaton$Sample))]
  
  Assie_levels <- c("CTNNB1","wt")
  Heaton_levels <- c("Nuclear","Membrane")
  tcga_levels <- c("CTNNB1_mut","wt")
  
}else{
  if(useZNRF3 == TRUE){
    in_tcga_select <- in_tcga
    in_Assie_select <- in_Assie
    in_Heaton_select <- in_Heaton
    
    Assie_levels <- c("CTNNB1","ZNRF3","wt")
    Heaton_levels <- c("Nuclear","Membrane")
    tcga_levels <- c("CTNNB1_mut", "ZNRF3_mut","wt")
  }
}

### Read metadata -----------------------------------------------------------
meta_Assie <- read.csv(file="./ACC_datasets_formatted/Assie/Metadata/Metadata_Assie.csv")
colnames(meta_Assie)[1] <- "Sample"
meta_Heaton <- read.csv(file="./ACC_datasets_formatted/Heaton/Metadata/Metadata_Heaton.csv")
colnames(meta_Heaton)[1] <- "Sample"
meta_Lefevre <- read.csv(file="./ACC_datasets_formatted/Lefevre/Metadata/Metadata_Lefevre.csv")
colnames(meta_Lefevre)[25] <- "Sample"
meta_tcga <- read.csv(file="./ACC_datasets_formatted/tcga/Metadata/Metadata_tcga.csv")
colnames(meta_tcga)[1] <- "Sample"

in_Assie_long  <- pivot_longer(in_Assie_select, cols=2:ncol(in_Assie_select), names_to="Sample", values_to="Value")
input_boxplot_Assie <- merge(in_Assie_long,meta_Assie,by="Sample",all=F)
input_boxplot_Assie <- input_boxplot_Assie[!(is.na(input_boxplot_Assie$CTNNB1_ZNRF3_Alteration)),]
input_boxplot_Assie$BcatStatus <- factor(input_boxplot_Assie$CTNNB1_ZNRF3_Alteration, levels=Assie_levels)
input_boxplot_Assie$Value <- as.numeric(as.vector(input_boxplot_Assie$Value))

in_tcga_long <- pivot_longer(in_tcga_select, cols=2:ncol(in_tcga_select), names_to="Sample", values_to="Value")
input_boxplot_tcga <- merge(in_tcga_long,meta_tcga,by="Sample",all=F)
input_boxplot_tcga <- input_boxplot_tcga[!(is.na(input_boxplot_tcga$BcatStatus)),]
input_boxplot_tcga$BcatStatus <- factor(input_boxplot_tcga$BcatStatus, levels=tcga_levels)
input_boxplot_tcga$Value <- as.numeric(as.vector(input_boxplot_tcga$Value))

in_Heaton_long  <- pivot_longer(in_Heaton_select, cols=2:ncol(in_Heaton_select), names_to="Sample", values_to="Value")
input_boxplot_Heaton <- merge(in_Heaton_long,meta_Heaton[meta_Heaton$Histotype %in% c("ACC"),],by="Sample",all=F)
input_boxplot_Heaton <- input_boxplot_Heaton[!(is.na(input_boxplot_Heaton$BetaCateninStaining)),]
input_boxplot_Heaton$Value <- as.numeric(as.vector(input_boxplot_Heaton$Value))
input_boxplot_Heaton$BcatStatus <- factor(input_boxplot_Heaton$BetaCateninStaining, levels=Heaton_levels)

in_Lefevre_long  <- pivot_longer(in_Lefevre, cols=2:ncol(in_Lefevre), names_to="Sample", values_to="Value")
input_boxplot_Lefevre <- merge(in_Lefevre_long, meta_Lefevre, by="Sample",all=F)
input_boxplot_Lefevre$Value <- as.numeric(as.vector(input_boxplot_Lefevre$Value))
input_boxplot_Lefevre$BcatStatus <- factor(input_boxplot_Lefevre$BcatStatus, levels=c("expressed","repressed","ctrl","ctrlDoxy"))

### Define function to plot TCGA boxplots ---------------------------------
plotTCGA <- function(input_boxplot, gene, analysis){
  
  n_CTNNB1 <- length(unique(input_boxplot[input_boxplot$BcatStatus=="CTNNB1_mut",]$Sample))
  n_ZNRF3 <- length(unique(input_boxplot[input_boxplot$BcatStatus=="ZNRF3_mut",]$Sample))
  n_wt <- length(unique(input_boxplot[input_boxplot$BcatStatus=="wt",]$Sample))

  df_sub <- input_boxplot[input_boxplot$GeneSymbol==gene, ]
  if(nrow(df_sub)>0){
    
    df_sub$BcatStatus <- stringr::str_replace(df_sub$BcatStatus, "CTNNB1_mut", paste0("CTNNB1 (n=",n_CTNNB1, ")"))
    df_sub$BcatStatus <- stringr::str_replace(df_sub$BcatStatus, "ZNRF3_mut", paste0("ZNRF3 (n=",n_ZNRF3, ")"))
    df_sub$BcatStatus <- stringr::str_replace(df_sub$BcatStatus, "wt", paste0("wt (n=", n_wt, ")"))
    df_sub$BcatStatus <- factor(df_sub$BcatStatus, levels=c(paste0("CTNNB1 (n=", n_CTNNB1, ")"),
                                                            paste0("ZNRF3 (n=", n_ZNRF3, ")"),
                                                            paste0("wt (n=", n_wt, ")")))
    
    my_comparisons <- list(c(levels(df_sub$BcatStatus)[1],levels(df_sub$BcatStatus)[3]))
    p <- ggplot(data=df_sub, aes(y=Value, x=BcatStatus)) + geom_boxplot(fill=c("indianred","rosybrown","lightblue")) +
      geom_jitter() + 
      labs(title="TCGA",
      y=paste(gene, "expression (log2 DESeq norm counts + 1)"),
           x="") +
      stat_compare_means(comparisons=my_comparisons, method="wilcox.test") +
      theme_bw()
    return(p)
  }
}

### Define function to plot Assie boxplots ---------------------------------
plotAssie <- function(input_boxplot, gene, analysis){
  
  n_wt <- length(unique(input_boxplot[input_boxplot$BcatStatus %in% c("wt"),]$Sample)) # 23
  n_CTNNB1 <- length(unique(input_boxplot[input_boxplot$BcatStatus %in% c("CTNNB1"),]$Sample)) # 8
  n_ZNRF3 <- length(unique(input_boxplot[input_boxplot$BcatStatus %in% c("ZNRF3"),]$Sample)) # 8
  
  df_sub <- input_boxplot[input_boxplot$GeneSymbol==gene, ]
  if(nrow(df_sub)>0){
    
    df_sub$BcatStatus <- stringr::str_replace(df_sub$BcatStatus, "CTNNB1", paste0("CTNNB1 (n=",n_CTNNB1, ")"))
    df_sub$BcatStatus <- stringr::str_replace(df_sub$BcatStatus, "ZNRF3", paste0("ZNRF3 (n=", n_ZNRF3, ")"))
    df_sub$BcatStatus <- stringr::str_replace(df_sub$BcatStatus, "wt", paste0("wt (n=", n_wt, ")"))
    df_sub$BcatStatus <- factor(df_sub$BcatStatus, levels=c(paste0("CTNNB1 (n=",n_CTNNB1, ")"),
                                                            paste0("ZNRF3 (n=",n_ZNRF3, ")"),
                                                            paste0("wt (n=", n_wt, ")")))
    
    my_comparisons <- list(c(levels(df_sub$BcatStatus)[1],levels(df_sub$BcatStatus)[3]))
    p <- ggplot(data=df_sub, aes(y=Value, x=BcatStatus)) + geom_boxplot(fill=c("indianred","rosybrown","lightblue")) +
      geom_jitter() + 
      stat_compare_means(comparisons=my_comparisons, method="wilcox.test") +
      labs(title="Assie",
           y=paste(gene, "expression (log2 norm counts)"),
           x="") +
      theme_bw()
    return(p)
  }
}


### Define function to plot Lefevre boxplots ---------------------------------
plotLefevre <- function(input_boxplot, gene, n_expr, n_repr, n_expr_ctrlshRNA_noDoxy, n_expr_ctrlshRNA_withDoxy){
  df_sub <- input_boxplot[input_boxplot$GeneSymbol==gene, ]
  if(nrow(df_sub)>0){
    
    df_sub$BcatStatus <- stringr::str_replace(df_sub$BcatStatus, "expressed", paste0("expressed (n=",n_expr, ")"))
    df_sub$BcatStatus <- stringr::str_replace(df_sub$BcatStatus, "repressed", paste0("repressed (n=", n_repr, ")"))
    df_sub$BcatStatus <- stringr::str_replace(df_sub$BcatStatus, "ctrl$", paste0("ctrl (n=", n_expr_ctrlshRNA_noDoxy, ")"))
    df_sub$BcatStatus <- stringr::str_replace(df_sub$BcatStatus, "ctrlDoxy", paste0("ctrlDoxy (n=", n_expr_ctrlshRNA_withDoxy, ")"))
    df_sub$BcatStatus <- factor(df_sub$BcatStatus, levels=c(paste0("expressed (n=",n_expr, ")"),
                                                            paste0("repressed (n=", n_repr, ")"),
                                                            paste0("ctrl (n=", n_expr_ctrlshRNA_noDoxy, ")"),
                                                            paste0("ctrlDoxy (n=", n_expr_ctrlshRNA_withDoxy, ")")))
    
    my_comparisons <- list(c(levels(df_sub$BcatStatus)[1],levels(df_sub$BcatStatus)[2]))
    p <- ggplot(data=df_sub, aes(y=Value, x=BcatStatus)) + geom_boxplot(fill=c("indianred","lightblue", "orange1", "orange4")) +
      geom_jitter() + 
      stat_compare_means(comparisons=my_comparisons, method="t.test") +
      labs(title="Lefevre (H295R cell line)",
           y=paste(gene, "expression (log2 norm counts)"),
           x="") +
      theme_bw()
    return(p)
  }
}

### Define function to plot Heaton boxplots ---------------------------------
plotHeaton <- function(input_boxplot, gene, analysis){
  
  n_nucl <- length(unique(input_boxplot[(input_boxplot$BetaCateninStaining=="Nuclear"),]$Sample))
  n_membr <- length(unique(input_boxplot[(input_boxplot$BetaCateninStaining=="Membrane"),]$Sample))

  df_sub <- input_boxplot[input_boxplot$GeneSymbol==gene, ]
  if(nrow(df_sub)>0){
     
    df_sub$BcatStatus <- stringr::str_replace(df_sub$BcatStatus, "Nuclear", paste0("Nuclear (n=", n_nucl, ")"))
    df_sub$BcatStatus <- stringr::str_replace(df_sub$BcatStatus, "Membrane", paste0("Membrane (n=", n_membr, ")"))
    
    df_sub$BcatStatus <- factor(df_sub$BcatStatus, levels=c(paste0("Nuclear (n=", n_nucl, ")"),
                                                            paste0("Membrane (n=", n_membr, ")")))
    
    my_comparisons <- list(c(levels(df_sub$BcatStatus)[1],levels(df_sub$BcatStatus)[2]))
    p <- ggplot(data=df_sub, aes(y=Value, x=BcatStatus)) + 
      geom_boxplot(fill=c("indianred","lightblue")) +
      geom_jitter() +
      stat_compare_means(comparisons=my_comparisons, method="wilcox.test") +
      labs(title="Heaton",
           y=paste(gene, "expression (log2 norm counts)"),
           x="") +
      theme_bw()
    return(p)
  }
}


# Plot expression known B-catenin targets ------------------------------
genesOfInterest <- c("ABCB1", "AFF3", "AXIN2", "BCL2L2", "BIRC5", "CCND1", "CDC25A", "CDKN2A", "CDX1", "CLDN1", "CTLA4", "DKK1", "EDN1", "ENAH", "ENC1", "FGF18", 
                 "FGF4", "FGFBP1", "FOSL1", "FSCN1", "FST", "FZD7", "GBX2", "HES1", "HNF1A", "ID2", "JAG1", "JUN", "KRT5", "L1CAM","LAMC2", "LEF1", "LGR5", "MMP14", 
                 "MMP7", "MYC", "MYCBP", "NEDD9", "NEUROD1", "NEUROG1", "NOS2", "NOTCH2", "NRCAM", "PDE2A", "PLAU", 
                 "PPARD", "PTGS2", "S100A4", "SGK1", "SMC3", "SP5", "SUZ12", "TBX3", "TBXT", "TCF4", "TERT", "TNC", "TNFRSF19", "VCAN", "VEGFA", 
                 "YY1AP1") # known B-catenin targets (those reported on Nusse lab website + Herbst et al 2014 paper + PDE2A + AFF3)
analysis <- "BcatMutWt"

pdf(file=paste0(outPath, "BoxplotExpr_allPublicACCDatasets_",analysis,"_AvgProbes_KnownBCateninTargets.pdf"), width=10, height = 10)
for(gene in genesOfInterest){
  print(gene)
  foundAssie <- nrow(input_boxplot_Assie[input_boxplot_Assie$GeneSymbol == gene,])
  foundHeaton <- nrow(input_boxplot_Heaton[input_boxplot_Heaton$GeneSymbol == gene,])
  foundTcga <- nrow(input_boxplot_tcga[input_boxplot_tcga$GeneSymbol == gene,])
  foundLefevre<- nrow(input_boxplot_Lefevre[input_boxplot_Lefevre$GeneSymbol == gene,])
  if( (foundAssie>0) & (foundHeaton>0) & (foundTcga>0) & (foundLefevre>0) ){
    tcga <- plotTCGA(input_boxplot=input_boxplot_tcga, gene=gene, analysis=analysis)
    Assie <- plotAssie(input_boxplot=input_boxplot_Assie, gene=gene, analysis=analysis)
    Heaton <- plotHeaton(input_boxplot=input_boxplot_Heaton, gene=gene, analysis=analysis)
    Lefevre <- plotLefevre(input_boxplot=input_boxplot_Lefevre, gene=gene, n_expr=3, n_repr=3, n_expr_ctrlshRNA_noDoxy=1, n_expr_ctrlshRNA_withDoxy=1)
    title1=text_grob(gene, size = 20, face = "bold") 
    print(grid.arrange(tcga,Assie,Heaton,Lefevre, nrow=2, ncol=2, top=title1))
  }
}
dev.off()

# Plot expression other genes of interest ------------------------------
genesOfInterest <- c("MYC","CDK6","DACH1","ABR","EFNA3","FAM169A","JARID2","LTBP1","SYTL2")

analysis <- "BcatMutWt"
for(gene in genesOfInterest){
  print(gene)
  foundAssie <- nrow(input_boxplot_Assie[input_boxplot_Assie$GeneSymbol == gene,])
  foundHeaton <- nrow(input_boxplot_Heaton[input_boxplot_Heaton$GeneSymbol == gene,])
  foundTcga <- nrow(input_boxplot_tcga[input_boxplot_tcga$GeneSymbol == gene,])
  foundLefevre<- nrow(input_boxplot_Lefevre[input_boxplot_Lefevre$GeneSymbol == gene,])
  if( (foundAssie>0) & (foundHeaton>0) & (foundTcga>0) ){
    png(file=paste0(outPath, gene, "_BoxplotExpr_allPublicCochinACCDatasets_",analysis,"_2025_04_25_DefForPaper_AvgMultipleProbes.png"),  width=1800, height=400)
    tcga <- plotTCGA(input_boxplot=input_boxplot_tcga, gene=gene, analysis=analysis)
    Assie <- plotAssie(input_boxplot=input_boxplot_Assie, gene=gene, analysis=analysis)
    Heaton <- plotHeaton(input_boxplot=input_boxplot_Heaton, gene=gene, analysis=analysis)
    Lefevre <- plotLefevre(input_boxplot=input_boxplot_Lefevre, gene=gene, n_expr=3, n_repr=3, n_expr_ctrlshRNA_noDoxy=1, n_expr_ctrlshRNA_withDoxy=1)
    title1=text_grob(gene, size = 20, face = "bold") 
    grid.arrange(Assie,Heaton,tcga,Lefevre, nrow=1, ncol=4, top=title1)
    dev.off()
  }else{print("Gene not detected in all datasets")}
}



# 2. GENES EXPRESSION ACROSS HYSTOTYPE ------------------------------------


### Settings ----------------------------------------------------------------
inPath <- "./ACC_datasets_formatted/"
outPath <- "./BoxplotsGeneExpression/"

signature <- c("AFF3","AXIN2","FSCN1","LEF1","MYC","PDE2A","SP5","TBX3","TNFRSF19")
genes_of_interest <- "MYC"


### Plot gene expression across histotypes, with signature levels -----------
for(dataset_name in c("Caramuta","JouinotFFPE","Demeure","Assie","Heaton","tcga")){
  print(dataset_name)
  
  # Read metadata and counts
  metadata <- read.csv(file=paste0(inPath, dataset_name, "/Metadata/Metadata_", dataset_name,".csv"))
  metadata$Histotype <- stringr::str_replace_all(metadata$Histotype, "ACTH-independant Macronodular Adrenal Hyperplasia", "Other")
  metadata$Histotype <- stringr::str_replace_all(metadata$Histotype, "Uncertain malignant potential ACT", "Other")
  
  if(dataset_name == "tcga"){
    dataset <- readRDS(file=paste0(inPath, dataset_name, "/Log2Pseudocount1_DESeqNormCounts_tcga.RDS"))
  }else{
    dataset <- readRDS(file=paste0(inPath, dataset_name, "/Log2NormData_all_SampleID_ok.RDS"))
  }
  
  # For genes with multiple probes use their mean value
  if(length(dataset$GeneSymbol) > length(unique(dataset$GeneSymbol))){
    dataset_nr <- dataset %>%
      group_by(GeneSymbol) %>%
      mutate(across(!(starts_with(c("GeneSymbol"))), ~ mean(.x, na.rm = TRUE), .names = "avg_{.col}")) %>%
      ungroup() %>% select(starts_with(c("GeneSymbol","avg"))) %>% unique()
    
    colnames(dataset_nr) <- stringr:: str_replace_all(colnames(dataset_nr), "avg_", "")
  }else{
    dataset_nr <- dataset
  }
  
  # Calculate B-Catenin signature  mean
  select <- dataset_nr[which(dataset_nr$GeneSymbol %in% signature),-1]
  select <- as.matrix(select)
  select <- apply(select, 2, as.numeric)
  mean <- apply(select, 2, mean, na.rm=T)
  mean <- as.data.frame(cbind(names(mean), mean))
  colnames(mean) <- c("TumorID","BcatSignatureMean")
  metadata <- merge(metadata, mean, by="TumorID", all=T)
  metadata$BcatSignatureMean <- as.numeric(as.vector(metadata$BcatSignatureMean))
  
  # Extract expression of gene of interest
  for(gene in genes_of_interest){
    print(gene)
    select <- dataset_nr[which(dataset_nr$GeneSymbol %in% c(gene)),-1]
    if(nrow(select)>0){
      select <- apply(select, 2, as.numeric)
      select <- as.data.frame(cbind(names(select),as.numeric(as.vector(select))))
      colnames(select) <- c("TumorID",gene)
      metadata_select <- merge(metadata, select, by="TumorID", all=T)
      metadata_select[,gene] <- as.numeric(as.vector(metadata_select[,gene]))
      
      png(paste0(outPath, gene, "_Signature_histotype_", dataset_name, ".png"), width=600)
      print(ggplot(data=metadata_select, aes(x=Histotype,y=get(gene))) +
        geom_boxplot() + geom_jitter() +
        scale_color_gradientn(colours=c("green4","yellow","red")) +
        ggtitle(dataset_name) +
        theme_bw())
      
      dev.off()
    }
  }
}


### Plot gene expression across histotypes (FOR FIGURE S8) -----------
genes_of_interest <- "MYC"
for(dataset_name in c("Caramuta","JouinotFFPE","Demeure","Heaton")){
  print(dataset_name)
  
  # Read metadata and counts
  metadata <- read.csv(file=paste0(inPath, dataset_name, "/Metadata/Metadata_", dataset_name,".csv"))
  metadata$Histotype <- stringr::str_replace_all(metadata$Histotype, "ACTH-independant Macronodular Adrenal Hyperplasia", "Other")
  metadata$Histotype <- stringr::str_replace_all(metadata$Histotype, "Uncertain malignant potential ACT", "Other")
  
  dataset <- readRDS(file=paste0(inPath, dataset_name, "/Log2NormData_all_SampleID_ok.RDS"))
  
  # For genes with multiple probes use their mean value
  if(length(dataset$GeneSymbol) > length(unique(dataset$GeneSymbol))){
    dataset_nr <- dataset %>%
      group_by(GeneSymbol) %>%
      mutate(across(!(starts_with(c("GeneSymbol"))), ~ mean(.x, na.rm = TRUE), .names = "avg_{.col}")) %>%
      ungroup() %>% select(starts_with(c("GeneSymbol","avg"))) %>% unique()
    
    colnames(dataset_nr) <- stringr:: str_replace_all(colnames(dataset_nr), "avg_", "")
  }else{
    dataset_nr <- dataset
  }
  
  # Extract expression of gene of interest
  for(gene in genes_of_interest){
    print(gene)
    select <- dataset_nr[which(dataset_nr$GeneSymbol %in% c(gene)),-1]
    if(nrow(select)>0){
      select <- apply(select, 2, as.numeric)
      select <- as.data.frame(cbind(names(select),as.numeric(as.vector(select))))
      colnames(select) <- c("TumorID",gene)
      metadata_select <- merge(metadata, select, by="TumorID", all=T)
      metadata_select[,gene] <- as.numeric(as.vector(metadata_select[,gene]))
      
      pdf(paste0(outPath, gene, "_byHistotype_", dataset_name, ".pdf"), width=10, height = 10, useDingbats = FALSE)
      print(ggplot(data=metadata_select, aes(x=Histotype, y=get(gene))) +
              geom_boxplot() + geom_jitter() +
              #scale_color_gradientn(colours=c("green4","yellow","red")) +
              ggtitle(dataset_name) +
              ylab("MYC mRNA expression") +
              theme_bw() +
              theme(axis.text = element_text(size = 17), axis.title = element_text(size = 20), plot.title = element_text(size = 30)))
      
      dev.off()
    }
  }
}


