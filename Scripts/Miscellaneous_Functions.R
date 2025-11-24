# VERSION 2: MSigDB hallmarks enrichment with fgsea -------------------------------------------------------
## Changes of this second version: i. the possibility to set type of gene identifier, ii. input is ranked list directly
## With signed p-values empty output, with log2FC too many ties
GSEA_MSigDBHallmark_v2 <- function(bulkType, analysis, RankedList, RankedListType, outPath, geneIdType){
  library(msigdbr)
  library(fgsea)
  library(gridExtra)
  library(ggplot2)
  options(digits=15)
  h_gene_sets = msigdbr(species = "Homo sapiens", category = "H")
  if(geneIdType=="gene_symbol"){msigdbr_list_H = split(x = h_gene_sets$gene_symbol, f = h_gene_sets$gs_name)}
  if(geneIdType=="ensembl_id"){msigdbr_list_H = split(x = h_gene_sets$ensembl_gene, f = h_gene_sets$gs_name)}
  
  msigdbr_list_H_nr <- lapply(msigdbr_list_H, unique)
  
  fgseaRes_FC_H <- fgsea(pathways = msigdbr_list_H_nr, 
                         stats    = RankedList)
  
  fgseaRes_FC_H$padj <- as.numeric(as.vector(fgseaRes_FC_H$padj))
  
  for(pval_cutoff in c(0.05, 0.1)){
    png(paste0(outPath, "fgsea_msigdbHallmark_", analysis, "_",bulkType,"_",RankedListType,"_RawP",pval_cutoff,".png"), width=1000, height = 400)
    p1 <- ggplot(data= fgseaRes_FC_H[which((fgseaRes_FC_H$pval<pval_cutoff)),], aes(x=NES, 
                                                                                    y=pathway, 
                                                                                    #colour=pval, 
                                                                                    size=size)) +
      geom_point() + geom_vline(xintercept=0) +
      expand_limits(x=(min(fgseaRes_FC_H$NES))) + theme_bw() +
      labs(x="Normalized Enrichment Score", y="pathway", colour="p value", size="N genes") +
      ggtitle(paste0("enriched pathways raw p-value < ",pval_cutoff))
    p2 <- ggplot(data= fgseaRes_FC_H[which((fgseaRes_FC_H$pval<pval_cutoff)),],aes(x=ES, 
                                                                                   y=pathway, 
                                                                                   colour=padj, 
                                                                                   size=size)) +
      geom_point() + geom_vline(xintercept=0) +
      expand_limits(x=(min(fgseaRes_FC_H$NES))) + theme_bw() +
      labs(x="Enrichment Score", y="pathway", colour="adjusted p value", size="N genes") 
    print(grid.arrange(p1,p2, ncol=2))
    dev.off()
    plotOut <- grid.arrange(p1,p2, ncol=2)
  }
  
  padj_cutoff <- 0.1
  png(paste0(outPath, "fgsea_msigdbHallmark_", analysis, "_",bulkType,"_",RankedListType,"_padj",padj_cutoff,".png"), width=1000, height = 400)
  print(ggplot(data= fgseaRes_FC_H[which((fgseaRes_FC_H$padj<padj_cutoff)),], aes(x=NES, 
                                                                                  y=pathway, 
                                                                                  colour=padj, 
                                                                                  size=size)) +
          geom_point() + geom_vline(xintercept=0) +
          expand_limits(x=(min(fgseaRes_FC_H$NES))) + theme_bw() +
          labs(x="Normalized Enrichment Score", y="pathway", colour="adj pvalue", size="N genes") +
          ggtitle(paste0("enriched pathways adj p-value < ",padj_cutoff)))
  dev.off()
  
  return(fgseaRes_FC_H)
  write.table(file=paste0(outPath, "fgsea_msigdbHallmark_", analysis, "_",bulkType,"_",RankedListType,"_RawP005.txt"), as.data.frame(fgseaRes_FC_H[fgseaRes_FC_H$pval <0.05,c(1:3,5:7)]), col.names = T, row.names = F, quote=F)
}


# Volcano plots public ACC data-----------------------------------------------------------
# Highligthing sign DE B-cat targets or highlighting independently of significance, 9 selected B-cat targets from B-cat activation
# signature (targets sign up in at least 3 of the 4 public ACC data in CTNNB1 mut vs wt comparisons).
# Need in input df with at least colupmn "logFC", "gene", "pFDR
volcanoPlots_BcatTargets <- function(outPath, res, dataset, analysis){
  NusseHerbstPDE2A <- scan(file="C:/Users/LF260934/Documents/WORKING_FOLDER/Results/ACC_Cochin_snRNAseq/MetadataAndCo/BCateninTargets_Nusse_Herbst2014_PDE2A_nr.txt", what=character())
  betacat_targets <- c(NusseHerbstPDE2A, "AFF3")
  betacat_targets <- betacat_targets[order(betacat_targets)]
  res$pFDR <- as.numeric(as.vector(res$pFDR))
  res$logFC <- as.numeric(as.vector(res$logFC))
  
  ### Plot labeling B-cat targets with p_adj<0.1
  targetsToLabel <- res[((res$pFDR<0.1)&(res$gene %in% betacat_targets)),]$gene
  res$label <- ""
  res[res$gene %in% targetsToLabel,]$label <- res[res$gene %in% targetsToLabel,]$gene
  res$col_label <- ifelse(res$label=="", "gray", "black")
  res$col_label <- factor(res$col_label, levels=c("gray","black"))
  plot1 <- ggplot(data = res, aes(x = logFC, y = -log10(pFDR), col = col_label, label=label)) +
    geom_vline(xintercept = c(-0.5, 0.5), col = "black", linetype = 'dashed') +
    geom_hline(yintercept = -log10(0.1), col = "black", linetype = 'dashed') +
    geom_point(size = 1, alpha=0.5) +
    geom_text_repel(max.overlaps = Inf) +
    scale_color_manual(values = c("gray","black")) +
    ggtitle(paste0(dataset, ", ", analysis, " (label sign DE BcatTargets)")) + theme_bw()
  
  # Show selected B-cat targets independently of significance
  res$label1 <- ""
  selectionOnPublicData9genes <- c("AXIN2","AFF3","LEF1","FSCN1", "MYC", "SP5","PDE2A","TBX3","TNFRSF19") # target genes found to be most consistently up in CTNNB1 mut vs wt in Assie,Heaton,TCGA,(Lefevre, only qualitative)
  targetsToLabel <- res[(res$gene %in% selectionOnPublicData9genes),]$gene
  res[res$gene %in% targetsToLabel,]$label1 <- res[res$gene %in% targetsToLabel,]$gene
  res$col_label <- ifelse(res$label1=="", "gray", "black")
  res$col_label <- factor(res$col_label, levels=c("gray","black"))
  plot2 <- ggplot(data = res, aes(x = logFC, y = -log10(pFDR), col = col_label, label=label1)) +
    geom_vline(xintercept = c(-0.5, 0.5), col = "black", linetype = 'dashed') +
    geom_hline(yintercept = -log10(0.1), col = "black", linetype = 'dashed') +
    geom_point(size = 1, alpha=0.5) +
    geom_text_repel(max.overlaps = Inf) +
    scale_color_manual(values = c("gray","black")) + # to set the colours of our variable
    ggtitle(paste0(dataset, ", ", analysis, " (label 9 BcatTargets signature)")) + theme_bw()
  
  return(list(BcatSignificant=plot1, NineSignatureTargets=plot2))
}


################################################################################
# Define function for DE analysis with microarray data ------------------------------
################################################################################
DE_microarray <- function(expr, metadata, condition, microarray_annotation){
  levelsCondition <- c("wt","altered")
  if(!(identical(levels(metadata$condition), levelsCondition))){
    print("condition must be a factor of metadata object with only two levels: wt and altered")
  }
  # Create the SummarizedExperiment object
  affyDataset <- SummarizedExperiment::SummarizedExperiment(assays = expr, colData = metadata)
  # Create a design matrix
  affyDesign <- model.matrix(~0 + condition, data = SummarizedExperiment::colData(affyDataset))
  # Avoid special characters in column names
  colnames(affyDesign) <- make.names(colnames(affyDesign))
  # Create a constrast matrix
  affyContrast <- limma::makeContrasts(conditionaltered-conditionwt, levels=affyDesign)
  # Run differential expression analysis
  affyDEExperiment <- RCPA::runDEAnalysis(affyDataset, method = "limma", design = affyDesign, contrast = affyContrast, annotation = microarray_annotation)
  # Extract the differential analysis result
  affyDEResults <- SummarizedExperiment::rowData(affyDEExperiment)
  return(list(affyDataset=affyDataset,affyDesign=affyDesign,affyContrast=affyContrast,affyDEExperiment=affyDEExperiment,affyDEResults=affyDEResults))
}

# Function heatmaps significant GSEA enrichments -----------------------------------
# To compare significant GSEA enrichments between different datasets/patients/cell types etc
# Arguments:
# geneSetsSuffix = suffix of GSEA output file indicating gene set DB used
# outpath = path to output
# grouping1 = first criterium to separate data (heatmap x axis)
# grouping2 = first criterium to separate data (heatmap facet in facet_grid/facet_wrap)

heatmapCommonGSEAEnrichments <- function(geneSetsSuffix, grouping1, grouping2, Correlation_type){
  ### Plot heatmaps
  index1 <- grep(grouping1, colnames(sub_pos))
  index2 <- grep(grouping2, colnames(sub_pos))
  sub_pos <- all_pos[all_pos$GeneSet == geneSetsSuffix,]
  pos <- ggplot(dat=sub_pos, aes(y=pathway, x=sub_pos[,index1])) + geom_tile() + facet_wrap(~sub_pos[,index2]) +
    theme_classic() +theme(axis.text.x = element_text(angle=45, hjust=1)) +
    xlab("dataset") +
    ggtitle(paste0("GSEA ",geneSet," significant (padj<0.1) positive enrichments\non ",Correlation_type,"correlation B-catenin signature- gene expression"))
  sub_neg <- all_neg[all_neg$GeneSet == geneSet,]
  
  neg <- ggplot(dat=sub_neg, aes(y=pathway, x=sub_neg[,index1])) + geom_tile() + facet_wrap(~sub_neg[,index2]) +
    theme_classic() +theme(axis.text.x = element_text(angle=45, hjust=1)) +
    xlab("dataset") +
    ggtitle(paste0("GSEA ",geneSet," significant (padj<0.1) negative enrichments\non ",Correlation_type,"correlation B-catenin signature- gene expression"))
  
  return(list(positiveEnrichmentsHeatmap=pos, negativeEnrichmentsHeatmap=neg))
}



# Density plot mut vs wt by gene/signature expression ---------------------
# Function to plot density of mean signature expression (or gene expression) by Bcat status (CTNNB1mut vs wt or nuclear vs membrane) --------------------
plot_Density_byBcatStatus <- function(dataset, dataset_name, signature, signatureName, Bcat){
  colnames(Bcat) <- c("TumorID","BcatStatus")
  
  select <- dataset[dataset$GeneSymbol %in% signature,]
  index <-  which(colnames(dataset) == "GeneSymbol")
  select <- as.matrix(select[,-index])
  select <- apply(select, 2, as.numeric)
  if(length(signature)>1){
    mean <- apply(select, 2, mean)
    select <- as.data.frame(cbind(names(mean), mean))
  }else{select <- as.data.frame(cbind(names(select), select))}
  colnames(select) <- c("TumorID","BcatSignature")
  
  all <- merge(select, Bcat, by="TumorID")
  all$BcatSignature <- as.numeric(as.vector(all$BcatSignature))
  all$BcatStatus <- as.factor(all$BcatStatus)
  
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
    colors_manual <- c('red','orange', 'gold1', 'green4')
  }
  if(length(intersect(unique(all$BcatStatus), c("Nuclear","Membrane")))==2){
    all$BcatStatus <- factor(all$BcatStatus, levels=c("Nuclear","Membrane"))
    colors_manual <- c('red', 'green4')
  }
  if(length(intersect(unique(all$BcatStatus), c("CTNNB1_mut","ZNRF3_mut","wt")))==2){
    all$BcatStatus <- factor(all$BcatStatus, levels=c("CTNNB1_mut","wt"))
    colors_manual <- c('red', 'green4')
  }
  if(length(intersect(unique(all$BcatStatus), c("mut","wt")))==2){
    all$BcatStatus <- factor(all$BcatStatus, levels=c("mut","wt"))
    colors_manual <- c('red', 'green4')
  }
  #Plot density distribution B-catenin activation signature separating CTNNB1 mutated vs wt samples
  ggplot(data=all, aes(x=BcatSignature, group=BcatStatus, fill=BcatStatus, color=BcatStatus)) +
    geom_density(alpha=0.5) +
    theme_bw() +
    xlab(signatureName) +
    ggtitle(dataset_name)
}


# Functions survival plots esthetical parameters ----------------------------
theme_juju <- function(){
  theme_classic()+
    theme(plot.title = element_text(hjust=0.5, face = "plain", size = 24),
          axis.text = element_text(color="black"),
          legend.title = element_text(hjust=0.5))
}

theme_cleantable_juju <- function(){
  theme_cleantable()+
    theme(plot.title = element_text(size = 16, hjust = 0.5),
          axis.text = element_text(size = 17))
  
}

# Function to generate fold change ranked list for GSEA -------------------
makeFCrankedList <- function(dataset_name, DE_res_file, outPath, analysis){
  RankedListType <- "FoldChange"
  
  ### Generate ranked list
  res <- readRDS(DE_res_file)
  res <- cbind(rownames(res), res)
  colnames(res) <- c("gene","baseMean","logFC","lfcSE","stat","p.value","pFDR")
  
  ### Generate ranked list
  colNameGenes <- "gene"
  RankedList <- res$logFC
  names(RankedList) <- res[,colNameGenes]
  if(length( which(is.na(RankedList)))>0){
    RankedList <- RankedList[-(which(is.na(RankedList)))]
  }
  RankedList <- sort(RankedList)
  RankedList <- as.data.frame(cbind(names(RankedList),RankedList))
  colnames(RankedList) <- c("GeneName","RankedList")
  write.table(file=paste0(outPath, "FoldChangeRanked_", analysis, "_",dataset_name,".rnk"), RankedList, col.names = T, row.names = T, sep="\t")
}

