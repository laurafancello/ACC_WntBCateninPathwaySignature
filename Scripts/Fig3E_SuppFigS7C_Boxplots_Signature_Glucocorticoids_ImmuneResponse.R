library(ggpubr)
library(smplot2)
library(ggpubr)
library(smplot2)
library(msigdbr)
library(dplyr)

# SETTINGS --------------
### Set nb groups to split into levels of signature expression
n_BcatLevels <- 2
#n_BcatLevels <- 3

norm <- "Zscore"

# READ INPUT --------------------------------------------------------------
inPath <- "./ACC_datasets_integration/"
outPath <- "./ACC_datasets_integration/ImmuneResponseAnalysis/"

metadata_all <- readRDS(file=paste0(inPath, "Metadata_allDatasets.RDS"))
metadata_all$CortisolYesNo <- factor(metadata_all$CortisolYesNo, levels=c("no","yes",NA))
metadata_all$HormoneExpressionYesNo <- factor(metadata_all$HormoneExpressionYesNo, levels=c("no","yes","NA"))
metadata_all$BcatStatusGeneral <- factor(metadata_all$BcatStatusGeneral, levels=c("Active","Inactive","ZNRF3_mut",NA)) # convert into factor for my_comparisons object to be used in stat_compare_means
metadata_all <- metadata_all[!((metadata_all$DatasetName=="tcga")&(is.na(metadata_all$BcatSignatureMeanZScore))),] # remove TCGAs for which I have metadata but no counts (otherwise quantile not working with NAs)

dataset_all <- readRDS(file=paste0(inPath, "Counts_", norm, "_allDatasets.RDS"))
df <- dataset_all[,-which(colnames(dataset_all) == "GeneSymbol")]
rownames(df) <- dataset_all$GeneSymbol

h_gene_sets = msigdbr(species = "Homo sapiens", category = "H")
msigdbr_list_H = split(x = h_gene_sets$gene_symbol, f = h_gene_sets$gs_name)
msigdbr_list_H_nr <- lapply(msigdbr_list_H, unique)

# Split into signature expression levels ----------------------------------
quantiles <- quantile(metadata_all$BcatSignatureMeanZScore)
if(n_BcatLevels == 2){
  metadata_all$BcatSignatureMeanZScore_group <- ifelse(metadata_all$BcatSignatureMeanZScore>=quantiles[3], "high", "low")
  metadata_all$BcatSignatureMeanZScore_group <- factor(metadata_all$BcatSignatureMeanZScore_group, levels=c("high","low"))
  my_comparisons_signature <- list(c(levels(metadata_all$BcatSignatureMeanZScore_group)[1],levels(metadata_all$BcatSignatureMeanZScore_group)[2]))
}else{
  if(n_BcatLevels == 3){
    metadata_all$BcatSignatureMeanZScore_group <- NA
    metadata_all[metadata_all$BcatSignatureMeanZScore<quantiles[2],]$BcatSignatureMeanZScore_group <- "low"
    metadata_all[((metadata_all$BcatSignatureMeanZScore>=quantiles[2])&(metadata_all$BcatSignatureMeanZScore<quantiles[3])),]$BcatSignatureMeanZScore_group <- "mid"
    metadata_all[metadata_all$BcatSignatureMeanZScore>=quantiles[3],]$BcatSignatureMeanZScore_group <- "high"
    metadata_all$BcatSignatureMeanZScore_group <- factor(metadata_all$BcatSignatureMeanZScore_group, levels=c("high","mid","low"))
    my_comparisons_signature <- list(c(levels(metadata_all$BcatSignatureMeanZScore_group)[1],levels(metadata_all$BcatSignatureMeanZScore_group)[3]))
  }
}


# Boxplots (or corrleation)  by Hallmark gene sets expression by signature levels, and/or by hormone expression ---------------------
gene_set <-  "LeadEdge"
hormones <- "cortisol"
suffix <- "Cortisol"
var<- "CortisolYesNo"

for(pathway in c("HALLMARK_INFLAMMATORY_RESPONSE","HALLMARK_IL2_STAT5_SIGNALING","HALLMARK_INTERFERON_ALPHA_RESPONSE","HALLMARK_INTERFERON_GAMMA_RESPONSE","HALLMARK_IL6_JAK_STAT3_SIGNALING","HALLMARK_ALLOGRAFT_REJECTION")){
  
  genes <- scan(file=paste0("./ACC_datasets_integration/GSEAs/LeadingEdge_gsea_ZscoreNormIntegrated_",pathway,".txt"), what=character())
  select <- df[which(rownames(df) %in% genes),]
  select <- as.matrix(select)
  select <- apply(select, 2, as.numeric)
  mean <- apply(select, 2, mean, na.rm=T)
  mean <- as.data.frame(cbind(names(mean), mean))
  colnames(mean) <- c("TumorID", paste0("LeadEdge", pathway))
  metadata_all <- merge(metadata_all, mean, by="TumorID", all.x=T, all.y=F)
  metadata_all[,paste0("LeadEdge", pathway)] <- as.numeric(as.vector(metadata_all[,paste0("LeadEdge", pathway)]))
  
  metadata_all$CortisolYesNo <- factor(metadata_all$CortisolYesNo, levels=c("no","yes",NA))
  my_comparisons_hormones <- list(c(levels(metadata_all$CortisolYesNo)[1],levels(metadata_all$CortisolYesNo)[2]))
  
  group_counts <- metadata_all %>%
    group_by(BcatSignatureMeanZScore_group, get(var)) %>%
    summarise(n = n(), .groups = "drop")
  
  print(paste0(gene_set, " ", pathway, ", ", hormones, ", n_BcatLevels=", n_BcatLevels, "--------------------"))
  print(group_counts)
  
  ### Compare signature groups within hormone groups
  png(file=paste0(outPath, "Boxplot_BCatSignature",n_BcatLevels,"Groups_", gene_set, pathway,"_",suffix,".png"), width=1000)
  print(ggplot(data=metadata_all, aes(x=BcatSignatureMeanZScore_group , y=get(paste0(gene_set, pathway)), fill=BcatSignatureMeanZScore_group)) +
          geom_boxplot() + geom_jitter(width=0.3) + theme_bw(base_size = 14) +  ggtitle(suffix) + 
          facet_wrap(.~get(var), scales = "free_x") + ylab(paste0(gene_set, pathway)) +
          stat_compare_means(comparisons=my_comparisons_signature, method="wilcox.test") +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)))
  dev.off()
  
  ### Compare hormone groups within signature groups
  # png(file=paste0(outPath, "Boxplot_",suffix,"_BCatSignature",n_BcatLevels,"Groups_", gene_set, pathway,".png"), width=1000)
  # print(ggplot(data=metadata_all, aes(x=get(var) , y=get(paste0(gene_set, pathway)))) +
  #         geom_boxplot() + geom_jitter(width=0.3) + theme_bw(base_size = 14) +
  #         ggtitle("Signature expression levels") + 
  #         facet_wrap(.~BcatSignatureMeanZScore_group) + ylab(paste0(gene_set, pathway)) +
  #         stat_compare_means(comparisons=my_comparisons_hormones, method="wilcox.test") +
  #         theme(axis.text.x = element_text(angle = 45, hjust = 1))) 
  # dev.off()
  
}


# Boxplots (or correlation) by MHC, immune checkpoint, GR, etc genes of interest by signature levels, and/or by hormone expression ---------------------
genesOfInterest <- list(majorHistocompatibilityComplex=c("HLA-G","HLA-F","HLA-E","HLA-DRB1","HLA-DRB5","HLA-DRA","HLA-DQA1","HLA-DQA2","HLA-DQB1","HLA-DQB2","HLA-DPA1","HLA-DPB1","HLA-DOA","HLA-DOB","HLA-DMA","HLA-DMB","HLA-A","HLA-B","HLA-C"),
                        TcellStimulator=c("TNFRSF9","TNFRSF8","TNFRSF4","TNFRSF25","TNFRSF18","TNFRSF14","ICOS","CD40LG","CD28","CD27","CD226","CD2"),
                        ImmuneCheckpointImmunotherMarker = c("CD274","CTLA4","HAVCR2","LAG3","PDCD1","PDCD1LG2","TIGIT","SIGLEC15","TGFB1","TGFBR2","ACVR1","MLH1","MSH2","MSH6","PMS2","CD47"),
                        GlucocorticoidReceptor = c("NR3C1"),
                        PDL1=c("CD274"),
                        GlucocorticoidActivitySignature_down = c("PTS","PCP4","FDX1","IFI27L2","NR4A3","CCN3","ALDOC","QPCT","ZNF711","FDPS","HSPD1","ATP4A","AMDHD1","EFNB1","FATE1","DPM3","CCDC107","IFI27L1","FADS1","MRPL41","ACAT2","TNFSF13B","PDE9A","MOB4","ACYP1","LAMC3","ERG28","BMP4","IFI27","RPP40","KCNK3"),
                        GlucocorticoidActivitySignature_up = c("SGPP2","ENPP5","COL14A1"))
gene_set <- ""
hormones <- "cortisol"
suffix <- "Cortisol"
var<- "CortisolYesNo"
metadata_all$CortisolYesNo <- factor(metadata_all$CortisolYesNo, levels=c("no","yes",NA))
my_comparisons_hormones <- list(c(levels(metadata_all$CortisolYesNo)[1],levels(metadata_all$CortisolYesNo)[2]))

for(pathway in names(genesOfInterest)){
   genes <- genesOfInterest[[pathway]]
  if(length(genes)>1){
    select <- df[which(rownames(df) %in% genes),]
    select <- as.matrix(select)
    select <- apply(select, 2, as.numeric)
    mean <- apply(select, 2, mean, na.rm=T)
    mean <- as.data.frame(cbind(names(mean), mean))
    colnames(mean) <- c("TumorID", pathway)
    metadata_all <- merge(metadata_all, mean, by="TumorID", all.x=T, all.y=F)
    metadata_all[,pathway] <- as.numeric(as.vector(metadata_all[,pathway]))
  }else{
    select <- t(df[which(rownames(df) %in% genes),])
    select <- as.data.frame(cbind(rownames(select), as.numeric(as.vector(select))))
    colnames(select) <- c("TumorID", pathway)
    select[,pathway] <- as.numeric(as.vector(select[,pathway]))
    metadata_all <- merge(metadata_all, select, by="TumorID", all.x=T, all.y=F)
  }
    
    group_counts <- metadata_all %>%
      group_by(BcatSignatureMeanZScore_group, get(var)) %>%
      summarise(n = n(), .groups = "drop")
    
    print(paste0(gene_set, " ", pathway, ", ", hormones, ", n_BcatLevels=", n_BcatLevels, "--------------------"))
    print(group_counts)
    
    ### Compare signature groups within hormone groups
    png(file=paste0(outPath, "Boxplot_BCatSignature",n_BcatLevels,"Groups_", gene_set, pathway,"_",suffix,".png"), width=1000)
    print(ggplot(data=metadata_all, aes(x=BcatSignatureMeanZScore_group , y=get(paste0(gene_set, pathway)), fill=BcatSignatureMeanZScore_group)) +
            geom_boxplot() + geom_jitter(width=0.3) + theme_bw(base_size = 14) +  ggtitle(suffix) + 
            facet_wrap(.~get(var), scales = "free_x") + ylab(paste0(gene_set, pathway)) +
            stat_compare_means(comparisons=my_comparisons_signature, method="wilcox.test") +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)))
    dev.off()
    
    ### Compare hormone groups within signature groups
    # png(file=paste0(outPath, "Boxplot_",suffix,"_BCatSignature",n_BcatLevels,"Groups_", gene_set, pathway,".png"), width=1000)
    # print(ggplot(data=metadata_all, aes(x=get(var) , y=get(paste0(gene_set, pathway)))) +
    #         geom_boxplot() + geom_jitter(width=0.3) + theme_bw(base_size = 14) +
    #         ggtitle("Signature expression levels") + 
    #         facet_wrap(.~BcatSignatureMeanZScore_group) + ylab(paste0(gene_set, pathway)) +
    #         stat_compare_means(comparisons=my_comparisons_hormones, method="wilcox.test") +
    #         theme(axis.text.x = element_text(angle = 45, hjust = 1))) 
    # dev.off()
}



# Nb CTNNB1 mutated in cortisol vs non-functioning vs other --------
# metadata_all$Cortisol_None_NA_Other <- NA
# metadata_all[metadata_all$HormoneExpressionSimplified %in% c("None"),]$Cortisol_None_NA_Other <- "non-functioning"
# metadata_all[metadata_all$HormoneExpressionSimplified %in% c("Cortisol","Cortisol+Androgen","Mineralocorticoids+Cortisol"),]$Cortisol_None_NA_Other <- "cortisol"
# metadata_all[metadata_all$HormoneExpressionSimplified %in% c("Mineralocorticoids", "Androgen_Estrogen"),]$Cortisol_None_NA_Other <- "other"
# 
# png(paste0(outPath, "Proportion_CTNNB1mut_Cortisol_NonFunctioning_Other_NA.png"))
# ggplot(data=metadata_all, aes(y=BcatStatusGeneral)) + geom_bar() +
#   theme_bw() + facet_wrap(.~DatasetName+Cortisol_None_NA_Other)
# dev.off()


#  Boxplot all MHC molecules in signature high vs low based on cortisol ---------------------
# gene_set <- ""
# majorHistocompatibilityComplex=c("HLA-G","HLA-E","HLA-DRB1","HLA-DRA","HLA-DQA1","HLA-DQB1","HLA-DPA1","HLA-DPB1","HLA-DOA","HLA-DMA")
# pathway <- "majorHistocompatibilityComplex"
# genes <- majorHistocompatibilityComplex
# 
# for(gene in genes){
#   print(gene)
#   select <- df[which(rownames(df) %in% gene),]
#   select <- as.matrix(select)
#   select <- apply(select, 2, as.numeric)
#   select <- as.data.frame(cbind(names(select), select))
#   colnames(select) <- c("TumorID", gene)
#   metadata_all <- merge(metadata_all, select, by="TumorID", all.x=T, all.y=F)
#   metadata_all[,gene] <- as.numeric(as.vector(metadata_all[,gene]))
# }
# 
# hormones <- "cortisol"
# suffix <- "Cortisol"
# var<- "CortisolYesNo"
# metadata_all$CortisolYesNo <- factor(metadata_all$CortisolYesNo, levels=c("no","yes",NA))
# subset <- metadata_all[,c("TumorID",genes,"CortisolYesNo", "BcatSignatureMeanZScore_group")]
#   
# group_counts <- subset %>%
#     group_by(BcatSignatureMeanZScore_group, get(var)) %>%
#     summarise(n = n(), .groups = "drop")
#   
# print(paste0(gene_set, " ", pathway, ", ", hormones, ", n_BcatLevels=", n_BcatLevels, "--------------------"))
# print(group_counts)
#   
# subset_long <- tidyr::pivot_longer(subset, cols = 2:11, names_to = "Gene", values_to = "Expression")
# my_comparisons_signature <- list(c(levels(subset_long$BcatSignatureMeanZScore_group)[1],levels(subset_long$BcatSignatureMeanZScore_group)[2]))
# 
# png(paste0(outPath, "Boxplot_MHCs_bySignature_byCortisol.png"), width=1800, height = 900)
# print(ggplot(data=subset_long, aes(x=BcatSignatureMeanZScore_group , y=Expression,  fill=BcatSignatureMeanZScore_group)) +
#         geom_boxplot() + 
#         theme_bw(base_size = 14) +
#         ggtitle("Signature expression levels") + 
#         facet_wrap(.~get(var)+Gene, nrow=3, scales="free_x") +
#         stat_compare_means(comparisons=my_comparisons_signature, method="wilcox.test", label.y=5) +
#         theme(axis.text.x = element_text(angle = 45, hjust = 1))) 
# dev.off()



#  Boxplot all immune checkpoint molecules in signature high vs low based on cortisol ---------------------
# gene_set <- ""
# ImmuneCheckpoint <- c("CD274","CTLA4","HAVCR2","LAG3","PDCD1","PDCD1LG2","TIGIT","SIGLEC15")
# pathway <- "ImmuneCheckpoint"
# genes <- ImmuneCheckpoint
# 
# for(gene in genes){
#   print(gene)
#   select <- df[which(rownames(df) %in% gene),]
#   select <- as.matrix(select)
#   select <- apply(select, 2, as.numeric)
#   select <- as.data.frame(cbind(names(select), select))
#   colnames(select) <- c("TumorID", gene)
#   metadata_all <- merge(metadata_all, select, by="TumorID", all.x=T, all.y=F)
#   metadata_all[,gene] <- as.numeric(as.vector(metadata_all[,gene]))
# }
# 
# hormones <- "cortisol"
# suffix <- "Cortisol"
# var<- "CortisolYesNo"
# metadata_all$CortisolYesNo <- factor(metadata_all$CortisolYesNo, levels=c("no","yes",NA))
# subset <- metadata_all[,c("TumorID",genes,"CortisolYesNo", "BcatSignatureMeanZScore_group")]
# 
# group_counts <- subset %>%
#   group_by(BcatSignatureMeanZScore_group, get(var)) %>%
#   summarise(n = n(), .groups = "drop")
# 
# print(paste0(gene_set, " ", pathway, ", ", hormones, ", n_BcatLevels=", n_BcatLevels, "--------------------"))
# print(group_counts)
# 
# subset_long <- tidyr::pivot_longer(subset, cols = 2:9, names_to = "Gene", values_to = "Expression")
# my_comparisons_signature <- list(c(levels(subset_long$BcatSignatureMeanZScore_group)[1],levels(subset_long$BcatSignatureMeanZScore_group)[2]))
# 
# png(paste0(outPath, "Boxplot_ImmuneCheckpoint_bySignature_byCortisol.png"), width=1800, height = 900)
# print(ggplot(data=subset_long, aes(x=BcatSignatureMeanZScore_group , y=Expression,  fill=BcatSignatureMeanZScore_group)) +
#         geom_boxplot() + 
#         theme_bw(base_size = 14) +
#         ggtitle("Signature expression levels") + 
#         facet_wrap(.~get(var)+Gene, nrow=3, scales="free_x") +
#         stat_compare_means(comparisons=my_comparisons_signature, method="wilcox.test", label.y=3) +
#         theme(axis.text.x = element_text(angle = 45, hjust = 1))) 
# dev.off()
# 
# 
