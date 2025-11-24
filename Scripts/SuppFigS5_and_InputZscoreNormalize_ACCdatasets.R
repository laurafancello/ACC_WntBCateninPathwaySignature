library(stringr)
library(ggplot2)
library(magrittr)
library(dplyr)
library(msigdbr)
library(factoextra)
library(smplot2)
library(ggpubr)


# Read input --------------------------------------------------------------
inPath <- "./ACC_datasets_formatted/"
outPath <- "./ACC_datasets_integration/"

metadata_all <- as.data.frame(matrix(ncol=14))
colnames(metadata_all) <- c("TumorID","Age","Sex","OS_Status","OS_Months","WeissScore","BcatStatusGeneral","TP53Alteration","DatasetName","Technology","MicroarrayPlatform","BcatSignatureMean","HormoneExpression","CortisolSecretion")
dataset_all <- as.data.frame(matrix(nrow=1))
colnames(dataset_all) <- "GeneSymbol"

signature <- c("AFF3","AXIN2","FSCN1","LEF1","MYC","PDE2A","SP5","TBX3","TNFRSF19")

norm <- "Zscore"

#  Generate dataset and metadata objects containing all datasets ----------
for(dataset_name in c("Caramuta","JouinotFFPE","Demeure","Assie", "Heaton", "tcga")){
  print(dataset_name)
  
  # Read metadata and counts
  metadata <- read.csv(file=paste0(inPath, dataset_name, "/Metadata/Metadata_", dataset_name,".csv"))
  
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
  
  # Select only ACC samples
  ACCsamples <- metadata[metadata$Histotype == "ACC",]$TumorID
  ACC <- as.matrix(dataset_nr[,which(colnames(dataset_nr) %in% ACCsamples)])
  rownames(ACC) <- dataset_nr$GeneSymbol
  metadata <- metadata[metadata$Histotype == "ACC",]
  
  # Calculate signature mean expression
  select <- ACC[which(rownames(ACC) %in% signature),]
  select <- as.matrix(select)
  select <- apply(select, 2, as.numeric)
  
  mean <- apply(select, 2, mean, na.rm=T)
  mean <- as.data.frame(cbind(names(mean), mean))
  colnames(mean) <- c("TumorID","BcatSignatureMean")
  metadata <- merge(metadata, mean, by="TumorID", all=T)
  
  # Apply normalization
  # Need to remove genes with st dev = 0 (for which Z-scores is NAs) 
  check_sd <- apply(ACC, 1, sd)
  toRemove <- which(check_sd == 0)
  if(length(toRemove)>0){ # it will be always 0 if only most variable genes selected. Otherwise it depends on the dataset
    ACC <- ACC[-(toRemove),]
    geneNames <- rownames(ACC)
  }else{
    geneNames <- rownames(ACC)
  }
  # Z-score normalization on a per gene basis
  ACC_Z <- as.data.frame(t(scale(t(ACC))))
  colnames(ACC_Z) <- paste0(dataset_name, "_", colnames(ACC_Z))
  ACC_Z$GeneSymbol <- geneNames
  
  # Select only metadata columns in common to all datasets and format values with same terminology
  if(dataset_name == "Heaton"){
    metadata$BcatStatusGeneral <- metadata$BetaCateninStaining
    metadata$BcatStatusGeneral <- stringr::str_replace(metadata$BcatStatusGeneral, "Nuclear", "Active")
    metadata$BcatStatusGeneral <- stringr::str_replace(metadata$BcatStatusGeneral, "Membrane", "Inactive")
  }
  
  if(dataset_name == "Assie"){
    metadata$BcatStatusGeneral <- metadata$CTNNB1_ZNRF3_Alteration 
    metadata$BcatStatusGeneral <- stringr::str_replace_all(metadata$BcatStatusGeneral, "CTNNB1", "Active")
    metadata$BcatStatusGeneral <- stringr::str_replace_all(metadata$BcatStatusGeneral, "wt", "Inactive")
    metadata$BcatStatusGeneral <- stringr::str_replace_all(metadata$BcatStatusGeneral, "ZNRF3", "ZNRF3_mut")
    
    metadata$TP53Alteration <- metadata$p53.pathway
  }
  
  if(dataset_name == "tcga"){
    metadata$BcatStatusGeneral <- metadata$BcatStatus
    metadata$BcatStatusGeneral <- stringr::str_replace(metadata$BcatStatusGeneral, "CTNNB1_mut", "Active")
    metadata$BcatStatusGeneral <- stringr::str_replace(metadata$BcatStatusGeneral, "wt", "Inactive")
    
    metadata$TP53Alteration <- metadata$TP53Status
    metadata$TP53Alteration <- stringr::str_replace_all(metadata$TP53Alteration, "0", "No")
    metadata$TP53Alteration <- stringr::str_replace_all(metadata$TP53Alteration, "1", "Yes")
  }
  
  # Add NAs if metadata variable not existing for a dataset 
  for(variable in c("Age","Sex","OS_Status","OS_Months","WeissScore","TP53Alteration","BcatStatusGeneral","HormoneExpression","CortisolSecretion")){
    if(length(which(colnames(metadata) == variable)) == 0){
      metadata[,variable] <- NA
    }
  }
  metadata$DatasetName <- dataset_name
  if(dataset_name %in% c("JouinotFFPE","tcga")){
    metadata$Technology <- "RNAseq"
    metadata$MicroarrayPlatform <- "no"
  }else{
    metadata$Technology <- "microarray"
  }
  metadata <- metadata[,c("TumorID","Age","Sex","OS_Status","OS_Months","WeissScore","BcatStatusGeneral","TP53Alteration","DatasetName","Technology","MicroarrayPlatform","BcatSignatureMean","HormoneExpression","CortisolSecretion")]
  metadata$TumorID <- paste0(dataset_name, "_", metadata$TumorID)
  
  metadata_all <- rbind(metadata_all, metadata)
  metadata_all$BcatSignatureMean <- as.numeric(as.vector(metadata_all$BcatSignatureMean))
  
  # Remove genes with name equal to "" or NA, due to Ensembl to Gene Symbol conversion
  if(length(which(ACC_Z$GeneSymbol==""))>0){
    ACC_Z <- ACC_Z[-which(ACC_Z$GeneSymbol==""),]
  }
  ACC_Z <- ACC_Z[!(is.na(ACC_Z$GeneSymbol)),]
  
  if(nrow(dataset_all)==1){
    dataset_all <- ACC_Z
  }else{
    dataset_all <- dplyr::inner_join(dataset_all, ACC_Z, by="GeneSymbol")
  }
  print(paste0("metadata:",dim(metadata)))
  print(paste0("dataset:",dim(ACC_Z)))
  print(paste0("dataset_all:",dim(dataset_all)))
}
metadata_all <- metadata_all[-1,]

### Simplify hormone expression representation
metadata_all$HormoneExpressionSimplified <- metadata_all$HormoneExpression
metadata_all$HormoneExpressionSimplified <- stringr::str_replace_all(metadata_all$HormoneExpressionSimplified, "^Androgen$", "Androgen_Estrogen")
metadata_all$HormoneExpressionSimplified <- stringr::str_replace_all(metadata_all$HormoneExpressionSimplified, "DHEAS", "Androgen_Estrogen")
metadata_all$HormoneExpressionSimplified <- stringr::str_replace_all(metadata_all$HormoneExpressionSimplified, "^Estrogen$", "Androgen_Estrogen")
metadata_all[metadata_all$HormoneExpressionSimplified %in% c("Androgen+Estrogen"),]$HormoneExpressionSimplified <- "Androgen_Estrogen"

metadata_all[(metadata_all$DatasetName=="JouinotFFPE")&(metadata_all$CortisolSecretion %in% c("yes")),]$HormoneExpressionSimplified <- "Cortisol"
metadata_all[(metadata_all$DatasetName=="JouinotFFPE")&(metadata_all$CortisolSecretion %in% c("no")),]$HormoneExpressionSimplified <- "NA"
metadata_all[(metadata_all$DatasetName=="JouinotFFPE")&(metadata_all$CortisolSecretion %in% c("yes")),]$HormoneExpression <- "Cortisol"
metadata_all[(metadata_all$DatasetName=="JouinotFFPE")&(metadata_all$CortisolSecretion %in% c("no")),]$HormoneExpression <- "NA"

metadata_all[(metadata_all$DatasetName=="Assie")&(metadata_all$HormoneExpression %in% c("yes")),]$HormoneExpressionSimplified <- "NA"
metadata_all[(metadata_all$DatasetName=="Assie")&(metadata_all$HormoneExpression %in% c("None")),]$HormoneExpressionSimplified <- "None"

#### HormoneExpresionYesNo
metadata_all$HormoneExpressionYesNo <- ifelse(metadata_all$HormoneExpression == "None", "no", "yes")
metadata_all[is.na(metadata_all$HormoneExpression),]$HormoneExpressionYesNo <- "NA"

### Get column only for cortisol secretion yes/no
metadata_all$CortisolYesNo <- "NA"
metadata_all[(metadata_all$DatasetName %in% c("Caramuta","Demeure","Assie")),]$CortisolYesNo <- "NA"
           
metadata_all[((metadata_all$DatasetName == "JouinotFFPE")&(metadata_all$CortisolSecretion %in% c("yes"))),]$CortisolYesNo <- "yes"
metadata_all[((metadata_all$DatasetName == "JouinotFFPE")&(metadata_all$CortisolSecretion %in% c("no"))),]$CortisolYesNo <- "no"

metadata_all[((metadata_all$DatasetName %in% c("Heaton","tcga"))&(metadata_all$HormoneExpressionSimplified %in% c("Cortisol","Cortisol+Androgen","Mineralocorticoids+Cortisol"))),]$CortisolYesNo <- "yes"
metadata_all[((metadata_all$DatasetName %in% c("Heaton","tcga"))&(metadata_all$HormoneExpressionSimplified %in% c("None","Androgen_Estrogen","Mineralocorticoids","Androgen+Mineralocorticoids"))),]$CortisolYesNo <- "no"
metadata_all[((metadata_all$DatasetName %in% c("Heaton","tcga"))&(is.na(metadata_all$HormoneExpressionSimplified))),]$CortisolYesNo <- "NA"


# Remove duplication TCGA sample ------------------------------------------
# I realize that in metadata a tumor is reported twice, thus I remove it (it's "tcga_TCGA-OR-A5JA")
metadata_all <- metadata_all[-which(duplicated(metadata_all$TumorID)),]

## B-catenin activation signature
df <- dataset_all[,-which(colnames(dataset_all) == "GeneSymbol")]
rownames(df) <- dataset_all$GeneSymbol
select <- df[which(rownames(df) %in% signature),]
select <- as.matrix(select)
select <- apply(select, 2, as.numeric)
mean <- apply(select, 2, mean, na.rm=T)
mean <- as.data.frame(cbind(names(mean), mean))
colnames(mean) <- c("TumorID","BcatSignatureMeanZScore")
metadata_all <- merge(metadata_all, mean, by="TumorID", all.x=T, all.y=F)
metadata_all$BcatSignatureMeanZScore <- as.numeric(as.vector(metadata_all$BcatSignatureMeanZScore))

## LEF1
select <- df[which(rownames(df) %in% c("LEF1")),]
select <- as.numeric(as.vector(select))
select <- as.data.frame(cbind(colnames(df), select))
colnames(select) <- c("TumorID","LEF1ZScore")
metadata_all <- merge(metadata_all, select, by="TumorID", all.x=T, all.y=F)
metadata_all$LEF1ZScore <- as.numeric(as.vector(metadata_all$LEF1ZScore))

saveRDS(file=paste0(outPath, "Metadata_allDatasets.RDS"), metadata_all)
saveRDS(file=paste0(outPath, "Counts_", norm, "_allDatasets.RDS"), dataset_all)

# Calculate B-cat signature or MYC targets mean expression on Z-score normalized values
metadata_all <- readRDS(file=paste0(outPath, "Metadata_allDatasets.RDS"))
dataset_all <- readRDS(file=paste0(outPath, "Counts_", norm, "_allDatasets.RDS"))


# PCA on z-score normalized counts all datasets together------------------------------------------------------------
metadata_all <- readRDS(file=paste0(outPath, "Metadata_allDatasets.RDS"))
dataset_all <- readRDS(file=paste0(outPath, "Counts_", norm, "_allDatasets.RDS"))

df <- dataset_all[,-which(colnames(dataset_all) == "GeneSymbol")]
rownames(df) <- dataset_all$GeneSymbol
meta_df <- metadata_all[match(colnames(df), metadata_all$TumorID), ]

meta_df <- as.data.frame(t(meta_df))
colnames(meta_df) <- meta_df[1,]
meta_df <- meta_df[-1,]
identical(colnames(df), colnames(meta_df))

df <- rbind(df, meta_df)
df <-as.data.frame(t(df))
df$BcatSignatureMeanZScore <- as.numeric(as.vector(df$BcatSignatureMeanZScore))
df$BcatSignatureMeanZScore_Discrete <- "medium"
df[as.numeric(as.vector(df$BcatSignatureMeanZScore)) <= quantile(as.numeric(as.vector(df$BcatSignatureMeanZScore)))[2],]$BcatSignatureMeanZScore_Discrete <- "low"
df[as.numeric(as.vector(df$BcatSignatureMeanZScore)) >= quantile(as.numeric(as.vector(df$BcatSignatureMeanZScore)))[4],]$BcatSignatureMeanZScore_Discrete <- "high"

res.pca <- prcomp(as.data.frame(apply(as.matrix(df[,1:12334]),2,as.numeric), scale=F))

### Plot percent variance explained by each principal component
png(paste0(outPath, "PCAs/PercVar_eachPC_",norm,"_integration.png"))
fviz_screeplot(res.pca, addlabels = TRUE, ylim = c(0, 50))
dev.off()
summary(res.pca)

### PNGs of PCAs colored by different categorical variables
for(var in c("BcatStatusGeneral","BcatSignatureMeanZScore_Discrete","Sex","DatasetName","Technology","MicroarrayPlatform")){
  png(paste0(outPath, "PCAs/Improved_PCAs_",norm,"_integration_",var,".png"))
  print(var)
  print(fviz_pca_ind(res.pca,
                     label = "none", # hide individual labels
                     habillage = df[,var], # color by groups
                     #palette = c("#00AFBB", "#E7B800", "#FC4E07"),
                     addEllipses = TRUE, # Concentration ellipses
                     title=var))
  dev.off()
}

### PNGs of PCAs colored by different BcatSignatureMeanZScore
color_var <- df$BcatSignatureMeanZScore
png(paste0(outPath, "PCAs/Improved_PCAs_",norm,"_integration_BcatSignatureMeanZScore.png"))
fviz_pca_ind(
  res.pca,
  geom.ind = "point",
  col.ind = color_var,              # color by the continuous variable
  gradient.cols = c("green4", "yellow", "red"),  # color gradient
  addEllipses = FALSE,              # no group ellipses 
  legend.title = "Signature\nMean Z-Score"
)
dev.off()
