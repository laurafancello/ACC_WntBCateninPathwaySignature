library(ggplot2)
library(ggpubr)
library(gridExtra)
library(ggExtra)
library(gridExtra)
library(cutpointr)
library(ggrepel)
library(survminer)
library(ggsurvfit)
source("./Scripts/Miscellaneous_Functions.R")

inpath <- "./ACC_datasets_integration/"
outpath <- paste0(inpath, "Survival/")
if (!(file.exists(outpath))){
  dir.create(file.path(outpath))
}

list_signatures <- list(ACC_BcatSignature = c("AFF3","AXIN2","FSCN1","LEF1","MYC","PDE2A","SP5","TBX3","TNFRSF19"), # Read signature genes
                        LEF1="LEF1")

data <- readRDS(file=paste0(inpath, "Counts_Zscore_allDatasets.RDS"))
metadata <- readRDS(file=paste0(inpath,"Metadata_allDatasets.RDS"))

samples_data <- colnames(data)[-which(colnames(data) == "GeneSymbol")]
metadata <- metadata[which(metadata$TumorID %in% samples_data),] # Keep only samples with metadata and expression data
metadata$OS_Status <- stringr::str_replace_all(metadata$OS_Status, "dead", "1")
metadata$OS_Status <- stringr::str_replace_all(metadata$OS_Status, "alive", "0")
metadata$OS_Status <- stringr::str_replace_all(metadata$OS_Status, "unknown", "NA")
metadata$OS_Status <- as.numeric(as.vector(metadata$OS_Status))
metadata <- metadata[!(is.na(metadata$OS_Months)),] # Keep only samples with overall survival data
genes <- data$GeneSymbol
data <- data[,which(colnames(data) %in% metadata$TumorID)]
data$GeneSymbol <- genes


# On merged dataset of Z-score transfoirmed expression values -------------
for(s in 1:length(list_signatures)){
  signature_name <- names(list_signatures)[s] 
  signature_genes <- list_signatures[[s]]
  signature <- data[data$GeneSymbol %in% signature_genes,-which(colnames(data) == "GeneSymbol")]
  if(length(signature_genes)>1){
    mean <- apply(as.matrix(signature), 2, mean)
    mean <- as.data.frame(cbind(names(mean),mean))
  }else{
    mean <- as.data.frame(t(rbind(colnames(signature),signature)))
  }
  colnames(mean) <- c("patientId","Signature")
  mean$Signature <- as.numeric(as.vector(mean$Signature))
  
  cutoff <- round(median(mean$Signature), digits=2)
  low <- mean[mean$Signature<cutoff,]
  high <- mean[mean$Signature>=cutoff,]
  
  low <- as.data.frame(cbind(low, rep("low", length(low$patientId))))
  colnames(low)[ncol(low)] <- "SignatureLevel"
  high <- as.data.frame(cbind(high, rep("high", length(high$patientId))))
  colnames(high)[ncol(high)] <- "SignatureLevel"
  SignatureLevels <- rbind(low, high)
  
  SignatureLevels$TumorID <- SignatureLevels$patientId
  
  input_surv <- merge(SignatureLevels, metadata, by="TumorID", all=F)
  input_surv$OS_Status <- as.numeric(as.vector(input_surv$OS_Status))
  input_surv$SignatureLevel <- factor(input_surv$SignatureLevel, levels=c("low","high"))
  
  fit <- survfit(Surv(input_surv$OS_Months, input_surv$OS_Status) ~
                   input_surv$SignatureLevel, data = input_surv)
  values <- surv_pvalue(fit)
  p <- signif(values$pval, digits=3)
  
  # PNG
  png(file = paste0(outpath, "Survival_curve_median", signature_name, "_publi.png"), width=600, height = 600)
  print(ggsurvplot(fit=fit,
                   linetype = c(1,3), censor.shape = 124, censor.size=2,
                   legend.title = "", font.x = 22, font.y = 22, title = paste0(signature_name, ", median\n p=", p),
                   font.tickslab = 18, font.legend = 18, pval.size = 0,
                   legend.labs = c("high","low"),
                   palette = c("red", "blue"), surv.scale = "percent",
                   ylab = "% of survival", xlab = "Time (in months)",
                   ggtheme = theme_juju(), risk.table = T, pval = T,
                   tables.theme = theme_cleantable_juju(), fontsize = 5.5,
                   risk.table.y.text.col = F))
  dev.off()
  # PDF
  pdf(file = paste0(outpath, "Survival_curve_median", signature_name, "_publi.pdf"), useDingbats = F)
  print(ggsurvplot(fit=fit,
                   linetype = c(1,3), censor.shape = 124, censor.size=2,
                   legend.title = "", font.x = 22, font.y = 22, title = paste0(signature_name, ", median\n p=", p),
                   font.tickslab = 18, font.legend = 18, pval.size = 0,
                   legend.labs = c("high","low"),
                   palette = c("red", "blue"), surv.scale = "percent",
                   ylab = "% of survival", xlab = "Time (in months)",
                   ggtheme = theme_juju(), risk.table = T, pval = T,
                   tables.theme = theme_cleantable_juju(), fontsize = 5.5,
                   risk.table.y.text.col = F))
  dev.off()
}


# On individual datasets of Z-score transfoirmed expression values -------------
for(dataset_name in c("Assie", "Heaton", "tcga", "JouinotFFPE", "Demeure")){
  
  genes <- data$GeneSymbol
  data_sub <- data[, grep(dataset_name, colnames(data))]
  data_sub$GeneSymbol <- genes
  
  samples_data <- colnames(data_sub)[-which(colnames(data_sub) == "GeneSymbol")]
  metadata_sub <- metadata[which(metadata$TumorID %in% samples_data),] # Keep only samples with metadata and expression data
  
  for(s in 1:length(list_signatures)){
    signature_name <- names(list_signatures)[s] 
    signature_genes <- list_signatures[[s]]
    signature <- data_sub[data_sub$GeneSymbol %in% signature_genes,-which(colnames(data_sub) == "GeneSymbol")]
    if(length(signature_genes)>1){
      mean <- apply(as.matrix(signature), 2, mean)
      mean <- as.data.frame(cbind(names(mean),mean))
    }else{
      mean <- as.data.frame(t(rbind(colnames(signature),signature)))
    }
    colnames(mean) <- c("patientId","Signature")
    mean$Signature <- as.numeric(as.vector(mean$Signature))
    
    cutoff <- round(median(mean$Signature), digits=2)
    low <- mean[mean$Signature<cutoff,]
    high <- mean[mean$Signature>=cutoff,]
    
    low <- as.data.frame(cbind(low, rep("low", length(low$patientId))))
    colnames(low)[ncol(low)] <- "SignatureLevel"
    high <- as.data.frame(cbind(high, rep("high", length(high$patientId))))
    colnames(high)[ncol(high)] <- "SignatureLevel"
    SignatureLevels <- rbind(low, high)
    
    SignatureLevels$TumorID <- SignatureLevels$patientId
    
    input_surv <- merge(SignatureLevels, metadata_sub, by="TumorID", all=F)
    input_surv$OS_Status <- as.numeric(as.vector(input_surv$OS_Status))
    input_surv$SignatureLevel <- factor(input_surv$SignatureLevel, levels=c("low","high"))
    
    fit <- survfit(Surv(input_surv$OS_Months, input_surv$OS_Status) ~
                     input_surv$SignatureLevel, data = input_surv)
    values <- surv_pvalue(fit)
    p <- signif(values$pval, digits=3)
    
    # PNG
    png(file = paste0(outpath, "Survival_curve_median", signature_name, "_", dataset_name, "_publi.png"), width=600, height = 600)
    print(ggsurvplot(fit=fit,
                     linetype = c(1,3), censor.shape = 124, censor.size=2,
                     legend.title = "", font.x = 22, font.y = 22, title = paste0(signature_name, ", ", dataset_name, ", median\n p=", p),
                     font.tickslab = 18, font.legend = 18, pval.size = 0,
                     legend.labs = c("high","low"),
                     palette = c("red", "blue"), surv.scale = "percent",
                     ylab = "% of survival", xlab = "Time (in months)",
                     ggtheme = theme_juju(), risk.table = T, pval = T,
                     tables.theme = theme_cleantable_juju(), fontsize = 5.5,
                     risk.table.y.text.col = F))
    dev.off()
    # PDF
    pdf(file = paste0(outpath, "Survival_curve_median", signature_name, "_", dataset_name, "_publi.pdf"), useDingbats = F)
    print(ggsurvplot(fit=fit,
                     linetype = c(1,3), censor.shape = 124, censor.size=2,
                     legend.title = "", font.x = 22, font.y = 22, title = paste0(signature_name, ", ", dataset_name, ", median\n p=", p),
                     font.tickslab = 18, font.legend = 18, pval.size = 0,
                     legend.labs = c("high","low"),
                     palette = c("red", "blue"), surv.scale = "percent",
                     ylab = "% of survival", xlab = "Time (in months)",
                     ggtheme = theme_juju(), risk.table = T, pval = T,
                     tables.theme = theme_cleantable_juju(), fontsize = 5.5,
                     risk.table.y.text.col = F))
    dev.off()
  }
}