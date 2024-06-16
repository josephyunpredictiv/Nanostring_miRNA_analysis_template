#v12: 2nd set
# Analyses

# HepG2 liver cell line
# - Analysis13_HepG2_HFD_vs_CD
# - Analysis14_HepG2_HFD_vs_C
# - Analysis15_HepG2_CD_vs_C

# HIEC-6 small intestinal cell line
# - Analysis16_HIEC6_HFD_vs_CD
# - Analysis17_HIEC6_HFD_vs_C
# - Analysis18_HIEC6_CD_vs_C

#v10: Based on the template code of Delirium_miRNA_integrated_20220909_v9.R
#     RandomForest integrated
#     Normalization by spikein control

rm(list=ls())

#################
#MODULES TO LOAD#
#################
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DESeq2")
#BiocManager::install("IHW")
library("DESeq2")
library(foreign)
#install.packages("caret")
library(caret)
library(Rtsne)
library("IHW")
library(plotly)
#install.packages("pheatmap")
library("pheatmap")
library(dplyr)
library(tidyr)
library(pheatmap)
#install.packages("randomForest")
library(randomForest)
#options(digits=2, scipen = 999)

###########
#Functions#
###########
get_summary_v2 <- function(resOrdered, adjusted_p_value_cutoff, title) {
  summary_result_file_name = paste0("summary_results_with_adjusted_p_value_",adjusted_p_value_cutoff,"_", title, ".txt")
  summary(resOrdered, alpha = adjusted_p_value_cutoff)
  sink(file=summary_result_file_name)
  summary(resOrdered, alpha = adjusted_p_value_cutoff)
  sink()
}


##########################################
#Processing commonly used master datasets#
##########################################
#Read the normalized raw data from Nanostring
setwd("/Users/josephyun/Desktop/Data_analysis_for_adipocyte_induced_exosome_project_2023_04")
dataset <- read.table("Normalized_housekeeping_cleaned_v2.csv", header=T, sep=",", row.names = 1, check.names = F, skip = 2)
dataset <- dataset[,7:ncol(dataset)]
View(dataset)

#dataset_previous <- read.table("Normalized_housekeeping_cleaned_v1.csv", header=T, sep=",", row.names = 1, check.names = F, skip = 2)
#dataset_previous <- dataset_previous[,7:ncol(dataset_previous)]
#View(dataset_previous)


dataset2 <- dataset
#Multiply by 10 (later, we will divide by 10) and make integer
for (j in 1:ncol(dataset)) {
  for (i in 1:nrow(dataset)) {
    dataset2[i,j] <- as.integer(dataset[i,j]*10)
  }
}
#Remove ^NEG & ^POS
dataset2 <- dataset2[-grep("NEG",rownames(dataset2)),]
dataset2 <- dataset2[-grep("POS",rownames(dataset2)),]
#Rowname (genename: Kir3dl1/2 => Kir3dl1_2)
rownames(dataset2) <- sub('/','_',rownames(dataset2))

#Prepare coldata, which contains sample_ID, batch, and condition
sample_ID <- colnames(dataset) #HepG2_C1, HepG2_C2
cell_line <- sub('_.*','',sample_ID) #HepG2, HepG2
sample_type_replicate_number <- sub('.*_','',sample_ID) #C1, C2,
condition <- substring(sample_type_replicate_number,1, nchar(sample_type_replicate_number)-1) #C, C
replicate_number <- substr(sample_type_replicate_number,(nchar(sample_type_replicate_number)+1)-1,nchar(sample_type_replicate_number)) #1, 2
coldata <- data.frame(sample_ID, cell_line, sample_type_replicate_number, replicate_number, condition)
View(coldata)

############################
#Analysis13_HepG2_HFD_vs_CD#
############################
name="Analysis13_HepG2_HFD_vs_CD"
this_cell_line = "HepG2"
condition1 = "HFD"
condition2 = "CD"

title=name
directory=paste0("/Users/josephyun/Desktop/Data_analysis_for_adipocyte_induced_exosome_project_2023_04/",name)
if (!(file.exists(directory)))
{
  dir.create(directory)
}
setwd(directory)
getwd()

#Get coldata2
coldata2 <- filter(coldata, cell_line == this_cell_line) %>% filter((condition == condition1)|(condition == condition2))
#Get dataset3
dataset3 <- dataset2 %>% select(coldata2$sample_ID)
View(dataset3)

#############################
#Prepare ddseq2 object (dds)#
#############################
ddsHTSeq <- DESeqDataSetFromMatrix(countData = dataset3, colData = coldata2, design = ~ condition)
dds <- ddsHTSeq[rowSums(counts(ddsHTSeq)) > 1, ]
dds$condition <- relevel(dds$condition, ref="CD") #For fold_change calculation, reference is CD
dds <- DESeq(dds)

#############################################################################################
save(dds, file = "dds.rds")
#############################################################################################
#dds <- load("dds.rds")
#dds <- updateObject(dds)
#############################################################################################


#############
#Get results#
#############
#Get results with adjusted p-value = 0.05
res <- results(dds, filterFun=ihw, alpha=0.05)
resOrdered <- res[order(res$padj),]
#Get results with adjusted p-value = 0.1
res01 <- results(dds, filterFun=ihw, alpha=0.1)
resOrdered01 <- res01[order(res$padj),]

#Get the summary
get_summary_v2(resOrdered, 0.05, title)
get_summary_v2(resOrdered01, 0.1, title)

#Get the CSV
#1: Fold_change_and_p_values_all_genes.csv
my_filename = paste0("results_Fold_change_and_p_values_", title, "_all_genes.csv")
write.csv(as.data.frame(resOrdered), file=my_filename)
#2: Fold_change_and_p_values_significant_genes.csv
resSig <- subset(resOrdered, padj < 0.05)
resSigOrdered <- resSig[order(resSig$padj),]
my_filename = paste0("results_Fold_change_and_p_values_", title, "_significant_genes_by_pvalue_0.05.csv")
write.csv(as.data.frame(resSigOrdered), file=my_filename)
resSig01 <- subset(resOrdered, padj < 0.1)
resSigOrdered01 <- resSig01[order(resSig01$padj),]
my_filename = paste0("results_Fold_change_and_p_values_", title, "_significant_genes_by_pvalue_0.1.csv")
write.csv(as.data.frame(resSigOrdered01), file=my_filename)

#3. Getting csv of normalized expression values for all samples for all genes
my_filename <- paste0("Normalized_data_by_DESEQ2_", title, ".csv")
write.csv(as.data.frame(assays(dds)[["mu"]])/10.0, file=my_filename)  #Since we multipled by 3761, here we divide by 3500 for count
#4. Getting csv of raw counts for all samples for all genes
my_filename <- paste0("Raw_counts_data_by_DESEQ2_", title, ".csv")
write.csv(as.data.frame(assay(dds[,]))/10.0, file=my_filename)  #Since we multipled by 3761, here we divide by 3500 for count

#############
#Count plots#
#############
#Count plot for genes padj < 0.05
genes=rownames(res[!is.na(res$padj) & res$padj<0.05,])
if (length(genes)>=1)
{
  if (!(file.exists("count_plots_padj_0.05")))
  {
    dir.create("count_plots_padj_0.05")
  }
  for (gene_name in genes) {
    my_filename = paste0("count_plots_padj_0.05/count_plot_", title, "_", gene_name, ".png")
    png(filename = my_filename)
    d <- plotCounts(dds, gene_name, intgroup="condition", returnData=TRUE)
    
    d2 <- ggplot(d, aes(x=condition, y=count/10.0)) + geom_point(position=position_jitter(w=0.1,h=0)) + scale_y_log10(breaks=c(1,as.integer(mean(d$count)/10.0),as.integer(max(d$count)/10.0)*1.1), limits=c(-1,as.integer(max(d$count)/10.0)*1.1+1))
    d2$layers <- c(geom_boxplot(), d2$layers)
    d3 <- d2 + xlab(gene_name) + ylab("Normalized Expression")
    d3 <- d3 + theme(text = element_text(size = 30))
    print(d3)
    
    dev.off()
  }
}

#Count plot for genes padj < 0.1
genes=rownames(res[!is.na(res$padj) & res$padj<0.1,])
if (length(genes)>=1)
{
  if (!(file.exists("count_plots_padj_0.1")))
  {
    dir.create("count_plots_padj_0.1")
  }
  for (gene_name in genes) {
    my_filename = paste0("count_plots_padj_0.1/count_plot_", title, "_", gene_name, ".png")
    png(filename = my_filename)
    d <- plotCounts(dds, gene_name, intgroup="condition", returnData=TRUE)
    
    d2 <- ggplot(d, aes(x=condition, y=count/10.0)) + geom_point(position=position_jitter(w=0.1,h=0)) + scale_y_log10(breaks=c(1,as.integer(mean(d$count)/10.0),as.integer(max(d$count)/10.0)*1.1), limits=c(-1,as.integer(max(d$count)/10.0)*1.1+1))
    d2$layers <- c(geom_boxplot(), d2$layers)
    d3 <- d2 + xlab(gene_name) + ylab("Normalized Expression")
    d3 <- d3 + theme(text = element_text(size = 30))
    print(d3)
    
    dev.off()
  }
}


#Count plot for all genes
genes=rownames(res[!is.na(res$padj),])
if (length(genes)>=1)
{
  if (!(file.exists("count_plots_all_genes")))
  {
    dir.create("count_plots_all_genes")
  }
  for (gene_name in genes) {
    my_filename = paste0("count_plots_all_genes/count_plot_", title, "_", gene_name, ".png")
    png(filename = my_filename)
    d <- plotCounts(dds, gene_name, intgroup="condition", returnData=TRUE)
    
    d2 <- ggplot(d, aes(x=condition, y=count/10.0)) + geom_point(position=position_jitter(w=0.1,h=0)) + scale_y_log10(breaks=c(1,as.integer(mean(d$count)/10.0),as.integer(max(d$count)/10.0)*1.1), limits=c(-1,as.integer(max(d$count)/10.0)*1.1+1))
    d2$layers <- c(geom_boxplot(), d2$layers)
    d3 <- d2 + xlab(gene_name) + ylab("Normalized Expression")
    d3 <- d3 + theme(text = element_text(size = 30))
    print(d3)
    
    dev.off()
  }
}


###########
#Heat maps#
###########
rld <- rlog(dds) 

#Heatmap for all significant genes, padj < 0.05
genes <- rownames(res[!is.na(res$padj) & res$padj<0.05,])
my_filename <- "heatmap_padj_0.05.png"
#Draw a heatmap
mat <- assay(rld)[genes,]
mat <- mat - rowMeans(mat)
#df <- as.data.frame(colData(rld)[,c("condition")])
png(filename = my_filename, units="in", width=4, height=6, res=300)
pheatmap(mat, show_rownames = T, show_colnames = T, cluster_cols = F, cluster_rows = T, treeheight_row = 0)
dev.off()

#Heatmap for all significant genes, padj < 0.1
genes <- rownames(res[!is.na(res$padj) & res$padj<0.1,])
my_filename <- "heatmap_padj_0.1.png"
#Draw a heatmap
mat <- assay(rld)[genes,]
mat <- mat - rowMeans(mat)
#df <- as.data.frame(colData(rld)[,c("condition")])
png(filename = my_filename, units="in", width=6, height=8, res=300)
pheatmap(mat, show_rownames = T, show_colnames = T, cluster_cols = F, cluster_rows = T, treeheight_row = 0)
dev.off()


##########
#PCA plot#
##########
rld <- varianceStabilizingTransformation(dds, blind=FALSE)
my_filename = paste0("PCA_plot_", title, ".png")
png(filename = my_filename, units="in", width=11, height=11, res=300)
z <- plotPCA(rld, intgroup="condition")
nudge <- position_nudge(y = 0.4)
print(z + geom_text(aes(label = name), position = nudge, size=5))  #ggplot2 requires print() if it is within a function
dev.off()

PCA_data <- plotPCA(rld, intgroup="condition", returnData = TRUE)
PCA_data <- PCA_data %>% arrange(PC1) 
write.table(PCA_data, file="PCA_results.csv",quote=F, sep=",")

#############################################################################################
save(dds, file = "dds.rds")
#############################################################################################
#dds <- load("dds.rds")
#dds <- updateObject(dds)
#############################################################################################
#############################################################################################
#############################################################################################

################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
###########################################
#Will not be delivered due to limited data#
###########################################
################################
#Prioritization by RandomForest#
################################
set.seed(1)

#Prepare a master dataframe with all the columns by merging resOrdered (padj, FC ...) and dataset3 (expression count)
#resOrdered is a DESEQ2 object
resOrdered_df <- as.data.frame(resOrdered)
#rowname as 'gene' column
resOrdered_df <- resOrdered_df %>% mutate(gene = row.names(resOrdered_df))
#Replace '-' as '_' for resOrdered_df
resOrdered_df$gene <- sub('-','_',resOrdered_df$gene)
dataset4 <- dataset3 %>% mutate(gene = row.names(dataset3))
#Replace '-' as '_' for dataset4
dataset4$gene <- sub('-', '_', dataset4$gene)

master_df <- inner_join(resOrdered_df, dataset4, by="gene")
#Create significance column 
master_df <- master_df %>% mutate(significance = case_when(padj < 0.05 ~ 1, padj >= 0.05 ~ 0))
#gene as rownames
row.names(master_df) <- master_df$gene
View(master_df)
transposed_master_df <- master_df %>% select(-gene) %>% t()
View(transposed_master_df)
#Run randomForest and build randomForest_model on master_df
randomForest_model <- randomForest(significance ~ . -weight -gene, data = master_df) #ntree=2000
print(randomForest_model)
summary(randomForest_model)
plot(randomForest_model)
getTree(randomForest_model,k = 1, labelVar = T)
importance(randomForest_model)
results_on_importance <- data.frame(importance(randomForest_model))
#View(results_on_importance)
varImpPlot(randomForest_model)

#Run randomForest and build randomForest_model on transposed_master_df against the most significant gene "Tfrc"
randomForest_model <- randomForest(Tfrc ~ . , data = transposed_master_df) #ntree=2000
print(randomForest_model)
summary(randomForest_model)
plot(randomForest_model)
getTree(randomForest_model,k = 1, labelVar = T)
importance(randomForest_model)
results_on_importance <- data.frame(importance(randomForest_model))
View(results_on_importance)
png(file="testRF.png")
varImpPlot(randomForest_model)
dev.off()

################################################################################################################################################
################################################################################################################################################
################################################################################################################################################