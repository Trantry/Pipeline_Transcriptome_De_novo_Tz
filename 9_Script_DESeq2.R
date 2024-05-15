#####################################################################
#####################################################################
########################## DE-seq 2  ################################
#####################################################################
#####################################################################
# DESeq() function
# Perform the median of ratios method of normalization
# Compare genes level between sample
# Create new columns : log2(FC), pval, padj...
# output : these normalized count are usefull for visualization of results
# output : do not use as input for DE analyses based on negative binomial model

BiocManager::install("DESeq2")
BiocManager::install("tximport")
install.packages("readr")
library(readr)
library(DESeq2)
library(tximport)
citation()
citation("DESeq2")
citation("tximport")
citation("readr")
citation("apeglm")
citation("pheatmap")
citation("ggplot2")
citation("dplyr")

setwd("/Users/trystan.bessonnier/Desktop/DESeq2_Transcriptome")
##########################################
## Import data & sample preparation ######
##########################################
#préparation des données:
#Création de 3 fichiers: file_list.txt, samples.txt, sampleTable.txt
files<-as.character(read.table("file_list.txt", header=FALSE)$V1)

samples<-as.character(read.table("samples.txt",header=FALSE)$V1)

names(files)<-samples

sampleTable <-read.table("sampleTable.txt",header=TRUE, row.names=1)

#IMPORTANT: MAKE SURE THE SAMPLES AND FILE_LIST ARE IN THE SAME ORDER- SAMPLES SHOULD MATCH UP WITH FILES
samples==rownames(sampleTable) #should return TRUE for all

#Import a file called tx2gene.csv which a csv file that contains the transcript id to gene id mapping
#For Tetraripis, this is located at: tzet_transcriptome_trinityV2_kraken2.Trinity.fasta.gene_trans_map
tx2gene <- read.csv("/Users/trystan.bessonnier/Desktop/DESeq2_Transcriptome/trinity_tzet_complete_transcriptome.Trinity.fasta.gene_trans_map",header=F, sep = "")
colnames(tx2gene) <- c("gene", "transcript")
tx2gene <- tx2gene[, c("transcript", "gene")]

#look at this data structure
#read in kallisto abundance files, summarizing by gene
txi <- tximport(files, type = "kallisto", tx2gene = tx2gene)
names(txi)
#txi <- na.omit(txi) 
#txi

#Optionally, if you want to save this count matrix as a file
#write.csv(assay(ddsMatrix), file="genecounts.raw.csv", quote=FALSE)

#Optionally, if you want to save normalized counts matrix as a file
#ddsMatrix<- estimateSizeFactors(ddsMatrix)
#normalized_counts <- counts(ddsMatrix, normalized=TRUE)
#write.csv(normalized_counts, file="genecounts.normalized.csv", quote=FALSE)

#Optionally, if you want to do variance stabilizing transformation or regularized log transformation on this count matrix and then save as a file: This can become input to things like wgcna, pca
#vsd<- vst(ddsMatrix, blind=TRUE)
#write.csv(assay(vsd), file="genecounts.variancestabilized.csv", quote=FALSE)

#####################################################################
#####################################################################
###################### Quality assessment  ##########################
#####################################################################
#####################################################################

#make a deseq2 object from the kallisto summarized counts
ddsMatrix <- DESeqDataSetFromTximport(txi, sampleTable, 
                                      design= ~condition) #class: DESeqDataSet dim: 476176 9 
#ddsMatrix <- ddsMatrix[rowSums(counts(ddsMatrix)) != 0, ] #class: DESeqDataSet dim: 407829 9 
#ddsMatrix<- na.omit(ddsMatrix) 

ddsMatrix<- estimateSizeFactors(ddsMatrix)
# Calculer le count moyen pour chaque gène
meanCounts <- rowMeans(counts(ddsMatrix))

# Filtrer les gènes qui ont un count moyen supérieur à 10
ddsMatrix <- ddsMatrix[meanCounts > 5, ]


dds<-DESeq(ddsMatrix)
dds #class: DESeqDataSet dim: 476176 9  
dds$condition <- factor(dds$condition) 
dds$condition <- relevel(dds$condition, ref = "A")
#dds <- DESeq(dds)

# this gives log2(n + 1)
ntd <- normTransform(dds)
#Extracting transformed values
# regularized lof transform (rlog): shrunk values for genes with low counts (because they show the largest absolute difference between samples)
# variance-stabilizing transformation (vst): faster and similar properties than rlog
# vsd transformation of normalized counts is only necessary for visualization during quality assessment
vsd <- vst(dds, blind=TRUE)
#Optionally, if you want to do variance stabilizing transformation or regularized log transformation on this count matrix and then save as a file: This can become input to things like wgcna, pca
#vsd<- vst(ddsMatrix, blind=TRUE)
rld <- rlog(dds, blind=TRUE)    #blind=TRUE : do the transformation in an unbiased manner
head(assay(vsd), 50)

BiocManager::install("vsn")
library("vsn")

meanSdPlot(assay(ntd)) # this gives log2(n + 1)
meanSdPlot(assay(vsd)) # variance-stabilizing transformation (vst): faster and similar properties than rlog
meanSdPlot(assay(rld)) # regularized lof transform (rlog): shrunk values for genes with low counts (because they show the largest absolute difference between samples)



#check the size factors for each samples
normalizationFactors(dds)
#total number of raw counts per sample
colSums(counts(dds))

#check the dispersion plot
plotDispEsts(dds)
#ylim=c(1e-06,1e+02)

#specify contrast? : if one factor has 3 levels; indicate which two levels classes we are interested in comparing
#if not specifying, results() returns the comparison of the last level of last variable in Design formula over the 1st level of this variable (ordered alphabetically if levels are not specified)
resultsNames(dds)

res <- results(dds, contrast=c("condition","F","A"), alpha = 0.05) #same sense
res #DataFrame with 476176 rows and 6 columns
#res <- subset(res, baseMean != 0) #DataFrame with 407829 rows and 6 columns
#res <- na.omit(res) #DataFrame with 323768 rows and 6 columns


#check res data
library(dplyr)
res %>%
  data.frame() %>%
  View()

#Shrinkage of log2 Fold Change
#It won't change the total number of genes significantly DE
#It helps with downstream assessment of results : to subset significant genes based on LFC for further evaluation
BiocManager::install("apeglm")
resultsNames(dds)
resLFC <- lfcShrink(dds, coef="condition_F_vs_A", type="apeglm")
resLFC #DataFrame with 407829 rows and 5 columns
#resLFC<- na.omit(resLFC) #DataFrame with 31646 rows and 5 columns

#order our results table by the smallest p value
resOrdered <- res[order(res$pvalue),]
summary(res)
#out of 323768 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)       : 49, 0.015%
#LFC < 0 (down)     : 0, 0%
#outliers [1]       : 0, 0%
#low counts [2]     : 0, 0%
#(mean count < 0)

#change the adjusted pval cutoff (automatically = 0.1)
#res10 <- results(dds, alpha=0.1)
#summary(res10)
#out of 120335 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 136, 0.11%
#LFC < 0 (down)     : 346, 0.29%
#outliers [1]       : 37695, 31%
#low counts [2]     : 42066, 35%
#(mean count < 15)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results


#How many adjusted adjusted p-values were less than 0.1 or 0.05?
#sum(res$padj < 0.05, na.rm=TRUE) #49
#sum(res10$padj < 0.1, na.rm=TRUE) #50


##################
#### MA-PLOT #####
##################
#more useful to visualize the MA-plot for the shrunken log2 fold changes, which remove the noise associated with log2 fold changes from low count genes without requiring arbitrary filtering thresholds.
DESeq2::plotMA(resLFC,ylim=c(-10,10),main="condition_Fan_vs_Arolium",alpha = 0.05,
               colLine = "grey40",  
               colNonSig = "gray60",
               colSig = "blue")
abline(h=c(-0.05,0.05), col="dodgerblue", lwd=2)
#saved as "process_DEseq2_M.GoS_vs_F.GoS_plotMA_resLFC.tiff"

DESeq2::plotMA(res,ylim=c(-25,25),main="condition_Fan_vs_Arolium", alpha = 0.05,
               colLine = "grey40",  
               colNonSig = "gray60",
               colSig = "blue")
abline(h=c(-0.05,0.05), col="dodgerblue", lwd=2)
#saved as "process_DEseq2_M.GoS_vs_F.GoS_plotMA_res.tiff"
#650x550

# Filtrer les résultats pour obtenir uniquement les gènes avec une valeur p ajustée inférieure à 0.1
resSig <- subset(res, padj < 0.05)
#resSig2 <- subset(resLFC, padj < 0.05)

# Afficher les identifiants des gènes differntiellemnt exprimées
row.names(resSig)
#row.names(resSig2)
# Sous-ensemble de 'res' où 'padj' est inférieur à 0.05
resSig <- subset(res, padj < 0.05)
# Créer un tableau à une colonne avec les identifiants des gènes
gene_ids <- data.frame(Gene_ID = row.names(resSig))
# Afficher le tableau
print(gene_ids)
# Écrire le tableau dans un fichier .txt 
write.table(gene_ids, file = "gene_ids.txt", sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)

# Créer un tableau avec les gènes surexprimés (LFC > 0)
overexpressed_genes <- subset(resSig, log2FoldChange > 0)
overexpressed_gene_ids <- data.frame(Gene_ID = row.names(overexpressed_genes))
write.table(overexpressed_gene_ids, file = "overexpressed_gene_ids.txt", sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)

# Créer un tableau avec les gènes sous-exprimés (LFC < 0)
underexpressed_genes <- subset(resSig, log2FoldChange < 0)
underexpressed_gene_ids <- data.frame(Gene_ID = row.names(underexpressed_genes))
write.table(underexpressed_gene_ids, file = "underexpressed_gene_ids.txt", sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)



####################################
######### Significant DEGs #########
####################################

#Create a tibble of shrunked results
resLFC_tb<- resLFC %>%
  data.frame() %>%
  tibble::rownames_to_column(var="gene") %>%
  as_tibble()
#write.table(resLFC_tb, file = "/Users/trystan.bessonnier/Desktop/DESeq2_Transcriptome/DEseq2_results_condition_F_vs_A.txt", sep = "\t", quote=F, row.names = F)

#Subset the tibble to keep only significant genes
sign.DEGs<-resLFC_tb %>%
  dplyr::filter(padj < 0.05 )#& (log2FoldChange < -0.58 | log2FoldChange > 0.58 ))  
sign.DEGs
summary(sign.DEGs)
#4073 DEGs (pval 0.05 L2FC 0.58)
sum(sign.DEGs$log2FoldChange>0)
#1896
sum(sign.DEGs$log2FoldChange<0)
#2177

######################################
## Heat Map all significant genes ####
######################################
normalized_counts<-counts(dds, normalized=T) %>%
  data.frame() %>%
  tibble::rownames_to_column(var="gene")

### Extract normalized expression for significant genes
norm_sig <- normalized_counts[,] %>% 
  dplyr::filter(gene %in% sign.DEGs$gene)  
str(norm_sig)
### Run pheatmap using the metadata data frame for the annotation
library(pheatmap)
# Assurez-vous que les noms des gènes sont les noms des lignes
rownames(norm_sig) <- norm_sig$gene
# Enlevez la colonne des noms de gènes des données de la heatmap
heatmap_data <- norm_sig[,2:10]
# Créez la heatmap
pheatmap(heatmap_data, 
         color=colorRampPalette(c("blue", "white", "red"))(10000),
         cluster_rows = T, 
         show_rownames = T, # Affiche les noms des gènes
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 10, 
         height = 20)
#saved as "WJ-3820_DEseq2_L3.M.P_vs_R_heatmap_padj_0.05_L2FC_0.58.tiff"
#1500x1000
##################
#### HEATMAP #####
## count matrix ##
##################
library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:160] #Play with the value (gene number)
df <- as.data.frame(colData(dds)[,c("condition")])
rownames(df) <- colnames(dds)
pheatmap(assay(ntd)[select,],color=colorRampPalette(c("blue", "white", "red"))(100000), 
         cluster_rows=TRUE, show_rownames=T,scale = "row",
         cluster_cols=F, annotation_col=df)

###########################
########## HEATMAP ########
## intra sample distance ##
###########################
library("pheatmap")
##Extract the vsd matrix from the object
vsd_mat<-assay(vsd) #transform the object in a matrix (2D data structure)
##Compute pairwisse correlation value
vsd_cor<-cor(vsd_mat)
head(vsd_cor)
pheatmap(vsd_cor)

ntd_mat<-assay(ntd) #transform the object in a matrix (2D data structure)
##Compute pairwisse correlation value
ntd_cor<-cor(ntd_mat)
head(ntd_cor)
pheatmap(ntd_cor)

#####################
## Volcanot plot ####
#####################
library(ggplot2)
## Obtain logical vector where TRUE values denote padj values < 0.05 and fold change > 1.5 in either direction
resLFC_tb <- resLFC_tb %>% 
  dplyr::mutate(threshold = padj < 0.05) #& abs(log2FoldChange) >= 0.58) #because log2(1.5)=0.58

ggplot(resLFC_tb) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold)) +
  ggtitle("Leg Fan VS Arolium") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  #scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))


## Create an empty column to indicate which genes to label
resLFC_tb <- resLFC_tb %>% dplyr::mutate(genelabels = "")

## Sort by padj values 
resLFC_tb <- resLFC_tb %>% dplyr::arrange(padj)

## Populate the genelabels column with contents of the gene symbols column for the first 10 rows, i.e. the top 10 most significantly expressed genes
resLFC_tb$genelabels[1:49] <- as.character(resLFC_tb$gene[1:49])

View(resLFC_tb)
library(ggrepel)

ggplot(resLFC_tb, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(colour = threshold)) +
  geom_text_repel(aes(label = genelabels)) +
  ggtitle("Leg Fan vs Arolium") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) 
#saved as "process_DEseq2_M.GoS_vs_F.GoS_volcanoplot_significant_genes.tiff"
#1500x1000

######################
######### PCA ########
######################
#by default, plotPCA uses the top 500 most variable genes, can change with ntop=argument

pcaData <- plotPCA(ntd, intgroup=c("condition"), returnData=TRUE,ntop=50)
pcaData["leg"]<- c("L1","L1","L1","L2","L2","L2","L3","L3","L3")
percentVar <- round(100 * attr(pcaData, "percentVar"))
str(pcaData)
#PCA 1 vs 2
pca<-ggplot(pcaData, aes(PC1, PC2, color=leg)) +
  geom_point(size=0) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  geom_hline(yintercept = 0, color = "black") +  
  geom_vline(xintercept = 0, color = "black") +
  geom_text(aes(label=name))
coord_fixed()
pca

pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE,ntop=50000)
pcaData["leg"]<- c("L1","L1","L1","L2","L2","L2","L3","L3","L3")
percentVar <- round(100 * attr(pcaData, "percentVar"))
str(pcaData)
#PCA 1 vs 2
pca<-ggplot(pcaData, aes(PC1, PC2, color=leg)) +
  geom_point(size=0) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  geom_hline(yintercept = 0, color = "black") +  
  geom_vline(xintercept = 0, color = "black") +
  geom_text(aes(label=name))
coord_fixed()
pca

pcaData <- plotPCA(rld, intgroup=c("condition"), returnData=TRUE,ntop=50000)
pcaData["leg"]<- c("L1","L1","L1","L2","L2","L2","L3","L3","L3")
percentVar <- round(100 * attr(pcaData, "percentVar"))
str(pcaData)
#PCA 1 vs 2
pca<-ggplot(pcaData, aes(PC1, PC2, color=leg)) +
  geom_point(size=0) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  geom_hline(yintercept = 0, color = "black") +  
  geom_vline(xintercept = 0, color = "black") +
  geom_text(aes(label=name))
coord_fixed()
pca



##################
## Count Plot ####
##################

# examine the counts of reads for a single gene across the groups
#gene with minimal pval
plotCounts(dds, gene=which.min(res$padj), intgroup="sex")
#plot the higher LFC
plotCounts(dds, gene=which.max(res$log2FoldChange), intgroup="sex")
plotCounts(dds, gene="g13223", intgroup="Condition")
#saved as 
#600x624

#Use ggplot
# Save plotcounts to a data frame object
d <- plotCounts(dds, gene="g13223", intgroup="Condition", returnData=TRUE)

# What is the data output of plotCounts()?
d %>% View()

# Plot the MOV10 normalized counts, using the samplenames (rownames(d) as labels)
library(ggplot2)
install.packages("ggrepel")
library(ggrepel)

ggplot(d, aes(x = sex, y = count, color = sex)) + 
  geom_point(position=position_jitter(w = 0.1,h = 0)) +
  geom_text_repel(aes(label = rownames(d))) + 
  theme_bw() +
  ggtitle("g17228") +
  theme(plot.title = element_text(hjust = 0.5))






BiocManager::install('EnhancedVolcano')
BiocManager::install("ashr")
library(EnhancedVolcano)


EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue')

EnhancedVolcano(resLFC,
                lab = rownames(resLFC),
                x = 'log2FoldChange',
                y = 'pvalue')

par(mfrow=c(2,2),mar=c(2,2,1,1))
ylim <- c(-28,28)
resGA <- results(dds, lfcThreshold=.5, altHypothesis="greaterAbs")
resLA <- results(dds, lfcThreshold=.5, altHypothesis="lessAbs")
resG <- results(dds, lfcThreshold=.5, altHypothesis="greater")
resL <- results(dds, lfcThreshold=.5, altHypothesis="less")
drawLines <- function() abline(h=c(-.5,.5),col="dodgerblue",lwd=2)
plotMA(resGA, ylim=ylim); drawLines()
plotMA(resLA, ylim=ylim); drawLines()
plotMA(resG, ylim=ylim); drawLines()
plotMA(resL, ylim=ylim); drawLines()
dev.off()


