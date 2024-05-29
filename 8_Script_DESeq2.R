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
citation("EnhancedVolcano")
citation("vsn")
citation("BiocManager")


setwd("/Users/trystan.bessonnier/Desktop/DESeq2_Transcriptome/Control_Deseq")
##########################################
## Import data & sample preparation ######
##########################################
#préparation des données:
#Création de 3 fichiers: file_list.txt, samples.txt, sampleTable.txt
files<-as.character(read.table("file_list_genes.txt", header=FALSE)$V1)

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
                                      design= ~ replicat + condition   ) #class: DESeqDataSet dim: 476176 9 



# Calculer le count moyen pour chaque gène
meanCounts <- rowMeans(counts(ddsMatrix))

# Filtrer les gènes qui ont un count moyen supérieur à 5
ddsMatrix <- ddsMatrix[meanCounts > 5, ]

#normalization
ddsMatrix<- estimateSizeFactors(ddsMatrix)

dds<-DESeq(ddsMatrix)
dds #class: DESeqDataSet 84107 12   
dds$condition <- factor(dds$condition) 
dds$condition <- relevel(dds$condition, ref = "A")


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

#specify contrast? : if one factor has 3 levels; indicate which two levels classes we are interested in comparing
#if not specifying, results() returns the comparison of the last level of last variable in Design formula over the 1st level of this variable (ordered alphabetically if levels are not specified)
resultsNames(dds)
res <- results(dds, contrast=c("condition","F","A"), alpha = 0.01) #same sense
res 



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
resLFC #DataFrame with 84107 rows and 5 columns


#order our results table by the smallest p value
resOrdered <- res[order(res$pvalue),]
summary(res)
#out of 84107 with nonzero total read count
#adjusted p-value < 0.01
#LFC > 0 (up)       : 2343, 2.8%
#LFC < 0 (down)     : 4595, 5.5%
#outliers [1]       : 0, 0%
#low counts [2]     : 45658, 54%
#(mean count < 15)


#change the adjusted pval cutoff (automatically = 0.05)
#res05 <- results(dds, alpha=0.05)
#summary(res05)

####################################
######### Significant DEGs #########
####################################

#Create a tibble of shrunked results
resLFC_tb<- resLFC %>%
  data.frame() %>%
  tibble::rownames_to_column(var="gene") %>%
  as_tibble()
#write.table(resLFC_tb, file = "/Users/trystan.bessonnier/Desktop/DESeq2_Transcriptome/Control_Deseq/DEseq2_results_condition_F_vs_A.txt", sep = "\t", quote=F, row.names = F)

#Subset the tibble to keep only significant genes
sign.DEGs<-resLFC_tb %>%
  dplyr::filter((padj < 0.01 ) & (log2FoldChange < -1 | log2FoldChange > 1 ))
sign.DEGs
summary(sign.DEGs)
#3743
# Filtrer les résultats pour obtenir uniquement les gènes avec une valeur p ajustée inférieure à 0.01
resLFCSig <- subset(resLFC, padj < 0.01 & (log2FoldChange < -1 | log2FoldChange > 1 ))
# Afficher les identifiants des gènes differntiellemnt exprimées
row.names(resLFCSig)
# Sous-ensemble de 'res' où 'padj' est inférieur à 0.01
resLFCSig <- subset(resLFC, padj < 0.01 & (log2FoldChange < -1 | log2FoldChange > 1 ))
# Créer un tableau à une colonne avec les identifiants des gènes
gene_ids <- data.frame(Gene_ID = row.names(resLFCSig))
# Afficher le tableau
print(gene_ids)
# Écrire le tableau dans un fichier .txt 
#write.table(gene_ids, file = "Candidate_shrinkLFC_genes_3743_pvalue_0.01.txt", sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)

sum(resLFCSig$log2FoldChange>1)
#926
# Créer un tableau avec les gènes surexprimés (LFC > 0)
overexpressed_genes <- subset(resLFCSig, log2FoldChange > 1)
overexpressed_gene_ids <- data.frame(Gene_ID = row.names(overexpressed_genes))
#write.table(overexpressed_gene_ids, file = "overexpressed_shrink_LFC_gene_926_ids.txt", sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)


sum(resLFCSig$log2FoldChange<(-1))
#2817
# Sous-ensemble de 'res' où 'padj' est inférieur à 0.01
resLFCSig <- subset(resLFC, padj < 0.01)
# Créer un tableau avec les gènes sous-exprimés (LFC < 0)
underexpressed_genes <- subset(resLFCSig, log2FoldChange < (-1))
underexpressed_gene_ids <- data.frame(Gene_ID = row.names(underexpressed_genes))
#write.table(underexpressed_gene_ids, file = "underexpressed_shrink_LFC_gene_2817_ids.txt", sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)


######################
######### PCA ########
######################
#by default, plotPCA uses the top 500 most variable genes, can change with ntop=argument
library(ggplot2)
pcaData <- plotPCA(rld, intgroup=c("condition","replicat"), returnData=TRUE,ntop=80000)
pcaData["leg"]<- c("L1","L1","L1","L1","L2","L2","L2","L2","L3","L3","L3","L3")
percentVar <- round(100 * attr(pcaData, "percentVar"))
str(pcaData)
#PCA 1 vs 2
tkt<-pca<-ggplot(pcaData, aes(PC1, PC2, color=leg)) +
  ggtitle("PCA des échantillons")+
  geom_point(size=0) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  geom_hline(yintercept = 0, color = "black") +  
  geom_vline(xintercept = 0, color = "black") +
  geom_text(aes(label=name))
coord_fixed()
pca

ggsave(filename = "Beautiful_PCA.png", plot = tkt, width = 10, height = 10, dpi = 300)

##################
#### HEATMAP #####
## count matrix ##
##################
library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:35000] #Play with the value (gene number)
df <- as.data.frame(colData(dds)[,c("condition","replicat")])
rownames(df) <- colnames(dds)
gg<-pheatmap(assay(ntd)[select,],color=colorRampPalette(c("blue", "white", "red"))(35000), 
         cluster_rows=TRUE, show_rownames=F,scale = "row",main="Heatmap des counts normalisés pour 35000 gènes entres les différents réplicats",
         cluster_cols=F, annotation_col=df)

ggsave(filename = "Count_matrix.png", plot = gg, width = 10, height = 10, dpi = 300)
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

ntd_mat<-assay(ntd) 
##Compute pairwisse correlation value
ntd_cor<-cor(ntd_mat)
head(ntd_cor)
pheatmap(ntd_cor)

##################
#### MA-PLOT #####
##################
DESeq2::plotMA(res,ylim=c(-25,25),main="condition_Fan_vs_Arolium", alpha = 0.01,
               colLine = "grey40",  
               colNonSig = "gray60",
               colSig = "blue")


#more useful to visualize the MA-plot for the shrunken log2 fold changes, which remove the noise associated with log2 fold changes from low count genes without requiring arbitrary filtering thresholds.
png(filename = "PlotMA_apelgm.png", width = 10, height = 10, units = "in", res = 300)
DESeq2::plotMA(resLFC,ylim=c(-15,15),main="condition_Fan_vs_Arolium",alpha = 0.01, #(resLFC$log2FoldChange < -1 | resLFC$log2FoldChange > 1 ),
               colLine = "grey40",  
               colNonSig = "gray60",
               colSig = "blue")

dev.off()

resLFCashr <- lfcShrink(dds, coef="condition_F_vs_A", type="ashr")
DESeq2::plotMA(resLFCashr,ylim=c(-15,15),main="condition_Fan_vs_Arolium",alpha = 0.01,
               colLine = "grey40",  
               colNonSig = "gray60",
               colSig = "blue")

resnormal <- lfcShrink(dds, coef="condition_F_vs_A", type="normal")
DESeq2::plotMA(resnormal,ylim=c(-15,15),main="condition_Fan_vs_Arolium",alpha = 0.01,
               colLine = "grey40",  
               colNonSig = "gray60",
               colSig = "blue")
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
heatmap_data <- norm_sig[,2:13]
# Créez la heatmap
beautiful<-pheatmap(heatmap_data, 
         color=colorRampPalette(c("blue", "white", "red"))(10000),
         cluster_rows = T, 
         show_rownames = F, # Affiche les noms des gènes
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 10, 
         height = 20,
         main="Heatmap des 3743 gènes significativement différentiellement exprimés")

ggsave(filename = "Heatmap_all_signi_genes_3743.png", plot = beautiful, width = 10, height = 10, dpi = 300)


#####################
## Volcanot plot ####
#####################
BiocManager::install('EnhancedVolcano')
BiocManager::install("ashr")
library(EnhancedVolcano)

res <- DESeq(dds)

# Effectuer le shrinkage des log fold changes
res <- lfcShrink(dds, coef="condition_F_vs_A", type='apeglm')

# Vérifiez les colonnes disponibles dans res pour confirmer la présence de padj
head(res)

# Renommer le gène Trinity
rownames(res)[rownames(res) == "TRINITY_DN158453_c0_g1"] <- "UBX"
rownames(res)[rownames(res) == "TRINITY_DN4045_c0_g1"] <- "mogsha"
rownames(res)[rownames(res) == "TRINITY_DN9986_c0_g2"] <- "SCR"
#rownames(res)[rownames(res) == "TRINITY_DN1658_c0_g1"] <- "R-ant-yellow"
rownames(res)[rownames(res) == "TRINITY_DN32148_c0_g1"] <- "CP19/7"
rownames(res)[rownames(res) == "TRINITY_DN68375_c0_g1"] <- "endocuticle structural glycoprotein SgAbd-1-like"
rownames(res)[rownames(res) == "TRINITY_DN51005_c0_g1"] <- "transposon TNT 1-94, partial"
rownames(res)[rownames(res) == "TRINITY_DN35_c0_g1"] <- "retrotransposable element R2DM"
rownames(res)[rownames(res) == "TRINITY_DN201171_c0_g1"] <- "transposon TNT 1-94, partial"
rownames(res)[rownames(res) == "TRINITY_DN18770_c0_g1"] <- "piggyBac transposable element"
rownames(res)[rownames(res) == "TRINITY_DN4669_c1_g1"] <- "piggyBac transposable element"
rownames(res)[rownames(res) == "TRINITY_DN12309_c0_g1"] <- "cell division cycle protein 23"

# Générer le volcano plot avec le gène renommé
super <- EnhancedVolcano(res,
                         lab = rownames(res),
                         x = 'log2FoldChange',
                         y = 'padj',
                         title = 'VolcanoPlot: Condition Fan vs Arolium',
                         subtitle = "Highlighting significant genes with log2FC < -1 | > 1  and p-adj < 0.01",
                         pCutoff = 0.01,
                         FCcutoff = 1,
                         pointSize = 2.0,
                         labSize = 4.0,
                         col = c("grey30", "forestgreen", "royalblue", "red2"),
                         legendLabels=c('NS','Log2FC','p-adj','p-adj & Log2FC'),
                         legendPosition = 'top',
                         legendLabSize = 10,
                         legendIconSize = 3.0,
                         drawConnectors = TRUE,
                         widthConnectors = 0.5,
                         colConnectors = 'grey30',
                         max.overlaps = 100,
                         selectLab = c('UBX', 'mogsha', 'SCR', 'CP19/7','endocuticle structural glycoprotein SgAbd-1-like','transposon TNT 1-94, partial','retrotransposable element R2DM', 'piggyBac transposable element', 'piggyBac transposable element','cell division cycle protein 23'), # Utilisation du nouveau nom
                         labCol = 'black',
                         labFace = 'bold',
                         boxedLabels = TRUE)

# Afficher le volcano plot
print(super)

# Sauvegarder le volcano plot
ggsave(filename = "VolcanoPlot_F_vs_A.png", plot = super, width = 10, height = 10, dpi = 300)

