#Vyom Shah
#Cold Spring Harbor Laboratory
#A metabolic switch that boosts immune cell fitness against cancer
#Naive vs Activated Bulk RNA sequencing

# Load In Libraries
library("tidyverse")
library("AnnotationHub")
library("DESeq2")
library('apeglm')
library('reshape')
library('EnhancedVolcano')
library('fgsea')
library('biomaRt')


load("~/analysis/Naive_Act_Immune_Bulk/deseq2_qc/deseq2.dds.RData")

dds
foo <- counts(dds, normalized = FALSE)
write.csv(foo, file="All_counts.csv")
rld <- rlog(dds, blind=TRUE)
pcaData <- plotPCA(rld, intgroup=c("Intersect"), returnData = TRUE) # vsd and plotPCA are part of DESeq2 package, nothing with my example below. 
percentVar <- round(100 * attr(pcaData, "percentVar")) 

ggplot(pcaData, aes(x = PC1, y = PC2, color = factor(Intersect))) + 
  geom_point(size =3) + 
  xlab(paste0("PC1: ", percentVar[1], "% variance")) + 
  ylab(paste0("PC2: ", percentVar[2], "% variance")) + theme_vyom
ggsave(file="VP16_Steady_Bulk_PCA.pdf",  width=6, height=3, units="in")

dds@colData$Intersect <- dds@colData$sample
dds@colData$Intersect <- factor(c("Activated_Control", "Activated_Control", "Activated_Control", "Activated_Control",
                                  "Activated_Experimental", "Activated_Experimental", "Activated_Experimental", "Activated_Experimental",
                                  "Naive_Control", "Naive_Control", "Naive_Control", "Naive_Control",
                                  "Naive_Experimental", "Naive_Experimental", "Naive_Experimental", "Naive_Experimental"))

dds@colData$Type <- dds@colData$sample
dds@colData$Type <- factor(c("Activated", "Activated", "Activated", "Activated",
                             "Activated", "Activated", "Activated", "Activated",
                             "Naive", "Naive", "Naive", "Naive",
                             "Naive", "Naive", "Naive", "Naive"))

dds@colData$Treatment <- dds@colData$sample
dds@colData$Treatment <- factor(c("Control", "Control", "Control", "Control",
                                  "Experimental", "Experimental", "Experimental", "Experimental",
                                  "Control", "Control", "Control", "Control",
                                  "Experimental", "Experimental", "Experimental", "Experimental"))

dds@colData
dds$sample
levels(dds@colData$Intersect)
my_levels <- c("Naive_Control", "Activated_Control",   "Activated_Experimental",   "Naive_Experimental")
dds@colData$Intersect <- factor(dds@colData$Intersect,levels = my_levels)
design(dds) <- formula(~Intersect)
dds <- DESeq(dds)
resultsNames(dds)
{
  results(dds)
  
  rownames(dds) <- gsub("\\.\\d*", "", rownames(dds))
  hub <- AnnotationHub()
  hubid <- "AH7799"
  anno <- hub[[hubid]]
  genemap <- tibble(gene_id=anno$gene_id,
                    symbol=anno$gene_name) %>%
    distinct()
  
  featureData <- tibble(gene_id=rownames(dds)) %>%
    left_join(genemap, by="gene_id") %>%
    mutate(symbol=case_when(is.na(symbol) ~ gene_id,
                            TRUE ~ symbol)) %>%
    dplyr::select(symbol)
  mcols(dds) <- DataFrame(mcols(dds), featureData)
  
  dds <- nbinomWaldTest(dds)
  res <- results(dds, name='Intersect_Naive_Experimental_vs_Naive_Control')
  res <- lfcShrink(dds, coef="Intersect_Naive_Experimental_vs_Naive_Control", res=res, type='apeglm')
  
  res$gene_name <- mcols(dds)$symbol
  res <- res[order(res$log2FoldChange),]
  res_format <- res %>%
    as.data.frame() %>%
    rownames_to_column(var="gene_id") %>%
    as_tibble() %>%
    rename_at(vars(-gene_id, -gene_name), ~ paste0(., ""))
  res_format[complete.cases(res_format),]
  res_format$log2FoldChange
  #res_format <- res_format[order(-res_format$log2FoldChange),]
  write.csv(res_format, 'Intersect_Naive_Experimental_vs_Naive_Control.csv')
  
  
  Deseq_results <- mcols(dds) 
  Normalized_mean_counts <- assays(dds)[['counts']]
  
  Normalized_mean_counts <- assays(dds)[['rlog']]
  Normalized_mean_counts <- data.frame(Normalized_mean_counts)
  Normalized_mean_counts$Gene_Name <- mcols(dds)$symbol
  Normalized_mean_counts$Zscore <- mcols(dds)$Intersect_Naive_Experimental_vs_Naive_Control
  Normalized_mean_counts[complete.cases(Normalized_mean_counts),]
  Normalized_mean_counts <- Normalized_mean_counts[order(-Normalized_mean_counts$Zscore),]
  
  write.csv(Normalized_mean_counts, 'Intersect_Naive_Experimental_vs_Naive_Control.csv')
}

#rld <- rlog(dds, blind=FALSE)

#Normalized_mean_counts <- assays(rld)[['rlog']]
EnhancedVolcano(res_format, lab = res_format$gene_name, x = 'log2FoldChange', y = 'padj', title = 'Naive_Experimental_vs_Naive_Control', pCutoff = .01, FCcutoff = .5, xlim = c(-5,11), ylim = c(0,60) )
tail(res_format)

Heatmap_data <- data.frame(res_format)
MHC_Genes <- c('H2-q1','H2-T3','H2-DMb2','H2-T10','H2-Ea-ps','H2-Q10','H2-D1')
Klr_Genes <-  c('Klra1','Klra4','Klrk1')
Metabolism_Genes <-  c('Ctsh','Ctse','Lyz2')
Hist_Genes <- c('Hist2h3c2','Hist1h2ag','Hist4h4','Hist1h2ak','Hist1h2ap','Hist2h2aa2')
All_genes <- c('H2-q1','H2-T3','H2-DMb2','H2-T10','H2-Ea-ps','H2-Q10','H2-D1','Klra1','Klra4','Klrk1', 'Ctsh','Ctse','Lyz2','Hist2h3c2','Hist1h2ag','Hist4h4','Hist1h2ak','Hist1h2ap','Hist2h2aa2')
rownames(Heatmap_data) <-  Heatmap_data$gene_name
Heatmap_data <- Heatmap_data[Heatmap_data$gene_name %in% All_genes,]
Heatmap_data <- Heatmap_data[,c('log2FoldChange', 'gene_name')]

GOBP_CHRONIC_INFLAMMATORY_RESPONSE

Heatmap_data <-melt(Heatmap_data)
pathwayColors =rev(viridis::magma(10))
All_genes <- c('H2-q1','H2-T3','H2-DMb2','H2-T10','H2-Ea-ps','H2-Q10','H2-D1','Klra1','Klra4','Klrk1', 'Ctsh','Ctse','Lyz2','Hist2h3c2','Hist1h2ag','Hist4h4','Hist1h2ak','Hist1h2ap','Hist2h2aa2')
Heatmap_data$gene_name <- factor(Heatmap_data$gene_name,levels = All_genes)
ggplot(Heatmap_data, aes(gene_name, variable, fill= value)) + geom_tile() + scale_fill_gradientn(name = "Log2FC", colors = pathwayColors) + theme(panel.background = element_blank())
ggsave(file="Bulk_logfc_scale.pdf",  width=15, height=2, units="in")


# heatmaps for genes of intrest
gene.list.activated <- c('Cdk5rap1','Lta','H2-D1','Plcb4','Ephb6')
gene.list.naive <- c('Ddr1','H2-T10','H2-D1','Rgs11','Dnahc8')
VP16_bulk <- res_format[res_format$gene_name %in% gene.list.activated,]
VP16_bulk <- VP16_bulk[c('log2FoldChange', 'gene_name' )]
VP16_bulk <-melt(VP16_bulk)
pathwayColorsDiff = rev(brewer.pal(9, "RdBu"))
ggplot(VP16_bulk, aes(gene_name, variable, fill= value)) + geom_tile() + 
  scale_fill_gradientn(name = "Log2FC", colors = pathwayColorsDiff, limits = c(-1.49,1.49)) + coord_flip() +
  theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45), axis.text.y =   element_text( vjust = 0.1, hjust = 1), panel.background = element_blank(), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())

ggsave(file="steady_activated_CVP_gene_heatmap.pdf",  width=2, height=3, units="in")
