#Vyom Shah
#Cold Spring Harbor Laboratory
#A metabolic switch that boosts immune cell fitness against cancer
#VP64 CAR-T Bulk RNA Sequencing

# Load In Libraries
library("tidyverse")
library("AnnotationHub")
library("DESeq2")
library('apeglm')
library('reshape')
library('EnhancedVolcano')
library('fgsea')
library('biomaRt')


load("~/data/VP16_Cart_Bulk/deseq2.dds.RData")
Vp64_counts <- counts(dds, normalized=FALSE)
write.csv(Vp64_counts, 'VP64_CART_Counts.csv')

coldata <- as.data.frame(colnames(Vp64_counts))
colnames(coldata) <- 'sample'
dds1 <- DESeqDataSetFromMatrix(countData = Vp64_counts,
                               colData = coldata,
                               design = ~ sample)
dds
rld <- rlog(dds, blind=TRUE)
pcaData <- plotPCA(rld, intgroup=c("sample"), returnData = TRUE) # vsd and plotPCA are part of DESeq2 package, nothing with my example below. 
percentVar <- round(100 * attr(pcaData, "percentVar")) 

ggplot(pcaData, aes(x = PC1, y = PC2, color = factor(sample))) + 
  geom_point(size =3) + 
  xlab(paste0("PC1: ", percentVar[1], "% variance")) + 
  ylab(paste0("PC2: ", percentVar[2], "% variance")) + theme_vyom
ggsave(file="ppar_agonist_Bulk_PCA_per_sample.pdf",  width=6, height=3, units="in")

dds@colData$Treatment <- dds@colData$sample
dds@colData$Treatment <- factor(c("CPT1AM", "CPT1AM", "CPT1AM", "CPT1AM", "CPT1AM", 
                                  "WT", "WT", "WT", "WT", "WT", 
                                  "Control", "Control", "Control", "Control", "Control", 
                                  "VP64", "VP64", "VP64", "VP64", "VP64"))
dds@colData$Treatment <- factor(c("WT", "WT", "WT", "WT", "VP64", "VP64", "VP64", "VP64"))
#dds <- dds[,-c(4, 8, 12, 16, 20)]
colnames(dds)
view(colData(dds))

levels(dds@colData$Treatment)
my_levels <- c( "WT", "VP64", "CPT1AM", "Control")
dds@colData$Treatment <- factor(dds@colData$Treatment,levels = my_levels)


design(dds) <- formula(~Treatment)
dds <- DESeq(dds)
resultsNames(dds)
result_levels <- c("Treatment_GW_vs_Control","Treatment_PPAR1_vs_Control","Treatment_PPAR2_vs_Control","Treatment_PPAR3_vs_Control")
result_levels <- c("Treatment_VP64_vs_WT")
for (i in result_levels){
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
  res <- results(dds, name=i)
  res <- lfcShrink(dds, coef=i, res=res, type='apeglm')
  
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
  
  write.csv(res_format, paste0(i,'_deseq2.csv'))
  res_format_filtered <- res_format[complete.cases(res_format),]
  write.csv(res_format_filtered, paste0(i,'_filter_deseq2.csv'))
  #Deseq_results <- mcols(dds) 
  #Normalized_mean_counts <- assays(dds)[['counts']]
  
  #Normalized_mean_counts <- assays(dds)[['rlog']]
  #Normalized_mean_counts <- data.frame(Normalized_mean_counts)
  #Normalized_mean_counts$Gene_Name <- mcols(dds)$symbol
  #Normalized_mean_counts$Zscore <- mcols(dds)[[i]]
  #Normalized_mean_counts[complete.cases(Normalized_mean_counts),]
  #Normalized_mean_counts <- Normalized_mean_counts[order(-Normalized_mean_counts$Zscore),]
  
  #write.csv(Normalized_mean_counts, paste0(i,'_deseq2_counts.csv'))
}

#rld <- rlog(dds, blind=FALSE)
rld <- rlog(dds, blind=TRUE)
pcaData <- plotPCA(rld, intgroup=c("Treatment"), returnData = TRUE) # vsd and plotPCA are part of DESeq2 package, nothing with my example below. 
percentVar <- round(100 * attr(pcaData, "percentVar")) 

ggplot(pcaData, aes(x = PC1, y = PC2, color = factor(Treatment))) + 
  geom_point(size =3) + 
  xlab(paste0("PC1: ", percentVar[1], "% variance")) + 
  ylab(paste0("PC2: ", percentVar[2], "% variance")) + theme_vyom
ggsave(file="VP64_Bulk_PCA_per_treatment.pdf",  width=4.5, height=3, units="in")
#Normalized_mean_counts <- assays(rld)[['rlog']]
EnhancedVolcano(res_format_filtered, lab = res_format_filtered$gene_name, x = 'log2FoldChange', y = 'padj', title = 'VP64 vs Control', pCutoff = .01, FCcutoff = .5, xlim = c(-3.5,3.5), ylim = c(0,175) )
tail(res_format)
res_format <- read.csv('Treatment_PPAR3_vs_Control_deseq2.csv')

Heatmap_data <- data.frame(res_format)
genes <- c("PDK4", "FCER2", "BCL6", "PLXND1", "CXCL5",  "ANGPTL4", "NR4A2","XCL1",'PLIN2','ASS1','IGF1','GSDMB', 'KDM6B','KLRK',"IL22", "IL9", "IL10")

Heatmap_data <- Heatmap_data[Heatmap_data$gene_name %in% genes,]
Heatmap_data <- Heatmap_data[,c('log2FoldChange', 'gene_name','padj')]


#Heatmap_data <-melt(Heatmap_data)
genes <- c("PDK4", "FCER2", "BCL6", "PLXND1", "CXCL5",  "ANGPTL4", "NR4A2","XCL1",'PLIN2','ASS1','IGF1','GSDMB', 'KDM6B','KLRK',"IL22", "IL9", "IL10")
pathwayColorsDiff = rev(brewer.pal(20, "RdBu"))
Heatmap_data$variable <- ''
Heatmap_data$gene_name <- factor(Heatmap_data$gene_name,levels = genes)
ggplot(Heatmap_data, aes(gene_name, variable, fill= log2FoldChange)) + geom_tile() + 
  scale_fill_gradientn(name = "Log2FC", colors = pathwayColorsDiff, limits = c(-5.54,5.54),trans = scales::transform_modulus(.3)) + coord_flip() +
  theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45), axis.text.y =   element_text( vjust = 0.1, hjust = 1), panel.background = element_blank(), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())

ggsave(file="VP64_CAR_T_Bulk_heatmap.pdf",  width=2, height=2.5, units="in")

Heatmap_data$padj <- -log(Heatmap_data$padj)
ggplot(Heatmap_data, aes(x = gene_name, y = variable, size = padj ,color= log2FoldChange)) + geom_point() + theme_vyom + 
  scale_color_gradientn(name = "Log2FC", colors = pathwayColorsDiff, limits = c(-5.54,5.54),trans = scales::transform_modulus(.3)) + coord_flip() + labs(size="-log(p)") +
  theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45), axis.text.y =  element_text( vjust = 0.1, hjust = 1), panel.background = element_blank(), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())       
ggsave(file="VP16Control_score_dotplot.pdf",  width=2.4, height=3.25, units="in")

library(enrichR)
dbs <- listEnrichrDbs()
dbs <- c("GO_Biological_Process_2023", "GO_Cellular_Component_2023", "GO_Molecular_Function_2023", "MSigDB_Hallmark_2020")
enriched <- enrichr(unique(res_format[(res_format$padj < .05 & res_format$log2FoldChange > .5 ),]$gene_id), dbs)
Enriched1 <- rbind(enriched[[1]],enriched[[2]],enriched[[3]],enriched[[4]] )
Enriched_filter <- Enriched1[Enriched1$Adjusted.P.value < .05,]
plotEnrich(Enriched_filter, showTerms = 100, numChar = 100, y = "Count", orderBy = "Adjusted.P.value")
write.csv(Enriched_filter, 'EnrichR_CART_VP64.csv')
view(enriched[[1]])
Enriched_filter
Enriched_filter$Term
pathways <- c('Positive Regulation Of Natural Killer Cell Mediated Cytotoxicity (GO:0045954)', 'Positive Regulation Of Leukocyte Mediated Cytotoxicity (GO:0001912)', 'Endothelial Cell Chemotaxis (GO:0035767)', 'Nuclear Glucocorticoid Receptor Binding (GO:0035259)', 'Hypoxia', 'TNF-alpha Signaling via NF-kB', 'IL-2/STAT5 Signaling', 'Epithelial Mesenchymal Transition')      

pathwayColors <- rev(c( "#FDC926FF", "#FA9E3BFF", "#ED7953FF", "#D8576BFF", "#BD3786FF", "#9C179EFF", "#7301A8FF", "#47039FFF"))
pathwayColors <- c("#440154", "#482475", "#414487", "#355f8d", "#2a788e", "#21918c", "#22a884", "#44bf70", "#7ad151", "#bddf26", "#fde725")
pathwayColors <- c("#ff0000", "#ff9700", "#d1ff00", "#3aff00", "#00ff5c", "#00fff3", "#0074ff", "#2200ff", "#b900ff", "#47039FFF")
Enriched_filter <- Enriched_filter[order(-Enriched_filter$Combined.Score),]
Enriched_filter <- Enriched_filter[ Enriched_filter$Term %in% pathways,]
Enriched_filter$Term <- factor(Enriched_filter$Term, levels = Enriched_filter$Term)
Enriched_filter$p_adj <- -log(Enriched_filter$Adjusted.P.value)
Enriched_filter$Overlap
Enriched_filter$num_genes <- c(4,4,7,27,10,23,17,16)
Enriched_filter <- Enriched_filter[order(Enriched_filter$num_genes),]
Enriched_filter$Term <- factor(Enriched_filter$Term, levels = Enriched_filter$Term)
ggplot(Enriched_filter, aes(x = Term, y = num_genes, size = p_adj ,color= Odds.Ratio)) + geom_point() + theme_vyom + 
  scale_color_gradientn(name = "Odds Ratio", colors = pathwayColors, limits = c(-0,30)) + coord_flip() + labs(size="-log(p)") +
  theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45), axis.text.y =  element_text( vjust = 0.1, hjust = 1), panel.background = element_blank(), axis.ticks.y=element_blank()) + ylab('Number of Genes')       
ggsave(file="VP64_Cart_GSEA_dotplot.pdf",  width=9.5, height=3, units="in")


