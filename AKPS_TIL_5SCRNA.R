#Vyom Shah
#Cold Spring Harbor Laboratory
#A metabolic switch that boosts immune cell fitness against cancer
#5' scRNA AKPS TIL analysis

#Load in libraries
library(Rmagic)
library(plyr)
library(Seurat)
library(tidyverse)
library(cowplot)
library(ggplot2)
library(MAST)
library(DESeq2)
library(EnhancedVolcano)
library(limma)
library(scales)
library(metR)
library(ggpubr)
library(rstatix)
library(svglite)
library(viridis)
library(reshape)
library(harmony)
library(nichenetr)
library(RColorBrewer)
library(Libra)
library(Nebulosa)

#Load in the data
Control1 <- Read10X_h5("/Users/vyom/cloud_computing/data_norepl/sequencing_backup/VP16_PPAR_AKPS/RNA/count/Beyaz_TM01_Control_1/per_sample_outs/Beyaz_TM01_Control_1/count/sample_filtered_feature_bc_matrix.h5")
Control2 <- Read10X_h5("/Users/vyom/cloud_computing/data_norepl/sequencing_backup/VP16_PPAR_AKPS/RNA/count/Beyaz_TM01_Control_2/per_sample_outs/Beyaz_TM01_Control_2/count/sample_filtered_feature_bc_matrix.h5")

VP161 <- Read10X_h5( "/Users/vyom/cloud_computing/data_norepl/sequencing_backup/VP16_PPAR_AKPS/RNA/count/Beyaz_TM01_Exp_1/per_sample_outs/Beyaz_TM01_Exp_1/count/sample_filtered_feature_bc_matrix.h5")
VP162 <- Read10X_h5( "/Users/vyom/cloud_computing/data_norepl/sequencing_backup/VP16_PPAR_AKPS/RNA/count/Beyaz_TM01_Exp_2/per_sample_outs/Beyaz_TM01_Exp_2/count/sample_filtered_feature_bc_matrix.h5")


Control1 <- CreateSeuratObject(Control1, project = "Control1")
Control2 <- CreateSeuratObject(Control2, project = "Control2")

VP161 <- CreateSeuratObject(VP161, project = "VP161")
VP162 <- CreateSeuratObject(VP162, project = "VP162")

d <- merge(Control1, y = c(Control2, VP161, VP162), add.cell.ids = c("Control1", "Control2","VP161", "VP162"), project = "Immune")

#Filter the data
d[["percent.mt"]] <- PercentageFeatureSet(d, pattern = "mt-")
VlnPlot(d, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),pt.size = 0 ,ncol = 3)
d <- subset(d, subset = nCount_RNA > 200 & nCount_RNA < 10000 & nFeature_RNA > 500 & nFeature_RNA < 3000 & percent.mt < 5)
d

#Normalization
Data.list <- SplitObject(d, split.by = "ident")
Data.list <- Data.list[c("Control1", "Control2","VP161", "VP162")]
for (i in 1:length(Data.list)) {
  
  Data.list[[i]] <- SCTransform(Data.list[[i]], verbose = FALSE)
}

#select highly variable genes 
Data.features <- SelectIntegrationFeatures(object.list = Data.list, nfeatures = 5000)
Data.list <- PrepSCTIntegration(object.list = Data.list, anchor.features = Data.features, 
                                verbose = FALSE)
Data.anchors <- FindIntegrationAnchors(object.list = Data.list, normalization.method = "SCT", 
                                       anchor.features = Data.features, verbose = FALSE)
vp16_obj <- IntegrateData(anchorset = Data.anchors, normalization.method = "SCT", 
                          verbose = TRUE)

# Visulization and Clustering
vp16_obj <- RunPCA(vp16_obj)
VizDimLoadings(vp16_obj, dims = 1:2, reduction = "pca")

#select appropriate PCs
DimPlot(vp16_obj, reduction = "pca")
ElbowPlot(vp16_obj, ndims = 50, reduction = "pca")

#calculate UMAP and Clustering
vp16_obj <- RunUMAP(vp16_obj, dims = 1:25, n.epochs = 500)

vp16_obj <- FindNeighbors(vp16_obj, dims = 1:25)

vp16_obj <- FindClusters(vp16_obj, resolution = 1)
DimPlot(vp16_obj, reduction = "umap", label = TRUE)

#identify cluster identifers
markers <- FindAllMarkers(vp16_obj, only.pos = TRUE, min.pct = 0.1, logfc.threshold = .5)
library(readxl)
write.csv(markers,'Gene_de_Sigs_Per_Clust_immune.csv')

DotPlot_Sig <- unique(c('Cd3g','Cd3e','Cd8a','Cd8b1','Cd4','Trac','Tcrg-C1','Lag3','Pdcd1','Havcr2','Tox','Tcf7','Gzmb','Tbx21','Ifng','Gata3','Il5','Il13','Rorc','Il17a','Il17f','Foxp3','Il10','Il2rb','Il2ra','Klrd1','Cd19','Jchain','Ighm','Ighg1','S100a8','S100a9','Ly6a','Cd14','Cd68','Cd74','Ciita','Nrc1','Klre1','Itgam','Itgax','H2-Eb1','H2-Ab1','Arg1','Mrc1','Tgfbi','Ccr2','Vegfa','Prdx1','Clec4d','Ccl5','Cd83','Ccr7','Fcn1','Msrb1','Ly6g','Col3a1','Sparc'))
DotPlot_Sig <- unique(c('Cd19','Ighm', 'Ighd','Sell', 'Cd69', 'Cd83','Jun', 'Irf4','Myc', 'Xbp1','Cd79b', 'Cd79a','Ms4a1','Ssr4', 'Derl3', 'Fkbp11', 'Prdx4','Igha','Iglc1','Iglc2')) 
DotPlot_Sig <- unique(c('Jchain','Cd5',' Oct2','Pax5','Havcr')) 

DotPlot(vp16_obj, features = DotPlot_Sig, assay = 'SCT') + labs(y= "Cell Treatment", x="") + scale_colour_distiller( palette ="RdYlBu") + scale_size(range = c(0, 1)) +
  theme(text = element_text(size=5), axis.text.x = element_text(size = 5, angle = 90, hjust = 1, vjust= .01), axis.text.y = element_text(size = 6), axis.title.x = element_blank(), axis.title.y = element_blank())


#Annotate clusters
new.cluster.ids <- c("B", "B", "B", "B", "Mac.", "B", "Neut.", "CD8 T", "CD4 T", "B", "B", "B", "Mono", "B", "Mac.", "B", "Plasma", "Mac.", "CD4 T", "Mac.", "Plasma", "Plasma", "NK", "Plasma", "Plasma", "B", "B", "Mac.")

vp16_obj[["Cell_Type"]] <- Idents(vp16_obj)
names(new.cluster.ids) <- levels(vp16_obj)
vp16_obj <- RenameIdents(vp16_obj, new.cluster.ids)
vp16_obj[["Cell_Type"]] <- Idents(vp16_obj)
DimPlot(vp16_obj, reduction = "umap", group.by= 'Cell_Type')

unique(vp16_obj$Cell_Type)
my_levels <- c("CD8 T","CD4 T",'NK',"B","Plasma", 'Mono', 'Mac.' ,  'Neut.')
vp16_obj$Cell_Type <- factor(vp16_obj$Cell_Type, levels = my_levels)

#Differential Expression
vp16_obj <- PrepSCTFindMarkers(vp16_obj)
Idents(vp16_obj) <- vp16_obj$Treatment
DE_all <- FindMarkers(vp16_obj, ident.1 = "VP16", ident.2 = "Control", test.use = "MAST", logfc.threshold = .01, min.pct = .1, assay = 'SCT')
write.csv(DE_all, paste0('All_AKPS_VP16_scrna_DE.csv')) 
if (max(-log(DE_all$p_val_adj)) < 321){
  ylim_num <- max(-log(DE_all$p_val_adj))} else {
    ylim_num <- 320}
EnhancedVolcano(DE_all, lab = rownames(DE_all), x = 'avg_log2FC', y = 'p_val_adj', title = 'VP16 vs Control',
                pCutoff = .05, FCcutoff = 0.25, xlim = c(min(DE_all$avg_log2FC)-.05,max(DE_all$avg_log2FC)+.05),ylim = c(0,ylim_num+5),
                subtitle = 'ALL CELLS', gridlines.minor = FALSE, gridlines.major = FALSE)
ggsave(file = paste0('All_DE_peak_volcano_SCRNA.pdf'), width=6, height=6, units="in")


Cell_Types <- c(levels(as.factor(vp16_obj$Cell_Type)))
Idents(vp16_obj) <- vp16_obj$Cell_Type
for(i in Cell_Types){
  subset_cell <- subset(vp16_obj,  idents = i)
  subset_cell <- PrepSCTFindMarkers(subset_cell, verbose = TRUE)
  Idents(subset_cell) <- subset_cell$Treatment
  DE_subset <- FindMarkers(subset_cell, ident.1 = "VP16", ident.2 = "Control",  test.use = "MAST", logfc.threshold = .1, min.pct = .01, assay = 'SCT')
  write.csv(DE_subset, paste0(i,'_AKPS_VP16_scrna_DE.csv')) 
  if (max(-log(DE_subset$p_val_adj)) < 321){
    ylim_num <- max(-log(DE_subset$p_val_adj))} else {
      ylim_num <- 320}
  EnhancedVolcano(DE_subset, lab = rownames(DE_subset), x = 'avg_log2FC', y = 'p_val_adj', title = 'VP16 vs Control',
                  pCutoff = .05, FCcutoff = 0.1, xlim = c(min(DE_subset$avg_log2FC)-.05,max(DE_subset$avg_log2FC)+.05),ylim = c(0,ylim_num+5),
                  subtitle = i, gridlines.minor = FALSE, gridlines.major = FALSE)
  ggsave(file = paste0(i,'_DE_peak_volcano_SCRNA.pdf'), width=6, height=6, units="in")
}

#violin plots:
use_python("/Users/vyom/miniconda3/bin/python")
vp16_obj <- magic(vp16_obj)

source("~/analysis/split_violin.R")
DefaultAssay(vp16_obj) <- "MAGIC_SCT"
gene.list <- unique(c('Cd8a','Ifng','Prf1',  'Gzmk','Cxcr3', 'Ifi208',"Irf7",'Irf1'))

gene.list <- unique(c('Cxcr3','Xcr1','Notch1',  'Ccr1','Il10ra'))

DefaultAssay(vp16_obj) <- "MAGIC_SCT"
selected_cells <- names(vp16_obj$Cell_Type[vp16_obj$Cell_Type %in% c("CD8 T")])
vln_data <- FetchData(vp16_obj,
                      vars = c(gene.list,"Treatment"),
                      cells = selected_cells,
                      slot = "data")

vln_data <- melt(vln_data)

All <- ggplot(vln_data, aes(x = variable, y = value, fill= Treatment)) + geom_split_violin(scale="width", trim = TRUE, size = .1) + theme_vyom + theme(legend.position = 'none') + 
  geom_boxplot(width=0.2,outlier.shape = NA, coef = 0, lwd=.2) + xlab('') + ylab('Expression Level')+ scale_fill_manual(values= c('#1b9e77' ,'#d95f02')) +
  theme(text = element_text(size=7), axis.text.x = element_text(size = 6, angle = 45, hjust = 1, vjust = 1), axis.text.y = element_text(size = 6), axis.title.y = element_text(size = 7)) +
  stat_compare_means( method = "wilcox.test",  label = "p.signif", size = 2) + scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
All

#score calculation
library(readxl)
IFN_gene_list <- read_excel("Documents/IFN gene list for Vyom.xlsx", 
                            sheet = "Vyom_Formating")


mean.exp <- zscore(colMeans(x = vp16_obj@assays$MAGIC_SCT@data[IFNG_Response, ], na.rm = TRUE), dist = 'norm')
if (all(names(x = mean.exp) == rownames(x = vp16_obj@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  vp16_obj@meta.data$IFNG_Response <- mean.exp
}

mean.exp <- zscore(colMeans(x = vp16_obj@assays$MAGIC_SCT@data[Ifna_Response, ], na.rm = TRUE), dist = 'norm')
if (all(names(x = mean.exp) == rownames(x = vp16_obj@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  vp16_obj@meta.data$Ifna_Response <- mean.exp
}

mean.exp <- zscore(colMeans(x = vp16_obj@assays$MAGIC_SCT@data[Inflammatory_Response, ], na.rm = TRUE), dist = 'norm')
if (all(names(x = mean.exp) == rownames(x = vp16_obj@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  vp16_obj@meta.data$Inflammatory_Response <- mean.exp
}

mean.exp <- zscore(colMeans(x = vp16_obj@assays$MAGIC_SCT$data[PPAR_targets, ], na.rm = TRUE), dist = 'norm')
if (all(names(x = mean.exp) == rownames(x = vp16_obj@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  vp16_obj@meta.data$PPAR_targets <- mean.exp
}


score.list <- c('Ifn_all', 'Ifna_Response', 'IFNG_Response', 'MGIe_cleanup')

for(i in gene.list) {
  DefaultAssay(vp16_obj) <- 'MAGIC_SCT'
  selected_cells <- names(vp16_obj$Cell_Type)
  vln_data <- FetchData(vp16_obj,
                        vars = c(i,"Treatment", "Cell_Type"),
                        cells = selected_cells,
                        layer = "data")
  vln_data$Treatment
  vln_data <- melt(vln_data)
  
  All <- VlnPlot(vp16_obj, split.by = "Treatment",group.by = 'Cell_Type', features = i, pt.size = 0, assay = "MAGIC_SCT",  cols = c('#1b9e77' ,'#d95f02'), log = TRUE, split.plot = TRUE) + 
    theme(legend.position = 'none') + 
    geom_boxplot(width=0.3,outlier.shape = NA, coef = 0, lwd= .2) + xlab('') + 
    theme(text = element_text(size=9), axis.text.x = element_text(size = 9), axis.text.y = element_text(size =9)) +
    stat_compare_means(data= vln_data, aes(x = Cell_Type, y = value, fill = Treatment), method = "wilcox.test",  label = "p.signif", size = 4, hide.ns = TRUE) #+ scale_y_continuous(expand = expansion(mult = c(0, .1)))
  All$layers[[1]]$aes_params$size = .15
  All
  ggsave(file = paste0('AKPS_VP16_vln_', i, '.pdf'), plot=All, width=2.5, height=2.5, units="in")
}

DefaultAssay(vp16_obj) <- 'MAGIC_SCT'

{
  MHC.score <- c('H2-D1', 'H2-K1', 'H2-Q1', 'H2-Aa', 'Cd74', 'H2-Eb1', 'H2-Ab1')
  survival.score <- c('Jun',	'Junb','Bcl2', 'Birc5', 'Pim2', 'Socs3')
  Chemokine.score <- c('Cxcr3', 'Ccl6', 'Ctsb','Ctsw', 'Ccl12', 'Ccl5', 'Ccl4', 'Ccl2')
  Interferon.score <- c('Irf7','Ifitm10', 'Ifi205', 'Ifit2', 'Ifit1bl1', 'Ifitm3', 'Ifit1', 'Ifi208', 'Irf1')
  Metabolism.score <- c('Acbd7', 'Acbd6', "Acot9" , "Mecr" , "Acot12" , "Acot13", "Mcat" , "Acsf2" , "Hadhb", "Acot9" , "Acot3" , "Acot7" , "Acadvl" , "Eci1", "Dbi" , "Pctp" , "Echs1" , "Mcee" , "Acadl" )
  PPAR.score <- c("Lpl", "Sema3c", "C1qa", "Pltp", "Lgals1", "Prelid2", "Anxa2", "Tpi1", 
                  "Anxa5", "Impa2", "Glrx", "Ddit4", "Cited2", "Ech1", "Tagln2", "Lgals4", 
                  "Pgk1", "Prdx6", "Gapdh", "Plbd1", "Ucp2", "Acot7", "Acat1", "Acaa1a", 
                  "Cyb5r3", "Ppa1", "Ndufv2", "Idh2", "Ifrd1", "Timm17a", "Hsd17b10", 
                  "Fdps", "Atxn10", "Jun", "Ndufa5", "Ndufa8", "Ier2", "Slc15a2", "Atf3", 
                  "Pnpla2", "Cd82", "Grpel1", "Dbi", "Cycs", "Ndufb6", "Ndufb5", "Uqcrc2", 
                  "Hadhb", "Pebp1", "Tomm40", "Etfb", "Mafb", "Ndufa11", "Ghitm", "Scp2", 
                  "Grhpr", "Cmpk1", "Uqcr10", "Twf2", "Mrps36", "Sdhd", "Timm23", "Ndufab1", 
                  "Ndufb10", "Uqcrq", "Cyc1", "Krtcap2", "Ndufb7", "Snrpb", "Pex13")

  All_Genes <- rownames(vp16_obj@assays$MAGIC_SCT)
  MHC.score <- intersect(All_Genes, MHC.score)
  survival.score <- intersect(All_Genes, survival.score)   
  Chemokine.score <- intersect(All_Genes, Chemokine.score)  
  Interferon.score <- intersect(All_Genes, Interferon.score)
  Metabolism.score <- intersect(All_Genes, Metabolism.score)
  PPAR.score <- intersect(All_Genes, PPAR.score)
  
  vp16_obj <- AddModuleScore(object = vp16_obj, features = list(MHC.score), name = 'MHC', assay = 'MAGIC_SCT')
  vp16_obj <- AddModuleScore(object = vp16_obj, features = list(survival.score), name = 'survival', assay = 'MAGIC_SCT')
  vp16_obj <- AddModuleScore(object = vp16_obj, features = list(Chemokine.score), name = 'Chemokine', assay = 'MAGIC_SCT')
  vp16_obj <- AddModuleScore(object = vp16_obj, features = list(Interferon.score), name = 'Interferon', assay = 'MAGIC_SCT')
  vp16_obj <- AddModuleScore(object = vp16_obj, features = list(Metabolism.score), name = 'Metabolism', assay = 'MAGIC_SCT')
  vp16_obj <- AddModuleScore(object = vp16_obj, features = list(PPAR.score), name = 'PPAR', assay = 'MAGIC_SCT')
  
  Idents(vp16_obj) <- vp16_obj$Cell_Type
  Cell_Types <- levels(vp16_obj$Cell_Type)
  MHC_FC <- list()
  survival_FC <- list()
  Meta_FC <- list()
  Chemokine_FC <- list()
  Interferon_FC <- list() 
  PPAR_FC <- list() 
  
  for(i in Cell_Types){
    MHC_FC[[i]] <- FoldChange(vp16_obj,ident.1 = "VP16", ident.2 = "Control",features = MHC.score, group.by = 'Treatment', subset.ident = i)
  }
  MHC_FC_aggr <- c()
  for(i in Cell_Types){
    MHC_FC_aggr <- c(MHC_FC_aggr, mean(as.numeric(MHC_FC[[i]]$avg_log2FC)))
  }
  
  for(i in Cell_Types){
    survival_FC[[i]] <- FoldChange(vp16_obj,ident.1 = "VP16", ident.2 = "Control",features = survival.score, group.by = 'Treatment', subset.ident = i)
  }
  survival_FC_aggr <- c()
  for(i in Cell_Types){
    survival_FC_aggr <- c(survival_FC_aggr, mean(as.numeric(survival_FC[[i]]$avg_log2FC)))
  }
  
  for(i in Cell_Types){
    Meta_FC[[i]] <- FoldChange(vp16_obj,ident.1 = "VP16", ident.2 = "Control",features = Metabolism.score, group.by = 'Treatment', subset.ident = i)
  }
  Meta_FC_aggr <- c()
  for(i in Cell_Types){
    Meta_FC_aggr <- c(Meta_FC_aggr, mean(as.numeric(Meta_FC[[i]]$avg_log2FC)))
  }
  
  for(i in Cell_Types){
    Chemokine_FC[[i]] <- FoldChange(vp16_obj,ident.1 = "VP16", ident.2 = "Control",features = Chemokine.score, group.by = 'Treatment', subset.ident = i)
  }
  Chemokine_FC_aggr <- c()
  for(i in Cell_Types){
    Chemokine_FC_aggr <- c(Chemokine_FC_aggr, mean(as.numeric(Chemokine_FC[[i]]$avg_log2FC)))
  }
  for(i in Cell_Types){
    Interferon_FC[[i]] <- FoldChange(vp16_obj,ident.1 = "VP16", ident.2 = "Control",features = Interferon.score, group.by = 'Treatment', subset.ident = i)
  }
  Interferon_FC_aggr <- c()
  for(i in Cell_Types){
    Interferon_FC_aggr <- c(Interferon_FC_aggr, mean(as.numeric(Interferon_FC[[i]]$avg_log2FC)))
  }
  
  for(i in Cell_Types){
    PPAR_FC[[i]] <- FoldChange(vp16_obj,ident.1 = "VP16", ident.2 = "Control",features = PPAR.score, group.by = 'Treatment', subset.ident = i)
  }
  PPAR_FC_aggr <- c()
  for(i in Cell_Types){
    PPAR_FC_aggr <- c(PPAR_FC_aggr, mean(as.numeric(PPAR_FC[[i]]$avg_log2FC)))
  }
  
  Module_Heatmap <- data.frame(MHC_FC_aggr, survival_FC_aggr, Chemokine_FC_aggr, Interferon_FC_aggr, PPAR_FC_aggr)
  rownames(Module_Heatmap) <- Cell_Types
  colnames(Module_Heatmap) <- c('MHC','survival','Chemokine','Interferon', 'PPARD')
  Module_Heatmap$Cell_Type <- rownames(Module_Heatmap)
  Module_Heatmap1 <-melt(Module_Heatmap)
  
  scores <- c('MHC1','survival1','Chemokine1','Interferon1', 'PPAR1' )
  Idents(vp16_obj) <- vp16_obj$Cell_Type
  pvals <- c()
  for(j in scores){
    for(i in Cell_Types){
      selected_cells <- WhichCells(object = vp16_obj ,idents = c(i))
      vln_data <- FetchData(vp16_obj,
                            vars = c(j,"Treatment", "Cell_Type"),
                            cells = selected_cells,
                            slot = "data")
      vln_data <- melt(vln_data)
      vln_data <- vln_data[c('Treatment', 'Cell_Type', 'value')]
      
      pvals <- c(pvals, t.test(vln_data[which(vln_data$Treatment == 'VP16'),]$value, vln_data[which(vln_data$Treatment == 'Control'),]$value,  alternative = c("greater"))$p.value)
      print(i)
    }
  }
  pvals <- p.adjust(pvals, method = "BH", n = length(pvals))
  Module_Heatmap1$pval <- -log(pvals) 
  Cell_Types
  my_levels <- c("CD8 T","CD4 T",'NK',"B","Plasma", 'Mono', 'Mac.' , 'Neut.')
  Module_Heatmap1 <- Module_Heatmap1[Module_Heatmap1$Cell_Type %in% my_levels,]
  #my_levels <- c("CD8+ T","CD4+ T","TReg","CD8+ TReg","DPT","B","Myeloid")
  Module_Heatmap1$Cell_Type <- factor(Module_Heatmap1$Cell_Type, levels = my_levels)
  min(Module_Heatmap1$pval)
  pathwayColorsDiff = rev(brewer.pal(9, "RdBu"))
  ggplot(Module_Heatmap1, aes(Cell_Type, variable, fill= value)) + geom_tile() + 
    scale_fill_gradientn(name = "Log2FC", colors = pathwayColorsDiff, limits = c(-1.7,1.05)) + coord_flip() +
    theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45), axis.text.y =   element_text( vjust = 0.1, hjust = 1), panel.background = element_blank(), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())
  ggsave(file="VP16Control_score_heatmap.pdf",  width=3, height=3, units="in")
  
  ggplot(Module_Heatmap1, aes(x = Cell_Type, y = variable, size = pval ,color= value)) + geom_point() + theme_vyom + 
    scale_color_gradientn(name = "Log2FC", colors = pathwayColorsDiff, limits = c(-1.9,1.81), trans = scales::transform_modulus(-2.5)) + coord_flip() + labs(size="-log(p)") +
    theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45), axis.text.y =  element_text( vjust = 0.1, hjust = 1), panel.background = element_blank(), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())       
  ggsave(file="VP16Control_score_dotplot.pdf",  width=2.85, height=2.5, units="in")
}

#write metadata
write.csv(as.data.frame(vp16_obj@meta.data), 'AKPS_metadata.csv')

#T Cell and B Cell Clonality

library(scRepertoire)
#add Vdj_T and Vdj_B information in T Cell clonality

# Load VDJ data (one csv per run)
Control_1<- read.csv("/Users/vyom/cloud_computing/data_norepl/sequencing_backup/VP16_PPAR_AKPS/RNA/count/Beyaz_TM01_Control_1/per_sample_outs/Beyaz_TM01_Control_1/vdj_t/filtered_contig_annotations.csv")
Control_2 <- read.csv("/Users/vyom/cloud_computing/data_norepl/sequencing_backup/VP16_PPAR_AKPS/RNA/count/Beyaz_TM01_Control_2/per_sample_outs/Beyaz_TM01_Control_2/vdj_t/filtered_contig_annotations.csv")
VP16_1 <- read.csv("/Users/vyom/cloud_computing/data_norepl/sequencing_backup/VP16_PPAR_AKPS/RNA/count/Beyaz_TM01_Exp_1/per_sample_outs/Beyaz_TM01_Exp_1/vdj_t/filtered_contig_annotations.csv")
VP16_2 <- read.csv("/Users/vyom/cloud_computing/data_norepl/sequencing_backup/VP16_PPAR_AKPS/RNA/count/Beyaz_TM01_Exp_2/per_sample_outs/Beyaz_TM01_Exp_2/vdj_t/filtered_contig_annotations.csv")

contig_list <- list(Control_1, Control_2, VP16_1, VP16_2)

combined.TCR <- combineTCR(contig_list,
                           ID = c("Control1", "Control2", "VP161", "VP162"),
                           samples = c("Control", "Control", "VP16", "VP16"),
                           removeNA = FALSE, 
                           removeMulti = FALSE, 
                           filterMulti = FALSE)


#basic Analysis
clonalQuant(combined.TCR, group.by = 'sample',
            cloneCall="strict", 
            chain = "both", 
            scale = TRUE, palette = 'Pastel 1') + scale_color_manual(values = c('#1b9e77' ,'#d95f02')) + scale_fill_manual(values = c('#1b9e77' ,'#d95f02'))
ggsave(file="AKPS_vdj_unique_clones.pdf",  width=3, height=4)

clonalAbundance(combined.TCR, group.by = 'sample',
                cloneCall = "gene", 
                scale = FALSE, palette = 'Pastel 1')+ scale_color_manual(values = c('#1b9e77' ,'#d95f02')) + scale_fill_manual(values = c('#1b9e77' ,'#d95f02'))
ggsave(file="AKPS_vdj_clone_abundances_num.pdf",  width=4, height=3)

clonalAbundance(combined.TCR, cloneCall = "gene",group.by = 'sample', scale = TRUE, palette = 'Pastel 1')

#visualizing clonal dynamics:
clonalRarefaction(combined.TCR, group.by = 'sample',
                  plot.type = 1,
                  hill.numbers = 1,
                  n.boots = 2) + scale_color_manual(values = c('#1b9e77' ,'#d95f02')) + scale_fill_manual(values = c('#1b9e77' ,'#d95f02'))
ggsave(file="AKPS_vdj_diversity_num.pdf",  width=5, height=4)


clonalRarefaction(combined.TCR, group.by = 'sample',
                  plot.type = 2,
                  hill.numbers = 1,
                  n.boots = 2, palette = 'Pastel 1') + scale_color_manual(values = c('#1b9e77' ,'#d95f02')) + scale_fill_manual(values = c('#1b9e77' ,'#d95f02'))
ggsave(file="AKPS_vdj_coverage_num.pdf",  width=5, height=4)

clonalRarefaction(combined.TCR,group.by = 'sample',
                  plot.type = 3,
                  hill.numbers = 1,
                  n.boots = 2, palette = 'Pastel 1') + scale_color_manual(values = c('#1b9e77' ,'#d95f02')) + scale_fill_manual(values = c('#1b9e77' ,'#d95f02'))
ggsave(file="AKPS_vdj_diversity_coverage.pdf",  width=5, height=4)

#Integrate clonal information into seurat object
vp16_obj <- combineExpression(combined.TCR, 
                              vp16_obj, 
                              cloneCall="gene", 
                              group.by = "sample", 
                              proportion = FALSE,
                              cloneSize=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))

DimPlot(vp16_obj, group.by = "cloneSize") 
FeaturePlot(vp16_obj, features = 'cloneSize')
vp16_obj$CTaa
vp16_obj$cloneSize
Prop_table<- prop.table(x = table(vp16_obj$cloneSize, vp16_obj$Treatment), margin = 2)
Prop_Table <- as.data.frame(Prop_table, row.names = NULL, optional = FALSE,
                            make.names = TRUE, stringsAsFactors = FALSE)
Prop_Table1 <- Prop_Table
Prop_Table1$Var1 <- factor(Prop_Table1$Var1)
plot <- ggplot(data = Prop_Table1, aes(Var1, Freq, fill=Var2)) + geom_bar(position="dodge", stat="identity", na.rm = TRUE) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), panel.background = element_rect(fill = "white", colour = "Black")) + xlab("Cell Treatment") + ylab("Fraction of Cells") + labs(fill = "Sample") +scale_y_continuous(expand = expansion(mult = c(0, .1)))
plot
pastel_colors <- c("#CC99FF","#FFA07A", "#CC99FF","#374259","#008080", "#374259", "#CCFF99" ,'#CBB279',  "#008080", "#FFA07A", "#D3D3D3")

plot <- ggplot(data = Prop_Table1, aes(Var2, Freq, fill=Var1)) + scale_fill_manual(values = pastel_colors)+ geom_bar(position="stack", stat="identity", na.rm = TRUE) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), panel.background = element_rect(fill = "white", colour = "Black")) + xlab("Cell Treatment") + ylab("Fraction of Cells") + labs(fill = "Sample") +scale_y_continuous(expand = expansion(mult = c(0, .1)))
plot
ggsave(file="AKPS_vdj_clonal_stack.pdf",  width=3.8, height=3)


DefaultAssay(clone_obj) <- "MAGIC_SCT"
gene.list <- c('Ciita', 'Cd74','H2-T3', 'H2-D1', 'H2-Q1', 'H2-Q10', 'H2-Q2', 'H2-Ab1','H2-Aa')
gene.list <- unique(c('Cd8a','Cd8b1','Ifng','Prf1', 'Bcl2', 'Gzmk','Fasl','Tnf','Il17ra','Irf4'))
gene.list <- unique(c('Cd8b1','Ifng','Prf1',  'Gzmk','Cxcr3', 'Ifi208',"Irf7",'Irf1'))

gene.list <- unique(c('Cxcr3','Ifng','Gzmk', 'Gzmb', 'Prf1', 'Foxp3', 'Lag3', 'Havcr2', 'Pdcd1', 'Irf1', 'Irf7','Ifi208'))
gene.list <- unique(c('Ifng','Gzmk', 'Gzmb','Irf1', 'Irf7','Ifi208'))


DefaultAssay(vp16_obj) <- "MAGIC_SCT"
selected_cells <- names(vp16_obj$Cell_Type[vp16_obj$Cell_Type %in% c('CD8 T')])
vln_data <- FetchData(vp16_obj,
                      vars = c(gene.list,"Clonality"),
                      cells = selected_cells,
                      slot = "data")

vln_data <- melt(vln_data)

All <- ggplot(vln_data, aes(x = variable, y = value, fill= Clonality)) + geom_split_violin(scale="width", trim = TRUE, size = .1) + theme_vyom + theme(legend.position = 'none') + 
  geom_boxplot(width=0.2,outlier.shape = NA, coef = 0, lwd=.2) + xlab('') + ylab('Expression Level')+ scale_fill_manual(values= c('#7CA1CC' ,'#FF4902')) +
  theme(text = element_text(size=7), axis.text.x = element_text(size = 6, angle = 45, hjust = 1, vjust = 1), axis.text.y = element_text(size = 6), axis.title.y = element_text(size = 7)) +
  stat_compare_means( method = "wilcox.test",  label = "p.signif", size = 2) + scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +ylim(0,1.5)
All
#ggsave(file = 'clone_obj_suppression_TReg_vln.pdf', plot=All, width=3, height=2, units="in")
ggsave(file = 'clone_Effector_CD8_on_clonal_populations.pdf', plot=All, width=1, height=1.75, units="in")


#BCR
Control_1<- read.csv("/Users/vyom/cloud_computing/data_norepl/sequencing_backup/VP16_PPAR_AKPS/RNA/count/Beyaz_TM01_Control_1/per_sample_outs/Beyaz_TM01_Control_1/vdj_b/filtered_contig_annotations.csv")
Control_2 <- read.csv("/Users/vyom/cloud_computing/data_norepl/sequencing_backup/VP16_PPAR_AKPS/RNA/count/Beyaz_TM01_Control_2/per_sample_outs/Beyaz_TM01_Control_2/vdj_b/filtered_contig_annotations.csv")
VP16_1 <- read.csv("/Users/vyom/cloud_computing/data_norepl/sequencing_backup/VP16_PPAR_AKPS/RNA/count/Beyaz_TM01_Exp_1/per_sample_outs/Beyaz_TM01_Exp_1/vdj_b/filtered_contig_annotations.csv")
VP16_2 <- read.csv("/Users/vyom/cloud_computing/data_norepl/sequencing_backup/VP16_PPAR_AKPS/RNA/count/Beyaz_TM01_Exp_2/per_sample_outs/Beyaz_TM01_Exp_2/vdj_b/filtered_contig_annotations.csv")
contig_list <- list(Control_1, Control_2, VP16_1, VP16_2)

combined.BCR <- combineBCR(contig_list,
                           ID = c("Control1", "Control2", "VP161", "VP162"),
                           samples = c("Control", "Control", "VP16", "VP16"),
                           threshold = 0.85)

clonalDiversity(combined.BCR, 
                cloneCall = "gene")

clonalRarefaction(combined.BCR,group.by = 'sample',
                  plot.type = 1,
                  hill.numbers = 1,
                  n.boots = 2, palette = 'Pastel 1') + scale_color_manual(values = c('#1b9e77' ,'#d95f02')) + scale_fill_manual(values = c('#1b9e77' ,'#d95f02'))
ggsave(file="AKPS_vdj_B_diversity_coverage.pdf",  width=5, height=4)


#merge to seurat obj
combined.BCR <- combineBCR(contig_list,
                           samples = c("Control1", "Control2", "VP161", "VP162"),
                           threshold = 0.85)
vp16_obj <- combineExpression(combined.BCR, 
                              vp16_obj, 
                              cloneCall="gene", 
                              group.by = "sample", 
                              proportion = FALSE,
                              cloneSize=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))

DimPlot(vp16_obj, group.by = "cloneSize")
B_Clone_Meta <- vp16_obj@meta.data[vp16_obj@meta.data$Cell_Type %in% c('B', 'Plasma'),]
B_Clone_Meta <- B_Clone_Meta[complete.cases(B_Clone_Meta),]
B_Clone_Meta$Cell_Type <- factor(B_Clone_Meta$Cell_Type, levels = c('B', 'Plasma'))
B_Clone_Meta$cloneSize[is.na(B_Clone_Meta$cloneSize)] <- 'None ( < X <= 0)'
B_Clone_Meta$Treatment_celltype <- paste0(B_Clone_Meta$Cell_Type, "_", B_Clone_Meta$Treatment)

Prop_table<- prop.table(x = table(B_Clone_Meta$cloneSize, B_Clone_Meta$Treatment), margin = 2)

Prop_Table <- as.data.frame(Prop_table, row.names = NULL, optional = FALSE,
                            make.names = TRUE, stringsAsFactors = FALSE)
Prop_Table1 <- Prop_Table
Prop_Table1$Var1 <- factor(Prop_Table1$Var1)
Prop_Table1$Var2 <- factor(Prop_Table1$Var2)
plot <- ggplot(data = Prop_Table1, aes(Var2, Freq, fill=Var1)) + geom_bar(position="dodge", stat="identity", na.rm = TRUE) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), panel.background = element_rect(fill = "white", colour = "Black")) + xlab("Clonal Populations") + ylab("Fraction of Cells") + labs(fill = "Sample") +scale_y_continuous(expand = expansion(mult = c(0, .1)))
plot
pastel_colors <- c("#CC99FF","#FFA07A", "#CBB279","#D3D3D3", "#374259","#008080")

plot <- ggplot(data = Prop_Table1, aes(Var2, Freq, fill=Var1)) + scale_fill_manual(values = pastel_colors)+ geom_bar(position="stack", stat="identity", na.rm = TRUE) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), panel.background = element_rect(fill = "white", colour = "Black")) + xlab("Clonal Populations") + ylab("Fraction of Cells") + labs(fill = "Sample") +scale_y_continuous(expand = expansion(mult = c(0, .1)))
plot
ggsave(file="AKPS_vdj_B_clonal_celltype_stack_by_treatment.pdf",  width=4, height=3)


#expression
Idents(vp16_obj) <- vp16_obj$cloneSize
new.cluster.ids <- c('Clonal','Diverse',  'Diverse', 'Diverse')
vp16_obj[["Clonality"]] <- Idents(vp16_obj)
names(new.cluster.ids) <- levels(vp16_obj)
vp16_obj <- RenameIdents(vp16_obj, new.cluster.ids)
vp16_obj[["Clonality"]] <- Idents(vp16_obj)
DimPlot(vp16_obj, group.by = 'Clonality')

Idents(vp16_obj) <- vp16_obj$cloneSize
clone_obj <- subset(vp16_obj,  idents = c("Single (0 < X <= 1)","Small (1 < X <= 5)",'Medium (5 < X <= 20)',"Large (20 < X <= 100)"))
clone_obj@meta.data$cloneSize <- factor(clone_obj@meta.data$cloneSize, levels = c("Single (0 < X <= 1)","Small (1 < X <= 5)",'Medium (5 < X <= 20)',"Large (20 < X <= 100)"))


DefaultAssay(clone_obj) <- "MAGIC_SCT"
gene.list <- unique(c('Ighg1', "Bcl6",'Bcl2','Myc','Cxcl9','Cxcr4','Ctla4','Ifng', 'Gzmb' ))

selected_cells <- names(clone_obj$Cell_Type[clone_obj$Cell_Type %in% c('Plasma')])

vln_data <- FetchData(clone_obj,
                      vars = c(gene.list,"Clonality"),
                      cells = selected_cells,
                      slot = "data")

vln_data <- melt(vln_data)

All <- ggplot(vln_data, aes(x = variable, y = value, fill= Clonality)) + geom_split_violin(scale="width", trim = TRUE, size = .1) + theme_vyom + theme(legend.position = 'none') + 
  geom_boxplot(width=0.2,outlier.shape = NA, coef = 0, lwd=.2) + xlab('') + ylab('Expression Level')+ scale_fill_manual(values= c('#7CA1CC' ,'#FF4902')) +
  theme(text = element_text(size=7), axis.text.x = element_text(size = 6, angle = 45, hjust = 1, vjust = 1), axis.text.y = element_text(size = 6), axis.title.y = element_text(size = 7)) +
  stat_compare_means( method = "wilcox.test",  label = "p.signif", size = 2) + scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +ylim(0,1.5)
All
#ggsave(file = 'clone_obj_suppression_TReg_vln.pdf', plot=All, width=3, height=2, units="in")
ggsave(file = 'CLonal_B_Plasma_resposne_on_clonal_populations.pdf', plot=All, width=1.6, height=1.75, units="in")


