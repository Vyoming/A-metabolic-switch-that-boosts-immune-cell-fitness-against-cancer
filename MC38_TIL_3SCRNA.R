#Vyom Shah
#Cold Spring Harbor Laboratory
#A metabolic switch that boosts immune cell fitness against cancer
#3' scRNA MC38 TIL analysis

# Load In Libraries
library(Rmagic)
library(plyr)
library(nichenetr)
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

# Load in the Data
Control <- Read10X(data.dir = "/Users/vyom/data/tumor_infiltrating_lymphocytes/IE01_control/outs/filtered_feature_bc_matrix/")
CPT1aKO <- Read10X(data.dir = "/Users/vyom/data/tumor_infiltrating_lymphocytes/IE01_CPT1aKO/outs/filtered_feature_bc_matrix/")
VP16PPARd <- Read10X(data.dir = "/Users/vyom/data/tumor_infiltrating_lymphocytes/IE01_VP16PPARd/outs/filtered_feature_bc_matrix/")

Control <- CreateSeuratObject(Control, project = "Control")
CPT1aKO <- CreateSeuratObject(CPT1aKO, project = "CPT1aKO")
VP16PPARd <- CreateSeuratObject(VP16PPARd, project = "VP16PPARd")

d <- merge(Control, y = c(CPT1aKO,VP16PPARd), add.cell.ids = c("Control", "CPT1aKO","VP16PPARd"), project = "Immune")

d[["percent.mt"]] <- PercentageFeatureSet(d, pattern = "mt-")
VlnPlot(d, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),pt.size = 0 ,ncol = 3)
d <- subset(d, subset = nCount_RNA > 100 & nCount_RNA < 20000 & nFeature_RNA > 100 & nFeature_RNA < 5000 & percent.mt < 8)

Data.list <- SplitObject(d, split.by = "ident")
Data.list <- Data.list[c("Control", "CPT1aKO","VP16PPARd")]
for (i in 1:length(Data.list)) {
  
  Data.list[[i]] <- SCTransform(Data.list[[i]], verbose = FALSE)
}

# Normilization
#select highly variable genes 
Data.features <- SelectIntegrationFeatures(object.list = Data.list, nfeatures = 5000)
Data.list <- PrepSCTIntegration(object.list = Data.list, anchor.features = Data.features, 
                                verbose = FALSE)
Data.anchors <- FindIntegrationAnchors(object.list = Data.list, normalization.method = "SCT", 
                                       anchor.features = Data.features, verbose = FALSE)
Immune_obj <- IntegrateData(anchorset = Data.anchors, normalization.method = "SCT", 
                            verbose = TRUE)
# Visulization and Clustering
Immune_obj <- RunPCA(Immune_obj)
VizDimLoadings(Immune_obj, dims = 1:2, reduction = "pca")

DimPlot(Immune_obj, reduction = "pca")
ElbowPlot(Immune_obj, ndims = 50, reduction = "pca")

Immune_obj <- RunUMAP(Immune_obj, dims = 1:10)
Immune_obj <- FindNeighbors(Immune_obj, dims = 1:10)
Immune_obj <- FindClusters(Immune_obj, resolution = 1)
DimPlot(Immune_obj, reduction = "umap", label = TRUE)

Idents(Immune_obj) <- Immune_obj$Cell_Type
markers <- FindAllMarkers(Immune_obj, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 1)
markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(markers,'Immutumor_Sigs_Per_Clust_anno.csv')

allmarkers <- FindAllMarkers(Immune_obj, min.pct = 0.5)

#Cluster Annotation
Immune_obj <- FindClusters(Immune_obj, resolution = 1)

new.cluster.ids <- c('CD8 T','CD8 T','GD T','CD4 T','CD4 T','CD4 T','NK Cell','B Cell','NK Cell','Mac.','Mac.','Mono','Mono','DC','Neut.')
Immune_obj[["Cell_type"]] <- Idents(Immune_obj)
names(new.cluster.ids) <- levels(Immune_obj)
Immune_obj <- RenameIdents(Immune_obj, new.cluster.ids)
Immune_obj[["Cell_type"]] <- Idents(Immune_obj)
DimPlot(Immune_obj, reduction = "umap", group.by= 'Cell_type')
DimPlot(Immune_obj, reduction = "umap", group.by= 'Type')

#proportions Control vs All
Prop_table<- prop.table(x = table(Immune_obj$Cell_type, Immune_obj$Type), margin = 2)
Prop_Table <- as.data.frame(Prop_table, row.names = NULL, optional = FALSE,
                            make.names = TRUE, stringsAsFactors = default.stringsAsFactors())
Prop_Table1 <- Prop_Table
my_levels <- c('T Cell','B Cell','NK Cell','Neutrophil','Macrophage','Monocyte','DC','Mast Cell','ILC')
Prop_Table1$Var1 <- factor(Prop_Table1$Var1,levels = my_levels)
Prop_Table1$Var2 <- factor(Prop_Table1$Var2,levels = c('Control','CPT1aKO','VP16PPARd'))
plot <- ggplot(data = Prop_Table1, aes(Var1, Freq, fill=Var2)) + geom_bar(position="dodge", stat="identity", na.rm = TRUE) +  scale_fill_manual(values = c('#1b9e77' ,'#d95f02', '#7570b3')) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), panel.background = element_rect(fill = "white", colour = "Black")) + xlab("Cell Type") + ylab("Fraction of Cells") + labs(fill = "Sample") +scale_y_continuous(expand = expansion(mult = c(0, .1)))
plot


#cluster sanity check
Idents(Immune_obj) <- Immune_obj$Cell_type
Immune_obj$Cell_Type <- Immune_obj$Cell_type
new.cluster.ids <- c('CD8 T','GD T','CD4 T','NK Cell','B Cell','Mac.','Mono','DC','Neut.')
Immune_obj$Cell_Type <- factor(x = Immune_obj$Cell_type, levels = my_levels)
DimPlot(Immune_obj, group.by = "Cell_Type", label = FALSE, pt.size=1, label.size = 1)

DotPlot_Sig <- unique(c('Cd3g','Cd3e','Cd8a','Cd4','Trac','Tcrg-C1','Lag3','Pdcd1','Havcr2','Tox','Tcf7','Gzmb','Tbx21','Ifng','Gata3','Il4','Rorc','Il17a','Il17f','Foxp3','Il10','Il2rb','Il2ra','Klrd1','Cd19','Ighm','Ighg1','Cd74','Ciita','Nrc1','Klre1','Itgam','Itgax','H2-Eb1','H2-Ab1','Arg1','Mrc1','Tgfbi','Ccr2','Vegfa','Prdx1','Clec4d','Ccl5','Cd83','Ccr7','Fcn1','Msrb1','Ly6g','Col3a1','Sparc'))
DotPlot(Immune_obj, features = DotPlot_Sig, assay = 'SCT') + labs(y= "Cell Treatment", x="") + scale_colour_distiller( palette ="RdYlBu") + scale_size(range = c(0, 1)) +
  theme(text = element_text(size=5), axis.text.x = element_text(size = 5, angle = 90, hjust = 1, vjust= .01), axis.text.y = element_text(size = 6), axis.title.x = element_blank(), axis.title.y = element_blank())

#Differential Expression
Immune_obj <- PrepSCTFindMarkers(Immune_obj)
Idents(Immune_obj) <- Immune_obj$Treatment
DE_all <- FindMarkers(Immune_obj, ident.1 = "VP16PPARd", ident.2 = "Control", test.use = "MAST", logfc.threshold = .01, min.pct = .1, assay = 'SCT')
write.csv(DE_all, paste0('All_MC38_VP16_scrna_DE.csv')) 
if (max(-log(DE_all$p_val_adj)) < 321){
  ylim_num <- max(-log(DE_all$p_val_adj))} else {
    ylim_num <- 320}
EnhancedVolcano(DE_all, lab = rownames(DE_all), x = 'avg_log2FC', y = 'p_val_adj', title = 'VP16PPARd vs Control',
                pCutoff = .05, FCcutoff = 0.25, xlim = c(min(DE_all$avg_log2FC)-.05,max(DE_all$avg_log2FC)+.05),ylim = c(0,ylim_num+5),
                subtitle = 'ALL CELLS', gridlines.minor = FALSE, gridlines.major = FALSE)
ggsave(file = paste0('All_DE_peak_volcano_SCRNA.pdf'), width=6, height=6, units="in")


Cell_Types <- c(levels(as.factor(Immune_obj$Cell_Type)))
Idents(Immune_obj) <- Immune_obj$Cell_Type
for(i in Cell_Types){
  subset_cell <- subset(Immune_obj,  idents = i)
  subset_cell <- PrepSCTFindMarkers(subset_cell, verbose = TRUE)
  Idents(subset_cell) <- subset_cell$Treatment
  DE_subset <- FindMarkers(subset_cell, ident.1 = "VP16PPARd", ident.2 = "Control",  test.use = "MAST", logfc.threshold = .1, min.pct = .01, assay = 'SCT')
  write.csv(DE_subset, paste0(i,'_MC38_VP16_scrna_DE.csv')) 
  if (max(-log(DE_subset$p_val_adj)) < 321){
    ylim_num <- max(-log(DE_subset$p_val_adj))} else {
      ylim_num <- 320}
  EnhancedVolcano(DE_subset, lab = rownames(DE_subset), x = 'avg_log2FC', y = 'p_val_adj', title = 'VP16PPARd vs Control',
                  pCutoff = .05, FCcutoff = 0.1, xlim = c(min(DE_subset$avg_log2FC)-.05,max(DE_subset$avg_log2FC)+.05),ylim = c(0,ylim_num+5),
                  subtitle = i, gridlines.minor = FALSE, gridlines.major = FALSE)
  ggsave(file = paste0(i,'_DE_peak_volcano_SCRNA.pdf'), width=6, height=6, units="in")
}

#violin plots:
use_python("/Users/vyom/miniconda3/bin/python")
Immune_obj <- magic(Immune_obj)

gene.list <- unique(c('Cd8a','Ifng','Prf1',  'Gzmk','Cxcr3', 'Ifi208',"Irf7",'Irf1'))

gene.list <- unique(c('Cxcr3','Xcr1','Notch1',  'Ccr1','Il10ra'))

DefaultAssay(Immune_obj) <- "MAGIC_SCT"
selected_cells <- names(Immune_obj$Cell_Type[Immune_obj$Cell_Type %in% c("CD8 T")])
vln_data <- FetchData(Immune_obj,
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


mean.exp <- zscore(colMeans(x = Immune_obj@assays$MAGIC_SCT@data[IFNG_Response, ], na.rm = TRUE), dist = 'norm')
if (all(names(x = mean.exp) == rownames(x = Immune_obj@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  Immune_obj@meta.data$IFNG_Response <- mean.exp
}

mean.exp <- zscore(colMeans(x = Immune_obj@assays$MAGIC_SCT@data[Ifna_Response, ], na.rm = TRUE), dist = 'norm')
if (all(names(x = mean.exp) == rownames(x = Immune_obj@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  Immune_obj@meta.data$Ifna_Response <- mean.exp
}

mean.exp <- zscore(colMeans(x = Immune_obj@assays$MAGIC_SCT@data[Inflammatory_Response, ], na.rm = TRUE), dist = 'norm')
if (all(names(x = mean.exp) == rownames(x = Immune_obj@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  Immune_obj@meta.data$Inflammatory_Response <- mean.exp
}

mean.exp <- zscore(colMeans(x = Immune_obj@assays$MAGIC_SCT$data[PPAR_targets, ], na.rm = TRUE), dist = 'norm')
if (all(names(x = mean.exp) == rownames(x = Immune_obj@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  Immune_obj@meta.data$PPAR_targets <- mean.exp
}


score.list <- c('Ifn_all', 'Ifna_Response', 'IFNG_Response', 'MGIe_cleanup')

for(i in gene.list) {
  DefaultAssay(Immune_obj) <- 'MAGIC_SCT'
  selected_cells <- names(Immune_obj$Cell_Type)
  vln_data <- FetchData(Immune_obj,
                        vars = c(i,"Treatment", "Cell_Type"),
                        cells = selected_cells,
                        layer = "data")
  vln_data$Treatment
  vln_data <- melt(vln_data)
  
  All <- VlnPlot(Immune_obj, split.by = "Treatment",group.by = 'Cell_Type', features = i, pt.size = 0, assay = "MAGIC_SCT",  cols = c('#1b9e77' ,'#d95f02'), log = TRUE, split.plot = TRUE) + 
    theme(legend.position = 'none') + 
    geom_boxplot(width=0.3,outlier.shape = NA, coef = 0, lwd= .2) + xlab('') + 
    theme(text = element_text(size=9), axis.text.x = element_text(size = 9), axis.text.y = element_text(size =9)) +
    stat_compare_means(data= vln_data, aes(x = Cell_Type, y = value, fill = Treatment), method = "wilcox.test",  label = "p.signif", size = 4, hide.ns = TRUE) #+ scale_y_continuous(expand = expansion(mult = c(0, .1)))
  All$layers[[1]]$aes_params$size = .15
  All
  ggsave(file = paste0('MC38_VP16_vln_', i, '.pdf'), plot=All, width=2.5, height=2.5, units="in")
}

DefaultAssay(Immune_obj) <- 'MAGIC_SCT'

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
  
  All_Genes <- rownames(Immune_obj@assays$MAGIC_SCT)
  MHC.score <- intersect(All_Genes, MHC.score)
  survival.score <- intersect(All_Genes, survival.score)   
  Chemokine.score <- intersect(All_Genes, Chemokine.score)  
  Interferon.score <- intersect(All_Genes, Interferon.score)
  Metabolism.score <- intersect(All_Genes, Metabolism.score)
  PPAR.score <- intersect(All_Genes, PPAR.score)
  
  Immune_obj <- AddModuleScore(object = Immune_obj, features = list(MHC.score), name = 'MHC', assay = 'MAGIC_SCT')
  Immune_obj <- AddModuleScore(object = Immune_obj, features = list(survival.score), name = 'survival', assay = 'MAGIC_SCT')
  Immune_obj <- AddModuleScore(object = Immune_obj, features = list(Chemokine.score), name = 'Chemokine', assay = 'MAGIC_SCT')
  Immune_obj <- AddModuleScore(object = Immune_obj, features = list(Interferon.score), name = 'Interferon', assay = 'MAGIC_SCT')
  Immune_obj <- AddModuleScore(object = Immune_obj, features = list(Metabolism.score), name = 'Metabolism', assay = 'MAGIC_SCT')
  Immune_obj <- AddModuleScore(object = Immune_obj, features = list(PPAR.score), name = 'PPAR', assay = 'MAGIC_SCT')
  
  Idents(Immune_obj) <- Immune_obj$Cell_Type
  Cell_Types <- levels(Immune_obj$Cell_Type)
  MHC_FC <- list()
  survival_FC <- list()
  Meta_FC <- list()
  Chemokine_FC <- list()
  Interferon_FC <- list() 
  PPAR_FC <- list() 
  
  for(i in Cell_Types){
    MHC_FC[[i]] <- FoldChange(Immune_obj,ident.1 = "VP16PPARd", ident.2 = "Control",features = MHC.score, group.by = 'Treatment', subset.ident = i)
  }
  MHC_FC_aggr <- c()
  for(i in Cell_Types){
    MHC_FC_aggr <- c(MHC_FC_aggr, mean(as.numeric(MHC_FC[[i]]$avg_log2FC)))
  }
  
  for(i in Cell_Types){
    survival_FC[[i]] <- FoldChange(Immune_obj,ident.1 = "VP16PPARd", ident.2 = "Control",features = survival.score, group.by = 'Treatment', subset.ident = i)
  }
  survival_FC_aggr <- c()
  for(i in Cell_Types){
    survival_FC_aggr <- c(survival_FC_aggr, mean(as.numeric(survival_FC[[i]]$avg_log2FC)))
  }
  
  for(i in Cell_Types){
    Meta_FC[[i]] <- FoldChange(Immune_obj,ident.1 = "VP16PPARd", ident.2 = "Control",features = Metabolism.score, group.by = 'Treatment', subset.ident = i)
  }
  Meta_FC_aggr <- c()
  for(i in Cell_Types){
    Meta_FC_aggr <- c(Meta_FC_aggr, mean(as.numeric(Meta_FC[[i]]$avg_log2FC)))
  }
  
  for(i in Cell_Types){
    Chemokine_FC[[i]] <- FoldChange(Immune_obj,ident.1 = "VP16PPARd", ident.2 = "Control",features = Chemokine.score, group.by = 'Treatment', subset.ident = i)
  }
  Chemokine_FC_aggr <- c()
  for(i in Cell_Types){
    Chemokine_FC_aggr <- c(Chemokine_FC_aggr, mean(as.numeric(Chemokine_FC[[i]]$avg_log2FC)))
  }
  for(i in Cell_Types){
    Interferon_FC[[i]] <- FoldChange(Immune_obj,ident.1 = "VP16PPARd", ident.2 = "Control",features = Interferon.score, group.by = 'Treatment', subset.ident = i)
  }
  Interferon_FC_aggr <- c()
  for(i in Cell_Types){
    Interferon_FC_aggr <- c(Interferon_FC_aggr, mean(as.numeric(Interferon_FC[[i]]$avg_log2FC)))
  }
  
  for(i in Cell_Types){
    PPAR_FC[[i]] <- FoldChange(Immune_obj,ident.1 = "VP16PPARd", ident.2 = "Control",features = PPAR.score, group.by = 'Treatment', subset.ident = i)
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
  Idents(Immune_obj) <- Immune_obj$Cell_Type
  pvals <- c()
  for(j in scores){
    for(i in Cell_Types){
      selected_cells <- WhichCells(object = Immune_obj ,idents = c(i))
      vln_data <- FetchData(Immune_obj,
                            vars = c(j,"Treatment", "Cell_Type"),
                            cells = selected_cells,
                            slot = "data")
      vln_data <- melt(vln_data)
      vln_data <- vln_data[c('Treatment', 'Cell_Type', 'value')]
      
      pvals <- c(pvals, t.test(vln_data[which(vln_data$Treatment == 'VP16PPARd'),]$value, vln_data[which(vln_data$Treatment == 'Control'),]$value,  alternative = c("greater"))$p.value)
      print(i)
    }
  }
  pvals <- p.adjust(pvals, method = "BH", n = length(pvals))
  Module_Heatmap1$pval <- -log(pvals) 
  Cell_Types
  my_levels <- c('CD8 T','GD T','CD4 T','NK Cell','B Cell','Mac.','Mono','DC','Neut.')
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
write.csv(as.data.frame(Immune_obj@meta.data), 'MC38_metadata.csv')


