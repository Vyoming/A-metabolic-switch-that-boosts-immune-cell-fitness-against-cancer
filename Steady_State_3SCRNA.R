#Vyom Shah
#Cold Spring Harbor Laboratory
#A metabolic switch that boosts immune cell fitness against cancer
#3' scRNA Steady state spleen analysis

# Load In Libraries
library(Rmagic)
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
Control <- Read10X(data.dir = "/Users/vyom/data/steady_state_vp16/Beyaz_IE02_Control/outs/filtered_feature_bc_matrix/")
VP16PPARd <- Read10X(data.dir = "/Users/vyom/data/steady_state_vp16/Beyaz_IE02_VB16EXP/outs/filtered_feature_bc_matrix/")

Control <- CreateSeuratObject(Control, project = "Control")

VP16PPARd <- CreateSeuratObject(VP16PPARd, project = "VP16PPARd")

d <- merge(Control, y = c(VP16PPARd), add.cell.ids = c("Control", "VP16PPARd"), project = "Immune")
d[["percent.mt"]] <- PercentageFeatureSet(d, pattern = "mt-")
VlnPlot(d, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),pt.size = 0 ,ncol = 3)
d <- subset(d, subset = nCount_RNA > 500 & nCount_RNA < 20000 & nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)

d$orig.ident
Data.list <- SplitObject(d, split.by = "ident")
Data.list <- Data.list[c("Control", "VP16PPARd")]
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
steady_obj <- IntegrateData(anchorset = Data.anchors, normalization.method = "SCT", 
                            verbose = TRUE)
# Visulization and Clustering
steady_obj <- RunPCA(steady_obj)
VizDimLoadings(steady_obj, dims = 1:2, reduction = "pca")

DimPlot(steady_obj, reduction = "pca")
ElbowPlot(steady_obj, ndims = 50, reduction = "pca")

steady_obj <- RunUMAP(steady_obj, dims = 1:10)

Idents(steady_obj) <- levels(factor(steady_obj@meta.data$orig.ident))
steady_obj[["Type"]] <- steady_obj$orig.ident
steady_obj$Type
steady_obj@meta.data$Type <- factor(steady_obj@meta.data$Type, levels = c("Control",  "VP16PPARd")) 

steady_obj <- FindNeighbors(steady_obj, dims = 1:10)
steady_obj <- FindClusters(steady_obj, resolution = 2)
DimPlot(steady_obj, reduction = "umap", label = TRUE)
markers <- FindAllMarkers(steady_obj, only.pos = TRUE, min.pct = 0.1, logfc.threshold = .5)
write.csv(markers,'Gene_de_Sigs_Per_Clust_immune.csv')

Prop_table<- prop.table(x = table(steady_obj$Cell_Type, steady_obj$Type), margin = 2)
Prop_Table <- as.data.frame(Prop_table, row.names = NULL, optional = FALSE,
                            make.names = TRUE, stringsAsFactors = default.stringsAsFactors())
Prop_Table1 <- Prop_Table
Prop_Table1$Var1 <- factor(Prop_Table1$Var1)
Prop_Table1$Var2 <- factor(Prop_Table1$Var2,levels = c('Control','VP16PPARd'))
plot <- ggplot(data = Prop_Table1, aes(Var1, Freq, fill=Var2)) + geom_bar(position="dodge", stat="identity", na.rm = TRUE) +  scale_fill_manual(values = c('#1b9e77' ,'#d95f02')) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), panel.background = element_rect(fill = "white", colour = "Black")) + xlab("Cell Type") + ylab("Fraction of Cells") + labs(fill = "Sample") +scale_y_continuous(expand = expansion(mult = c(0, .1)))
plot
DimPlot(steady_obj)
DotPlot_Sig <- unique(c('Cd3g','Cd3e','Cd8a','Cd4','Trac','Tcrg-C1','Lag3','Pdcd1','Havcr2','Tox','Tcf7','Gzmb','Tbx21','Ifng','Gata3','Il5','Il13','Rorc','Il17a','Il17f','Foxp3','Il10','Il2rb','Il2ra','Klrd1','Cd19','Ighm','Ighg1','Cd74','Ciita','Nrc1','Klre1','Itgam','Itgax','H2-Eb1','H2-Ab1','Arg1','Mrc1','Tgfbi','Ccr2','Vegfa','Prdx1','Clec4d','Ccl5','Cd83','Ccr7','Fcn1','Msrb1','Ly6g','Col3a1','Sparc'))
DotPlot(steady_obj, features = DotPlot_Sig, assay = 'SCT') + labs(y= "Cell type", x="") + scale_colour_distiller( palette ="RdYlBu") + scale_size(range = c(0, 1)) +
  theme(text = element_text(size=5), axis.text.x = element_text(size = 5, angle = 90, hjust = 1, vjust= .01), axis.text.y = element_text(size = 6), axis.title.x = element_blank(), axis.title.y = element_blank())
ggsave(file="steady_dotplot.pdf",  width=4.25, height=2)

new.cluster.ids <- c('B Cell', 'CD4+ T Cell', 'CD8+ T Cell', 'CD4+ T Cell', 'CD8+ T Cell', 'Neutrophil', 'B Cell', 'CD4+ T Cell', 'B Cell', 'NK Cell', 'CD4+ T Cell', 'Neutrophil', 'Monocytes', 'NKT Cell', 'NKT Cell', 'Treg', 'AB T Lymphocytes', 'B Cell', 'NKT Cell', 'Neutrophil', 'B Cell', 'M1 Macrophage', 'M2 Macrophage', 'Active NKT Cell', 'Stromal', 'Stromal', 'Stromal', 'CD4+ T Cell', 'CD8+ T Cell')           
steady_obj[["Cell_type"]] <- Idents(steady_obj)
names(new.cluster.ids) <- levels(steady_obj)
steady_obj <- RenameIdents(steady_obj, new.cluster.ids)
steady_obj[["Cell_type"]] <- Idents(steady_obj)

DimPlot(steady_obj, reduction = "umap", label = TRUE)


my_levels <- c('CD8+ T Cell','AB T Lymphocytes','CD4+ T Cell','Treg','NKT Cell','Active NKT Cell','NK Cell','B Cell','M1 Macrophage','M2 Macrophage','Monocytes','Neutrophil','Stromal')
steady_obj$Cell_Type <- factor(x = steady_obj$Cell_type, levels = my_levels)
Idents(steady_obj) <- steady_obj$Cell_Type
DotPlot_Sig <- unique(c('Cd3g','Cd3e','Cd8a','Cd4','Trac','Tcrg-C1','Lag3','Pdcd1','Havcr2','Tox','Tcf7','Gzmb','Tbx21','Ifng','Gata3','Il5','Il13','Rorc','Il17a','Il17f','Foxp3','Il10','Il2rb','Il2ra','Klrd1','Cd19','Ighm','Ighg1','Cd74','Ciita','Nrc1','Klre1','Itgam','Itgax','H2-Eb1','H2-Ab1','Arg1','Mrc1','Tgfbi','Ccr2','Vegfa','Prdx1','Clec4d','Ccl5','Cd83','Ccr7','Fcn1','Msrb1','Ly6g','Col3a1','Sparc'))
DotPlot(steady_obj, features = DotPlot_Sig, assay = 'SCT') + labs(y= "Cell type", x="") + scale_colour_distiller( palette ="RdYlBu") + scale_size(range = c(0, 5)) +
  theme(text = element_text(size=5), axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust= .01), axis.text.y = element_text(size = 6), axis.title.x = element_blank(), axis.title.y = element_blank())

#Violin Plots
for(i in Inflammation_Score) {
  selected_cells <- names(steady_obj$Cell_Type)
  vln_data <- FetchData(steady_obj,
                        vars = c(i,"Type", "Cell_Type"),
                        cells = selected_cells,
                        slot = "data")
  vln_data$Type
  vln_data <- melt(vln_data)
  
  All <- VlnPlot(steady_obj, split.by = "Type", features = i, pt.size = 0, assay = "MAGIC_RNA",  cols = c('#1b9e77' ,'#d95f02'), log = FALSE, split.plot = TRUE) + 
    theme(legend.position = 'none') + 
    geom_boxplot(width=0.3,outlier.shape = NA, coef = 0, lwd= .2) + xlab('') + 
    theme(text = element_text(size=7), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6)) +
    stat_compare_means(data= vln_data, aes(x = Cell_Type, y = value, fill = Type), method = "wilcox.test",  label = "p.signif", size = 2, hide.ns = TRUE) #+ scale_y_continuous(expand = expansion(mult = c(0, .1)))
  All$layers[[1]]$aes_params$size = .15
  All
  ggsave(file = paste0('Immune_VP16_', i, '.pdf'), plot=All, width=3, height=3, units="in")
}

write.csv(as.data.frame(steady_obj@meta.data), 'steady_state_metadata.csv')