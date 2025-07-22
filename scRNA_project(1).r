rm(list=ls())
options(stringsAsFactors = F)
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)
# remotes::install_github('chris-mcginnis-ucsf/DoubletFinder') 
# remotes::install_github('immunogenomics/harmony')


dir='GSE147528_RAW' 
samples=list.files( dir )
samples 
sceList = lapply(samples,function(pro){ 
  # pro=samples[1] 
  print(pro)  
  
  sce =CreateSeuratObject(counts =  Read10X_h5(file.path(dir,pro)) ,
                          project =   gsub('_raw_gene_bc_matrices_h5.h5','',gsub('^GSM[0-9]*_','',pro) ) ,
                          min.cells = 5,
                          min.features = 300 )
  print(sce)
  return(sce)
})
names(sceList)  
# gsub('^GSM[0-9]*','',samples)
sce.all=merge(x=sceList[[1]],
              y=sceList[ -1 ],
              add.cell.ids =  gsub('_raw_gene_bc_matrices_h5.h5','',gsub('^GSM[0-9]*_','',samples) )     )

as.data.frame(sce.all@assays$RNA@counts[1:10, 1:2])
head(sce.all@meta.data, 10)
table(sce.all$orig.ident) 

sce.all$Tissue=NA
sce.all$Tissue[sce.all$orig.ident %in% c("EC1" , "EC10" , "EC2" ,"EC3", "EC4", "EC5" , "EC6" ,"EC7", "EC8" ,"EC9")]="EC"
sce.all$Tissue[sce.all$orig.ident %in% c("SFG1", "SFG10" , "SFG2",  "SFG3" , "SFG4" , "SFG5",  "SFG6" , "SFG7" , "SFG8",  "SFG9")]="SFG"
table(sce.all$Tissue) 



#rm(list=ls())
###### step2:QC质控 ######
dir.create("./1-QC")
setwd("./1-QC")

sce.all <- SetIdent(sce.all, value = "GSE147528")

# sce.all=readRDS("../sce.all_raw.rds")
#Calculate the proportion of mitochondrial genes
#The genetic names of humans and mice are slightly different
mito_genes=rownames(sce.all)[grep("^MT-", rownames(sce.all))] 
mito_genes
sce.all=PercentageFeatureSet(sce.all, "^MT-", col.name = "percent_mito")
fivenum(sce.all@meta.data$percent_mito)
#Calculate the proportion of ribosomal genes
ribo_genes=rownames(sce.all)[grep("^Rp[sl]", rownames(sce.all),ignore.case = T)]
ribo_genes
sce.all=PercentageFeatureSet(sce.all, "^RP[SL]", col.name = "percent_ribo")
fivenum(sce.all@meta.data$percent_ribo)
#Calculate the proportion of red blood cell genes
rownames(sce.all)[grep("^Hb[^(p)]", rownames(sce.all),ignore.case = T)]
sce.all=PercentageFeatureSet(sce.all, "^HB[^(P)]", col.name = "percent_hb")
fivenum(sce.all@meta.data$percent_hb)
#Visualize the above proportion of cells
feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "percent_hb")
feats <- c("nFeature_RNA", "nCount_RNA")
p1=VlnPlot(sce.all,  features = feats, pt.size = 0.01, ncol = 2) + 
    NoLegend()
p1
library(ggplot2)
ggsave(filename="Vlnplot1.pdf",plot=p1,height = 6,width = 14)
ggsave(filename="Vlnplot1.png",plot=p1,height = 6,width = 14)
feats <- c("percent_mito", "percent_ribo", "percent_hb")
p2=VlnPlot(sce.all, features = feats, pt.size = 0.01, ncol = 3, same.y.lims=T) + 
    scale_y_continuous(breaks=seq(0, 100, 5)) +
	NoLegend()
p2	
ggsave(filename="Vlnplot2.pdf",plot=p2,height = 6,width = 18)

p3=FeatureScatter(sce.all, "nCount_RNA", "nFeature_RNA",  pt.size = 0.5,raster=FALSE)
ggsave(filename="Scatterplot.pdf",plot=p3)
#Filter low-quality cells/genes based on the above indicators
#Filter indicator 1: Cells with the least number of expressed genes&Genes with the least number of expressed cells
selected_c <- WhichCells(sce.all, expression = nFeature_RNA < 5000)

selected_f <- rownames(sce.all)[Matrix::rowSums(sce.all@assays$RNA@counts > 0 ) > 3]

sce.all.filt <- subset(sce.all, features = selected_f, cells = selected_c)
dim(sce.all) 
dim(sce.all.filt) 

# par(mar = c(4, 8, 2, 1))
C=sce.all.filt@assays$RNA@counts
dim(C)
C=Matrix::t(Matrix::t(C)/Matrix::colSums(C)) * 100
C=C[,sample(1:ncol(C),1000)]
most_expressed <- order(apply(C, 1, median), decreasing = T)[50:1]
pdf("TOP50_most_expressed_gene.pdf",width=14)
boxplot(as.matrix(Matrix::t(C[most_expressed, ])),
        cex = 0.1, las = 1, 
        xlab = "% total count per cell", 
    col = (scales::hue_pal())(50)[50:1], 
    horizontal = TRUE)
dev.off()
rm(C)

# #Filter indicator 2: Ratio of mitochondrial/ribosomal genes (according to the violin chart above)
# selected_mito <- WhichCells(sce.all.filt, expression = percent_mito < 15)
# selected_ribo <- WhichCells(sce.all.filt, expression = percent_ribo > 3)
# selected_hb <- WhichCells(sce.all.filt, expression = percent_hb < 0.1)
# length(selected_hb)
# length(selected_ribo)
# length(selected_mito)
# 
# 
# sce.all.filt <- subset(sce.all.filt, cells = selected_mito)
# sce.all.filt <- subset(sce.all.filt, cells = selected_ribo)
# sce.all.filt <- subset(sce.all.filt, cells = selected_hb)
dim(sce.all.filt)
dim(sce.all)
table(sce.all.filt$orig.ident) 
feats <- c("nFeature_RNA", "nCount_RNA")
p1_filtered=VlnPlot(sce.all.filt,features = feats, pt.size = 0.1, ncol = 2) + 
    NoLegend()
ggsave(filename="Vlnplot1_filtered.pdf",plot=p1_filtered,height = 6,width = 14)

feats <- c("percent_mito", "percent_ribo", "percent_hb")
p2_filtered=VlnPlot(sce.all.filt,  features = feats, pt.size = 0.1, ncol = 3) + 
    NoLegend()
ggsave(filename="Vlnplot2_filtered.pdf",plot=p2_filtered,height = 6,width = 18)

#Filter metric 3: Filter specific genes
#Filter MALAT1 Butler Gene
sce.all.filt <- sce.all.filt[!grepl("MALAT1", rownames(sce.all.filt),ignore.case = T), ]
# Filter Mitocondrial
sce.all.filt <- sce.all.filt[!grepl("^MT-", rownames(sce.all.filt),ignore.case = T), ]

dim(sce.all.filt)
# [1] 29053 67202

#Cell cycle score
sce.all.filt = NormalizeData(sce.all.filt)
s.genes=Seurat::cc.genes.updated.2019$s.genes
g2m.genes=Seurat::cc.genes.updated.2019$g2m.genes
sce.all.filt=CellCycleScoring(object = sce.all.filt, 
                              s.features = s.genes, 
                              g2m.features = g2m.genes, 
                              set.ident = TRUE)
p4=VlnPlot(sce.all.filt, features = c("S.Score", "G2M.Score"), group.by = "orig.ident", 
    ncol = 2, pt.size = 0.1)
ggsave(filename="Vlnplot4_cycle.pdf",plot=p4)

sce.all.filt@meta.data  %>% ggplot(aes(S.Score,G2M.Score))+geom_point(aes(color=Phase))+
  theme_minimal()
ggsave(filename="cycle_details.pdf" )

dim(sce.all.filt) 
dim(sce.all) 

# saveRDS(sce.all.filt, "sce.all_qc.rds")


setwd('../')
# rm(list=ls()) 

dir.create("2-harmony")
getwd()
setwd("2-harmony")

# sce.all=readRDS("../1-QC/sce.all_qc.rds")
sce=sce.all.filt 
sce
sce <- NormalizeData(sce, 
                     normalization.method = "LogNormalize",
                     scale.factor = 1e4) 
sce <- FindVariableFeatures(sce)
sce <- ScaleData(sce)
sce <- RunPCA(sce, features = VariableFeatures(object = sce))

library(harmony)
seuratObj <- RunHarmony(sce, "orig.ident")
names(seuratObj@reductions)
seuratObj <- RunUMAP(seuratObj,  dims = 1:15, 
                     reduction = "harmony")
DimPlot(seuratObj,reduction = "umap",label=T,raster=FALSE ) 
DimPlot(seuratObj,reduction = "harmony",label=T ,raster=FALSE ) 
sce=seuratObj
sce <- FindNeighbors(sce, reduction = "harmony",
                     dims = 1:15) 
sce.all=sce
for (res in c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5,0.8,1)) {
  sce.all=FindClusters(sce.all, #graph.name = "CCA_snn", 
                       resolution = res, algorithm = 1)
}
colnames(sce.all@meta.data)
apply(sce.all@meta.data[,grep("RNA_snn",colnames(sce.all@meta.data))],2,table)

p1_dim=plot_grid(ncol = 3, DimPlot(sce.all, reduction = "umap", group.by = "RNA_snn_res.0.01",raster=FALSE ) + 
                   ggtitle("louvain_0.01"), DimPlot(sce.all, reduction = "umap", group.by = "RNA_snn_res.0.1",raster=FALSE ) + 
                   ggtitle("louvain_0.1"), DimPlot(sce.all, reduction = "umap", group.by = "RNA_snn_res.0.2",raster=FALSE ) + 
                   ggtitle("louvain_0.2"))
ggsave(plot=p1_dim, filename="Dimplot_diff_resolution_low.pdf",width = 14)

p1_dim=plot_grid(ncol = 3, DimPlot(sce.all, reduction = "umap", group.by = "RNA_snn_res.0.2",raster=FALSE ) + 
                   ggtitle("louvain_0.8"), DimPlot(sce.all, reduction = "umap", group.by = "RNA_snn_res.1",raster=FALSE ) + 
                   ggtitle("louvain_1"), DimPlot(sce.all, reduction = "umap", group.by = "RNA_snn_res.0.3",raster=FALSE ) + 
                   ggtitle("louvain_0.3"))
ggsave(plot=p1_dim, filename="Dimplot_diff_resolution_high.pdf",width = 18)


p2_tree=clustree(sce.all@meta.data, prefix = "RNA_snn_res.")
ggsave(plot=p2_tree, filename="Tree_diff_resolution.pdf",width = 10,height = 10)


sel.clust = "RNA_snn_res.0.2"
sce.all <- SetIdent(sce.all, value = sel.clust)
table(sce.all@active.ident) 


setwd('../')
dir.create("3-cell")
setwd("3-cell")  

DimPlot(sce.all, reduction = "umap", group.by = "seurat_clusters",label = T,raster=FALSE ) 
DimPlot(sce.all, reduction = "umap", group.by = "RNA_snn_res.0.2",label = T,raster=FALSE ) 
ggsave('1.umap_by_RNA_snn_res.0.2.pdf')



# T Cells (CD3D, CD3E, CD8A), 
# B cells (CD19, CD79A, MS4A1 [CD20]), 
# Plasma cells (IGHG1, MZB1, SDC1, CD79A), 
# Monocytes and macrophages (CD68, CD163, CD14),
# NK Cells (FGFBP2, FCG3RA, CX3CR1),  
# Photoreceptor cells (RCVRN), 
# Fibroblasts (FGF7, MME), 
# Endothelial cells (PECAM1, VWF). 
# epi or tumor (EPCAM, KRT19, PROM1, ALDH1A1, CD24).
# immune (CD45+,PTPRC), epithelial/cancer (EpCAM+,EPCAM), 
# stromal (CD10+,MME,fibo or CD31+,PECAM1,endo) 
# DC cells ("THBD","IL3RA"  ,"ITGAX" ,"CD1C" )
library(ggplot2)
genes_to_check = c("SLC17A7","CAMK2A", "GAD1", "GAD2","MBP","MOG" , "CLDN5","FLT1", "PDGFRA","SOX10" , "AQP4", "SLC1A2","CD74","CX3CR1" )
library(stringr)
p <- DotPlot(sce.all, features = genes_to_check,
             assay='RNA'  )  + coord_flip()

p
# ggsave(plot=p, filename="check_marker_by_seurat_cluster.pdf")


library(SingleR)
library(celldex)

DefaultAssay(sce.all) <- "RNA"
sce.all
DimPlot(sce.all, reduction = "umap", label = TRUE, repel = TRUE,raster=FALSE )

immune_singler <- GetAssayData(sce.all,slot = "data") 

clusters <- sce.all$RNA_snn_res.0.2
table(sce.all$RNA_snn_res.0.2)



load("pathway_to_HumanPrimaryCellAtlasData.Rdata")
pred.immune <- SingleR(test = immune_singler,
                       ref = HumanPrimaryCellAtlasData, 
                       labels = HumanPrimaryCellAtlasData$label.fine,
                       method = "cluster", 
                       clusters = clusters)
table(pred.immune$labels)
plotScoreHeatmap(pred.immune, clusters = pred.immune$labels)

new.clusterID <- pred.immune$labels
new.clusterID

names(new.clusterID) <- levels(sce.all)
new.clusterID
table()

new.clusterID[0+1]="Oligodendrocytes"
new.clusterID[1+1]="Excitatory_neurons"
new.clusterID[2+1]="Astrocytes"
new.clusterID[3+1]="Excitatory_neurons"
new.clusterID[4+1]="Inhibitory_neurons"
new.clusterID[5+1]="Microglia"
new.clusterID[6+1]="Excitatory_neurons"
new.clusterID[7+1]="Inhibitory_neurons"
new.clusterID[8+1]="OPCs"
new.clusterID[9+1]="Excitatory_neurons"
new.clusterID[10+1]="Endothelial_cells"
new.clusterID[11+1]="Microglia"
new.clusterID[12+1]="Microglia"




new.clusterID


sce.all <- RenameIdents(sce.all,new.clusterID)
table(sce.all$seurat_clusters)
DimPlot(sce.all, reduction = "umap", label = TRUE, repel = TRUE,raster = F)

# celltype=data.frame(ClusterID=0:17,
#                     celltype='na')
# celltype[celltype$ClusterID %in% c(0, 1, 2,3,4,5,6,7,8,9,10,12,15,16),2]='Malignant'
# celltype[celltype$ClusterID %in% c(17),2]='Endothelial cells'
# celltype[celltype$ClusterID %in% c(14),2]='T cells'
# celltype[celltype$ClusterID %in% c(11),2]='Mono.Macr'
# celltype[celltype$ClusterID %in% c( 13),2]='Fibroblast'
# # 
# # 
# head(celltype)
# celltype
# table(celltype$celltype)
# sce.all@meta.data$celltype = "NA"
# for(i in 1:nrow(celltype)){
#   sce.all@meta.data[which(sce.all@meta.data$RNA_snn_res.0.8 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
# table(sce.all@meta.data$celltype)

sce.all@meta.data$celltype=sce.all@active.ident
sce.all@meta.data$celltype[sce.all@meta.data$celltype=="Inhibitory_neurons,"]="Inhibitory_neurons"

tab.1=table(sce.all@meta.data$celltype,sce.all@meta.data$RNA_snn_res.0.2) 
library(gplots)
tab.1
pro='cluster'
pdf(file = paste0("2." , pro,'  celltype VS  res.0.2.pdf'),width = 14)
balloonplot(tab.1, main =" celltype VS  res.0.2 ", xlab ="", ylab="",
            label = T, show.margins = F)
dev.off()

library(Seurat)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db) 

#Confirm the gene name of TIM3 (possibly HAVCR2)
#Check if the gene name is in the dataset
gene_names <- rownames(sce.all)
tim3_gene <- "HAVCR2" 
if(!tim3_gene %in% gene_names) {
  potential_names <- c("TIM3", "HAVCR2", "CD366")
  for(name in potential_names) {
    if(name %in% gene_names) {
      tim3_gene <- name
      message("Find in other names: ", tim3_gene)
      break
    }
  }
}


# Extract Microglia cells
microglia_cells <- WhichCells(sce.all, idents = "Microglia")
microglia_obj <- subset(sce.all, cells = microglia_cells)

# Define TIM3 positive threshold (adjustable based on data distribution)
tim3_expression <- GetAssayData(microglia_obj, slot = "data")[tim3_gene,]
hist(tim3_expression, breaks = 50, main = "TIM3 Expression in Microglia", xlab = "Expression Level")

#Use expression distribution to determine the threshold, here using expression value>0 as an example
#In practical applications, stricter thresholds may be required, such as the mean plus one standard deviation
threshold <- 0
tim3_positive_cells <- names(tim3_expression[tim3_expression > threshold])

#Check the proportion of TIM3 positive cells
message("Microglia_cell_total: ", length(microglia_cells))
message("TIM3_Microglia_cellcount:", length(tim3_positive_cells))
message("TIM3_Microglia_cellrate:", round(length(tim3_positive_cells)/length(microglia_cells)*100, 2), "%")






new_cell_types <- sce.all@meta.data$celltype 
names(new_cell_types) <- colnames(sce.all)


new_cell_types <- as.character(sce.all@meta.data$cell_type)
names(new_cell_types) <- colnames(sce.all)

new_cell_types[tim3_positive_cells] <- "TIM3_positive_Microglia"


# Add to metadata
sce.all$cell_type_with_tim3 <- new_cell_types
table(sce.all$cell_type_with_tim3 )
sce.all$celltype1=sce.all$celltype
sce.all$celltype=sce.all$cell_type_with_tim3
# Set a new identification label
Idents(sce.all) <- "celltype"

table(sce.all$celltype)

if(!any(grepl("umap", names(sce.all@reductions)))) {
  sce.all <- RunPCA(sce.all, npcs = 30)
  sce.all <- RunUMAP(sce.all, dims = 1:30)
}

mycolors <-c("#B07AA1","#F89C74","#66C5CC","#75ACC3","#FE88B1","#80BA5A","#F6CF71","#EDC948","#B84D64",
             "#008695","#59A14F","#FF9DA7","#F28E2B","#DCB0F2","#EE7072","#E73F74","#f69896","#9B8E8C",
             "#E68310","#E15759")
p1 <- DimPlot(sce.all, reduction = "umap", group.by = "cell_type_with_tim3", raster=FALSE,
              label = TRUE, repel = TRUE,cols = mycolors) + 
  ggtitle("Cell Types with TIM3+ Microglia")

p2 <- FeaturePlot(sce.all, features = tim3_gene, reduction = "umap",raster=FALSE) + 
  ggtitle(paste(tim3_gene, "Expression"))

microglia_subset <- subset(sce.all, cell_type_with_tim3 %in% c("Microglia", "TIM3_positive_Microglia"))
p3 <- DimPlot(microglia_subset, reduction = "umap", group.by = "cell_type_with_tim3", 
              label = TRUE, repel = TRUE,raster=FALSE,cols = mycolors) + 
  ggtitle("Microglia by TIM3 Status")


p_combined <- p1 | p2
print(p_combined)
print(p3)

pdf("TIM3_Expression.pdf",height = 8,width = 10)
print(p_combined)
dev.off()

#Set the cell identifier to a new cell type
#Perform differential expression analysis
de_results <- FindMarkers(sce.all, 
                          ident.1 = "TIM3_positive_Microglia", 
                          ident.2 = "Microglia",
                          min.pct = 0.25,
                          logfc.threshold = 0)

#View results of differentially expressed genes
head(de_results, 20)

write.csv(de_results, "TIM3_positive_vs_negative_Microglia_DEGs.csv")



phe=sce.all@meta.data
save(phe,file = 'phe_by_markers.Rdata')

unique(sce.all$celltype)
library(dplyr)
library(Seurat)
library(ComplexHeatmap)
library(GSVA)
library(GSEABase)
library(limma)
library(ggplot2)
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 

library(magrittr)

table(sce.all@meta.data$celltype)
table(sce.all@active.ident)
length(sce.all@meta.data$celltype)
length(sce.all@active.ident)

sce.all <- SetIdent(sce.all, value = "celltype")
table(sce.all@active.ident)
library(future)
options(future.globals.maxSize = 89128960000)  
# check the current active plan
plan()
plan("multicore", workers = 8)
plan()
start = Sys.time()
sce.all.markers1 <- FindAllMarkers(sce.all,only.pos = FALSE,  min.pct = 0.25, logfc.threshold = 0.25)
end = Sys.time()
dur = end-start
dur

sce.all.markers=dplyr::filter(sce.all.markers1,avg_log2FC > 0.25)


library(scRNAtoolVis)


#cluster.order :
# ajust cluster orders
colour=c("#B07AA1","#F89C74","#66C5CC","#75ACC3","#FE88B1","#80BA5A","#F6CF71","#EDC948","#B84D64",
         "#008695","#59A14F","#FF9DA7","#F28E2B","#DCB0F2","#EE7072","#E73F74","#f69896","#9B8E8C",
         "#E68310","#E15759")
p=jjVolcano(diffData = sce.all.markers1,
            size  = 3.5,
            fontface = 'italic',
            legend.position = c(0.8,0.2),tile.col = colour,
            flip = T)
ggsave(filename = "3.Marker_gene_pointplot_celltype.pdf", plot = p, width = 10, height = 6)



p=markerVocalno(markers = sce.all.markers1,
                topn = 5,log2FC=0.25,
                labelCol = colour)

ggsave(filename = "3.Marker_gene_pointplot_celltype2.pdf", plot = p, width = 16, height = 10)

sce.all.markers1=read.csv("output_sce.all.markers_celltype.csv",row.names = 1)

write.csv(sce.all.markers1, "output_sce.all.markers_celltype.csv", quote = F)



col =c("#3176B7","#F78000","#3FA116","#CE2820","#9265C1",
       "#885649","#DD76C5","#7F7F7F","#BBBE00","#41BED1")
pdf("1 cell.anno.pdf",height = 8,width = 10)
DimPlot(sce.all,reduction = "umap",label=T ,group.by = "celltype",cols = mycolors,raster=FALSE,label.box = T) 
dev.off()

pdf("2 cell.anno.pdf",height = 8,width = 10)
DimPlot(sce.all,reduction = "umap",label=T ,group.by = "celltype1",cols = mycolors,raster=FALSE,label.box = T) 
dev.off()





library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(ggsci)
library(scater)
library(SingleR)
library(hdf5r)
library(tidyverse)
library(RColorBrewer)
library(dplyr)
library(Seurat)
library(ComplexHeatmap)
library(GSVA)
library(GSEABase)
library(limma)
library(ggplot2)



p1 <-DotPlot(sce.all, features = tim3_gene, assay='RNA' )+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 0.5,vjust=0.5))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))

ggsave(print(p1),filename = "TIM3_celltype.pdf",height = 6,width =5)
dev.off()










save(sce.all,file = "GSE147528.RData")
setwd('../')



options(stringsAsFactors = F)
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)


if(T){
  
  pro = 'cosg_celltype_'
  library(COSG)
  marker_cosg <- cosg(
    sce.all,
    groups='all',
    assay='RNA',
    slot='data',
    mu=1,
    n_genes_user=100)
  
  save(marker_cosg,file = paste0(pro,'_marker_cosg.Rdata'))
  
  sce = sce.all
  ## Top10 genes
  library(dplyr)  
  top_10 <- unique(as.character(apply(marker_cosg$names,2,head,10)))
  # width <-0.006*dim(sce)[2];width
  # height <- 0.25*length(top_10)+4.5;height
  
  width <- 15+0.5*length(unique(Idents(sce)));width
  height <- 8+0.1*length(top_10);height
  
  DoHeatmap( subset(sce,downsample=100), top_10 , 
             size=3)
  
  ggsave(filename=paste0(pro,'DoHeatmap_check_top10_markers_by_clusters.pdf') ,
         # limitsize = FALSE,
         units = "cm",width=width,height=height)
  width <- 15+0.6*length(unique(Idents(sce)));width
  height <- 8+0.2*length(top_10);height
  DotPlot(sce, features = top_10 ,
          assay='RNA'  )  + coord_flip() +FontSize(y.text = 6)
  ggsave(paste0(pro,'DotPlot_check_top10_markers_by_clusters.pdf'),
         units = "cm",width=width,height=height)
  
  
  ## Top3 genes
  top_3 <- unique(as.character(apply(marker_cosg$names,2,head,3)))
  
  width <- 15+0.2*length(unique(Idents(sce)));width
  height <- 8+0.1*length(top_3);height
  
  DoHeatmap( subset(sce,downsample=100), top_3 ,
             size=3)
  ggsave(filename=paste0(pro,'DoHeatmap_check_top3_markers_by_clusters.pdf') ,
         units = "cm",width=width,height=height)
  
  width <- 8+0.2*length(unique(Idents(sce)));width
  height <- 8+0.1*length(top_3);height
  DotPlot(sce, features = top_3 ,
          assay='RNA'  )  + coord_flip()
  ggsave(paste0(pro,'DotPlot_check_top3_markers_by_clusters.pdf'),width=width,height=height)
  
  
}


getwd()
pro = 'cosg_celltype_'
load(file = paste0(pro,'_marker_cosg.Rdata'))

top_10 <- unique(as.character(apply(marker_cosg$names,2,head,10))) 
sce.Scale <- ScaleData(subset(sce.all,downsample=100),features =  top_10  )  
table(sce.Scale$celltype)
library(paletteer) 
color <- c(paletteer_d("awtools::bpalette"),
           paletteer_d("awtools::a_palette"),
           paletteer_d("awtools::mpalette"))

table(Idents(sce.Scale)) 
ord=c("Oligodendrocytes" ,     "Excitatory_neurons" ,   "Microglia",  "Inhibitory_neurons" ,     "OPCs", 
      "Astrocytes" ,"TIM3_positive_Microglia"    ,   "Endothelial_cells" )


sce.Scale$celltype = factor(sce.Scale$celltype ,levels = ord)
ll =  as.list(as.data.frame(apply(marker_cosg$names,2,head,10)))
ll = ll[ord]
rmg=names(table(unlist(ll))[table(unlist(ll))>1])
ll = lapply(ll, function(x) x[!x %in% rmg])
ll
DoHeatmap(sce.Scale,
          features = unlist(ll),
          group.by = "celltype",
          assay = 'RNA', label = T)+
  scale_fill_gradientn(colors = c("white","grey","firebrick3"))

ggsave(filename = "fig4-top10-marker_pheatmap.pdf",units = "cm",width = 25,height = 25)

library(ggsci)


# load(file = '../3-cell/phe-by-markers.Rdata')
phe = sce.all@meta.data
head(phe)
table(phe$group,phe$orig.ident)
cal_table = function(x,y,prefix ){
  # x = phe$orig.ident
  # y = phe$celltype
  library(sur)
  library(reshape2)
  tbl =  table(x,y)
  pdf(paste0(prefix,'-table.pdf'),width = 10,height = 10)
  gplots::balloonplot( tbl )
  dev.off() 
  df = dcast(as.data.frame(tbl),x~y)
  head(df)
  write.csv(  df ,paste0(prefix,'-table.csv'))
  
  # ptb = round(sur::percent.table(x,y),2)
  ptb = round(100*tbl/rowSums(df[,-1]),2)
  
  pdf(paste0(prefix,'-percent-table.pdf'),width = 10,height = 10)
  gplots::balloonplot( ptb )
  dev.off()
  write.csv(  dcast(as.data.frame(ptb),x~y) ,paste0(prefix,'-percent-table.csv')) 
  
}
cal_table(phe$orig.ident,phe$celltype,prefix = 'celltype-vs-orig.ident')



phe=sce.all@meta.data
library(tidyr)
library(reshape2) 
colnames(phe)
tb=table(phe$celltype,
         phe$orig.ident)
head(tb)
library (gplots) 
balloonplot(tb)
bar_data <- as.data.frame(tb)


bar_per <- bar_data %>% 
  group_by(Var2) %>%
  mutate(sum(Freq)) %>%
  mutate(percent = Freq / `sum(Freq)`)
head(bar_per) 
write.csv(bar_per,file = "celltype_by_patient_percent.csv")
col =mycolors
ggplot(bar_per, aes(x = percent, y = Var2)) +
  geom_bar(aes(fill = Var1) , stat = "identity") + coord_flip() +
  theme(axis.ticks = element_line(linetype = "blank"),
        legend.position = "top",
        panel.grid.minor = element_line(colour = NA,linetype = "blank"), 
        panel.background = element_rect(fill = NA),
        plot.background = element_rect(colour = NA)) +
  labs(y = "% Relative cell source", fill = NULL)+labs(x = NULL)+
  scale_fill_manual(values=col)

ggsave("2. celltype_by_patient_percent.pdf",
       width = 8,height = 4)



rm(list=ls())
options(stringsAsFactors = F)
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)
library(stringr)
library(ggsci) 
library(patchwork) 
library(ggsci)
library(ggpubr)
library(RColorBrewer) 
getwd()
dir.create('com_go_kegg-figures/')
setwd('com_go_kegg-figures/') 

load('../cosg_celltype__marker_cosg.Rdata') 
head(marker_cosg)
symbols_list <-  as.list(as.data.frame(apply(marker_cosg$names,2,head,100)))
TIM3GENE=symbols_list$TIM3_positive_Microglia
save(TIM3GENE,file = "TIM3GENE.RDATA")
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(ggplot2)
library(stringr)

source('com_go_kegg_ReactomePA_human.R')
com_go_kegg_ReactomePA_human(symbols_list, pro='pbmc' )














