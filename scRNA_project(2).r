library(Seurat)
load("3-cell/GSE147528.RData")
dir.create("4-further")
setwd("4-further")
unique(sce.all$celltype)
mycolors <-c("#B07AA1","#F89C74","#66C5CC","#75ACC3","#FE88B1","#80BA5A","#F6CF71","#EDC948","#B84D64", 
             "#008695","#59A14F","#FF9DA7","#F28E2B","#DCB0F2","#EE7072","#E73F74","#f69896","#9B8E8C",
             "#E68310","#E15759")

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(RColorBrewer)
library(slingshot)
library(clusterProfiler)
library(org.Hs.eg.db) 
library(CellChat)
library(monocle3)
library(pheatmap)
library(viridis)
library(ComplexHeatmap)

output_dir <- "TIM3_Microglia_Analysis/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)


celltype_counts <- table(sce.all$celltype)
celltype_percentage <- prop.table(celltype_counts) * 100

p1 <- ggplot(data.frame(CellType = names(celltype_counts), 
                        Count = as.numeric(celltype_counts)),
             aes(x = reorder(CellType, -Count), y = Count, fill = CellType)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = mycolors[1:length(celltype_counts)]) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        legend.position = "right") +
  labs(title = "Cell Type Distribution", x = "Cell Type", y = "Cell Count")

ggsave(paste0(output_dir, "CellType_Distribution.pdf"), p1, width = 10, height = 6)

p2 <- DimPlot(sce.all, reduction = "umap", group.by = "celltype",raster=FALSE, 
              label = TRUE, label.size = 4, repel = TRUE) +
  scale_color_manual(values = mycolors[1:length(unique(sce.all$celltype))]) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14)) +
  labs(title = "UMAP Visualization of All Cell Types")

ggsave(paste0(output_dir, "UMAP_All_Cells.pdf"), p2, width = 10, height = 8)











