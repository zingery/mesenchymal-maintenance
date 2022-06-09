## Zinger Yang Loureiro
## <zinger.yang@umassmed.edu>
## code below was implemented in R version 4.0.2 (2020-06-22)

# [1] viridis_0.6.2        viridisLite_0.4.0    gplots_3.1.1         ggpubr_0.4.0         scales_1.1.1        
# [6] ggrepel_0.9.1        velocyto.R_0.6       Matrix_1.4-1         SeuratWrappers_0.3.0 Seurat_4.1.0        
# [11] forcats_0.5.1        stringr_1.4.0        dplyr_1.0.8          purrr_0.3.4          readr_2.1.2         
# [16] tidyr_1.2.0          tibble_3.1.6         ggplot2_3.3.5        tidyverse_1.3.1      SeuratObject_4.0.4  

library(tidyverse)
library(Seurat)
library(SeuratWrappers)
library(velocyto.R)
library(ggplot2)
library(ggrepel)
library(scales)
library(ggpubr)
library(gplots)
library(viridis)


## ------------- Preparing Seurat object -------------
## preparing the merged object of cell ranger generated result with velocyto velocity containing object

# loading the 10X Cellranger output
undiff_adi.data <- Read10X(data.dir = "CellRanger output/")
undiff_adi.seurat <- CreateSeuratObject(counts = undiff_adi.data, 
                                        project = "Hu57", min.cells = 250, min.features = 2000) 
undiff_adi.seurat@meta.data$orig.ident = ifelse(grepl("-1",rownames(undiff_adi.seurat@meta.data)),"Adipo3d","Undiff")
undiff_adi.seurat.count <- undiff_adi.seurat@assays$RNA@counts %>% as.matrix()
colnames(undiff_adi.seurat.count)  <- (data.frame( original_name = colnames(undiff_adi.seurat.count)) %>% 
                                         mutate(header = ifelse(grepl("-1",original_name), "2019NOV7_Adipo3d:","2019NOV7_Undiff:")) %>%
                                         mutate(new_name = paste0(header, str_sub(original_name,0,-3), "x"))) $new_name

# loading velocyto.R loom file
merge.velo <- ReadVelocity("RNA Velocity/2019NOV7_Merged.loom")
merge.velo.Seurat.nofilter <- as.Seurat(x = merge.velo)
merge.velo.Seurat <- subset(merge.velo.Seurat.nofilter,cells = colnames(undiff_adi.seurat.count))


# Combine the original Seurat object with object containing velocity information
undiff_adi.seurat.obj <- CreateAssayObject(counts=undiff_adi.seurat.count[,colnames(merge.velo.Seurat)] )
merge.velo.Seurat[["RNA"]] <- undiff_adi.seurat.obj
 

rm(undiff_adi.data, undiff_adi.seurat, merge.velo, merge.velo.Seurat.nofilter)
gc()


## ------------- Seurat analysis -------------  
## Normalization, PCA, UMAP
merge.velo.Seurat <- SCTransform(merge.velo.Seurat, assay="RNA")  
merge.velo.Seurat <- RunPCA(object = merge.velo.Seurat)
merge.velo.Seurat <- FindNeighbors(object = merge.velo.Seurat, dims = 1:10)
merge.velo.Seurat <- FindClusters(object = merge.velo.Seurat, resolution = 0.6)

# Annotate Clusters
merge.velo.Seurat@meta.data$my_cluster <- (merge.velo.Seurat@meta.data %>%
                                             as.data.frame() %>%
                                             mutate(my_cluster = case_when(
                                               SCT_snn_res.0.6 == 5 ~ "ADIPOQ+",    
                                               SCT_snn_res.0.6 == 3 ~ "MGP+",    
                                               SCT_snn_res.0.6 %in% c(2,6) ~ "others",    
                                               TRUE ~ "Non-induced progenitors"     
                                             )))$my_cluster

# Annotate samples
merge.velo.Seurat@meta.data$sample <-  
  (merge.velo.Seurat@meta.data %>%
     as.data.frame() %>%
     tibble::rownames_to_column(var="the_cols") %>%
     mutate(sample = ifelse(grepl("Adipo3d",the_cols), "Induced","Non-induced")) )$sample

## load precalculated RDS
merge.velo.Seurat <- readRDS("Merged_Velocity_Seruat.RDS")



## ----- QC (Extended Data Figure 1) -----

# Estimate mitochondrial reads percentage
merge.velo.Seurat[["percent.mt"]] <- PercentageFeatureSet(merge.velo.Seurat, pattern = "^MT-")


## Extended data figure 1b
VlnPlot(merge.velo.Seurat, assay = "unspliced", group.by = "sample",features = c("nFeature_RNA","nCount_RNA","percent.mt"), 
        pt.size = 0, ncol = 3) +
  NoLegend()


## Extended data figure 1c
FeaturePlot(merge.velo.Seurat, features = c("MALAT1","RPL4","ACTB"), ncol = 3, 
            cols=c("#cccccc","#d62728"), 
            reduction = "pca",
            pt.size = 0.5, order=T) & NoAxes()

## Extended data figure 1d
VlnPlot(merge.velo.Seurat, features = c("MALAT1","RPL4","ACTB"),pt.size = 0, group.by = "my_cluster")


 
## ----- PCA projections -----

## Figure 2b
DimPlot(merge.velo.Seurat, reduction = "pca", group.by = "sample", cols=c("#2ca02c","#9467bd"))
  


## ----- Velocity analysis -----

merge.velo.Seurat <- RunVelocity(object = merge.velo.Seurat, deltaT = 1, kCells = 25, fit.quantile = 0.02)

ident.colors <- c("#DEDEDE","#DEDEDE","#DEDEDE","#F8766D","#DEDEDE","#00BADE","#DEDEDE","#DEDEDE","#DEDEDE","#DEDEDE" )
names(x = ident.colors) <- levels(x = merge.velo.Seurat)
cell.colors <- ident.colors[Idents(object = merge.velo.Seurat)]
names(x = cell.colors) <- colnames(x = merge.velo.Seurat)

## Figure 2c without velocity
DimPlot(merge.velo.Seurat, reduction = "pca", cols=ident.colors)

## Figure 2c  
show.velocity.on.embedding.cor(emb = Embeddings(object = merge.velo.Seurat, reduction = "pca"), 
                               vel = Tool(object = merge.velo.Seurat, slot = "RunVelocity"), 
                               n = 250, scale = "sqrt", 
                               cell.colors = ac(x = cell.colors, alpha = 0.5), 
                               cex = 0.8, arrow.scale = 40, show.grid.flow = TRUE, 
                               min.grid.cell.mass = 0.5, grid.n = 25, arrow.lwd = 1, 
                               do.par = FALSE, cell.border.alpha = 0.1)


## ----- Gene View -----

## Figure 2d,f,h (PCA projection)
FeaturePlot(merge.velo.Seurat, features = c("ADIPOQ","PLIN1","LPL","MGP","CTHRC1","DCN","THY1","ENG","NT5E"), ncol = 3, 
            cols=c("#cccccc","#d62728"), 
            reduction = "pca",
            pt.size = 0.1, order=T) & NoAxes()

## Figure 2e,g,i (dot plots) 
p1 <- DotPlot(merge.velo.Seurat, features = c("ADIPOQ","PLIN1","LPL"),scale.min=0, scale.max = 100,scale=F,
        group.by = "my_cluster", idents = c(0,1,3,5,7,8,9))  + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("")+ylab("")+
  scale_color_gradientn(limits = c(0,2), colors=c("lightgrey","blue"), oob=squish)

p2 <- DotPlot(merge.velo.Seurat, features = c("MGP","CTHRC1","DCN"), scale.min=0, scale.max = 100,scale=F,
        group.by = "my_cluster", idents = c(0,1,3,5,7,8,9))  + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))   +
  xlab("")+ylab("")+
  scale_color_gradientn(limits = c(0,2), colors=c("lightgrey","blue"), oob=squish)

p3 <- DotPlot(merge.velo.Seurat, features = c("THY1","ENG","NT5E"), scale.min=0, scale.max = 100, scale=F,
        group.by = "my_cluster", idents = c(0,1,3,5,7,8,9))  + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))    +
  xlab("")+ylab("")+
  scale_color_gradientn(limits = c(0,2), colors=c("lightgrey","blue"), oob=squish)

ggarrange(p1,p2,p3,ncol=1,common.legend = T, legend = "right")
 

## Figure 3e
DotPlot(merge.velo.Seurat, features = c("SNAI2","SFRP2","DPP4","WISP2","DKK1"),scale=F,
        dot.scale = 11,  
        group.by = "my_cluster", idents = c(0,1,3,5,7,8,9))  + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("")+ylab("") +
  theme(legend.box = "horizontal",
        axis.text = element_text(size = 18), text = element_text(size = 18))

## Figure 3f
DotPlot(merge.velo.Seurat, features = c("WNT5A","WNT5B","CTNNB1","TCF7L2"), scale=F,
        group.by = "my_cluster", idents = c(0,1,3,5,7,8,9))  + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))   +
  xlab("")+ylab("") +
  theme(legend.box = "horizontal",
        axis.text = element_text(size = 18), text = element_text(size = 18))

## Figure 3i
DotPlot(merge.velo.Seurat, features = c("EFNB1","EPHB6"),scale=F,
        group.by = "my_cluster", idents = c(0,1,3,5,7,8,9))  + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  +
  xlab("")+ylab("") +
  theme(legend.box = "horizontal",
        axis.text = element_text(size = 18), text = element_text(size = 18))
 

## ----- Differential expression analysis -----

# Induced nondifferentiating vs ADIPOQ+, pre-filter features that are detected at <25% frequency in either cell population
de.no_adi.adi <- FindMarkers(merge.velo.Seurat, ident.1 = "5", ident.2 = "3", min.pct = 0.25, logfc.threshold=0)  

## Figure 3a (volcano plot)
fig3_gene <- c("DCN","CTSK","VCAN","COL1A1","COL3A1","CTHRC1","MGP","FABP4","FABP5","ADIPOQ","FASN","PNPLA2","LPL","DGAT2")
highlight <- data.frame(gene = fig3_gene,interesting = fig3_gene) 
  
de.no_adi.adi.plot <- de.no_adi.adi %>% 
  tibble::rownames_to_column(var = "gene") %>%
  mutate(significance = 
           case_when(
             avg_log2FC >1 & p_val_adj < 0.001 ~ "Up in ADIPOQ+ cells (121)",
             avg_log2FC <(-1) & p_val_adj < 0.001 ~ "Up in MGP+ cells (114)",
             TRUE ~ "not significant"
           )
  ) %>%
  mutate(neg_log10_pval = -log10(p_val_adj)) %>%
  mutate(neg_log10_pval = ifelse(is.infinite(neg_log10_pval), 320, neg_log10_pval)) %>%
  left_join(highlight)


ggplot(de.no_adi.adi.plot, 
       aes(x=avg_log2FC,y=neg_log10_pval, color=significance, label=interesting)) +
  geom_point(alpha=0.9) +
  geom_label_repel(size = 4, color="black",
                   box.padding = unit(0.5, "lines"),
                   segment.color = 'grey50', max.overlaps=25)+
  theme_classic(base_size=14) +
  scale_color_manual(values = c("#555555","#00BADE","#F8766D"))+
  scale_y_continuous(limits = c(0, 350))+
  scale_x_continuous(limits = c(-6, 7))+
  xlab("Average log2 Fold Change")+
  ylab("-log10 adjusted p-value")+
  theme(legend.title = element_blank()) + theme(legend.position="top")

 
## Figure 3d (volcano plot)
highlight <- data.frame(gene = c("SFRP2","DKK1","WISP2","DPP4","SNAI2"),
                        interesting = c("SFRP2","DKK1","WISP2","DPP4","SNAI2")) 

de.no_adi.adi.plot <- de.no_adi.adi %>% 
  tibble::rownames_to_column(var = "gene") %>%
  mutate(significance = 
           case_when(
             avg_log2FC >1 & p_val_adj < 0.001 ~ "Up in ADIPOQ+ cells (121)",
             avg_log2FC <(-1) & p_val_adj < 0.001 ~ "Up in MGP+ cells (114)",
             TRUE ~ "not significant"
           )
  ) %>%
  mutate(neg_log10_pval = -log10(p_val_adj)) %>%
  mutate(neg_log10_pval = ifelse(is.infinite(neg_log10_pval), 320, neg_log10_pval)) %>%
  left_join(highlight)


ggplot(de.no_adi.adi.plot, 
       aes(x=avg_log2FC,y=neg_log10_pval, color=significance, label=interesting)) +
  geom_point(alpha=0.9) +
  geom_label_repel(size = 4, color="black",
                   box.padding = unit(0.5, "lines"),
                   segment.color = 'grey50', max.overlaps=25)+
  theme_classic(base_size=14) +
  scale_color_manual(values = c("#555555","#00BADE","#F8766D"))+
  scale_y_continuous(limits = c(0, 350))+
  scale_x_continuous(limits = c(-6, 7))+
  xlab("Average log2 Fold Change")+
  ylab("-log10 adjusted p-value")+
  theme(legend.title = element_blank()) + theme(legend.position="top")

 
 
##----- Heatmap of top markers-----

heatmap_genes <- c( rownames(de.no_adi.adi %>% slice_max(order_by=avg_log2FC,n=40)) ,
                    rownames(de.no_adi.adi %>% slice_min(order_by=avg_log2FC,n=40)) )

merge.velo.Seurat.norm <- GetAssayData(
  object = subset(x = merge.velo.Seurat, idents = c("3","5"), 
                  features = heatmap_genes), slot = "data") 
  
merge.velo.Seurat.norm.df.t <- merge.velo.Seurat.norm %>% as.data.frame() %>% 
  tibble::rownames_to_column(var="gene") %>%
  pivot_longer(-gene,names_to="umi",values_to="count") 
  
merge.velo.Seurat.norm.scale <- merge.velo.Seurat.norm.df.t %>% left_join(
  merge.velo.Seurat.norm.df.t %>% group_by(gene) %>%
      summarize(maxCount = max(count)) %>%
      ungroup()
  ) %>%
  mutate(scaled_exp = count / maxCount) %>%
  select(-count,-maxCount) %>%
  pivot_wider(names_from = "umi", values_from="scaled_exp") %>%
  tibble::column_to_rownames(var="gene") %>%
  as.matrix()


## Figure 4e, left panel
heatmap.2(merge.velo.Seurat.norm.scale, Rowv = F, Colv = T, scale = "none", dendrogram = "column",
          density.info='none',
          lhei=c(1,4),
          col=viridis(100), trace="none", margins = c(1,5),labCol = FALSE) 

## Figure 4e, right panel
DotPlot(merge.velo.Seurat, features = heatmap_genes,  scale.min=0, scale.max = 100,scale=F,
        dot.scale = 3,
        cols = c("lightgrey", "blue"),
        group.by = "my_cluster", idents = c(3,5)) +coord_flip() +
        scale_x_discrete(limits=rev)
 


## ------------- Ligand Receptor Analysis -------------

# Ligand Receptor database from Ramilowski 2015 (total of 2422 pairs)
dat_lr <- data.table::fread(" fantom5_LR_all.txt") %>%
  filter(Pair.Evidence %in% c("literature supported","putative"))

# Average gene expression in the single cell clusters
merge.velo.Seurat.avgCount <- AverageExpression(merge.velo.Seurat, group.by = "my_cluster")[["SCT"]] 
merge.velo.Seurat.avgCount.df <- merge.velo.Seurat.avgCount %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var="gene") 

# Top 25 percentile expressed genes in each cluster are considered for matching with known ligand receptors
exp.in.adi <-  (merge.velo.Seurat.avgCount.df  %>% filter(`ADIPOQ+`>=quantile(merge.velo.Seurat.avgCount.df$`ADIPOQ+`)[["75%"]]))$gene
exp.in.mgp <- (merge.velo.Seurat.avgCount.df  %>% filter(`MGP+`>=quantile(merge.velo.Seurat.avgCount.df$`ADIPOQ+`)[["75%"]]))$gene

# ligands/receptors in each group
dat_lr %>% filter(Ligand.ApprovedSymbol %in% exp.in.adi) %>% select(Ligand.ApprovedSymbol) %>% distinct() %>% dim() #81
dat_lr %>% filter(Receptor.ApprovedSymbol %in% exp.in.adi) %>% select(Receptor.ApprovedSymbol) %>% distinct() %>% dim() #53
dat_lr %>% filter(Ligand.ApprovedSymbol %in% exp.in.mgp) %>% select(Ligand.ApprovedSymbol) %>% distinct() %>% dim() #90
dat_lr %>% filter(Receptor.ApprovedSymbol %in% exp.in.mgp) %>% select(Receptor.ApprovedSymbol) %>% distinct() %>% dim() #67

## Extended data table 1
lr.adi_to_mgp <- dat_lr %>% filter(Ligand.ApprovedSymbol  %in% exp.in.adi &
                                   Receptor.ApprovedSymbol  %in% exp.in.mgp )
lr.mgp_to_adi <- dat_lr %>% filter(Ligand.ApprovedSymbol  %in% exp.in.mgp &
                                   Receptor.ApprovedSymbol  %in% exp.in.adi )
lr.adi_to_adi <- dat_lr %>% filter(Ligand.ApprovedSymbol  %in% exp.in.adi &
                                     Receptor.ApprovedSymbol  %in% exp.in.adi )
lr.mgp_to_mgp <- dat_lr %>% filter(Ligand.ApprovedSymbol  %in% exp.in.mgp &
                                     Receptor.ApprovedSymbol  %in% exp.in.mgp )

 