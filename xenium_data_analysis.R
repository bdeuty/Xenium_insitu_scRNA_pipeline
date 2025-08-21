#Xenium in situ scRNA data analysis pipeline
#from: https://satijalab.org/seurat/articles/seurat5_spatial_vignette_2

library(Seurat)
#library(SeuratData)
library(future)
#ibrary("multisession", workers = 10)
library(ggplot2)
library(patchwork)
library(dplyr)
library(spacexr)
library(SummarizedExperiment)

path <- "./Xenium_V1_FF_Mouse_Brain_Coronal_Subset_CTX_HP_outs/"
#Load xenium data
xenium.obj <- LoadXenium(path, fov = "fov")
#remove cells with 0 counts
xenium.obj <- subset(xenium.obj, subset = nCount_Xenium > 0)

#vlnplot
VlnPlot(xenium.obj, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0)

#ImageDimPlot
ImageDimPlot(xenium.obj, fov = "fov", molecules = c("Gad1", "Sst", "Pvalb", "Gfap"), nmols = 20000)

#Image Feature Plot
ImageFeaturePlot(xenium.obj, features = c("Cux2", "Rorb", "Bcl11b", "Foxp2"), max.cutoff = c(25, 35, 12, 10), size = 0.75, cols = c("white", "red"))

#zoom functionality
options(future.globals.maxSize = 2 * 1024^3) #2GB
cropped.coords <- Crop(xenium.obj[['fov']], x = c(1200, 2900), y = c(3750, 4550), coords = "plot")
xenium.obj[["zoom"]] <- cropped.coords
#vizualize cropped area with cell segmentation and selected molecules
DefaultBoundary(xenium.obj[['zoom']]) <- "centroids"
ImageDimPlot(xenium.obj, fov = "zoom", axes = TRUE, border.color = "white", border.size = 0.1, cols = "polychrome", coord.fixed = FALSE, molecules = c("Gad1", "Sst", "Npy2r", "Pvalb", "Nrn1", nmols = 10000))

#normalization
xenium.obj <- SCTransform(xenium.obj, assay = "Xenium")
xenium.obj <- RunPCA(xenium.obj, npcs = 30, features = rownames(xenium.obj))
xenium.obj <- RunUMAP(xenium.obj, dims = 1:30)
xenium.obj <- FindNeighbors(xenium.obj, reduction = "pca", dims = 1:30)
xenium.obj <- FindClusters(xenium.obj, resolution = 0.3)

DimPlot(xenium.obj)

FeaturePlot(xenium.obj, features = c("Cux2", "Bcl11b", "Foxp2", "Gad1", "Sst", "Gfap"))

ImageDimPlot(xenium.obj, cols = "poltchrome", size = 0.75)

xenium.obj <- LoadXenium(path)
p1 <- ImageFeaturePlot(xenium.obj, feature = "Slc17a7", axes = TRUE, max.cutoff = "q90")
p1

crop <- Crop(xenium.obj[["fov"]], x = c(600, 2100), y = c(900, 4700))
xenium.obj[["crop"]] <- crop
p2 <- ImageFeaturePlot(xenium.obj, fov = "crop", features = "Slc17a7", size = 1, axes = TRUE, max.cutoff = "q90")
p2

#Robust Cell Type Decomposition (RCTD)

query.counts <- GetAssayData(xenium.obj, assay = "Xenium", slot = "counts")[, Cells(xenium.obj[["crop"]])]
coords <- GetTissueCoordinates(xenium.obj[["crop"]], which = "centroids")
rownames(coords) <- coords$cell
coords$cell <- NULL
query <- createSpatialRNA(coords, query.counts, colSums(query.counts))

allen.cortex.ref <- readRDS("./allen_cortex.rds")
allen.cortex.ref <- UpdateSeuratObject(allen.cortex.ref)

Idents(allen.cortex.ref) <- "subclass"
#remove CR cells because there aren't enough of them for annotation
allen.cortex.ref <- subset(allen.cortex.ref, subset = subclass != "CR")
counts <- GetAssayData(allen.cortex.ref, assay = "RNA", slot = "counts")
cluster <- as.factor(allen.cortex.ref$subclass)
names(cluster) <- colnames(allen.cortex.ref)
nUMI <- allen.cortex.ref$nCount_RNA
names(nUMI) <- colnames(allen.cortex.ref)
nUMI <- colSums(counts)
levels(cluster) <- gsub("/", "-", levels(cluster))
reference <- createReference(counts, cluster, nUMI)

#make experiment SummarizedExperiment objects
spatial_experiment <- SummarizedExperiment(
  assays = list(counts = query.counts),
  rowData = data.frame(gene = rownames(query.counts)),
  colData = data.frame(
    cell_id = colnames(query.counts),
    nUMI = colSums(query.counts),
    x = coords[,1],
    y = coords[,2]
  )
)

reference_experiment <- SummarizedExperiment(
  assays = list(counts = counts),
  rowData = data.frame(gene = rownames(counts)),
  colData = data.frame(
    cell_id = colnames(counts),
    cell_type = cluster,
    nUMI = colSums(counts)
  )
)

RCTD <- createRctd(spatial_experiment, reference_experiment)
RCTD <- runRctd(RCTD)

annotations.df <- as.data.frame(colData(RCTD))
#annotations.df <- RCTD@results$results_df
annotations <- annotations.df$first_type
names(annotations) <- rownames(annotations.df)
xenium.obj$predicted.celltype <- annotations
keep.cells <- Cells(xenium.obj)[!is.na(xenium.obj$predicted.celltype)]
xenium.obj <- subset(xenium.obj, cells = keep.cells)

xenium.obj <- BuildNicheAssay(object = xenium.obj, fov = "crop", group.by = "predicted.celltype", niches.k = 5, neighbors.k = 30)

celltype.plot <- ImageDimPlot(xenium.obj, group.by = "predicted.celltype", size = 1.5, cols = "polychrome", dark.background = F) + ggtitle("Cell type")
niche.plot <- ImageDimPlot(xenium.obj, group.by = "niches", size = 1.5, dark.background = F) + ggtitle("Niches") + scale_fill_manual(values = c("#442288", "#6CA2EA", "#B5D33D", "#FED23F", "#EB7D5B"))
celltype.plot | niche.plot

table(xenium.obj$predicted.celltype, xenium.obj$niches)
