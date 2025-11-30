suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(glmGamPoi))

# import data -----------------------------------------------------------------
seurat_object = readRDS(snakemake@input[[1]])
pca_components = snakemake@config$integrate$wnn$pca$rna_pcs
adt_components = snakemake@config$integrate$wnn$pca$adt_pcs

umap_neighbors = snakemake@config$integrate$wnn$umap$n_neighbors
umap_components = snakemake@config$integrate$wnn$umap$n_components
integrated_resolution = snakemake@config$integrate$FindClusters$resolution

# prepare to find multi modal neighbors ---------------------------------------
seurat_object = FindMultiModalNeighbors(
	seurat_object,
	reduction.list = list("pca", "apca"),
	dims.list = list(1:pca_components, 1:adt_components),
	modality.weight.name = "RNA.weight"
)

# run umap
seurat_object <- RunUMAP(
	seurat_object,
	nn.name = "weighted.nn",
	reduction.name = "wnn.umap",
	reduction.key = "wnnUMAP_",
	n.neighbors = umap_neighbors,
	n.components = umap_components
)

# find clusters
seurat_object <- FindClusters(
	seurat_object,
	graph.name = "wsnn",
	algorithm = 3,
	resolution = integrated_resolution,
	verbose = FALSE
)

# generate umap for wnn, rna-, and adt-only assays ----------------------------
seurat_object <- RunUMAP(seurat_object, reduction = 'pca', dims = 1:30, assay = 'RNA', reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')
seurat_object <- RunUMAP(seurat_object, reduction = 'apca', dims = 1:18, assay = 'ADT', reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')

p1 <- DimPlot(seurat_object, reduction = 'wnn.umap', group.by = 'orig.ident', label = TRUE, repel = TRUE, label.size = 2.5) + NoLegend()
p2 <- DimPlot(seurat_object, reduction = 'rna.umap', group.by = 'orig.ident', label = TRUE, repel = TRUE, label.size = 2.5) + NoLegend()
p3 <- DimPlot(seurat_object, reduction = 'adt.umap', group.by = 'orig.ident', label = TRUE, repel = TRUE, label.size = 2.5) + NoLegend()

ggsave(snakemake@output$wnn_umap, p1, width = 16, height = 9, dpi = 600)
ggsave(snakemake@output$rna_umap, p2, width = 16, height = 9, dpi = 600)
ggsave(snakemake@output$adt_umap, p3, width = 16, height = 9, dpi = 600)

# export results --------------------------------------------------------------
saveRDS(seurat_object, snakemake@output$rds)