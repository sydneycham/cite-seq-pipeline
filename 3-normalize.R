suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(Seurat))

# import data -----------------------------------------------------------------
seurat_object = readRDS(snakemake@input$sample)
samplesheet = fread(snakemake@input$samplesheet)
DefaultAssay(seurat_object) = "RNA"

# define constants
cluster_resolution = samplesheet %>%
  filter(sample == snakemake@wildcards$sample) %>%
  pull(resolution)
seed = snakemake@config$seed
pca_components = snakemake@config$PCA$npcs

# RNA normalization -----------------------------------------------------------
seurat_object = 

  # perform scTransform
  SCTransform(
    seurat_object,
    vst.flavor = "v2") %>%

  # pca
  RunPCA(
    npcs = pca_components,
    seed.use = seed) %>%

  # umap
  RunUMAP(
    reduction = "pca",
    n.neighbors = snakemake@config$UMAP$n_neighbors,
    dims = 1:pca_components,
    metric = snakemake@config$UMAP$metric,
    seed.use = seed) %>%

  # find neighbors
  FindNeighbors(
    reduction = "pca",
    dims = 1:pca_components) %>%

  # find clusters
  FindClusters(
    resolution = cluster_resolution,
    random.seed = seed
  )

# ADT normalization -----------------------------------------------------------
DefaultAssay(seurat_object) = "ADT"
VariableFeatures(seurat_object) = rownames(seurat_object[["ADT"]])
seurat_object = NormalizeData(seurat_object, normalization.method = "CLR", margin = 2) %>%
	ScaleData() %>%
	RunPCA(reduction.name = "apca")

# set default assay back to RNA
DefaultAssay(seurat_object) = "RNA"

# generate umap image
p1 = DimPlot(seurat_object, label = TRUE, repel = TRUE)
ggsave(snakemake@output$umap, p1, width = 16, height = 9, dpi = 600)

# generate scree plots for pca
p2 = ElbowPlot(seurat_object, ndims = pca_components, reduction = "pca")
ggsave(snakemake@output$scree, p2, width = 16, height = 9, dpi = 600)

saveRDS(seurat_object, snakemake@output$rds)