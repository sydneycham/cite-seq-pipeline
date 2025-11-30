suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(Seurat))

# import data -----------------------------------------------------------------

seurat_objects = lapply(unlist(snakemake@input$samples), readRDS)

# merge objects
merged_object = merge(
  seurat_objects[[1]],
  seurat_objects[2:length(seurat_objects)]
)

# RNA normalization -----------------------------------------------------------
merged_object = 

  # perform scTransform on all samples and cells
  SCTransform(
    merged_object,
    vst.flavor = "v2",
    variable.features.n = length(rownames(merged_object)),
    seed.use = snakemake@config$seed) %>%

  # pca
  RunPCA(
    npcs = snakemake@config$PCA$npcs,
    seed.use = snakemake@config$seed)

# ADT normalization -----------------------------------------------------------
DefaultAssay(merged_object) = "ADT"
VariableFeatures(merged_object) = rownames(merged_object[["ADT"]])
merged_object = NormalizeData(merged_object, normalization.method = "CLR", margin = 2) %>%
	ScaleData() %>%
	RunPCA(reduction.name = "apca")

# set default assay back to RNA
DefaultAssay(merged_object) = "RNA"

saveRDS(merged_object, snakemake@output$rds)

# generate scree plots per omic
rna_scree = ElbowPlot(merged_object, ndims = snakemake@config$PCA$npcs, reduction = "pca")
ggsave(snakemake@output$rna_scree, rna_scree, width = 16, height = 9, dpi = 600)

adt_scree = ElbowPlot(merged_object, ndims = snakemake@config$PCA$npcs, reduction = "apca")
ggsave(snakemake@output$adt_scree, adt_scree, width = 16, height = 9, dpi = 600)

# generate PCA colored by samples
p2 = DimPlot(merged_object, reduction = "pca")
ggsave(snakemake@output$pca, p2, width = 16, height = 9, dpi = 600)