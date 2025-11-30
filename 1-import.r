suppressPackageStartupMessages(library(Seurat))

# import data ------------------------------------------------------------------
gex_matrix = Read10X_h5(snakemake@input[[1]])$`Gene Expression`
adt_matrix = Read10X_h5(snakemake@input[[1]])$`Antibody Capture`

# create seurat object
seurat_object = CreateSeuratObject(gex_matrix, project = snakemake@wildcards$sample)
seurat_object[["ADT"]] = CreateAssay5Object(adt_matrix)

# include mitochondrial content into each cell
seurat_object[["percent_mt"]] = PercentageFeatureSet(seurat_object, pattern = "^MT-|^mt-")

# include genes per UMI for each cell
seurat_object[["log10GenesPerUMI"]] = log10(seurat_object[["nFeature_RNA"]]) / log10(seurat_object[["nCount_RNA"]])

# include antibodies per UMI for each cell
seurat_object[["log10ProteinPerUMI"]] = log10(seurat_object[["nFeature_ADT"]]) / log10(seurat_object[["nCount_ADT"]])

# export data -----------------------------------------------------------------
saveRDS(seurat_object, snakemake@output[[1]])