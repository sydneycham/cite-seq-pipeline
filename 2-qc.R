suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(tibble))

# import data -----------------------------------------------------------------
seurat_object = readRDS(snakemake@input$rds)
samplesheet = fread(snakemake@input$samplesheet)

# collect metadata ------------------------------------------------------------
cell_md_before = seurat_object@meta.data %>%
  rownames_to_column("barcode") %>%
  mutate(type = "BEFORE") # add before QC label

# filter based on samplesheet qc metrics
nCount_RNA_min_metric = samplesheet %>% filter(sample == snakemake@wildcards$sample) %>% pull("nCount_RNA_pct_filter_min")
nCount_RNA_max_metric = samplesheet %>% filter(sample == snakemake@wildcards$sample) %>% pull("nCount_RNA_pct_filter_max")

nFeature_RNA_min_metric = samplesheet %>% filter(sample == snakemake@wildcards$sample) %>% pull("nFeature_RNA_pct_filter_min")
nFeature_RNA_max_metric = samplesheet %>% filter(sample == snakemake@wildcards$sample) %>% pull("nFeature_RNA_pct_filter_max")

percent_mt_min_metric = samplesheet %>% filter(sample == snakemake@wildcards$sample) %>% pull("mt_pct_filter_min")
percent_mt_max_metric = samplesheet %>% filter(sample == snakemake@wildcards$sample) %>% pull("mt_pct_filter_max")

nCount_RNA_min_metric = quantile(seurat_object$nCount_RNA, nCount_RNA_min_metric)[[1]]
nCount_RNA_max_metric = quantile(seurat_object$nCount_RNA, nCount_RNA_max_metric)[[1]]

nFeature_RNA_min_metric = quantile(seurat_object$nFeature_RNA, nFeature_RNA_min_metric)[[1]]
nFeature_RNA_max_metric = quantile(seurat_object$nFeature_RNA, nFeature_RNA_max_metric)[[1]]

percent_mt_min_metric = quantile(seurat_object$percent_mt, percent_mt_min_metric)[[1]]
percent_mt_max_metric = quantile(seurat_object$percent_mt, percent_mt_max_metric)[[1]]

print(paste("Sample:", snakemake@wildcards$sample))
print(paste("nCount_RNA_min_metric:", nCount_RNA_min_metric))
print(paste("nCount_RNA_max_metric:", nCount_RNA_max_metric))
print(paste("nFeature_RNA_min_metric:", nFeature_RNA_min_metric))
print(paste("nFeature_RNA_max_metric:", nFeature_RNA_max_metric))
print(paste("percent_mt_min_metric:", percent_mt_min_metric))
print(paste("percent_mt_max_metric:", percent_mt_max_metric))

# filter cells
seurat_object = subset(

    seurat_object,

    subset = nFeature_RNA > nCount_RNA_min_metric &
      nFeature_RNA < nCount_RNA_max_metric &

      nCount_RNA > nFeature_RNA_min_metric &
      nCount_RNA < nFeature_RNA_max_metric &

      percent_mt > percent_mt_min_metric &
      percent_mt < percent_mt_max_metric
  )

# export resulst --------------------------------------------------------------

# write md file per sample
cell_md_after = seurat_object@meta.data %>%
  rownames_to_column("barcode") %>%
  mutate(type = "AFTER")

cell_md = bind_rows(cell_md_before, cell_md_after) %>%
  mutate(type = factor(type, levels = c("BEFORE", "AFTER")))
fwrite(cell_md, snakemake@output$md)

# write qc metrics used per sample
metrics_df = data.frame(
  "Sample" = snakemake@wildcards$sample,
  "nCount_RNA_min_metric" = nCount_RNA_min_metric,
  "nCount_RNA_max_metric" = nCount_RNA_max_metric,
  "nFeature_RNA_min_metric" = nFeature_RNA_min_metric,
  "nFeature_RNA_max_metric" = nFeature_RNA_max_metric,
  "percent_mt_min_metric" = percent_mt_min_metric,
  "percent_mt_max_metric" = percent_mt_max_metric
)
fwrite(metrics_df, snakemake@output$metrics)

# write filtered Seurat object to RDS
saveRDS(seurat_object, snakemake@output$rds)