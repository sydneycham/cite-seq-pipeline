---
title: "fgsea_data"
author: "Sydney Hamilton"
date: "2024-08-19"
output: html_document
---
```{r}
library(fgsea)
library(data.table)
library(ggplot2)
```


```{r setup, include=FALSE}
# Read in csvs of genesets from Van Galen
Normal_Derived_Combined_GMP <- as.vector(t(read.csv("/Users/hamiltsy/Documents/Pipeline/new_cite-seq/Gene_sets/Normal_Derived_Combined_GMP.csv", header = FALSE)))

Normal_Derived_Combined_HSC_Prog <- as.vector(t(read.csv("/Users/hamiltsy/Documents/Pipeline/new_cite-seq/Gene_sets/Normal_Derived_Combined_HSC_Prog.csv", header = FALSE)))

Normal_Derived_Combined_Myeloid <- as.vector(t(read.csv("/Users/hamiltsy/Documents/Pipeline/new_cite-seq/Gene_sets/Normal_Derived_Combined_Myeloid.csv", header = FALSE)))

Tumor_Derived_Combined_GMPlike<- as.vector(t(read.csv("/Users/hamiltsy/Documents/Pipeline/new_cite-seq/Gene_sets/Tumor_Derived_Combined_GMPlike.csv", header = FALSE)))

Tumor_Derived_Combined_HSC_Prog <- as.vector(t(read.csv("/Users/hamiltsy/Documents/Pipeline/new_cite-seq/Gene_sets/Tumor_Derived_Combined_HSC_Prog.csv", header = FALSE)))

Tumor_Derived_Combined_Myeloidlike <- as.vector(t(read.csv("/Users/hamiltsy/Documents/Pipeline/new_cite-seq/Gene_sets/Tumor_Derived_Combined_Myeloidlike.csv", header = FALSE)))

Tumor_Derived_percelltype_CDClike <- as.vector(t(read.csv("/Users/hamiltsy/Documents/Pipeline/new_cite-seq/Gene_sets/Tumor_Derived_percelltype_CDClike.csv", header = FALSE)))

Tumor_Derived_percelltype_GMPlike <- as.vector(t(read.csv("/Users/hamiltsy/Documents/Pipeline/new_cite-seq/Gene_sets/Tumor_Derived_percelltype_GMPlike.csv", header = FALSE)))

Tumor_Derived_percelltype_HSClike <- as.vector(t(read.csv("/Users/hamiltsy/Documents/Pipeline/new_cite-seq/Gene_sets/Tumor_Derived_percelltype_HSClike.csv", header = FALSE)))

Tumor_Derived_percelltype_Monocytelike <- as.vector(t(read.csv("/Users/hamiltsy/Documents/Pipeline/new_cite-seq/Gene_sets/Tumor_Derived_percelltype_Monocytelike.csv", header = FALSE)))

Tumor_Derived_percelltype_Progenitorlike <- as.vector(t(read.csv("/Users/hamiltsy/Documents/Pipeline/new_cite-seq/Gene_sets/Tumor_Derived_percelltype_Progenitorlike.csv", header = FALSE)))

Tumor_Derived_percelltype_Promonolike <- as.vector(t(read.csv("/Users/hamiltsy/Documents/Pipeline/new_cite-seq/Gene_sets/Tumor_Derived_percelltype_Promonolike.csv", header = FALSE)))

cell_types <- list(
  NDC_GMP = Normal_Derived_Combined_GMP,
  NDC_HSCPROG = Normal_Derived_Combined_HSC_Prog, 
  NDC_MY = Normal_Derived_Combined_Myeloid, 
  TDC_GMP = Tumor_Derived_Combined_GMPlike, 
  TDC_HSCPROG = Tumor_Derived_Combined_HSC_Prog, 
  TDC_MY = Tumor_Derived_Combined_Myeloidlike, 
  TDP_CDC = Tumor_Derived_percelltype_CDClike, 
  TDP_GMPL = Tumor_Derived_percelltype_GMPlike, 
  TDP_HSCL = Tumor_Derived_percelltype_HSClike, 
  TDP_MonoL = Tumor_Derived_percelltype_Monocytelike, 
  TDP_ProgL = Tumor_Derived_percelltype_Progenitorlike, 
  TDP_PromL = Tumor_Derived_percelltype_Promonolike
)
```

```{r}
#per cluster markers
table(DMSO_scRNAseq_markers$cluster)

# Filter out only significant genes
DMSO_scRNAseq_markers.filtered <- DMSO_scRNAseq_markers[DMSO_scRNAseq_markers$p_val_adj < 0.05,]

DMSO_scRNAseq_markers.filtered.cluster0 <- DMSO_scRNAseq_markers[DMSO_scRNAseq_markers$cluster == 0,]
DMSO_scRNAseq_markers.filtered.cluster1 <- DMSO_scRNAseq_markers[DMSO_scRNAseq_markers$cluster == 1,]
DMSO_scRNAseq_markers.filtered.cluster2 <- DMSO_scRNAseq_markers[DMSO_scRNAseq_markers$cluster == 2,]
DMSO_scRNAseq_markers.filtered.cluster3 <- DMSO_scRNAseq_markers[DMSO_scRNAseq_markers$cluster == 3,]
DMSO_scRNAseq_markers.filtered.cluster4 <- DMSO_scRNAseq_markers[DMSO_scRNAseq_markers$cluster == 4,]
DMSO_scRNAseq_markers.filtered.cluster5 <- DMSO_scRNAseq_markers[DMSO_scRNAseq_markers$cluster == 5,]
DMSO_scRNAseq_markers.filtered.cluster6 <- DMSO_scRNAseq_markers[DMSO_scRNAseq_markers$cluster == 6,]
DMSO_scRNAseq_markers.filtered.cluster7 <- DMSO_scRNAseq_markers[DMSO_scRNAseq_markers$cluster == 7,]
DMSO_scRNAseq_markers.filtered.cluster8 <- DMSO_scRNAseq_markers[DMSO_scRNAseq_markers$cluster == 8,]
DMSO_scRNAseq_markers.filtered.cluster9 <- DMSO_scRNAseq_markers[DMSO_scRNAseq_markers$cluster == 9,]


# Create a named vector for fgsea
fgsea_input <- setNames(object = DMSO_scRNAseq_markers.filtered$avg_log2FC, nm = rownames(DMSO_scRNAseq_markers.filtered))

fgseaRes <- fgsea(pathways = cell_types, 
                  stats    = fgsea_input,
                  eps      = 0.0,
                  minSize  = 15,
                  maxSize  = 500)

head(fgseaRes[order(pval), ])

plotNDC_GMP = plotEnrichment(cell_types[["NDC_GMP"]],
               fgsea_input) + labs(title="NDC_GMP")

plotNDC_HSCPROG = plotEnrichment(cell_types[["NDC_HSCPROG"]],
               fgsea_input) + labs(title="NDC_HSCPROG")

plotNDC_MY = plotEnrichment(cell_types[["NDC_MY"]],
               fgsea_input) + labs(title="NDC_MY")

plotTDC_GMP = plotEnrichment(cell_types[["TDC_GMP"]],
               fgsea_input) + labs(title="TDC_GMP")

plotTDC_HSCPROG = plotEnrichment(cell_types[["TDC_HSCPROG"]],
               fgsea_input) + labs(title="TDC_HSCPROG")

plotTDC_MY = plotEnrichment(cell_types[["TDC_MY"]],
               fgsea_input) + labs(title="TDC_MY")

plotTDP_CDC = plotEnrichment(cell_types[["TDP_CDC"]],
               fgsea_input) + labs(title="TDP_CDC")

plotTDP_GMPL = plotEnrichment(cell_types[["TDP_GMPL"]],
               fgsea_input) + labs(title="TDP_GMPL")

plotTDP_HSCL = plotEnrichment(cell_types[["TDP_HSCL"]],
               fgsea_input) + labs(title="TDP_HSCL")

plotTDP_MonoL = plotEnrichment(cell_types[["TDP_MonoL"]],
               fgsea_input) + labs(title="TDP_MonoL")

plotTDP_ProgL = plotEnrichment(cell_types[["TDP_ProgL"]],
               fgsea_input) + labs(title="TDP_ProgL")

plotTDP_PromL = plotEnrichment(cell_types[["TDP_PromL"]],
               fgsea_input) + labs(title="TDP_PromL")

#creating grid from plots
grid <- grid.arrange(plotNDC_GMP, plotNDC_HSCPROG, plotNDC_MY, plotTDC_GMP, plotTDC_HSCPROG, plotTDC_MY, plotTDP_CDC, plotTDP_GMPL, plotTDP_HSCL, plotTDP_MonoL, plotTDP_ProgL, plotTDP_PromL, ncol=3, nrow=4)

ggsave(filename = "DMSO_Enrichment_plots.jpg", plot = grid, path ="~/Documents/Pipeline/new_cite-seq/Olga/GSEA_results")
```

```{r}
topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
DMSO_GSEA = plotGseaTable(cell_types[topPathways], fgsea_input, fgseaRes, 
              gseaParam=0.5) + labs(title = "DMSO")
ggsave(filename = "DMSO_GSEA.png", path ="~/Documents/Pipeline/new_cite-seq/Olga/GSEA_results")
```

```{r}
DMSO_scRNAseq_markers.filtered.cluster0 <- DMSO_scRNAseq_markers[DMSO_scRNAseq_markers$cluster == 0,]
DMSO_scRNAseq_markers.filtered.cluster1 <- DMSO_scRNAseq_markers[DMSO_scRNAseq_markers$cluster == 1,]
DMSO_scRNAseq_markers.filtered.cluster2 <- DMSO_scRNAseq_markers[DMSO_scRNAseq_markers$cluster == 2,]
DMSO_scRNAseq_markers.filtered.cluster3 <- DMSO_scRNAseq_markers[DMSO_scRNAseq_markers$cluster == 3,]
DMSO_scRNAseq_markers.filtered.cluster4 <- DMSO_scRNAseq_markers[DMSO_scRNAseq_markers$cluster == 4,]
DMSO_scRNAseq_markers.filtered.cluster5 <- DMSO_scRNAseq_markers[DMSO_scRNAseq_markers$cluster == 5,]
DMSO_scRNAseq_markers.filtered.cluster6 <- DMSO_scRNAseq_markers[DMSO_scRNAseq_markers$cluster == 6,]
DMSO_scRNAseq_markers.filtered.cluster7 <- DMSO_scRNAseq_markers[DMSO_scRNAseq_markers$cluster == 7,]
DMSO_scRNAseq_markers.filtered.cluster8 <- DMSO_scRNAseq_markers[DMSO_scRNAseq_markers$cluster == 8,]
DMSO_scRNAseq_markers.filtered.cluster9 <- DMSO_scRNAseq_markers[DMSO_scRNAseq_markers$cluster == 9,]

DMSO_scRNAseq_markers.cluster0 <- DMSO_scRNAseq_markers[DMSO_scRNAseq_markers$cluster == 0,]
DMSO_scRNAseq_markers.cluster1 <- DMSO_scRNAseq_markers[DMSO_scRNAseq_markers$cluster == 1,]
DMSO_scRNAseq_markers.cluster2 <- DMSO_scRNAseq_markers[DMSO_scRNAseq_markers$cluster == 2,]
DMSO_scRNAseq_markers.cluster3 <- DMSO_scRNAseq_markers[DMSO_scRNAseq_markers$cluster == 3,]
DMSO_scRNAseq_markers.cluster4 <- DMSO_scRNAseq_markers[DMSO_scRNAseq_markers$cluster == 4,]
DMSO_scRNAseq_markers.cluster5 <- DMSO_scRNAseq_markers[DMSO_scRNAseq_markers$cluster == 5,]
DMSO_scRNAseq_markers.cluster6 <- DMSO_scRNAseq_markers[DMSO_scRNAseq_markers$cluster == 6,]
DMSO_scRNAseq_markers.cluster7 <- DMSO_scRNAseq_markers[DMSO_scRNAseq_markers$cluster == 7,]
DMSO_scRNAseq_markers.cluster8 <- DMSO_scRNAseq_markers[DMSO_scRNAseq_markers$cluster == 8,]
DMSO_scRNAseq_markers.cluster9 <- DMSO_scRNAseq_markers[DMSO_scRNAseq_markers$cluster == 9,]


# Create a named vector for fgsea
fgsea_input <- setNames(object = DMSO_scRNAseq_markers.filtered$avg_log2FC, nm = rownames(DMSO_scRNAseq_markers.filtered))

fgseaRes <- fgsea(pathways = cell_types, 
                  stats    = fgsea_input,
                  eps      = 0.0,
                  minSize  = 15,
                  maxSize  = 500)

head(fgseaRes[order(pval), ])

plotGMP = plotEnrichment(cell_types[["GMP"]],
               fgsea_input) + labs(title="GMP")

plotHSCPROG = plotEnrichment(cell_types[["HSC_PROG"]],
               fgsea_input) + labs(title="HSC_PROG")

plotMY = plotEnrichment(cell_types[["MY"]],
               fgsea_input) + labs(title="MY")

plotGMP = plotEnrichment(cell_types[["TDC_GMP"]],
               fgsea_input) + labs(title="TDC_GMP")

plotHSCPROG = plotEnrichment(cell_types[["TDC_HSCPROG"]],
               fgsea_input) + labs(title="TDC_HSCPROG")

plotMY = plotEnrichment(cell_types[["TDC_MY"]],
               fgsea_input) + labs(title="TDC_MY")

plotCDC = plotEnrichment(cell_types[["TDP_CDC"]],
               fgsea_input) + labs(title="TDP_CDC")

plotGMPL = plotEnrichment(cell_types[["TDP_GMPL"]],
               fgsea_input) + labs(title="TDP_GMPL")

plotHSCL = plotEnrichment(cell_types[["TDP_HSCL"]],
               fgsea_input) + labs(title="TDP_HSCL")

plotMonoL = plotEnrichment(cell_types[["TDP_MonoL"]],
               fgsea_input) + labs(title="TDP_MonoL")

plotProgL = plotEnrichment(cell_types[["TDP_ProgL"]],
               fgsea_input) + labs(title="TDP_ProgL")

plotPromL = plotEnrichment(cell_types[["TDP_PromL"]],
               fgsea_input) + labs(title="TDP_PromL")

#creating grid from plots
grid <- grid.arrange(plotNDC_GMP, plotNDC_HSCPROG, plotNDC_MY, plotTDC_GMP, plotTDC_HSCPROG, plotTDC_MY, plotTDP_CDC, plotTDP_GMPL, plotTDP_HSCL, plotTDP_MonoL, plotTDP_ProgL, plotTDP_PromL, ncol=3, nrow=4)

ggsave(filename = "DMSO_Enrichment_plots.jpg", plot = grid, path ="~/Documents/Pipeline/new_cite-seq/Olga/GSEA_results")

# Save the grid of plots
ggsave(filename = "DMSO_GSEA_clusters.png", plot = grid, path ="~/Documents/Pipeline/new_cite-seq/Olga/GSEA_results")
```
# Load necessary libraries
library(fgsea)
library(ggplot2)
library(gridExtra)

# Define the number of clusters
num_clusters <- 10

# Initialize a list to store the plots
plotList <- list()

# Read the GMT file into R
C2_pathways <- gmtPathways("~/Documents/Pipeline/new_cite-seq/Olga/GSEA_results/c2.all.v2024.1.Hs.symbols.gmt")

# Loop over the clusters
for (i in 0:(num_clusters-1)) {
  # Dynamically get the filtered data frame for the current cluster
  filtered_df <- get(paste0("DMSO_scRNAseq_markers.cluster", i))
  
  # Add a small amount of random noise to break ties
  filtered_df$avg_log2FC <- filtered_df$avg_log2FC + runif(nrow(filtered_df), min = -1e-10, max = 1e-10)
  
  # Run fgsea
  fgsea_input <- setNames(object = filtered_df$avg_log2FC, nm = rownames(filtered_df))
  fgseaRes <- fgsea(pathways = C2_pathways, 
                    stats    = fgsea_input,
                    eps      = 0.0,
                    minSize  = 15,
                    maxSize  = 500)
  
  # Create the plot
  topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
  topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
  plot <- plotGseaTable(C2_pathways[topPathways], fgsea_input, fgseaRes, 
                gseaParam=0.5) + labs(title = paste0("DMSO Cluster ", i))
  
  # Save the plot
  ggsave(filename = paste0("DMSO_Enrichment_plot_cluster", i, ".png"), plot = plot, path ="~/Documents/Pipeline/new_cite-seq/Olga/GSEA_results/C2_Pathway_Results")
}



