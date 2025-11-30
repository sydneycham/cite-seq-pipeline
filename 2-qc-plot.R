suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tibble))

# import data -----------------------------------------------------------------
md = lapply(unlist(snakemake@input), fread) %>%
	bind_rows() %>%
	mutate(type = factor(type, levels = c("BEFORE", "AFTER")))

# visualize cell counts before and after filtration ---------------------------
cell_counts = md %>%
  group_by(orig.ident, type) %>%
  summarize(n = n())

p = cell_counts %>%
  ggplot(aes(x = orig.ident, y = n, fill = type, label = n, group = type)) +
    geom_col(position = "dodge") +
    geom_text(position = position_dodge(width = .9), vjust = -0.9) +
    ggtitle("Cells before and after QC") +
    theme(
      text = element_text(size = 12), 
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

ggsave(snakemake@output$cell_counts, p, width = 16, height = 9, dpi = 600)

# visualize QC metrics with cutoffs
anno = c("orig.ident", "type")

p = md %>%
  select(-barcode) %>%
  tidyr::pivot_longer(!all_of(anno)) %>%
  ggplot(aes(x = orig.ident, y = value, fill = type, col = type)) +
    # geom_violin() +
		geom_boxplot() +
    facet_wrap(name ~ ., scales = "free") +
    theme(
      text = element_text(size = 12), 
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

ggsave(snakemake@output$qc_metrics, p, width = 16, height = 9, dpi = 600)