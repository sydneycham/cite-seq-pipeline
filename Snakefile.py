__author__ = "Garth Kong"
__license__ = "MIT"
__email__ = "kongg@ohsu.edu"

"""
Analyze single-cell CITE-Seq (RNA + Protein) data after CellRanger read alignment
"""

import os
import pandas as pd

configfile: 'config/config.yaml'
sheet = pd.read_table("config/samplesheet.tsv")
sheet.index = sheet["sample"]
SAMPLES = sheet["sample"]

rule all:
    input:
        # import cells
        expand("data/import/{sample}.rds", sample = SAMPLES),
        
        # quality control
        expand("data/qc/{sample}/{sample}.rds", sample = SAMPLES),
        "data/qc/cell-counts.png",
        "data/qc/qc-metrics.png",

        # SCTransform - individual and group
        expand("data/normalize/{sample}.rds", sample = SAMPLES),
        "data/merge/merged_object.rds",

        # weighted nearest neighbor
        "data/integrate/wnn/integrated.rds"

def get_CR(wildcards):
    return(sheet.loc[wildcards.sample, "cellranger_path"])

rule import:
    input:
        get_CR
    output:
        "data/import/{sample}.rds"
    conda:
        "workflow/envs/signac.yaml"
    script:
        "workflow/src/1-import.R"

rule qc:
    input:
        rds = "data/import/{sample}.rds",
        samplesheet = "config/samplesheet.tsv"
    output:
        rds = "data/qc/{sample}/{sample}.rds",
        md = "data/qc/{sample}/{sample}-md.txt",
        metrics = "data/qc/{sample}/{sample}-metrics.txt"
    conda:
        "workflow/envs/seurat.yaml"
    script:
        "workflow/src/2-qc.R"

rule plot_qc:
    """
    Plot all before/after filtering of all samples together.
    """
    input:
        expand("data/qc/{sample}/{sample}-md.txt", sample = SAMPLES)
    output:
        cell_counts = "data/qc/cell-counts.png",
        qc_metrics = "data/qc/qc-metrics.png"
    conda:
        "workflow/envs/seurat.yaml"
    script:
        "workflow/src/2-qc-plot.R"

rule individual_normalize:
    input:
        sample = "data/qc/{sample}/{sample}.rds",
        samplesheet = "config/samplesheet.tsv"
    output:
        rds = "data/normalize/{sample}.rds",
        umap = "data/normalize/{sample}-umap.png",
        scree = "data/normalize/{sample}-scree-plot.png"
    conda:
        "workflow/envs/seurat.yaml"
    script:
        "workflow/src/3-normalize.R"

rule group_normalize:
    input:
        samples = expand("data/qc/{sample}/{sample}.rds", sample = SAMPLES)
    output:
        rds = "data/merge/merged_object.rds",
        rna_scree = "data/merge/rna-scree-plot.png",
        adt_scree = "data/merge/adt-scree-plot.png",
        pca = "data/merge/pca.png"
    conda:
        "workflow/envs/seurat.yaml"
    script:
        "workflow/src/4-merge.R"

rule wnn:
    """
    Weighted Nearest Neighbors for sample integration
    """
    input:
        merged = "data/merge/merged_object.rds",
        samplesheet = "config/samplesheet.tsv"
    output:
        rds = "data/integrate/wnn/integrated.rds",
        wnn_umap = "data/integrate/wnn/wnn-umap.png",
        rna_umap = "data/integrate/wnn/rna-umap.png",
        adt_umap = "data/integrate/wnn/adt-umap.png"
    conda:
        "workflow/envs/seurat.yaml"
    script:
        "workflow/src/5-wnn.R"
