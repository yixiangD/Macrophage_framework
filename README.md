# Macrophage framework
Macrophages populate every organ during homeostasis and disease, displaying features of tissue imprinting and heterogeneous activation. The disjointed picture of macrophage biology that has emerged from these observations can be difficult to integrate across models or with in vitro macrophage activation paradigms. For these reasons, we set out to contextualize macrophage heterogeneity across tissues and inflammatory conditions, specifically aiming to define a common framework of macrophage activation. We built a predictive model with which we mapped the activation of macrophages across 12 tissues and 25 biological conditions, finding a strikingly common and finite number of transcriptional profiles which we modelled as stages along 4 conserved activation paths. We verified this model with adoptive cell transfer experiments and identified transient RELMɑ expression as a feature of macrophage tissue engraftment. We propose that this integrative approach of macrophage classification allows the establishment of a common predictive framework of macrophage activation that may serve to contextualize the future study of these cells.

An interactive shiny application to visualize critical aspects of the data in this publication can be accessed at https://t.jh.edu/macrophage-framework. scRNAseq datasets generated specifically for this publication may be retrieved from publicly available repositories with no restrictions on their use, under the following accession numbers: helminth infection of adipose tissue (GSE157313), bacterial infection of adipose tissue (GSE171328), High fat diet lamina propria (GSE171330) and skin wound (GSE183489).

This repository contains all relevant code to reproduce the findings and figures in the associated publication. Please cite us!

## Learning scRNA-seq using R
[Tutorial 1](https://nbisweden.github.io/excelerate-scRNAseq/session-qc/Quality_control.html)

### Quality control
1. filtering cells with lower and upper bounds on RNA expression and mitochondria expression


### Batch effect correction and integration
- [Seurat batch integration tutorial](https://satijalab.org/seurat/articles/integration_introduction.html)
- [Benchmark paper recommending Harmony](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1850-9)

## Reference papers
- Normalizing single-cell RNA sequencing data: challenges and opportunities, 2017 Nature Methods

