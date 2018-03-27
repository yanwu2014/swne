# swne
Similarity Weighted Nonnegative Embedding (SWNE), is a method for visualizing high dimensional datasets.
SWNE uses Nonnegative Matrix Factorization to decompose datasets into latent factors, projects
those factors onto 2 dimensions, and embeds samples and key features in 2 dimensions relative to the
factors. SWNE can capture both the local and global dataset structure, and allows
factors and relevant features to be embedded directly onto the visualization, facilitating interpretation
of the data.

You can find our bioRxiv preprint [here](https://www.biorxiv.org/content/early/2018/03/05/276261.1).

## Installation instructions

1. Install devtools with `install.packages("devtools")` if not already installed
2. Install swne with `devtools::install_github("yanwu2014/swne")`

### Optionally install:
1. liger with `devtools::install_github("JEFworks/liger")` for additional geneset enrichment functionality

## Seurat walkthrough
For a quick example using the pbmc dataset from 10x genomics, see this [walkthrough](https://yanwu2014.github.io/swne/Examples/pbmc3k_swne_seurat.html) which generates SWNE visualizations using data processed with the [Seurat](http://satijalab.org/seurat/) pipeline. The pre-computed [Seurat object](https://github.com/yanwu2014/swne/blob/master/Examples/pbmc3k_seurat.Robj) can be downloaded from the Examples directory.

## Pagoda2 walkthrough
See this [walkthrough](https://yanwu2014.github.io/swne/Examples/pbmc3k_swne_pagoda2.html) to generate SWNE visualizations using data processed with the [Pagoda2](https://github.com/hms-dbmi/pagoda2) pipeline. The pre-computed [Pagoda2 object](https://github.com/yanwu2014/swne/blob/master/Examples/pbmc3k_pagoda2.Robj) can also be downloaded from the Examples directory.


## Recreating Figures
To recreate the figures from our preprint, see the Scripts/ directory. 

The raw data for the hematopoietic cells can be found, courtesy of the monocle2 developers, [here](http://www.gs.washington.edu/~xqiu/proj2/RGE_analysis_data.tar.gz). After untarring the directory,
cd to the Notebook/ sub-directory. The debatched expression matrix and informative genes can be found in
the Paul_Cell_MARSseq_GSE72857.RData file, and the metdata can be found in the GSE72857_experimental_design.txt
file.

The raw data for the neural cells can be found at the GEO accession GSE97930. We only used the cerebellar
(GSE97930_CerebellarHem_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt.gz) and the visual cortex
(GSE97930_VisualCortex_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt.gz) datasets. The metadata file can be 
found in the Scripts/ directory.

The raw PBMC dataset can be found at the 10X genomics [website](https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/pbmc3k).
