# SWNE
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

## Walkthroughs and examples
Since SWNE is primarily meant for visualization and interpretation of the data, we typically use either [Seurat](http://satijalab.org/seurat/) or [Pagoda2](https://github.com/hms-dbmi/pagoda2) as a primary scRNA-seq pipeline. All the R markdown files and pre-computed objects used to generate the walkthroughs can be found under the Examples/ directory.

1. A basic [walkthrough](https://yanwu2014.github.io/swne/Examples/pbmc3k_swne_seurat.html) of 3k PBMC cells starting from a pre-computed [Seurat object](https://github.com/yanwu2014/swne/blob/master/Examples/pbmc3k_seurat.Robj).
2. A basic [walkthrough](https://yanwu2014.github.io/swne/Examples/pbmc3k_swne_pagoda2.html) of 3k PBMC cells starting from a pre-computed [Pagoda2 object](https://github.com/yanwu2014/swne/blob/master/Examples/pbmc3k_pagoda2.Robj)
3. An [walkthrough](https://yanwu2014.github.io/swne/Examples/dentate_gyrus_swne_velocyto.html) integrating [RNA velocity](https://www.biorxiv.org/content/early/2017/10/19/206052) with SWNE for a developmental mouse dentate gyrus dataset. The pre-computed Pagoda2 and Velocyto objects can be found [here](https://github.com/yanwu2014/swne/blob/master/Examples/dentate_gyrus.p2.velocyto.RData.gz)

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
