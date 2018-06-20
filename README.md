# SWNE
Similarity Weighted Nonnegative Embedding (SWNE), is a method for visualizing high dimensional datasets.
SWNE uses Nonnegative Matrix Factorization to decompose datasets into latent factors, projects
those factors onto 2 dimensions, and embeds samples and key features in 2 dimensions relative to the
factors. SWNE can capture both the local and global dataset structure, and allows
factors and relevant features to be embedded directly onto the visualization, facilitating interpretation
of the data.

You can find our bioRxiv preprint [here](https://www.biorxiv.org/content/early/2018/04/26/276261).

## Installation instructions

1. Install devtools with `install.packages("devtools")` if not already installed
2. liger with `devtools::install_github("JEFworks/liger")` for geneset enrichment functionality
3. Install swne with `devtools::install_github("yanwu2014/swne")`

## Walkthroughs and examples
Since SWNE is primarily meant for visualization and interpretation of the data, we typically use either [Seurat](http://satijalab.org/seurat/) or [Pagoda2](https://github.com/hms-dbmi/pagoda2) as a primary scRNA-seq pipeline. All the R markdown files used to generate the walkthroughs can be found under the Examples/ directory.

1. A basic [walkthrough](https://yanwu2014.github.io/swne/Examples/pbmc3k_swne_seurat.html) of 3k PBMC cells starting from a pre-computed Seurat object.
2. A basic [walkthrough](https://yanwu2014.github.io/swne/Examples/pbmc3k_swne_pagoda2.html) of 3k PBMC cells starting from a pre-computed Pagoda2 object.
3. A [walkthrough](https://yanwu2014.github.io/swne/Examples/multiple_pancreas_alignment_swne.html) using SWNE to visualize four pancreas datasets that have undergone batch alignment with Seurat's [manifold alignment](https://www.nature.com/articles/nbt.4096). The script used to generate the object can be found [here](https://yanwu2014.github.io/swne/Examples/multiple_pancreas_workflow.R) and the raw datasets can be found [here](http://bit.ly/IAexpmat).
4. A [walkthrough](https://yanwu2014.github.io/swne/Examples/dentate_gyrus_swne_velocyto.html) integrating [RNA velocity](https://www.biorxiv.org/content/early/2017/10/19/206052) with SWNE for a developmental mouse dentate gyrus dataset.

## Recreating Figures
To recreate the figures from our preprint, see the Scripts/ directory. 

To generate the simulated discrete and trajectory datasets, use `splatter_generate.R`. The simulated datasets we generated can be found [here](https://bit.ly/2JQDDNc)

To generate the visualizations and embedding evaluations, run `splatter_discrete_swne.R` and `splatter_trajectory_swne.R` for the discrete and trajectory simulations, respectively. To benchmark SWNE runtimes, use `splatter_runtime_analysis.R`.

The data needed to run `hemato_swne.R` can be found [here](https://bit.ly/2MFiByO) and the raw data for the hematopoietic cells can be found, courtesy of the monocle2 developers, [here](http://www.gs.washington.edu/~xqiu/proj2/RGE_analysis_data.tar.gz).

The data needed to run `snDropSeq_swne.R` on the cerebellar and visual cortex data can be found [here](https://bit.ly/2I6R5XL) and the raw data can be found at the GEO accession GSE97930.

The raw PBMC dataset can be found at the 10X genomics [website](https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/pbmc3k).
