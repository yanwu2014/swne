# swne
Similarity Weighted Nonnegative Embedding (SWNE), is a method for visualizing high dimensional datasets.
SWNE uses Nonnegative Matrix Factorization to decompose datasets into latent factors, projects
those factors onto 2 dimensions, and embeds samples and key features in 2 dimensions relative to the
factors. SWNE can capture both the local and global dataset structure, and allows
factors and relevant features to be embedded directly onto the visualization, facilitating interpretation
of the data.

You can find our bioRxiv preprint here: https://www.biorxiv.org/content/early/2018/03/05/276261.1

## Installation instructions

1. Install devtools with `install.packages("devtools")` if not already installed
2. Install swne with `devtools::install_github("yanwu2014/swne")`

### Optionally install:
1. liger with `devtools::install_github("JEFworks/liger")` for additional geneset enrichment functionality

## Usage
We highly recommend using SWNE with either Seurat (http://satijalab.org/seurat/) or Pagoda2 (https://github.com/hms-dbmi/pagoda2), two general single cell RNA-seq analysis pipelines. 

For a quick example using the pbmc dataset from 10x genomics, see pbmc3k_swne.R and the pre-computed
Seurat object under the Examples directory. The SWNE plot for the pbmc dataset with key marker genes embedded
is shown below.

<img src="Figures/pbmc3k_swne_plot.png" width="550" height="550" />

## Recreating Figures
To recreate the figures from our preprint, see the Scripts/ directory. 

The raw data for the hematopoietic cells can be found, courtesy of the monocle2 developers, at
http://www.gs.washington.edu/~xqiu/proj2/RGE_analysis_data.tar.gz. After untarring the directory,
cd to the Notebook/ sub-directory. The debatched expression matrix and informative genes can be found in
the Paul_Cell_MARSseq_GSE72857.RData file, and the metdata can be found in the GSE72857_experimental_design.txt
file.

The raw data for the neural cells can be found at the GEO accession GSE97930. We only used the cerebellar
(GSE97930_CerebellarHem_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt.gz) and the visual cortex
(GSE97930_VisualCortex_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt.gz) datasets. The metadata file can be 
found in the Scripts/ directory.

The raw PBMC dataset can be found at the 10X genomics website: https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/pbmc3k.
