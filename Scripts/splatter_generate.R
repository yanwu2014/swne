#### Generate discrete clusters using Splatter ####

library(splatter)
library(pagoda2)

setwd("/media/Scratch_SSD/Yan/R_Analysis/swne-dev/simulations")

## Fit params from pbmc3k dataset
pbmc3k.path <- "/media/Scratch_SSD/Yan/R_Analysis/swne-dev/pbmc3k_matrix/"
pbmc3k.data <- read.10x.matrices(pbmc3k.path)

splat.seed <- 3250879
init.params <- newSplatParams(group.prob = c(0.45, 0.15, 0.15, 0.15, 0.1),
                              de.prob = c(0.3, 0.15, 0.15, 0.15, 0.3),
                              de.facLoc = 0.3, de.facScale = 0.6,
                              seed = splat.seed)
est.params <- splatEstimate(as.matrix(pbmc3k.data), params = init.params)

## Generate synthetic data using fitted parameters
splat.groups <- splatSimulateGroups(est.params, verbose = F)
counts <- as(splat.groups@assays@.xData$data$counts, "dgCMatrix")

metadata <- splat.groups@colData@listData
colnames(counts) <- paste(colnames(counts), metadata$Group, sep = "_")
# write.table(counts, file = "splatter.discrete.counts.tsv", sep = "\t")

## Use pagoda2 for data processing
library(pagoda2)

r <- Pagoda2$new(counts, modelType = 'plain', trim = 3, log.scale = T)
r$adjustVariance(plot = F, do.par = F, gam.k = 10)
r$calculatePcaReduction(nPcs = 30, n.odgenes = 3e3, maxit = 1000)
r$getEmbedding(type = 'PCA', embeddingType = 'tSNE', perplexity = 30, verbose = F)

clusters <- factor(metadata$Group); names(clusters) <- colnames(counts);
r$plotEmbedding(type = 'PCA', embeddingType = 'tSNE', groups = clusters)

## Save results
save(list = c("r", "metadata"), file = "splatter.discrete.RData")


#### Generate trajectories using splatter ####

library(splatter)
setwd("/media/Scratch_SSD/Yan/R_Analysis/swne-dev/simulations")

splat.seed <- 3250879
splat.params <- newSplatParams(seed = splat.seed,
                               group.prob = c(0.3,0.3,0.2,0.2),
                               de.prob = 0.3, de.facLoc = 0.4, de.facScale = 0.6,
                               path.from = c(0,1,1,3), 
                               path.length = c(100,200,50,50),
                               path.skew = c(0.5,0.5,0.5,0.5))

## Fit params from hematopoeisis dataset
hemato.path <- "/media/Scratch_SSD/Yan/R_Analysis/swne-dev/hemato_expr_debatched.tsv"
hemato.data <- as.matrix(read.table(hemato.path, sep = "\t"))

splat.seed <- 3250879
splat.params <- splatEstimate(hemato.data, params = splat.params)

## Simulate trajectory dataset
splat.groups <- splatSimulatePaths(splat.params, verbose = F)

counts <- as.matrix(splat.groups@assays@.xData$data$counts)
metadata <- data.frame(splat.groups@colData@listData)

colnames(counts) <- paste(colnames(counts), metadata$Group, sep = "_")
rownames(metadata) <- colnames(counts)


## Use pagoda2 for data processing
library(pagoda2)

r <- Pagoda2$new(counts, modelType = 'plain', trim = 2, log.scale = T)
r$adjustVariance(plot = F, do.par = F, gam.k = 10)
r$calculatePcaReduction(nPcs = 20, n.odgenes = 1e3, maxit = 1000)
r$getEmbedding(type = 'PCA', embeddingType = 'tSNE', perplexity = 40, verbose = F)

clusters <- factor(metadata$Group); names(clusters) <- colnames(counts);
r$plotEmbedding(type = 'PCA', embeddingType = 'tSNE', groups = clusters)

## Save results
save(list = c("r", "metadata"), file = "splatter.trajectory.RData")
