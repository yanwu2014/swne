#' Embeds genes based on promoter accessibility relative to topic coordinates
#'
#' @param swne.embedding SWNE embedding computed using RunSWNE
#' @param cisTopicObject cisTopic object
#' @param genes.embed Genes to embed
#' @param peaks.use Subset of peaks to use
#' @param alpha.exp Increasing alpha.exp increases how much the factors "pull" the features
#' @param n_pull Number of factors pulling on each feature. Must be >= 3
#' @param scale.cols Whether or not to scale the input columns to 0 - 1
#' @param overwrite Whether or not to overwrite any existing feature embedding
#'
#' @return swne.embedding with feature coordinates (feature.coords)
#'
#' @import org.Hs.eg.db
#' @import TxDb.Hsapiens.UCSC.hg38.knownGene
#' @import cisTopic
#'
#' @export
#'
EmbedPromoters <- function(swne.embedding, cisTopicObject, genes.embed, peaks.use = NULL, alpha.exp = 1, n_pull = 3,
                           scale.cols = T, overwrite = T) {
  cisTopicObject <- getRegionsScores(cisTopicObject)
  cisTopicObject <- annotateRegions(cisTopicObject, txdb = TxDb.Hsapiens.UCSC.hg38.knownGene,
                                    annoDb = "org.Hs.eg.db")
  regions.anno <- cisTopicObject@region.data

  if (!is.null(peaks.use)) regions.anno <- regions.anno[peaks.use,]
  regions.anno <- subset(regions.anno, SYMBOL %in% genes.embed & grepl("Promoter", annotation))
  genes <- factor(regions.anno$SYMBOL)

  regions.scores <- as.matrix(regions.anno[,grepl("Scores", colnames(regions.anno))])
  regions.scores <- apply(regions.scores, 2, function(x) tapply(x, genes, mean))
  missing.genes <- genes.embed[!genes.embed %in% rownames(regions.scores)]
  if (length(missing.genes) > 0) {
    warning(paste(paste(missing.genes, collapse = ",", sep = ","),
                  "have no peaks in their promoter regions"))
  }

  swne.emb$H.coords$name <- ""
  rownames(regions.scores) <- paste0(rownames(regions.scores), "-pr")
  swne.emb <- EmbedFeatures(swne.emb, regions.scores, alpha.exp = alpha.exp, scale.cols = scale.cols,
                            n_pull = n_pull, overwrite = overwrite)

  swne.emb
}
