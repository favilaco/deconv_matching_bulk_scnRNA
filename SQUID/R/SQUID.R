#' Single-cell RNA Quantity Informed Deconvolution (SQUID)!
#'
#' This is an executable R function named 'SQUID'
#' which analyse cell-type composition predictions based on concurrent RNA-Seq and scRNA-Seq profiles.
#'
#' You can learn more about SQUID compare to some other approaches at:
#'
#'   https://github.com/favilaco/deconv_matching_bulk_scnRNA
#'
#' Citation:
#' Avila Cobos, F., Najaf Panah, M. J., Epps, J., Long, X., Man, T. K., Chiu, H. S., ..., Mestdagh, P & Sumazin, P. (2022). Effective
#' methods for bulk RNA-Seq deconvolution using scnRNA-Seq transcriptomes. bioRxiv, 2022-12.
#'
#' GitHub: https://github.com/favilaco/deconv_matching_bulk_scnRNA
#'
#' Licence: MIT License
#'
#' Copyright (c) [2023] [Francisco Avila Cobos, Mohammad Javad Najaf Panah, Pieter Mestdagh, Pavel Sumazin]
#'
#'
#' Permission is hereby granted, free of charge, to any person obtaining a copy
#' of this software and associated documentation files (the "Software"), to deal
#' in the Software without restriction, including without limitation the rights
#' to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#' copies of the Software, and to permit persons to whom the Software is
#' furnished to do so, subject to the following conditions:
#'
#' The above copyright notice and this permission notice shall be included in
#' all copies or substantial portions of the Software.
#'
#' THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#' IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#' FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#' AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#' LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#' OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
#' THE SOFTWARE.
#'
#'
#' For commercial use inquiries, please contact the corresponding authors [Pavel Sumazin at Pavel.Sumazin@bcm.edu; Pieter Mestdagh at pieter.mestdagh@ugent.be].
#'
#'

# Load the example dataset:
load("SQUID_Toy_Example_Data.RData", verbose = T)

#' Single-cell RNA Quantity Informed Deconvolution (SQUID).
#'
#' SQUID executes a combination of RNA-Seq transformation and dampened weighted least-squares deconvolution approaches in predicting the composition of
#' cell mixtures and tissue samples.
#'
#' SQUID is an improved deconvolution method based on concurrent RNA-Seq and scRNA-Seq profiles.
#'
#' NOTE: To systematically test the benefit of bulk transformation and deconvolution with SQUID, a leave-one-out cross-validation strategy can also be used:
#' iteratively, concurrent RNA-Seq and scnRNA-Seq profiles of all but one of the samples were used to predict the composition of the remaining sample based
#' on its bulk RNA-Seq profile.
#'
#' @param B           A bulk RNA-seq numeric matrix which could be either count or normalized values. Rows should be genes and columns should be sample.ids.
#' @param scC         A single-cell numeric matrix which could be either count or normalized values. Rows should be genes and columns should be cell.ids.
#' @param scMeta      A single-cell annotation data.frame which rows should be cell.ids and columns including cell.id, sample.id, cluster.id, cellType, etc.
#' @param pB          (Optional) A pseudo-bulk matrix generated based on cluster.id/cellType from single-cell analysis. Rows should be genes and columns should be cluster.id/cellType.
#' @param P           (Optional) A numeric matrix of expected cell fractions including the composition of cluster.id/cellTypes per samples can be calculated from scRNA-seq analysis. Rows should be unique cluster.ids/cellTypes and columns should be sample.ids.
#' @param LeaveOneOut A logical variable which allow users run SQUID with/without a leave-one-out cross validation strategy. Default is FALSE.
#'
#' @return            A table including cellType, sample.id, observed_fraction (predicted cluster/cell-type fraction), expected_fraction (P) for each bulk sample.
#'
#' @export RESULTS    Save the RESULTS.rds and RESULTS.csv files in your working directory.
#'
#' @usage             RESULTS <- SQUID(B=B, scC=scC , scMeta=scMeta, pB=NULL, P=NULL, LeaveOneOut=FALSE)
#'
#'
SQUID <- function(B=B, scC=scC , scMeta=scMeta, pB=NULL, P=NULL, LeaveOneOut=FALSE) {
  cat('



              ----------------------------------- Welcome To SQUID ---------------------------------------

                                  ########    ########    #      #    ########    #######
                                  #           #      #    #      #       ##       #      #
                                  #           #      #    #      #       ##       #       #
                                  ########    #  ##  #    #      #       ##       #        #
                                         #    #    # #    #      #       ##       #       #
                                         #    #     ####  #      #       ##       #      #
                                  ########    ########    ########    ########    #######


               Execute an Effective method for bulk RNA-seq deconvolution using scnRNA-seq transcriptomes

              --------------------------------------------------------------------------------------------



')
  cat('SQUID just began the analysis ... \n\n')

  # Make pseudo-Bulk (pB) based on clusters from single-cell analysis, if not provided
  if(is.null(pB)){
    cellType <- scMeta$cellType # could be 'scMeta$cluster.id' as well.
    group = list()
    for(i in unique(cellType)){
      group[[i]] <- which(cellType %in% i)
    }
    pB = lapply(group,function(x) base::rowMeans(scC[,x, drop = FALSE]))
    pB = do.call(cbind.data.frame, pB)
  } else {
    pB=pB
  }

  # Make proportional table (P), if not provided
  if(is.null(P)){
    P <- table(scMeta$sample.id, scMeta$cellType) %>% t()
    P <- P[,gtools::mixedsort(colnames(P))]
    P <- P[,colnames(P) %in% colnames(B)]
    P <- apply(P, 2, function(x) x / sum(x))
  } else {
    P=P
  }

  # Unify row names between B and pB
  perc.overlap <- (length(intersect(rownames(pB),rownames(B))) * 60)/min(nrow(B),nrow(pB))
  if(perc.overlap < 33.3){
    if(length(grep("ENSG000|ENSMUSG000",rownames(B))) > 100){
      rownames.B <- ENSG_to_HGNC$HGNC[match(rownames(B), ENSG_to_HGNC$ENSG_ID)]
      B = B[!is.na(rownames.B),]
      rownames.B = rownames.B[!is.na(rownames.B)]
    } else {
      rownames.B <- ENSG_to_HGNC$ENSG_ID[match(rownames(B), ENSG_to_HGNC$HGNC)]
      B = B[!is.na(rownames.B),]
      rownames.B = rownames.B[!is.na(rownames.B)]
    }
    B <- collapsing(B, rownames.B)
  }

  # To avoid issues with sample.ids that are fully integers
  if(length(grep("[A-Z|a-z]", colnames(B))) == 0){colnames(B) <- paste("X",colnames(B),sep="")}
  if(length(grep("[A-Z|a-z]", colnames(P))) == 0){colnames(P) <- paste("X",colnames(P),sep="")}

  # Remove NA values (if any) after row.name interconversion
  if(sum(is.na(rownames(B))) != 0){B = B[!is.na(rownames(B)), ]}

  # Transform the bulk data
  if(LeaveOneOut) {
    cat("SQUID is transforming the Bulk samples using Leave-One-Out cross validation method ... ")
  } else {
    cat("SQUID is transforming the Bulk samples using all informations ... ")
  }

  Z = as.matrix(pB)
  Y = Z %*% P
  B.new = transformation2(X = B, Y = Y, leave.one.out = LeaveOneOut)
  B.new[!is.finite(B.new)] <- 0

  cat(' DONE!\n\n')

  # take common genes
  Genes <- intersect(rownames(Z),rownames(B.new))
  B.new <- as.matrix(B.new[Genes,])
  Z <- as.matrix(Z[Genes,])

  # solve Dampened WLS
  cat('SQUID is solving the dampened weighted least squares model for Bulk samples ... \n\n')
  RESULTS <- apply(B.new, 2, function(x){
    b = setNames(x, rownames(B.new))
    tr <- trimData(Z, b)
    RES <- t(solveDampenedWLS(S = Z, B = x))
  })

  rownames(RESULTS) <- colnames(Z)
  RESULTS = apply(RESULTS,2,function(x) ifelse(x < 0, 0, x)) # explicit non-negativity constraint
  RESULTS = apply(RESULTS,2,function(x) x/sum(x)) # explicit STO constraint
  RESULTS = RESULTS[gtools::mixedsort(rownames(RESULTS)), , drop = FALSE]
  RESULTS = suppressMessages(reshape2::melt(RESULTS))
  colnames(RESULTS) <- c("cellType","sample.id","observed_fraction")
  RESULTS$cellType = as.character(RESULTS$cellType)
  RESULTS$sample.id = as.character(RESULTS$sample.id)
  if(!is.null(P)){
    P = P[gtools::mixedsort(rownames(P)),,drop = FALSE] %>% data.frame(., check.names = FALSE)
    P$cellType = rownames(P)
    P = suppressMessages(reshape2::melt(P))
    colnames(P) <- c("cellType","sample.id","expected_fraction")
    P$cellType = as.character(P$cellType)
    P$sample.id = as.character(P$sample.id)
    RESULTS = merge(RESULTS, P, by = c("cellType", "sample.id"), all = TRUE)
    RESULTS[is.na(RESULTS)] <- 0
    RESULTS$expected_fraction <- round(RESULTS$expected_fraction, 3)
    RESULTS$observed_fraction <- round(RESULTS$observed_fraction, 3)
  }
  saveRDS(RESULTS, "RESULTS.rds")
  write.csv(RESULTS, "RESULTS.csv")
  cat('\nSQUID has saved the results table in .rds and .csv formats.\n\n')
  cat('Congrats! SQUID has finished the analyses!\n')
  return(RESULTS)
} # The End
