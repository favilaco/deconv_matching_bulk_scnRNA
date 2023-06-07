\name{SQUID}

\alias{SQUID}

\title{Single-cell RNA Quantity Informed Deconvolution (SQUID)}

\description{SQUID predicts the composition of cell mixtures and tissue samples based on concurrent RNA-Seq and scnRNA-Seq profiles.}

\usage{SQUID(B=B, scC=scC , scMeta=scMeta, pB=NULL, P=NULL, LeaveOneOut=FALSE)}

\arguments{
  \item{B}{A bulk RNA-seq numeric matrix which could be either count or normalized values. Rows should be genes and columns should be samples.}
  \item{scC}{A single-cell numeric matrix which could be either count or normalized values. Rows should be genes and columns should be cell.id.}
  \item{scMeta}{A single-cell annotation data.frame which rows should be cell.ids and columns including cell.id, sample.id, cluster.id, cellType, etc.}
  \item{pB}{(Optional) A pseudo-bulk matrix generated based on cluster.id/cellType from single-cell analysis. Rows should be genes and columns should be cluster.id/cellType.}
  \item{P}{ (Optional) A numeric matrix of expected cell fractions including the composition of clusters/cellTypes per samples can be calculated from scRNA-seq analysis. Rows should be unique cluster.ids/cellTypes and columns should be samples.}
  \item{LeaveOneOut}{A logical variable which allow users run SQUID with/without a leave-one-out cross validation strategy. Default is FALSE. Note that the LeaveOneOut argument included here to give users an extra option to systematically test the benefit of bulk transformation and deconvolution accuracy by SQUID.}
}

\value{A matrix including cellType, sample.id, observed_fraction (predicted cluster/cell-type fraction), expected_fraction (P) for each bulk sample. By default, SQUID also saves the RESULTS.rds and RESULTS.csv files on the working directory.}


\examples{
RESULTS <- SQUID(B=B, scC=scC , scMeta=scMeta, pB=NULL, P=NULL, LeaveOneOut=FALSE)
}
