# SQUID: Single-cell RNA Quantity Informed Deconvolution

# Overview
SQUID is an R package for conducting (tumor) deconvolution analyses. SQUID executes a combination of RNA-Seq transformation and dampened weighted least-squares deconvolution to predict the composition of cell mixtures and tissue samples based on the concurrent RNA-Seq and scnRNA-Seq profiles. SQUID harnesses the power of concurrent RNA-Seq and scnRNA-Seq profiling, outperforming other methods in predicting the composition of cell mixtures and tissue samples.

# Installation
First load the requirements: 
```r
library(devtools)
devtools::source_url("https://github.com/favilaco/deconv_matching_bulk_scnRNA/blob/master/helper_functions.R?raw=TRUE")
library(dplyr)
```
Now, you can install the latest version of SQUID from GitHub using the devtools package:
```r
install_github("favilaco/deconv_matching_bulk_scnRNA/SQUID")
```
# Usage
To use SQUID, load the package into your R session using the library function:
```r
library(SQUID)
```
You can then use the provided SQUID function to conduct the experimantal analyses. You need to provide the required inputs as described in the package documentation in details. SQUID package contains a toy example dataset including bulk and single-cell RNA-seq simulated count data which you can use to run the program and check how it works: 
```r
RESULTS <- SQUID(B = B, scC = scC , scMeta = scMeta, pB = NULL, P = NULL, LeaveOneOut = FALSE)

#' @param B           A bulk RNA-seq numeric matrix which could be either count or normalized values. Rows should be genes and columns should be sample.ids.
#' @param scC         A single-cell numeric matrix which could be either count or normalized values. Rows should be genes and columns should be cell.ids.
#' @param scMeta      A single-cell annotation data.frame which rows should be cell.ids and columns including cell.id, sample.id, cluster.id, cellType, etc.
#' @param pB          (Optional) A pseudo-bulk matrix generated based on cluster.id/cellType from single-cell analysis. Rows should be genes and columns should be cluster.id/cellType.
#' @param P           (Optional) A numeric matrix of expected cell fractions including the composition of cluster.id/cellTypes per samples can be calculated from scRNA-seq analysis. Rows should be unique cluster.ids/cellTypes and columns should be sample.ids.
#' @param LeaveOneOut A logical variable which allow users run SQUID with/without a leave-one-out cross validation strategy. Default is FALSE.
#'
```
Notice: To improve the performance of SQUID prediction, you may need to set up your research-specific pipeline for the quality contol, normalization, imputation, and single-cell clustering procedures on your datasets prior to run SQUID. You are welcome to look at our preprint [paper](https://www.biorxiv.org/content/10.1101/2022.12.13.520241v2) as an example.

# License
This package is licensed under the MIT license. See the LICENSE file for details.

# Contact
For questions or comments about SQUID, please contact the package maintainers at Francisco Avila Cobos <francisco.avilacobos@ugent.be> and Mohammad Javad Najaf Panah <mohammadjavad.najafpanah@bcm.edu>. If you find a bug or have a feature request, please submit an issue on the GitHub repository.

# Citation
 Francisco Avila Cobos, Mohammad Javad Najaf Panah, Jessica Epps, Xiaochen Long, Tsz-Kwong Man, Hua-Sheng Chiu, Elad Chomsky, Evgeny Kiner, Michael J Krueger, Diego di Bernardo, Luis Voloch, Jan Molenaar, Sander R. van Hooff, Frank Westermann, Selina Jansky, Michele L. Redell, Pieter Mestdagh, Pavel Sumazin Effective methods for bulk RNA-Seq deconvolution using scnRNA-Seq transcriptomes (bioRxiv; https://www.biorxiv.org/content/10.1101/2022.12.13.520241v2)
