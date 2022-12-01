require(tidyverse)
require(data.table)
library(ggrastr)
require(viridis)
library(ggpointdensity)
library(egg)
require(Matrix)
require(matrixStats)
require(limma)
library(ggplot2)
library(RColorBrewer)
require(ggrepel)
require(ggpubr)
require(monocle3)
require(Seurat)
require(foreach)
require(R.utils)

args <- R.utils::commandArgs(trailingOnly=TRUE)

if(length(args) != 5){

  	print("Please check that all required parameters are indicated or are correct")
  	stop()
} 

set.seed(4)
options(digits = 3)
options(future.globals.maxSize = 2000 * 1024^2) #https://satijalab.org/seurat/archive/v3.0/future_vignette.html

ENSG_to_HGNC = readRDS("./ENSG_to_HGNC.rds")

source('./helper_functions.R')

##############################################################################
## Input files and parameters:

directory <- args[1] 
ref.normalization = args[2] # normalization choice for scC or C
mix.normalization = args[3] # normalization choice for T
method = args[4]

cell.type.info.present <- args[5] #if "yes", no need to apply unsupervised clustering  

##############################################################################
## Reading input into memory:

setwd(directory)

scC <- readRDS("./scC.rds") %>% as.matrix()
T = readRDS("./T.rds") %>% as.matrix()
phenoDataC = read.table("phenoDataC_clusters_wo_regrouping.txt", header = TRUE) # It must contain at least two columns: "cellID" & "SubjectName". "cellType" is optional (e.g. if cell.type.info.present = "no")

#Unless proportions are provided from IHC, it would be computed directly from scRNA-seq/snRNA-seq
if(!file.exists("./P.rds")){
	P = NULL
} else {
	P = readRDS("./P.rds")	
}

##############################################################################
## PRE-PROCESSING AND QUALITY CONTROL

# Filter to remove 1% of extreme values (top 0.5% and bottom 0.5%)
filterCells <- function(filterParam){
  cellsToRemove <- which(filterParam < quantile(filterParam, .005) | filterParam > quantile(filterParam, .995)) 
  cellsToRemove
}

libSizes <- colSums(scC)
gene_names <- rownames(scC)

if(length(grep("ENSG000|ENSMUSG000",rownames(scC))) > 100){

  mtID.Ensembl <- ENSG_to_HGNC$ENSG_ID[grep("^MT-|_MT-", ENSG_to_HGNC$HGNC, ignore.case = TRUE)]
  rbID.Ensembl <- ENSG_to_HGNC$ENSG_ID[grep("^RPL|^RPS|_RPL|_RPS", ENSG_to_HGNC$HGNC, ignore.case = TRUE)]
  mtID <- gene_names %in% mtID.Ensembl
  rbID <- gene_names %in% rbID.Ensembl

} else {

  mtID <- grepl("^MT-|_MT-", gene_names, ignore.case = TRUE)
  rbID <- grepl("^RPL|^RPS|_RPL|_RPS", gene_names, ignore.case = TRUE)

}

mtPercent <- colSums(scC[mtID, ])/libSizes
rbPercent <- colSums(scC[rbID, ])/libSizes

lapply(list(libSizes = libSizes, mtPercent = mtPercent, rbPercent = rbPercent), filterCells) %>% 
  unlist() %>% 
  unique() -> cellsToRemove

if(length(cellsToRemove) != 0){
  scC <- scC[,-cellsToRemove]
  phenoDataC <- phenoDataC[-cellsToRemove,]
}

# Keep only "detectable" genes: min(round(0.01 * ncol(scC), 10) cells (regardless of the group) have a read/UMI count different from 0
keep <- which(Matrix::rowSums(scC > 0) >= min(round(0.01 * ncol(scC)), 10)) #otherwise it removes many genes from Ellis dataset
scC = scC[keep,]

##############################################################################
## Unsupervised clustering (Monocle3) if metadata information not provided

if(cell.type.info.present == "no"){

	if(file.exists("phenoDataC_clusters_after_regrouping.txt")){

		phenoDataC = read.table("phenoDataC_clusters_after_regrouping.txt", header = TRUE)

	} else {

		cds <- monocle3::new_cell_data_set(scC)
		cds <- monocle3::preprocess_cds(cds, num_dim = 100, norm_method = "log", method = "PCA", scaling = TRUE)
		cds <- monocle3::reduce_dimension(cds, max_components = 2, umap.metric = "cosine", umap.fast_sgd = FALSE, preprocess_method = 'PCA') 
		cds <- monocle3::cluster_cells(cds, k = 6, resolution = NULL, partition_qval = 0.05, num_iter = 1) #with resolution = NULL, the parameter is determined automatically

		Monocle3_metadata <- data.frame(cellID = names(cds@clusters$UMAP$clusters), cellType = paste("cluster",  cds@clusters$UMAP$clusters, sep ="."))
		phenoDataC$clusterID <- Monocle3_metadata$cellType

		m1 <- monocle3::plot_cells(cds, color_cells_by = "cluster", group_label_size = 4)

		ggsave(paste("Monocle3_",basename(directory),".pdf",sep=""), gridExtra::grid.arrange(m1, nrow = 1), height = 6, width = 6, limitsize = FALSE)
		dev.off()

		## Make C based on clusters from unsupervised clustering
		cellType <- phenoDataC$clusterID
		group = list()

		for(i in unique(cellType)){ 
		  group[[i]] <- which(cellType %in% i)
		}

		C = lapply(group,function(x) Matrix::rowMeans(scC[,x])) #C should be made with the mean (not sum) to agree with the way markers were selected
		C = do.call(cbind.data.frame, C)

		###############
		# Regrouping if correlation > 0.95  

		correlation.matrix = cor(C)
		annotation.regrouped <- regroup.cor(correlation.matrix, correlation.threshold = 0.95)

		if(nrow(annotation.regrouped) > 0){

			pdf(paste("Monocle3_cor_full_", basename(getwd()), "_log.ls.norm.pdf", sep = ""), width = max(7, ncol(C)/2), height = max(7, ncol(C)/2))
				pheatmap::pheatmap(correlation.matrix, annotation = annotation.regrouped, scale = "none", display_numbers = TRUE, number_format = "%.2f")
			dev.off()

			counter = 1
			for(elem in unique(annotation.regrouped$new.group)){
				phenoDataC$clusterID[phenoDataC$clusterID %in% rownames(annotation.regrouped)[annotation.regrouped$new.group == elem]] <- paste("regrouped.cluster", counter, sep = ".")
				counter = counter + 1
			}

			phenoDataC$cellType <- phenoDataC$clusterID
			# re-order phenoDataC in function of scC
			phenoDataC = phenoDataC[phenoDataC$cellID %in% colnames(scC),]
			phenoDataC <- phenoDataC[match(colnames(scC),phenoDataC$cellID),]

		} 

		write.table(phenoDataC, file = "phenoDataC_clusters_after_regrouping.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

		###############
		# re-making UMAP plots with new information:

		cds@colData$new.cluster = phenoDataC$clusterID[match(cds@colData$cell, phenoDataC$cellID)]

		m2 <- monocle3::plot_cells(cds, color_cells_by = "new.cluster", group_label_size = 4, show_trajectory_graph = FALSE)

		ggsave(paste("MONOCLE.UMAP_REGROUPED_",basename(directory),".pdf",sep=""), gridExtra::grid.arrange(m2, nrow = 1), height = 6, width = 6, limitsize = FALSE)
		dev.off()

	}

} else {

	# re-order phenoDataC in function of scC
	phenoDataC = phenoDataC[phenoDataC$cellID %in% colnames(scC),]
	phenoDataC <- phenoDataC[match(colnames(scC),phenoDataC$cellID),]

}

# Create matrix C (made with the mean expression values for all cells across cell types)
cellType <- phenoDataC$cellType
group = list()

for(i in unique(cellType)){ 
  group[[i]] <- which(cellType %in% i)
}

# Generating C (raw) from the mean (raw) expression values across all cells from each cell type
C = lapply(group,function(x) Matrix::rowMeans(scC[,x])) 
C = do.call(cbind.data.frame, C) %>% as.matrix() 

##############################################################################
## Remove row duplicates

C = collapsing(C, rownames(C))
scC = collapsing(scC, rownames(scC))
scC.raw = scC # needed to run marker selection + deconvolution
T = collapsing(T, rownames(T))

##############################################################################
## Single-cell and bulk normalization 

C <- Scaling(matrix = C, option = ref.normalization)
scC <- Scaling(matrix = scC.raw, option = ref.normalization, phenoDataC = phenoDataC)

T <- Scaling(matrix = T, option = mix.normalization)

##############################################################################
## Seurat - FindAllMarkers (Wilcoxon) on TMM normalized scRNA-seq/snRNA-seq 
# MuSiC & DWLS: without markers

if(!method %in% c("MuSiC", "DWLS", "SQUID")){

	if(! file.exists("./markers_TMM_FC1.5_Seuratwilcox.rds")){

		scC <- Scaling(matrix = scC.raw, option = "TMM", phenoDataC = phenoDataC)

		scC <- Seurat::CreateSeuratObject(counts = scC)
		scC$groups <- as.character(phenoDataC$clusterID)
		Idents(object = scC) <- as.character(phenoDataC$clusterID) #Needed. See: https://github.com/satijalab/seurat/issues/1482

		markers <- Seurat::FindAllMarkers(scC, 
		                                 logfc.threshold = log2(1.5), 
		                                 grouping.var = "groups",
		                                 test.use = "wilcox") 

		markers <- markers[,c("gene","avg_log2FC","cluster")]
		markers$AveExpr <- Inf
		colnames(markers) <- c("gene","log2FC","cell_type","AveExpr")
		markers$HGNC <- ENSG_to_HGNC$HGNC[match(markers$gene, ENSG_to_HGNC$ENSG_ID)]
		markers$ENSG <- ENSG_to_HGNC$ENSG[match(markers$gene, ENSG_to_HGNC$HGNC)]
		if(sum(is.na(markers$HGNC)) == nrow(markers)){markers$HGNC <- NULL; markers$HGNC <- rownames(markers)}
		if(sum(is.na(markers$ENSG)) == nrow(markers)){markers$ENSG <- NULL; markers$ENSG <- rownames(markers)}

		#if markers$gene is both ENSG and HGNC, split:
		if(length(grep("--",markers$gene)) == nrow(markers)){ markers$gene <- strsplit(markers$gene, "--") %>% lapply(., function(x) x[1]) %>% unlist()}
		if(length(grep("__",markers$gene)) == nrow(markers)){ markers$gene <- strsplit(markers$gene, "__") %>% lapply(., function(x) x[1]) %>% unlist()}

		## if needed, put markers$gene in same gene format as rownames C & T:
		if(length(intersect(markers$gene, rownames(C))) == 0){

		    if(length(grep("ENSG000|ENSMUSG000",markers$gene)) > 100){

		        markers.gene <- ENSG_to_HGNC$HGNC[match(markers$gene, ENSG_to_HGNC$ENSG_ID)] 

		    } else {

		        markers.gene <- ENSG_to_HGNC$ENSG_ID[match(markers$gene, ENSG_to_HGNC$HGNC)] 

		    }

		    markers$gene <- markers.gene

		}

		no.dups <- which(table(markers$gene)[markers$gene] == 1) %>% as.numeric()
		dups <- which(table(markers$gene)[markers$gene] > 1) #keeps ALL occurrences in original table + original ordering
		dups <- sapply(unique(names(dups)), function(x) grep(x, markers$gene)[1]) %>% as.numeric() #keeps ONLY FIRST occurrence of duplicates (in original table) + original ordering
		markers <- markers[c(no.dups, dups),]

		saveRDS(markers, "./markers_TMM_FC1.5_Seuratwilcox.rds")

	} else {

		markers = readRDS("./markers_TMM_FC1.5_Seuratwilcox.rds")

	}

} else { 

    markers = ""

}

##############################################################################
## Deconvolution

STRING <- paste(basename(directory), method, ref.normalization, mix.normalization, sep = "_")

output.filename <- paste(STRING, ".rds", sep = "")
output_folder <- paste(directory, "OUTPUT_DECONVOLUTION/", sep = "/")
ifelse(!dir.exists(output_folder), dir.create(output_folder), FALSE)
output_name <- paste(output_folder, "RESULTS_", output.filename, sep = "")

if(is.null(P)){

	P <- table(phenoDataC$SubjectName,phenoDataC$cellType) %>% t() %>% as.data.table() %>% data.table::dcast.data.table(., V1 ~ V2, fill = "N")
	P <- data.frame(P, check.names = FALSE) 
	rownames(P) <- P$V1
	P$V1 <- NULL
	P <- P[,colnames(P) %in% colnames(T)]
	P <- apply(P, 2, function(x) x / sum(x))

}

## Unify row names between T and C (after normalization!)
perc.overlap <- (length(intersect(rownames(C),rownames(T))) * 100)/min(nrow(T),nrow(C))
if(perc.overlap < 33.3){

    if(length(grep("ENSG000|ENSMUSG000",rownames(T))) > 100){

        rownames.T <- ENSG_to_HGNC$HGNC[match(rownames(T), ENSG_to_HGNC$ENSG_ID)] 
        T = T[!is.na(rownames.T),]
        rownames.T = rownames.T[!is.na(rownames.T)]

    } else {

        rownames.T <- ENSG_to_HGNC$ENSG_ID[match(rownames(T), ENSG_to_HGNC$HGNC)] 
        T = T[!is.na(rownames.T),]
        rownames.T = rownames.T[!is.na(rownames.T)]

    }

    T <- collapsing(T, rownames.T)

}

# To avoid issues with sample names that are fully integers
if(length(grep("[A-Z|a-z]", colnames(T))) == 0){colnames(T) <- paste("X",colnames(T),sep="")}
if(length(grep("[A-Z|a-z]", colnames(P))) == 0){colnames(P) <- paste("X",colnames(P),sep="")}

## Remove NA values (if any) after row.name interconversion
if(sum(is.na(rownames(T))) != 0){T = T[!is.na(rownames(T)), ]}

#To avoid MAST-DWLS to be repeated!
if(method == "DWLS"){STRING <- paste("DWLS", ref.normalization, sep = "_")} 

if(method %in% c("CIBERSORT", "RLR", "FARDEEP", "nnls")){

	RESULTS = Deconvolution(T = T, 
	                        C = C, 
	                        method = method,
	                        phenoDataC = phenoDataC,
	                        P = P, 
							STRING = STRING,
	                        marker_distrib = markers, 
	                        refProfiles.var = NULL)
	
} else {

	RESULTS = Deconvolution(T = T, 
                        C = scC, 
                        method = method,
                        phenoDataC = phenoDataC,
                        P = P, 
						STRING = STRING,
                        marker_distrib = markers, 
                        refProfiles.var = NULL)

}

RESULTS$dataset <- basename(directory)
RESULTS$method <- method
RESULTS$ref.normalization <- ref.normalization
RESULTS$mix.normalization <- mix.normalization

ann_text.global = RESULTS %>% dplyr::summarise(RMSE = sqrt(mean((observed_fraction - expected_fraction)^2)) %>% round(.,3), 
									   		   Pearson = cor(observed_fraction, expected_fraction) %>% round(.,3))

label = paste("Pearson = ", ann_text.global$Pearson, "\n RMSE = ", ann_text.global$RMSE, sep = "")

scaleFUN <- function(x) sprintf("%.3f", x)

P <- ggplot(RESULTS, aes(x = observed_fraction , y = expected_fraction, colour = cell_type)) +
      geom_point() +
      geom_abline(intercept = 0, slope = 1, lty = 2, col = "blue") +
      facet_grid( ~ method) +
      geom_text(vjust = "inward", hjust = "inward", data = ann_text.global, aes(x = 0.02, y = 1.04, label = label), size = 4, inherit.aes = FALSE)+
      ggtitle("all_normalizations_per_cell") +
      theme_bw() +
      xlim(0,1) +
      ylim(0,1.05) +
      xlab("observed fraction") +
      ylab("expected fraction") +
      ggtitle(STRING) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
           legend.position = "bottom")

ggsave(P, filename = paste(output_folder, paste(STRING, ".pdf", sep = ""), sep = "/"), width = 6, height = 6)
saveRDS(RESULTS, output_name)
