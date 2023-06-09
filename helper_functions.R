regroup.cor <- function(correlation.matrix, correlation.threshold = 0.95){

      top.cor <- which(abs(correlation.matrix) >= correlation.threshold & row(correlation.matrix) < col(correlation.matrix), arr.ind = TRUE)
      
      if(nrow(top.cor) != 0){

        ## reconstruct names from positions
        high_cor <- matrix(colnames(correlation.matrix)[top.cor], ncol = 2)

        count = 0; groups = list()

        for(elem in unique(high_cor[,1])){
          count = count + 1
          items = high_cor[high_cor[,1] %in% elem,]
          if(!is.null(dim(items))){items <- items %>% melt() %>% .[3] %>% .$value %>% unique()}
          groups[[paste("group",count,sep="_")]] <- items
        }

        for(elem in unique(high_cor[,2])){
          count = count + 1
          items = high_cor[high_cor[,2] %in% elem,]
          if(!is.null(dim(items))){items <- items %>% melt() %>% .[3] %>% .$value %>% unique()}
          groups[[paste("group",count,sep="_")]] <- items
        }

        # If one element present in multiple groups, keep largest
        groups <- groups[sapply(groups, length) %>% order() %>% rev()] #re-order NEEDED to ensure biggest groups and no repetition
        new.groups <- list()
        used <- c()
        count = 0
        for(i in names(groups)){
          count = count + 1
          to.regroup <- sapply(groups, function(x) groups[[i]] %in% x) %>% colSums() != 0
          to.regroup <- names(to.regroup)[to.regroup]
          if(sum(to.regroup %in% used) == 0){
            new.groups[[paste("group",count,sep="_")]] <- groups[to.regroup] %>% unlist() %>% unique()
            used <- c(used, to.regroup)
          }
          
        }

        annotation = melt(new.groups)
        rownames(annotation) <- annotation$value
        annotation$value <- NULL
        colnames(annotation) <- "new.group"

      } else {

        annotation = data.frame(new.group = "")

      }

      return(annotation)

}



collapsing <- function(input, duplicated.names){
    
    if(sum(duplicated(duplicated.names)) != 0){

      exp.var <- apply(input, 1, mean) %>% melt()
      exp.var$gene  <- duplicated.names
      exp.var$IQR <- apply(input, 1, IQR) %>% as.numeric()
      exp.var$row.number <- 1:nrow(exp.var)
      max.values = exp.var %>% 
        group_by(gene) %>% 
        summarise(across(c("value", "IQR"), ~ max(.x)))
      exp.var <- merge(max.values, exp.var, by = "gene")
      to.keep <- which(exp.var[,2] == exp.var[,4] & exp.var[,3] == exp.var[,5]) %>% sort()

      #if after this criteria, still duplicates, choose the first one (smallest row.number)!
      exp.var <- exp.var[to.keep,]
      min.row = exp.var %>% 
        group_by(gene) %>% 
        summarise(across(c("row.number"), ~ min(.x)))
      exp.var <- merge(min.row, exp.var, by = "gene")
      to.keep <- which(exp.var[,2] == exp.var[,7]) %>% sort()
      exp.var <- exp.var[to.keep,] %>% na.omit()
      input <- input[exp.var[,2], ]

      if(length(grep("ENS.*\\.[0-9]+", rownames(input))) > 1000){ # For datasets with ENSG names with format "ENSG00...003.14": keep genes with highest mean expression AND highest variability (first criterium insufficient, there are cases with same mean)
      
        rownames(input) <- strsplit(rownames(input),"\\.") %>% lapply(., function(x) x[1]) %>% unlist()

      } else if(length(grep("ENS.*\\__", rownames(input))) > 1000){  # For datasets with gene names as "ENSG00...__HGNC": remove the HGNC part

        rownames(input) <- strsplit(rownames(input),"__") %>% lapply(., function(x) x[1]) %>% unlist()

      } else {

        rownames(input) <- exp.var$gene

      }

    }

  return(input)

}



transformation2 <- function(X, Y, leave.one.out = TRUE) {

    # use the same genes for all input datasets
    Genes <- intersect(row.names(Y), row.names(X))

    X <- as.matrix(X[Genes,])
    Y <- as.matrix(Y[Genes,])

    if(leave.one.out){

        pred <- intersect(colnames(X), colnames(Y)) # matching samples with sc and bulk data
        X.new <- matrix(0, nrow=dim(X)[1], ncol=length(pred))

        for(j in 1:length(pred)){

            X.train <- as.matrix(X[,pred[-j]])
            X.test <- as.matrix(X[,pred[j]])
            Y.train <- as.matrix(Y[,pred])

            # track the mean and sd after leaving one out:
            X.mean <- rowMeans(X.train)
            X.sd <- apply(X.train,1,sd) 
            Y.mean <- rowMeans(Y)
            Y.sd <- apply(Y,1,sd)

            for (i in 1:dim(X)[1]){

                # transforming over genes
                x <- X.test[i]
                sigma_j <- Y.sd[i]*sqrt((length(Y[i,])-1)/(length(Y[i,])+1))
                x.new <- (x-X.mean[i])/X.sd[i]
                x.new <- x.new*sigma_j + Y.mean[i]

                X.new[i,j] <- x.new

            }

        }

    } else {

        # transforming over genes
        X.new <- matrix(0, nrow=dim(X)[1], ncol=dim(X)[2])

        for (i in 1:dim(X)[1]){

            y <- Y[i,]
            x <- X[i,]
            l <- length(y)
            sigma_j <- sd(y)*sqrt((l-1)/(l+1))
            x.new <- (x-mean(x))/sd(x)
            x.new <- x.new*sigma_j + mean(y)
            X.new[i,] <- x.new

        }

    }

    # explicit non-negativity constraint:
    X.new = apply(X.new,2,function(x) ifelse(x < 0, 0, x))

    rownames(X.new) = rownames(X); colnames(X.new) = colnames(X)

    return(X.new)

}



bulkC.fromSC <- function(scC, phenoDataC){

    phenoDataC = phenoDataC[colnames(scC),]
    cellType <- phenoDataC$cellType
    group = list()

    for(i in unique(cellType)){ 

        group[[i]] <- which(cellType %in% i)

    }

    C = lapply(group,function(x) Matrix::rowMeans(scC[,x])) 
    C = do.call(cbind.data.frame, C)

    return(as.matrix(C))

}



Scaling <- function(matrix, option, phenoDataC=NULL){

    #Avoid Error: Input matrix x contains at least one null or NA-filled row.
    matrix = matrix[rowSums(matrix) != 0,]
    #Avoid error if all elements within a row are equal (e.g. all 0, or all a common value after log/sqrt/vst transformation)
    matrix = matrix[!apply(matrix, 1, function(x) var(x) == 0),]
    matrix = matrix[,colSums(matrix) != 0]
 
    if (option == "LogNormalize"){
        
       matrix = expm1(Seurat::LogNormalize(matrix, verbose = FALSE)) %>% as.matrix()

    } else if (option == "TMM"){# CPM counts coming from TMM-normalized library sizes; https://support.bioconductor.org/p/114798/

        if(!is.null(phenoDataC)){

            Celltype = as.character(phenoDataC$cellType[phenoDataC$cellID %in% colnames(matrix)])

        } else {

            Celltype = colnames(matrix)

        }

        matrix <- edgeR::DGEList(counts=matrix, group=Celltype)
        CtrlGenes <- grep("ERCC-",rownames(data))
  
        if(length(CtrlGenes) > 1){
            
            spikes <- data[CtrlGenes,]
            spikes <- edgeR::calcNormFactors(spikes, method = "TMM") 
            matrix$samples$norm.factors <- spikes$samples$norm.factors
        
        } else {
        
            matrix <- edgeR::calcNormFactors(matrix, method = "TMM")  
        
        }

        matrix <- edgeR::cpm(matrix)

    } else if (option == "TPM"){

        # MULTI-CORE VERSION: adapted from #devtools::source_url('https://raw.githubusercontent.com/dviraran/SingleR/master/R/HelperFunctions.R', sha1 = "df5560b4ebb28295349aadcbe5b0dc77d847fd9e")

        TPM <- function(counts,lengths=NULL){

            require(foreach)
            require(Matrix)

            if (is.null(lengths)) {
                data('gene_lengths')
            }

            A = intersect(rownames(counts),names(lengths))
            counts = as.matrix(counts[A,])
            lengths = lengths[A]
            rate = counts / lengths 

            num_cores <- min(4, parallel::detectCores()) #slurm will automatically restrict this to parameter --cpus-per-task
            print(paste("num_cores = ", num_cores, sep = ""))

            doMC::registerDoMC(cores = num_cores)
            jump = 100
            tpm <- foreach::foreach(elem = seq(1,ncol(counts), jump), .combine='cbind.data.frame') %dopar% {   
                    
                    lsizes <- colSums(rate[, elem:min(ncol(counts),(elem + jump - 1))])
                    lo <- 1e6*rate[, elem:min(ncol(counts),(elem + jump - 1))]
                    tpm <- lo %*% diag(1/lsizes)
                    tpm

                }

            colnames(tpm) <- colnames(counts)
            return(tpm)

        }

        if(! file.exists("human_lengths.rda")){
          download.file("https://github.com/dviraran/SingleR/blob/master/data/human_lengths.rda?raw=true","./human_lengths.rda")
        }
        load("./human_lengths.rda")

        # Doesn't work with Ensembl IDs:
        if(length(grep("ENSG000",rownames(matrix))) > 50){
            
            suppressMessages(library("AnnotationDbi"))
            suppressMessages(library("org.Hs.eg.db"))
            temp = mapIds(org.Hs.eg.db,keys = names(human_lengths), column = "ENSEMBL", keytype = "SYMBOL", multiVals = "first")
            names(human_lengths) = as.character(temp)

        }
        
        matrix = TPM(counts = matrix, lengths = human_lengths)
        rownames(matrix) = toupper(rownames(matrix))

    } else if (option == "TPM.murine"){

        # MULTI-CORE VERSION: adapted from #devtools::source_url('https://raw.githubusercontent.com/dviraran/SingleR/master/R/HelperFunctions.R', sha1 = "df5560b4ebb28295349aadcbe5b0dc77d847fd9e")

        TPM <- function(counts,lengths=NULL){

            require(foreach)
            require(Matrix)

            if (is.null(lengths)) {
                data('gene_lengths')
            }


            A = intersect(rownames(counts),names(lengths))
            counts = as.matrix(counts[A,])
            lengths = lengths[A]
            rate = counts / lengths 

            num_cores <- min(4, parallel::detectCores()) #slurm will automatically restrict this to parameter --cpus-per-task
            print(paste("num_cores = ", num_cores, sep = ""))

            doMC::registerDoMC(cores = num_cores)
            jump = 100
            tpm <- foreach::foreach(elem = seq(1,ncol(counts), jump), .combine='cbind.data.frame') %dopar% {   
                    
                    lsizes <- colSums(rate[, elem:min(ncol(counts),(elem + jump - 1))])
                    lo <- 1e6*rate[, elem:min(ncol(counts),(elem + jump - 1))]
                    tpm <- lo %*% diag(1/lsizes)
                    tpm

                }

            colnames(tpm) <- colnames(counts)
            return(tpm)

        }

        if(! file.exists("mouse_lengths.rda")){
          download.file("https://github.com/dviraran/SingleR/blob/master/data/mouse_lengths.rda?raw=true","./mouse_lengths.rda")
        }

        load("./mouse_lengths.rda")

        ## Append version with ENSMUSG000 names
        mouse_lengths2 <- mouse_lengths

        # Doesn't work with Ensembl IDs:
        suppressMessages(library("AnnotationDbi"))
        suppressMessages(library("org.Mm.eg.db"))
        temp = mapIds(org.Mm.eg.db,keys = names(mouse_lengths), column = "ENSEMBL", keytype = "SYMBOL", multiVals = "first")
        names(mouse_lengths) = as.character(temp)

        mouse_lengths = c(mouse_lengths,mouse_lengths2)
        matrix = TPM(counts = matrix, lengths = mouse_lengths)

    ####################################################################################
    ## scRNA-seq specific  

    } else if (option == "SCTransform"){

        matrix = as(matrix, "dgCMatrix")
        matrix = sctransform::vst(matrix, return_corrected_umi = TRUE, show_progress = FALSE)$umi_corrected

    } else if (option == "scran"){
        
        sce = SingleCellExperiment::SingleCellExperiment(assays = list(counts=as.matrix(matrix)))
        sce = scran::computeSumFactors(sce, clusters=NULL)
        sce = scater::logNormCounts(sce, log = FALSE)
        matrix = normcounts(sce)
        
    } else if (option == "scater"){  

        sce = SingleCellExperiment::SingleCellExperiment(assays = list(counts=as.matrix(matrix)))
        size_factors = scater::librarySizeFactors(sce)
        sce = scater::logNormCounts(sce, log = FALSE)
        matrix = normcounts(sce)

    }

    return(matrix)

}



Deconvolution <- function(T, C, method, phenoDataC, P = NULL, elem = NULL, STRING = NULL, marker_distrib, refProfiles.var){ 

    bulk_methods = c("CIBERSORT","nnls","FARDEEP","RLR")
    sc_methods = c("MuSiC","MuSiC_with_markers","DWLS","SQUID", "Bisque", "Bisque_with_markers")

    ########## Using marker information for bulk_methods

    if(method %in% bulk_methods){

        C = C[rownames(C) %in% unique(marker_distrib$gene), , drop = FALSE]
        T = T[rownames(T) %in% unique(marker_distrib$gene), , drop = FALSE]
        refProfiles.var = refProfiles.var[rownames(refProfiles.var) %in% unique(marker_distrib$gene), , drop = FALSE]

    } else { ### For scRNA-seq methods 

        if(length(grep("[N-n]ame",colnames(phenoDataC))) > 0){
            sample_column = grep("[N-n]ame",colnames(phenoDataC))
        } else {
            sample_column = grep("[S-s]ample|[S-s]ubject",colnames(phenoDataC))
        }

        colnames(phenoDataC)[sample_column] = "SubjectName"
        rownames(phenoDataC) = phenoDataC$cellID
        # establish same order in (sc)C and phenoDataC:
        phenoDataC <- phenoDataC[match(colnames(C),phenoDataC$cellID),]

        if(method %in% c("MuSiC","MuSiC_with_markers", "Bisque", "Bisque_with_markers")){

            require(xbioc)
            C.eset <- Biobase::ExpressionSet(assayData = as.matrix(C), phenoData = Biobase::AnnotatedDataFrame(phenoDataC))
            T.eset <- Biobase::ExpressionSet(assayData = as.matrix(T))
        
        }
        
    }

    ##########    MATRIX DIMENSION APPROPRIATENESS    ##########

    keep = intersect(rownames(C),rownames(T)) 
    C = C[keep, , drop = FALSE]
    T = T[keep, , drop = FALSE]

    ###################################

    if(method == "CIBERSORT"){ 

        #source("./CIBERSORT.R")
        RESULTS = CIBERSORT(sig_matrix = C, mixture_file = T, QN = FALSE) 
        RESULTS = t(RESULTS[,1:(ncol(RESULTS)-3),drop=FALSE]) 

    } else if (method == "nnls"){

        require(nnls)
        RESULTS = do.call(cbind.data.frame,lapply(apply(T,2,function(x) nnls::nnls(as.matrix(C),x)), function(y) y$x))
        RESULTS = apply(RESULTS,2,function(x) x/sum(x)) #explicit STO constraint
        rownames(RESULTS) <- colnames(C)

    } else if (method == "FARDEEP"){

        require(FARDEEP)
        RESULTS = t(FARDEEP::fardeep(C, T, nn = TRUE, intercept = TRUE, permn = 10, QN = FALSE)$abs.beta)
        RESULTS = apply(RESULTS,2,function(x) x/sum(x)) #explicit STO constraint
        # These 2 lines are similar as retrieving $relative.beta instead of abs.beta + re-scaling

    } else if (method == "RLR"){ #RLR = robust linear regression

        require(MASS)
        RESULTS = do.call(cbind.data.frame,lapply(apply(T,2,function(x) MASS::rlm(x ~ as.matrix(C), maxit=100)), function(y) y$coefficients[-1]))
        RESULTS = apply(RESULTS,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
        RESULTS = apply(RESULTS,2,function(x) x/sum(x)) #explicit STO constraint
        rownames(RESULTS) <- unlist(lapply(strsplit(rownames(RESULTS),")"),function(x) x[2]))

    ###################################
    ###################################
    
    } else if (method == "MuSiC"){

        require(MuSiC)
        RESULTS = t(MuSiC::music_prop(bulk.eset = T.eset, sc.eset = C.eset, clusters = 'cellType',
                                            markers = NULL, normalize = FALSE, samples = 'SubjectName', 
                                            verbose = F)$Est.prop.weighted)

    } else if (method == "MuSiC_with_markers"){

        require(MuSiC)
        RESULTS = t(MuSiC::music_prop(bulk.eset = T.eset, sc.eset = C.eset, clusters = 'cellType',
                                            markers = unique(marker_distrib$gene), normalize = FALSE, samples = 'SubjectName', 
                                           verbose = F)$Est.prop.weighted)
    
    } else if (method == "Bisque"){#By default, Bisque uses all genes for decomposition. However, you may supply a list of genes (such as marker genes) to be used with the markers parameter

        require(BisqueRNA)
        RESULTS <- BisqueRNA::ReferenceBasedDecomposition(T.eset, C.eset, markers = NULL, use.overlap = FALSE)$bulk.props 

    } else if (method == "Bisque_with_markers"){#By default, Bisque uses all genes for decomposition. However, you may supply a list of genes (such as marker genes) to be used with the markers parameter

        require(BisqueRNA)
        RESULTS <- BisqueRNA::ReferenceBasedDecomposition(T.eset, C.eset, markers = marker_distrib$gene, use.overlap = FALSE)$bulk.props 

    } else if (method == "DWLS"){
        
        require(DWLS)
        path = paste(getwd(),"/results_",STRING,sep="")

        if(! dir.exists(path)){ dir.create(path) } #to avoid repeating marker_selection step when removing cell types; Sig.RData automatically created

        if(!file.exists(paste(path,"Sig.RData",sep="/"))){

            Signature <- DWLS::buildSignatureMatrixMAST(scdata = C, id = as.character(phenoDataC$cellType), path = path, diff.cutoff = 0.5, pval.cutoff = 0.01)

        } else {#re-load signature and remove CT column + its correspondent markers

            load(paste(path,"Sig.RData",sep="/"))
            Signature <- Sig
            
        }

        Signature = as.matrix(Signature)
        
        RESULTS <- apply(T, 2, function(x){
            b = setNames(x, rownames(T))
            tr <- DWLS::trimData(Signature, b)
            RES <- t(DWLS::solveDampenedWLS(tr$sig, tr$bulk))
        })

        rownames(RESULTS) <- as.character(unique(phenoDataC$cellType))
        RESULTS = apply(RESULTS,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
        RESULTS = apply(RESULTS,2,function(x) x/sum(x)) #explicit STO constraint
        print(head(RESULTS))

    } else if (method == "SQUID"){ 
            
        # Transforming T :
        Z = bulkC.fromSC(scC = C, phenoDataC = phenoDataC)
        X = T
        P <- as.matrix(P[base::colnames(Z),,drop = FALSE])
        Y = Z %*% P

        X.new = transformation2(X = X, Y = Y, leave.one.out = FALSE)
        X.new[!is.finite(X.new)] <- 0 
        
        # take common genes
        Genes <- intersect(rownames(Z),rownames(X.new))

        X.new <- as.matrix.Vector(X.new[Genes,])
        Z <- as.matrix.Vector(Z[Genes,])

        RESULTS <- apply(X.new, 2, function(x){
            RES <- t(DWLS::solveDampenedWLS(S = Z, B = x))
        })

        rownames(RESULTS) <- colnames(Z)
        RESULTS = apply(RESULTS,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
        RESULTS = apply(RESULTS,2,function(x) x/sum(x)) #explicit STO constraint

    }

    RESULTS = RESULTS[gtools::mixedsort(rownames(RESULTS)), , drop = FALSE] 
    RESULTS = data.table::melt(RESULTS)
    colnames(RESULTS) <- c("cell_type","tissue","observed_fraction")
    RESULTS$cell_type = as.character(RESULTS$cell_type)
    RESULTS$tissue = as.character(RESULTS$tissue)

    if(!is.null(P)){

        P = P[gtools::mixedsort(rownames(P)),,drop = FALSE] %>% data.frame(., check.names = FALSE)
        P$cell_type = rownames(P)
        P = data.table::melt(P, id.vars="cell_type")
        colnames(P) <-c("cell_type","tissue","expected_fraction")
        P$cell_type = as.character(P$cell_type)
        P$tissue = as.character(P$tissue)

        RESULTS = merge(RESULTS, P, by = c("cell_type", "tissue"), all = TRUE)
        RESULTS[is.na(RESULTS)] <- 0
        RESULTS$expected_fraction <- round(RESULTS$expected_fraction, 3)
        RESULTS$observed_fraction <- round(RESULTS$observed_fraction, 3)

    }

    return(RESULTS) 

}
