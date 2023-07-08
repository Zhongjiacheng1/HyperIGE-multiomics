# The following function was modified by ArchR author 
# https://github.com/GreenleafLab/ArchR/commit/f98354ce07f72c86213982e5fbcd9a5bd70d743b

.MarkersSC_fixed <- function(
  ArchRProj = NULL,
  groupBy = "Clusters",
  useGroups = NULL,
  bgdGroups = NULL,
  normBy = NULL,
  maxCells = 500,
  scaleTo = 10^4,
  bufferRatio = 0.8,
  bias = NULL,
  k = 100,
  threads = 1,
  binarize = FALSE,
  useSeqnames = NULL,
  testMethod = "wilcoxon",
  useMatrix = "GeneScoreMatrix",
  markerParams = list(),
  verbose = TRUE,
  logFile = NULL
){
  tstart <- Sys.time()
  
  #####################################################
  # Feature Info
  #####################################################
  ArrowFiles <- getArrowFiles(ArchRProj)
  featureDF <- .getFeatureDF(head(ArrowFiles, 2), useMatrix)
  matrixClass <- as.character(h5read(getArrowFiles(ArchRProj)[1], paste0(useMatrix, "/Info/Class")))
  .logThis(range(as.vector(table(paste0(featureDF$seqnames)))), "FeaturesPerSeqnames", logFile = logFile)
  isDeviations <- FALSE
  if(all(unique(paste0(featureDF$seqnames)) %in% c("z", "dev"))){
    isDeviations <- TRUE
  }
  .logThis(featureDF, "FeatureDF", logFile=logFile)
  .logMessage(paste0("MatrixClass = ", matrixClass), logFile=logFile)
  seqnames <- unique(as.vector(featureDF$seqnames))
  useSeqnames <- useSeqnames[useSeqnames %in% seqnames]
  if(length(useSeqnames)==0){
    useSeqnames <- NULL
  }
  if(!is.null(useSeqnames)){
    if(matrixClass == "Sparse.Assays.Matrix"){
      if(length(useSeqnames) == 1){
        featureDF <- featureDF[BiocGenerics::which(featureDF$seqnames %bcin% useSeqnames),]
      }else{
        .logMessage("When accessing features from a matrix of class Sparse.Assays.Matrix it requires 1 seqname!\n",
                    "Continuing with first seqname '", seqnames[1], "'!\n",
                    "If confused, try getSeqnames(ArchRProj, '", useMatrix,"'') to list out available seqnames for input!", verbose = verbose, logFile = logFile)
        useSeqnames <- seqnames[1]
        featureDF <- featureDF[BiocGenerics::which(featureDF$seqnames %bcin% useSeqnames),]
      }
    }else{
      featureDF <- featureDF[BiocGenerics::which(featureDF$seqnames %bcin% useSeqnames),]
    }
  }else{
    if(matrixClass == "Sparse.Assays.Matrix"){
      .logMessage("When accessing features from a matrix of class Sparse.Assays.Matrix it requires 1 seqname!\n",
                  "Continuing with first seqname '", seqnames[1], "'!\n",
                  "If confused, try getSeqnames(ArchRProj, '", useMatrix,"'') to list out available seqnames for input!", verbose = verbose, logFile = logFile)
      useSeqnames <- seqnames[1]
      featureDF <- featureDF[BiocGenerics::which(featureDF$seqnames %bcin% useSeqnames),]
    }
  }
  if(!(nrow(featureDF) > 1)){
    .logStop("Less than 1 feature is remaining in featureDF please check input!", logFile = logFile)
  }
  #####################################################
  # Match Bias Groups
  #####################################################
  .logDiffTime("Matching Known Biases", t1 = tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
  groups <- getCellColData(ArchRProj, groupBy, drop = TRUE)
  colDat <- getCellColData(ArchRProj)
  matchObj <- .matchBiasCellGroups(
    input = colDat, 
    groups = groups,
    useGroups = useGroups,
    bgdGroups = bgdGroups,
    bias = bias,
    k = k,
    n = maxCells,
    bufferRatio = bufferRatio,
    logFile = logFile
  )
  
#####################################################
# Pairwise Test Per Seqnames
#####################################################
#ColSums
mColSums <- tryCatch({
  suppressMessages(tmpColSum <- .getColSums(ArrowFiles, seqnames = featureDF$seqnames@values, useMatrix = useMatrix, threads = threads))
  tmpColSum[ArchRProj$cellNames]
}, error = function(x){
  rep(1, nCells(ArchRProj))
})

if(all(mColSums==1) & is.null(normBy)){
  normBy <- "none"
}

if(is.null(normBy)){
  if(tolower(testMethod) == "binomial"){
    normFactors <- NULL
  }else{
    if(tolower(useMatrix) %in% c("tilematrix", "peakmatrix")){
      normBy <- "ReadsInTSS"
      normFactors <- getCellColData(ArchRProj, normBy, drop=FALSE)
      normFactors[,1] <- median(normFactors[,1]) / normFactors[,1]
    }else{
      normFactors <- scaleTo / mColSums
      normFactors <- DataFrame(normFactors)
    }
  }
}else{
  if(tolower(normBy) == "none"){
    normFactors <- NULL
  }else{
    normFactors <- scaleTo / mColSums
    normFactors <- DataFrame(normFactors)
  }
}
if(!is.null(normFactors)){
  normFactors[,1] <- normFactors[,1] * (scaleTo / median(normFactors[names(mColSums), 1] * mColSums))
}

diffList <- .safelapply(seq_along(matchObj[[1]]), function(x){
  .logDiffTime(sprintf("Computing Pairwise Tests (%s of %s)", x, length(matchObj[[1]])), tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
  .testMarkerSC(
    ArrowFiles = ArrowFiles,
    matchObj = matchObj, 
    group = names(matchObj[[1]])[x], 
    testMethod = testMethod, 
    threads = 1, 
    useMatrix = useMatrix,
    featureDF = featureDF,
    normFactors = normFactors,
    binarize = binarize,
    logFile = logFile
  )
}, threads = threads)
.logDiffTime("Completed Pairwise Tests", tstart, addHeader = TRUE, verbose = verbose, logFile = logFile)
#####################################################
# Summarize Output
#####################################################
if(tolower(testMethod) == "wilcoxon"){
  pse <- SummarizedExperiment::SummarizedExperiment(
    assays = 
      SimpleList(
        Log2FC = lapply(seq_along(diffList), function(x) data.frame(x = diffList[[x]]$log2FC)) %>% Reduce("cbind",.),
        Mean = lapply(seq_along(diffList), function(x) data.frame(x = diffList[[x]]$mean1)) %>% Reduce("cbind",.),
        FDR = lapply(seq_along(diffList), function(x) data.frame(x = diffList[[x]]$fdr)) %>% Reduce("cbind",.),
        Pval = lapply(seq_along(diffList), function(x) data.frame(x = diffList[[x]]$pval)) %>% Reduce("cbind",.),
        MeanDiff = lapply(seq_along(diffList), function(x) data.frame(x = diffList[[x]]$mean1 - diffList[[x]]$mean2)) %>% Reduce("cbind",.),
        AUC = lapply(seq_along(diffList), function(x) data.frame(x = diffList[[x]]$auc)) %>% Reduce("cbind",.),
        MeanBGD = lapply(seq_along(diffList), function(x) data.frame(x = diffList[[x]]$mean2)) %>% Reduce("cbind",.)
      ),
    rowData = featureDF
  )
}else if(tolower(testMethod) == "ttest"){
  pse <- SummarizedExperiment::SummarizedExperiment(
    assays = 
      SimpleList(
        Log2FC = lapply(seq_along(diffList), function(x) data.frame(x = diffList[[x]]$log2FC)) %>% Reduce("cbind",.),
        Mean = lapply(seq_along(diffList), function(x) data.frame(x = diffList[[x]]$mean1)) %>% Reduce("cbind",.),
        Variance = lapply(seq_along(diffList), function(x) data.frame(x = diffList[[x]]$var1)) %>% Reduce("cbind",.),
        FDR = lapply(seq_along(diffList), function(x) data.frame(x = diffList[[x]]$fdr)) %>% Reduce("cbind",.),
        Pval = lapply(seq_along(diffList), function(x) data.frame(x = diffList[[x]]$pval)) %>% Reduce("cbind",.),
        MeanDiff = lapply(seq_along(diffList), function(x) data.frame(x = diffList[[x]]$mean1 - diffList[[x]]$mean2)) %>% Reduce("cbind",.),
        MeanBGD = lapply(seq_along(diffList), function(x) data.frame(x = diffList[[x]]$mean2)) %>% Reduce("cbind",.),
        VarianceBGD = lapply(seq_along(diffList), function(x) data.frame(x = diffList[[x]]$var2)) %>% Reduce("cbind",.)
      ),
    rowData = featureDF
  )
}else if(tolower(testMethod) == "binomial"){
  pse <- SummarizedExperiment::SummarizedExperiment(
    assays = 
      SimpleList(
        Log2FC = lapply(seq_along(diffList), function(x) data.frame(x = diffList[[x]]$log2FC)) %>% Reduce("cbind",.),
        Mean = lapply(seq_along(diffList), function(x) data.frame(x = diffList[[x]]$mean1)) %>% Reduce("cbind",.),
        FDR = lapply(seq_along(diffList), function(x) data.frame(x = diffList[[x]]$fdr)) %>% Reduce("cbind",.),
        Pval = lapply(seq_along(diffList), function(x) data.frame(x = diffList[[x]]$pval)) %>% Reduce("cbind",.),
        MeanDiff = lapply(seq_along(diffList), function(x) data.frame(x = diffList[[x]]$mean1 - diffList[[x]]$mean2)) %>% Reduce("cbind",.),
        MeanBGD = lapply(seq_along(diffList), function(x) data.frame(x = diffList[[x]]$mean2)) %>% Reduce("cbind",.)
      ),
    rowData = featureDF
  )
}else{
  stop("Error Unrecognized Method!")
}
colnames(pse) <- names(matchObj[[1]])
metadata(pse)$MatchInfo <- matchObj
if(isDeviations){
  assays(pse)[["Log2FC"]] <- NULL #This measure does not make sense with deviations matrices better to just remove
}
return(pse)
}

# Override the original function in ArchR 1.0.1

environment(.MarkersSC_fixed) <- asNamespace('ArchR')
assignInNamespace(".MarkersSC", .MarkersSC_fixed, ns = 'ArchR')
