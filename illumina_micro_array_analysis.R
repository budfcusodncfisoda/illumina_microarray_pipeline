#libraries required to run pipeline
library("affy")
# lumi is a package specific for illumina bead arrays. contains many functions for formatting
# GenomeStudio output as well as normalization
library("lumi")
# used for differential gene expression (building linear models, topTable)
library('limma')
library("gplots")
# used for comparing unsupervised differences between samples
library("entropy")
library("biomaRt") 
library("annotate")
library('GOstats')
library("org.Mm.eg.db")
library("lumiMouseIDMapping")
library("topGO")

# normFilt Function
# uploads, normalizes, filters gene expression data
# takes genome studio output, makes an expression set based upon treatment vector

#normFilt arguments:
#treatmentVector  = simple vector containing string variables that describe the specific
#                   experimental conditions. (Must match conditions in target file)
#                   example: treatmentVector <- c("Epi_10", "Aza_10")
#exprData       = genome studio output containing: probeIDs, signal intensities, detection p-values,
#                   and annotations for the 37 samples of the micro-array experiment 
#controlFile      = used for general QC analysis, control probes contain beads matched to common 
#                   house keeping genes, negative controls, etc. used in normalization process
#                   found at the bottom of genomeStudio output. "Control Probe Profile"
#targetFile     = tab file with sample IDs (from Illumina genome studio output) and treatment
#                   description (i.e. "Cis_10" for Cisplatin 10wk treatment). Make sure sample
#                   sample IDs in the Target File are in the same order as sample IDs in the ExprData
normFilt <- function(treatmentVector, exprData, controlFile, targetFile, path){
  
  # sets working directory to specified "path." contains necessary files to run pipeline.
  setwd(path)
  
  # lumiR.batch: reads in illumina GenomeStudio output.
  # expreData contains a tab file with 8 header lines, and control data.
  # this specific upload function does not perform log2 transform
  # the detection threshold is .05 (p-value)
  mouse.lumi.data37<- suppressMessages( lumiR.batch(exprData, detectionTh=.05,  inputAnnotation=TRUE,
                                        annotationColumn=c('ACCESSION', 'SYMBOL', 'PROBE_START','CHROMOSOME',
                                        'PROBE_CHR_ORIENTATION', 'PROBE_COORDINATES', 'DEFINITION'), verbose=FALSE))
  
  # column 1: Sample-ID; Column 2: Treatment (Norm, Veh, Epi, etc.)
  Target37 <- read.delim(targetFile, header=TRUE)
  
  # create new data frame with just treatments specified by the (input) treatment vector
 
  extrct <- Target37[which(as.character(Target37$Treatment) %in% treatmentVector),]
  
  # create file name variable for new targetFile based upon treatment vector 
  filenameT <- paste("target",paste(treatmentVector, collapse=""),".txt", sep="")
  
  # create new target file from treatment vector and "master" target file 
  write.table(extrct, filenameT, row.names=FALSE, quote=FALSE, sep="\t")
  
  # extract new expression set object using specified sample IDs that match treatments 
  # specified in the treatment vector and matched using the target file
  exprSet <- suppressWarnings(mouse.lumi.data37[,mouse.lumi.data37@phenoData@data$sampleID %in% extrct$Sample_ID])
  
  # When the new batch object is created, the control data is lost. luckily the "addControlData2lumi"
  # function is able to add the control file onto the new extracted expression set, despite them having
  # a different number of Sample IDs
  # the control data must be in a different file than the GenomeStudio output. (copy paste it into new speadsheet)
  print("fixing control data")
  exprSet <- addControlData2lumi(controlFile, exprSet)
  
  # create a density plot of the samples being compared.
  # this illustrates the difference in distibution of expression values
  #windows()
  #plot(exprSet, what='density')
  #windows()
  #boxplot(exprSet,xlab="", main="Gene expression levels (RSN)", margins=c(10,25),cex.names=1.5, cex.lab=1.3, cex.main=2)
  # Variance stablizing transform (VST), takes advantage of the technical replicates available on an Illumina microarray.
  lumi.T37 <- lumiT(exprSet, verbose=TRUE)
  
  # robust spline normalization (RSN) algorithm combines the features of quantile and loess normalization. It is designed
  # to normalize the variance-stabilized data.
  exprNorm <- lumiN(lumi.T37, method="rsn", verbose=FALSE)
  
  # display change in density distribution after normalization
  #windows()
  #plot(exprNorm, what='density')
  
  #windows()
  #boxplot(exprNorm,xlab="", main="Gene expression levels (RSN)", margins=c(10,25),cex.names=1.5, cex.lab=1.3, cex.main=2)
 
  #windows()
  # outlier detection is based on the distance from the sample to the center (average of all samples), euclidean
  #plot(exprNorm, what = c("outlier"), main="Detect Outliers (Euclidean)")
  
  # remove uneeded variables
  rm(extrct, exprSet, mouse.lumi.data37, targetFile, treatmentVector, lumi.T37, controlFile)
  
  # get number of probes in expression set
  unfilt.len <- length(rownames(exprNorm))
  cat("\n\ntotal number of probes in array: ", unfilt.len, "\n")
  # filters probe-sets with at least 1 probe that is significantly detected (i.e. p-val < .05) 
  filt.Norm.ExprSet <-exprNorm[which(rowSums(exprNorm@assayData$detection < 0.05) >= 1),]
  
  # get length of probesets after filtering
  filt.len <- length(rownames(filt.Norm.ExprSet))
  cat("number of filtered probes with a detection p-val < .05: ", filt.len,"\n")
  # display change in number of probesets, before and after filtering
  cat("percentage of probes remaining: ", (round((filt.len/unfilt.len), 4)*100),"\n")
  
  # remove uneeded variables
  rm(unfilt.len, filt.len, exprNorm)
  
  # return normalized and filtered expression set 
  return(c(filt.Norm.ExprSet, filenameT))
}

# diffExpr()
# performs differenetial expression calculations on normalized filtered data
# filt.Norm.ExprSet = Expression Set without non-detected genes, with treatment vector samples
# and normalized data. (example: vehNorm.normFilt[[1]])
# filenameT         = target file name
# treatmentVector   = types of treatments being compared
# :RETURNS: entire topTable (f.top) and eBayes output

diffExpr <- function(filt.Norm.ExprSet, filenameT, treatmentVector){
  # get target file containing cample ids and treatment types
  target <- readTargets(file=filenameT)
  # vector containing different treatments (i.e. "Norm" "Vehicle")
  lev <- levels(factor(target$Treatment))
  # each column corresponds to a treatmentment type
  # the values of the matrix are 0 or 1, based upon whether the sample corresponds to the treatment 
  design<-model.matrix(~0+target$Treatment)
  #rename column names
  colnames(design)<-lev
  # if both sample sizes are greater than three then the arrayWeights function will assign differing
  # weights to each slide within the arrays, based upon quality control
  if(sum(design[,1]) > 3 && sum(design[,2]) > 3){
     arrayw <- arrayWeights(filt.Norm.ExprSet, design)
     imageplot(arrayw)
     # models systematic parts of the data
     # each row corresponds to an array in the experiment
     # each column is a coefficient (Epi, Cis, etc.)
     fit<-lmFit(filt.Norm.ExprSet, design, weights=arrayw)
  }
  # if one of the sample sizes is less than three the arrays are not weighted
  else if(sum(design[,1]) <= 3 || sum(design[,2]) <= 3){
    # build linear model without weighted arrays
    fit<-lmFit(filt.Norm.ExprSet, design)
  }
  # Else if a sample has 0 arrays in it
  else{
    print("ERROR: check sample size of treatment vector groups")
  }
  # vector to hold strings comparisons ("Vehicle-Norm") comparisons
  # used in makeContrasts function (contr = contrasts)
  contr.str<-c()
  # get length of the samples being compared
  len<-length(lev)
  # for-loop for taking levels (or targets) and combining them to string and putting them into vector 
  for(i in 1:(len-1)){
    contr.str<-c(contr.str, paste(lev[(i+1):len], lev[i], sep="-"))
    }
  # expresses contrasts between a set of parameters as a numeric matrix
  # contasts = c("Vehicle-Norm", "Vehicle-Epi_10", etc.)
  contr.mat<-makeContrasts(contrasts=contr.str, levels=lev)
  # given a linear model fit to microarray data, compute estimated coefficients and standard errors
  # given set of contrasts
  fit2<-contrasts.fit(fit, contr.mat)
  # compute moderated t-statistics, moderated F-statistic, log-odds of differential expression by
  # empirical Bayes moderation of standard errors towards common value
  fit2<-eBayes(fit2)
  # extract table of top ranked genes
  # important values include logFC, adj.P.val
  f.top<-topTable(fit2, number=nrow(filt.Norm.ExprSet))
  # return topTable (f.top) and eBayes output...
  return(list(f.top, fit2))
}

# principle component analysis:
# performed as the first step after low-level data processing to obtain a "big picture" of the
# data.  capture the variance in a dataset in terms of principal components. Can identify the sample 
# outliers in the dataset by reducing the dimensionality of the data
pCompAnalysis <- function(filt.Norm.ExprSet, targetFile){
  # transpose expression matrix, "scale." is explicitly turned on
  # so that all the genes contribute equally to the analysis regardless of
  # the magnitude of the change.
  target <- read.delim(targetFile, header=TRUE)
  pca.res <- prcomp(t(exprs(filt.Norm.ExprSet)), scale.=TRUE, retx=TRUE)
  #plot(pca.res, las=1)
  # percentage of variation from each principle componet
  pc.var<-pca.res$sdev^2
  pc.per<-round(pc.var/sum(pc.var)*100, 1)
  # plot the samples on to the first two PCs that carry the most variance
  plot.pch <- (1:length(levels(factor(target$Treatment))))[as.numeric(factor(target$Treatment))]
  #windows()
  plot(pca.res$x, col=1, pch=plot.pch, las=1, cex=2, xlab=paste("PC1 (",pc.per[1], "%)", sep=""),ylab=paste("PC2 (", pc.per [2], "%)", sep=""))
  legend(60, 125, levels(factor(target$Treatment)), pch=1:length(levels(factor(target$Treatment))),pt.cex=1.5)
  rm(pc.var,pca.res, pc.per, plot.pch)
}

# retProbeID.cv()
# this function has 5 options (i.e. answers) c("entropy", "CV","StDev", "ryba", "meanSub")
# it will return the top 5% of significant probes based upon one of the 5 options 
# mouse.exprSet = expression set with filtered and normalized data (i.e. vehNorm.normFilt[[1]], 
# filt.Norm.ExprSet)
# targetFileName = filenameT -- second output from normFilt 
# answer = type of distance formula used to compare different experiments
retProbeID.cv <- function(mouse.exprSet, targetFileName, answer) {
  cat("calculating differential expression using:", answer, "\n")
  # calculate the coefficient of variance
  co.var <- function(x) ( 100*sd(x)/mean(x) )
  # function for standard deviation
  SD <- function(x) sd(x) 
  # ::: this is used for user input :::
  # ::: the code has been modifide b/c the wrapper function will pass the answer parameter :::
  #options <- c("entropy", "CV","StDev", "ryba", "meanSub")
  #answer <- readline("return top probe IDs for (entropy, CV, StDev, ryba, meanSub)? ")
  #while(!(answer %in% options)){
  #  answer <- readline("return top probe IDs for (entropy, CV, StDev, ryba, meanSub)? ") 
  #}
  # read target file
  target <- readTargets(file=targetFileName)
  # extract pure expression data
  expr <- exprs(mouse.exprSet)
  # re-name colnames to Treatment names
  colnames(expr) <- target$Treatment
  # create a blank matrix with a column for each sample
  mat <- matrix(, nrow= nrow(expr), ncol = sum(!duplicated(colnames(expr))))
  # create extra blank matrix with a column for each sample
  mat.WS <- matrix(, nrow= nrow(expr), ncol = sum(!duplicated(colnames(expr))))
  # treatments meant for column names
  treatments <- colnames(expr)[!duplicated(colnames(expr))]
  # assign column names
  colnames(mat) <- treatments
  # row names same as expression set
  rownames(mat) <- rownames(expr)
  # same for mat.WS
  colnames(mat.WS) <- treatments
  rownames(mat.WS) <- rownames(expr)
  # the for-loop iterates over the columns of the new matrixes (or the number of treatments 
  # (which is   the same))
  for(i in 1:length(treatments)){
    # get mean of expression values for each treatment then store those in "mat" matrix
    # each sample gets its expression value mean for each probe stored in "mat" 
    mat[,i] <- apply(expr[, colnames(expr) %in% treatments[i]],1, mean)
    # does the same thing, except calculates the standard deviation and stores it inside 
    # mat.WS
    mat.WS[,i] <- apply(expr[, colnames(expr) %in% treatments[i]],1, sd)
  }
  # answer passed as a parameter
  if(answer == "entropy"){
    # calculate the entropy between each sample
    CV <- apply(mat, 1, entropy)
    # sort probe IDs by which have the highest entropy (high to low)
    CV <- sort(CV, decreasing=FALSE)
  }
  # calculate coefficient of variation between
  else if(answer == "CV"){
    CV <- apply(mat, 1, co.var)
    CV <- sort(CV, decreasing=TRUE)
  }
  # calculate the standard deviation between the two means...
  # this doesn't work very well...
  else if(answer == "StDev"){ 
    CV <- apply(mat, 1, sd)
    CV <- sort(CV, decreasing=TRUE)
  }
  # subtract the means by each other
  else if(answer == "meanSub"){
    CV <- apply(mat, 1, '-')
    # abs used because subtraction creates negatives
    CV <- sort(abs(CV), decreasing=TRUE)
  }
  # Ryba's anova like function
  else if(answer == "ryba"){
    mat.WS <- rowMeans(mat.WS)
    mat.SD <- sd(mat)
    mat.BSWS <- mat.SD / mat.WS
    CV <- sort(mat.BSWS, decreasing=FALSE)
   }
  else{
    print("Error: idk what method of top probeIDs...")
  }
  
  if(length(treatments) == 2){
    upDown.mat <- matrix(, nrow= nrow(expr), ncol = sum(!duplicated(colnames(expr))))
    for(i in 1:nrow(upDown.mat)){
        if(mat[i,1] > mat[i,2] ){
          upDown.mat[i,1] <- '+' 
          upDown.mat[i,2] <- '-'
        }
        else if(mat[i,1] < mat[i,2] ){
          upDown.mat[i,1] <- '-' 
          upDown.mat[i,2] <- '+'
        }
        else if(mat[i,1] == mat[i,2]){ 
          upDown.mat[i,1] <- '=' 
          upDown.mat[i,2] <- '='
        }
        else{ 
          print("unforeseen value")
          upDown.mat[i,1] <- 'NULL' 
          upDown.mat[i,2] <- 'NULL'
        }
    }  
    print("making 'up-down' regulation table")
    upDown.frame <-as.data.frame(upDown.mat)
    rownames(upDown.frame)<-rownames(mouse.exprSet)
    colnames(upDown.frame)<-treatments
    
    topCVnum <- .05*nrow(expr)
    # get expression data for top 10% of probes
    topCVprobID <- head(CV, topCVnum)
    print("extract top 10% of differentially expressed genes")
    cat("number of genes remaining", length(topCVprobID), "\n\n")
    upDown.frame <- upDown.frame[rownames(upDown.frame) %in% names(topCVprobID),]
    return(list(topCVprobID, upDown.frame))
  }
  
  # get 5% of number of rows
  topCVnum <- .05*nrow(expr)
  # get expression data for top 5% of probes
  topCVprobID <- head(CV, topCVnum)
  print("extract top 10% of differentially expressed genes")
  cat("number of genes remaining", length(topCVprobID), "\n\n")
  # return top 10%
  return(topCVprobID)
}

# analysis() Function
# this function runs the normalization and filtering function 
# along with the differential expression function
#
# treatmentVector = simple vector containing string variables that describe the specific
#                   experimental conditions. This gives the user the ability to specify
#                   which treatments they would like to compare
# exprData      = name of text file containing genome studio output
# controlFile     = name of text file containing control data. found at the bottom of 
#                   genomeStudio output, it should be labelled "Control Probe Profile." 
#                   The control file must be created by cutting and pasting all the data
#                   into a new tab document (I used excel). 
# targetFile    = tab file with sample IDs (from Illumina genome studio output) and treatment
#                   description (i.e. "Cis_10" for Cisplatin 10wk treatment). Make sure sample
#                   sample IDs in the Target File are in the same order as sample IDs in the ExprData
# path            = directory path to where the data is being stored
# lfc             = specified log fold change necessary for cut-off (positive)
# adj.pval        = specified adjusted pvalue for cut-off
# answer          = types of differential expression equations (i.e "entropy", "CV","StDev", "ryba", "meanSub") distance metrics
#                   you use to compare differences in expression levels between two groups of genes
# ::RETURNS::
# [[1]] = ExprSet supervised, 
# [[2]] = ExprSet unsupervised,") 
# [[3]] = top table of signifcantly diff. expr. genes,
# [[4]] = normalized and filtered expression set,
# [[5]] = target file
# [[6]] = data from linear model 
# [[7]] = TopTable from all norm filtered data
# [[8]] = data frame with + or minus for up down regulation

analysis <- function(treatmentVector,exprData, controlFile, targetFile, path, answer, lfc, adj.pval){
  # uploads, normalizes, filters gene expression data
  # :RETURNS: [[1]] normalized and filtered expression set, [[2]] name of a
  # target file for specified treatment vector
  norm.filt <- normFilt(treatmentVector,exprData, controlFile, targetFile, path)
  #performs quality control, plots principle componets
  pCompAnalysis(norm.filt[[1]], norm.filt[[2]])
  if(length(treatmentVector) == 2){ 
    cat("comparing gene expression between 2 treatments:", treatmentVector[1], ":", treatmentVector[2],"\n")
  # retProbeID.cv() performs unsupervised differential expression
  # this function has 5 options (i.e. answers) c("entropy", "CV","StDev", "ryba", "meanSub")
  # alternative to lfc and Adj.P.Val
  # top 5% of differentially expressed genes between two groups    
    UnSuperv.ID <- retProbeID.cv(norm.filt[[1]], norm.filt[[2]], answer)
    upDown.frame <- UnSuperv.ID[[2]]
    UnSuperv.ID <- UnSuperv.ID[[1]]
    expr.Unsup <- norm.filt[[1]][rownames(norm.filt[[1]]) %in% names(UnSuperv.ID),]
  }
  # creates expression set with top 5% of differentially exressed genes based upon 'answer' equation selected
  # sum of diff expr genes
  # supervised data
  # performs differenetial expression calculations on normalized filtered data
  # :RETURNS: entire topTable (f.top) and eBayes output
  top <- diffExpr(norm.filt[[1]], norm.filt[[2]])
  # if the user does not specify a log fold change parameter 
  if(missing(lfc)){
     print("without lfc")
     cat("Adjusted P Val: ", adj.pval, "\n" )
     top.lfc <- top[[1]][top[[1]]$adj.P.Val<adj.pval,]
     # sort table based on adj.P.Val
     top.lfc <- top.lfc[order(top.lfc$adj.P.Val),]
     # create expression set with pvalues that meet the reqirement
     expr.Sup <- norm.filt[[1]][rownames(norm.filt[[1]]) %in% rownames(top.lfc),]
     cat("number of differentiall expressed genes: ",nrow(expr.Sup),"\n")
  }
  else if(missing(adj.pval)){
    print("without adj.p.val")
    cat("log fold change greater than: ", lfc,"\n")
    neg.lfc <- lfc-(2*lfc)
    cat("log fold change less than: ", neg.lfc, "\n")
    # The binary logarithm of n is the power to which the number 2 must
    # be raised to obtain the value n.
    top.lfc <- top[[1]][top[[1]]$logFC>lfc | top[[1]]$logFC<neg.lfc,]
    top.lfc <- top.lfc[order(-top.lfc$logFC),]
    expr.Sup <- norm.filt[[1]][rownames(norm.filt[[1]]) %in% rownames(top.lfc),]
    cat("differentially expressed genes: ",nrow(expr.Sup), "\n")
  }
  else if(!missing(adj.pval) & !missing(lfc)){
    cat("Adjusted P Val less than: ", adj.pval, "\n" )
    cat("log fold change greater than: ", lfc, "\n")
    # calculates a negative log fold change for down regulation 
    neg.lfc <- lfc-(2*lfc)
    cat("log fold change less than: ", neg.lfc, "\n")
    # creates new table with probes that have appropriate adj.pval and lfc
    top.lfc <- top[[1]][(top[[1]]$logFC>lfc | top[[1]]$logFC<neg.lfc) & top[[1]]$adj.P.Val<adj.pval,]
    # sort table
    top.lfc <- top.lfc[order(-top.lfc$logFC, top.lfc$adj.P.Val),]
    # create expression set with top differentiall expressed genes
    expr.Sup <- norm.filt[[1]][rownames(norm.filt[[1]]) %in% rownames(top.lfc),]
    cat("differentially expressed genes: ",nrow(expr.Sup), "\n")
    nrow(expr.Sup)
  }
  else{
    print("You did not input an 'lfc' or 'adj.pval' argument(s) to the analysis function") 
    return(c(expr.Unsup, norm.filt[[1]], norm.filt[[2]], top[[1]], top[[2]]))
  }
  
  # display number of genes that were in both supervised and 
  # unsupervised differential expression
  #print("returns:  ")
  #print("[[1]] = ExprSet supervised, [[2]] = ExprSet unsupervised,") 
  #print("[[3]] = top table of signifcantly diff. expr. genes,")
  #print("[[4]] = normalized and filtered expression set,[[5]] = target file")
  if(length(treatmentVector) == 2){
    print("number of genes in both Superv. and Unsuperv: ")
    cat(sum(rownames(expr.Unsup) %in% rownames(expr.Sup)))
    return(list(expr.Sup, expr.Unsup, top.lfc, norm.filt[[1]], norm.filt[[2]], top[[1]], top[[2]], upDown.frame))
  }
  else{
    return(list(expr.Sup, top.lfc, norm.filt[[1]], norm.filt[[2]], top[[1]], top[[2]]))
  }
}

# cluster() Function:
# creates heatmap of differentially expressed genes
# exp.cl = expression set (i.e. analysis() output [[1]])
# targetFile = target file name (i.e. analysis() output [[5]])
# Do ::auto make color ar for different experiments::
cluster<- function(exp.cl, targetFile){
  # target file used to match sample IDs with experiment names
  target <- read.delim(targetFile, header=TRUE)
  # returns distance mat
  dist.eu <- dist(exp.cl, method="minkowski")
  # performs hierarchical clustering on the distance matrix
  hc.res<-hclust(dist.eu, method= "ward.D")
  #plot(hc.res, labels=hc.res$labels)
  numTreat<-length(unique(target$Treatment))
  clus.res<-cutree(hc.res, k=numTreat)
  hclust.ward<-function(d){hclust(d, method="ward.D")}
  windows()
  # produce heat map 
  # lables are automatic, as well as column side bar colors
   if(length(exp.cl@featureData@data$SYMBOL) <= 76){
  heat.res<-heatmap.2(exprs(exp.cl), scale="none", labRow=exp.cl@featureData@data$SYMBOL,labCol=as.character(target$Treatment),hclustfun = hclust.ward, col=topo.colors(100), cexCol=1.2, margins=c(8,10), trace="none", main="Differentially Expressed Genes Between Samples", dendrogram = c("column"), ColSideColors=topo.colors(sum(!duplicated(target$Treatment)))[as.numeric(target$Treatment)])
  }
  else{
   heat.res<-heatmap.2(exprs(exp.cl), scale="none", labRow="",labCol=as.character(target$Treatment),hclustfun = hclust.ward, col=topo.colors(100), cexCol=1.2, margins=c(8,10), trace="none", main="Differentially Expressed Genes Between Samples", dendrogram = c("column"), ColSideColors=topo.colors(sum(!duplicated(target$Treatment)))[as.numeric(target$Treatment)]) 
  }
  return(heat.res)
}

# annotation() function
# converts probeIDs to nuIDs to EntrezIDs using 'lumiMouseIDMapping' library
# requires 1 arguement: expression set
# returns expression set with entrez IDs and nuIDS
annotation <- function(exp.cl){
  # adds nuIDs to expression set
  exp.cl <- addNuID2lumi(exp.cl, lib.mapping='lumiMouseIDMapping')
  if(require(lumiMouseIDMapping)) {
    # makes the nuIDs feature names
    nuIDs <- featureNames(exp.cl)
    # converts the nuIDs to entrez IDs
    mappingInfo <- nuID2EntrezID(nuIDs, lib.mapping='lumiMouseIDMapping')
    head(mappingInfo)
  }
  # adds a new column in the feature data for entrez IDs
  featureData(exp.cl)$Entrez <- mappingInfo
  # remove blanks (nuIDs that didn't match to entrez IDs)
  exp.cl2 <- exp.cl[exp.cl@featureData@data$Entrez!="",]
  cat("duplicates = ")
  print(anyDuplicated(exp.cl@featureData@data$Entrez))
  exp.cl <- exp.cl2[!duplicated(exp.cl2@featureData@data$Entrez),]
  cat("duplicates = ")
  print(anyDuplicated(exp.cl@featureData@data$Entrez))
  rm(nuIDs)
  return(exp.cl)
}

# hyperG.comp() function
# performs hyper geometric tests on expression sets, retireves GO annotation
# this is important for pathway analysis
# requires 2 arguements:
# vehNorm.ann = annotation output containing the smaller group of genes
# vehX.ann = annotation for the gene universe (I use top 10% differentially expressed genes (after filtering))
#            calculated by the coefficient of variation 
# testDir = "over" or "under" expression
hyperG.comp<-function(vehNorm.ann, vehX.ann, testDir){ 
  # gets entrez ids to form gene universe
  universe.entrez<-vehX.ann@featureData@data$Entrez
  # gets smaller set to build sample
  clus.entrz<-vehNorm.ann@featureData@data$Entrez
  # create a new class that specifies the cluster, the universe, the annotation database, the type of ontology
  # (Biological Process), p-value cut off, etc 
  # over representation
  params<- new("GOHyperGParams",geneIds=clus.entrz, universeGeneIds=universe.entrez, annotation="org.Mm.eg.db", ontology='BP',pvalueCutoff=.05, conditional=FALSE,testDirection=testDir)
  # perform test
  hyperG <- hyperGTest(params)
  #returns: [[1]] hyperG output, and [[2]] entrez IDs to 
  return(c(hyperG, clus.entrz))
}

# vehNorm.ann = differential expressed eset annotation output (i.e. annotation(vehNormComp[[1]]))
# vehX.ann = top 10% diff. expressed eset, annotation out (i.e. annotation(vehNormComp[[2]]))
# topTab.all = top table for all filtered genes (i.e. vehNormComp[[6]])
# gsea.comp(vehNorm.ann, vehNorm.ann.all, vehNormComp[[6]])
gsea.comp<-function(vehTop.ann, vehX.ann, topTab.all){ 
  topTab.all<- head(topTab.all, .1*nrow(topTab.all))
  vehEpi.probe.adjpval.all <- data.frame(ProbeID=topTab.all$ProbeID, P.Value=topTab.all$adj.P.Val)
  vehEpi.probe.entrez.all<-data.frame(ProbeID=vehX.ann@featureData@data$ProbeID,Entrez=vehX.ann@featureData@data$Entrez)
  vehEpi.topTab.all <-merge(vehEpi.probe.adjpval.all,vehEpi.probe.entrez.all, by="ProbeID")
  vehEpi.vect.all<- c(vehEpi.topTab.all$P.Value)
  names(vehEpi.vect.all)<- vehEpi.topTab.all$Entrez
  
  topdiff.tab <-c(rownames(topTab.all) %in% vehTop.ann@featureData@data$ProbeID) 
  names(topdiff.tab)<- rownames(topdiff.tab)
  
  topDiffGenes <- function(allScore) {
  return(allScore < 0.2)
  }
  
  params<- new("topGOdata", description = "Simple session", ontology='BP', allGenes = vehEpi.vect.all, geneSel = topDiffGenes, nodeSize=5, annot=annFUN.org, mapping="org.Mm.eg.db", ID="entrez")
  test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
  resultFisher <- getSigGroups(params, test.stat)
  test.stat <- new("classicScore", testStatistic = GOKSTest, name = "KS tests")
  resultKS <- getSigGroups(params, test.stat)
  test.stat <- new("elimScore", testStatistic = GOKSTest, name = "Fisher test", cutOff = 0.01)
  resultElim <- getSigGroups(params, test.stat)
  test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
  resultWeight <- getSigGroups(params, test.stat)
  resultFis <- runTest(params, algorithm = "classic", statistic = "fisher")
  
  allRes <- GenTable(params, classic = resultFis, KS = resultKS, weight = resultWeight, orderBy = "weight", ranksOf = "classic", topNodes = 250)

  
   
  return(allRes)
}

# getGO() function
# uses biomaRt data base to convert entrez IDs to gene ontologies
# required arguements:
# exp.cl.both = output from annotation()
# posneg = output from analysis() contains table of probeIDs with plus and negative symbols
# returns:
# [[1]] = expression set that was mapped to GO ids
# [[2]] = a merged table containing information about each probe
getGO <- function(exp.cl.both, posneg){
  entrez <- exp.cl.both@featureData@data$Entrez
  ensembl <- useMart("ensembl")
  # specify mouse ensemble gene IDs
  mouse <- useDataset("mmusculus_gene_ensembl", mart=ensembl)
  # get BioMarker with entrez gene-ids, go-ids, and a biologucal description
  entrez.GO <- getBM(attributes=c('entrezgene','go_id', "description"), filters='entrezgene', values=entrez, mart=mouse)
  # write.table(entrez.GO$go_id, sep="\t", row.names=F, quote=F,file="Gorrilla_GO.txt")
  entrez.GO <- entrez.GO[entrez.GO$go_id != "",]
  # add probe IDs to the entrez GO annotation table
  probe.entrez<-data.frame(ProbeID=exp.cl.both@featureData@data$ProbeID,entrezgene=exp.cl.both@featureData@data$Entrez)
  entrez.GO <- merge(probe.entrez, entrez.GO, by="entrezgene")
  # creates expression set with probes that were matched to GO ontologies
  exp.cl.both <- exp.cl.both[exp.cl.both@featureData@data$Entrez %in% entrez.GO[,1]]
  # add positive negative expression to each probe in the GO annotation table
  posneg$ProbeID <- rownames(posneg) 
  entrez.GO<-merge(posneg, entrez.GO, by="ProbeID")
  return(list(exp.cl.both, entrez.GO))
}

# upDown()
# looks for matching negative and positive regulation patterns within the Riken transciption Factor database
upDown <- function(goset, transciptFact.db){
  colnames(goset)[5] <- "GO.term"
  transcrFact.DB<- read.delim(transciptFact.db, header=FALSE,skip=10)
  colnames(transcrFact.DB) <- c("GO.term", "Function", "Where")
  head(transcrFact.DB)
  sum(transcrFact.DB$GO.term %in% goset$GO.term)
  transcrFact.DB <- transcrFact.DB[transcrFact.DB$GO.term %in% goset$GO.term,]
  
  # clust1 <- clust.ALL[clust.ALL$cluster == 1,]
  clust1 <- goset[goset[,2]== "+",]
  trans.DB.cl1 <- transcrFact.DB[transcrFact.DB$GO.term %in% clust1$GO.term,]
  trans.DB.cl1.neg<-trans.DB.cl1[grepl("negative", trans.DB.cl1$Function),]
  
  clust2 <- goset[goset[,2]== "-",]
  trans.DB.cl2 <- transcrFact.DB[transcrFact.DB$GO.term %in% clust2$GO.term,]
  trans.DB.cl2.pos<-trans.DB.cl2[grepl("positive", trans.DB.cl2$Function),]
  up.down.neg <- data.frame(GO.term=character(), Function=character(), Where=character())
  for(i in 1:nrow(trans.DB.cl1.neg)){
    up.down.neg<-rbind(up.down.neg,trans.DB.cl1.neg[agrep(trans.DB.cl1.neg[i,2], trans.DB.cl2.pos$Function),])
    up.down.neg<-rbind(up.down.neg,trans.DB.cl2.pos[agrep(trans.DB.cl1.neg[i,2], trans.DB.cl2.pos$Function),])
  }
  up.down.neg<-na.omit(up.down.neg)
  nrow(up.down.neg)
  
  up.down.pos <- data.frame(GO.term=character(), Function=character(), Where=character())
  for(i in 1:nrow(trans.DB.cl1.neg)){
    up.down.pos<-rbind(up.down.pos,trans.DB.cl1.neg[agrep(trans.DB.cl2.pos[i,2], trans.DB.cl1.neg$Function),])
    up.down.pos<-rbind(up.down.pos,trans.DB.cl2.pos[agrep(trans.DB.cl2.pos[i,2], trans.DB.cl1.neg$Function),])
  }
  up.down.pos<-na.omit(up.down.pos)
  nrow(up.down.pos)
  
  up.down<-up.down.pos[up.down.pos$GO.term %in% up.down.neg$GO.term,]
  goset.trans<-merge(goset, up.down, by="GO.term")
  print(goset.trans[,c(3,4,7)])
  
  write.table(goset.trans, sep="\t", row.names=F, quote=F,file="up.down.reg.ALL.txt")
  return(goset.trans)
}

# wordSearch.GO() searchers Gene Ontologies for a specific word (i.e. oncogene)
# [[1]] = goset is the returned list() from the getGO function
# [[2]] = term is a string containing the specific word you are searching for in the GO ontologies
wordSearch.GO<-function(goset, term){
  goMatch.table <- goset[[2]][0,]
  for(i in 1:length(goset[[2]]$description)){
    stringVect.descr<-unlist(strsplit(goset[[2]]$description[i], " "))
    if(term %in% stringVect.descr){
      goMatch.table<-rbind(goMatch.table,goset[[2]][i,])
    }  
  }
  goMatch.expr <- goset[[1]][goset[[1]]@featureData@data$ProbeID %in% goMatch.table$ProbeID]
  return(list(goMatch.expr, goMatch.table))
}

# wordSearch.eset()
# searches for a specific word within the annotation of an expression set
wordSearch.eset <- function(exprSet, term){
  exprSet <- exprSet[exprSet@featureData@data$DEFINITION != "",]
  probeID <- c()
  for(i in 1:length(exprSet@featureData@data$DEFINITION)){
    stringVect <- unlist(strsplit(gsub("[[:punct:]]", " ", exprSet@featureData@data$DEFINITION[i]), " "))
    if(term %in% stringVect){
      probeID <- append(probeID, exprSet@featureData@data$ProbeID[i])
    }
  }
  matchExpr <- exprSet[exprSet@featureData@data$ProbeID %in% probeID,]
  return(matchExpr)
}

# compSupUnsup() 
# returns expression set 1 with probes that were matched to the 10% of differentiall expressed genes
# should contain all of exprset 1
compSupUnsup <- function(exprSet1_2){
  both <- exprSet1_2[[1]][exprSet1_2[[1]]@featureData@data$ProbeID %in% exprSet1_2[[2]]@featureData@data$ProbeID]
  return(both)
}

# compTreatment() 
# compares expression sets and returns a new expression containing data from exprset1 in exprset2 
compTreatment <- function(exprSet1, exprSet2){
  cat("number of probes in expression set 1: ", nrow(exprSet1), "\n")
  cat("number of probes in expression set 2: ", nrow(exprSet2), "\n")
  both <- exprSet1[exprSet1@featureData@data$ProbeID %in% exprSet2@featureData@data$ProbeID,]
  cat("number of probes in both: ", nrow(both), "\n")
  return(both)
}

# tumorMatch()
# matches gene symbols to known tumor suppressors found in the tumor_suppressor_gene_symbol.txt file
tumorMatch<-function(expr){
  tumor_suppressor_gene_symbol <- read.delim("tumor_suppressor_gene_symbol.txt", header=TRUE)
  exprTumor<-expr[expr@featureData@data$SYMBOL %in% tumor_suppressor_gene_symbol$Gene_symbol, ]
  cat("number of tumor suppressors found:", length(featureNames(exprTumor)) ,"\n")
  return(exprTumor)
}







