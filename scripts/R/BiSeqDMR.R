################################################################################
# Copyright (C) 2016 Johanna Sandoval
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# If not, see <http://www.gnu.org/licenses/>.
################################################################################


# using the BiSeq library to test for differentially methylated regions
# source("http://bioconductor.org/biocLite.R")
# biocLite("BiSeq")
# library("BiSeq")



usage=function(errM) {
     cat("\nUsage : Rscript BsSeqDMR.R [option] \n")
     cat("       -i        : input file(s), a comma separated string for path to input files (methylated profiles containting chromosome, position, coverage and methyl counts, tab separated, header=T) \n")
     cat("       -s        : sampleNames, a comma separated string (used if multiple files and each file correspond to a specific sample) \n")
     cat("       -o        : output prefix\n")
     cat("       -c        : coverage counts columns in input files\n")
     cat("       -m        : methylated counts in input files\n\n")
     cat("       -d        : path to the design file\n")
     cat("       -p        : number of cores to use\n")
     cat("       -t        : path to the targeted capture bed file\n")
     cat("       -h        : this help\n\n")     
     stop(errM)
}

path_target=NULL
targets.gr=NULL
outFile="BiSeqDMR"
set.seed(123456789)
mc.cores=8

ARG = commandArgs(trailingOnly = T)
 ## get arg variables
for (i in 1:length(ARG)) {
 if (ARG[i] == "-i") {
 input_files=ARG[i+1]
 } else if (ARG[i] == "-s") {
 sampleNames=ARG[i+1]
 } else if (ARG[i] == "-o") {
 outFile=ARG[i+1]
 } else if (ARG[i] == "-c") {
   coverage.counts.cols=as.numeric(split(ARG[i+1],","))
 } else if (ARG[i] == "-m") {
   methylated.counts.cols=as.numeric(split(ARG[i+1],","))
 } else if (ARG[i] == "-d") {
   design_file=ARG[i+1]
} else if (ARG[i] == "-p") {
    mc.cores=as.integer(ARG[i+1])
} else if (ARG[i] == "-t") {    
    path_target=ARG[i+1]       
} else if (ARG[i] == "-h") {        
    usage("Help:")
 }
}
#     list_files=list.files(path="hg19/methylation", pattern="profile.cg_strand_combined.csv$", full.names=T, include.dirs=T);
#      list_files=list.files(path="test", pattern="profile.cg_strand_combined.csv$", full.names=T, include.dirs=T);
#     input_files=paste(list_files,collapse=",")
#      sampleNames=paste(unlist(lapply(basename(list_files), function(x) gsub(".profile.cg_strand_combined.csv","", x, fixed=T))),collapse=",")
#      outFile="DMR"
#      coverage.counts.cols=11
#      methylated.counts.cols=10
#      design_file="design.tsv"
#      path_target="SeqCap_EPI_4M_CPGiant.targets.bed"

library("BiSeq")
library(data.table)

sample_names=unlist(strsplit(sampleNames,","))
input_files=unlist(strsplit(input_files,","))
# Minimun coverage for methylated sites
minCoverage=2
# targets
if(file.exists(path_target)){
    targets <- read.table(path_target,header=F,stringsAsFactors=F)
    targets.gr <- GRanges(targets[,1], IRanges(targets[,2] + 1, targets[,3]))
    print(paste("Notice: added a target bed file"))
}
print(targets.gr)

if(length(input_files)<length(sample_names)) { 
    if(length(input_files)==1){
        warn("ERROR: Sample names and input files are not from the same length, only one input file? "); 
        sampleNames=lapply(sample_names, function(x) sampleNames[[x]]=sample_names)
    }            
}else{
    sampleNames=sample_names
}

# Parse NxtGenUtils output, tab separated file containing chromosome, position, and strand merged methylation/total counts 
readMyOutput <- function(files, colData, methylated.counts.cols=10,coverage.counts.cols=11, target.gr=NULL) {
    if (nrow(colData) != length(files)) {
        stop("Row number of colData must equal length of files.")
        }
    methData = list()
    for (i in 1:length(files)) {
        cat(paste("Processing sample ", rownames(colData)[i], " ... \n", sep=""))
        #         bismark <- scan(files[i], skip=0, sep="\t",
        #                         what=list("character", integer(), NULL, NULL, integer(), integer()))
        #                 
        df0=fread(files[i], header=T)
        M = as.integer(unlist(df0[,methylated.counts.cols,with=FALSE]))
        Cov = as.integer(unlist(df0[,coverage.counts.cols,with=FALSE]))
        rowRanges=GRanges(
            seqnames=unlist(df0[,1,with=FALSE]),
            ranges=IRanges(start=as.numeric(unlist(df0[,2,with=FALSE])), width=1),
            methylated=M,
            reads=Cov)
        
        if(!is.null(target.gr)){
                methData[[i]]=subsetByOverlaps(rowRanges, target.gr)            
            }else{
                methData[[i]]=rowRanges                
            }
        print(methData[[i]])
	rm(df0)
    }
    cat("Building BSraw object.\n")
    fData <- methData[[1]]
    if(length(methData) > 1){
        for(i in seq(along=methData)[-1]){
            fData <- unique(c(fData, methData[[i]]))
            }
        }
        elementMetadata(fData) <- NULL
        names(fData) <- as.character(1:length(fData))
        tReads <- matrix(integer(length = length(fData) * length(methData)), nrow=length(fData))
        mReads <- matrix(integer(length = length(fData) * length(methData)), nrow=length(fData))
        for(i in seq(along=methData)){
            mtch <- findOverlaps(fData, methData[[i]])
            ind1 <- queryHits(mtch)
            ind2 <- subjectHits(mtch)
            tReads[ind1, i] <- mcols(methData[[i]])[ind2,"reads"]
            mReads[ind1, i] <- mcols(methData[[i]])[ind2,"methylated"]
            }
        # 
        colnames(tReads) <- rownames(colData)
        colnames(mReads) <- rownames(colData)
        rownames(tReads) <- names(fData)
        rownames(mReads) <- names(fData)
        rrbs = BSraw(
            colData = colData,
            rowRanges = fData,
            totalReads = tReads,
            methReads = mReads)
        return(rrbs)
}

# Perform DMR 
perform_DRM<-function(rrbs, group , cutoff.quantile=0.975, minCoverage, output="BiSeqDMR.tsv", plot.top.regions=50, FDR.cluster=0.1, FDR.loc=0.05, mc.cores=1, regression='logistic')
{
    # Cluster sites must be done without groups
#     # rrbs.clust.unlim <- clusterSites(object = rrbs, groups = colData(rrbs)[,group], perc.samples = 4/5, min.sites = 20, max.dist = 100)
    rrbs.clust.unlim <- clusterSites(object = rrbs, perc.samples = 4/5, min.sites = 20, max.dist = 100)
    print(rowRanges(rrbs.clust.unlim))
    clusterSitesToGR(rrbs.clust.unlim)
    # Smooth methylation data
    ind.cov <- totalReads(rrbs.clust.unlim) > minCoverage
    quant <- quantile(totalReads(rrbs.clust.unlim)[ind.cov], cutoff.quantile)
    rrbs.clust.lim <- limitCov(rrbs.clust.unlim, maxCov = quant )
    predictedMeth <- predictMeth( object = rrbs.clust.lim , mc.cores=mc.cores )
    print(predictedMeth)
    # Model and test group effect
    form=as.formula(paste0("~",group))
    if ( regression == "logistic"){
        betaResults <- logisticRegression(formula = form, link = "logit", object = predictedMeth, 
                                          mc.cores=mc.cores)
    }else{
        betaResults <- betaRegression(formula = form,
                                  link = "probit",
                                  object = predictedMeth,
                                  type = "BR", mc.cores=mc.cores)
    }
    # sort by p-val and print
    betaResults <- betaResults[order(abs(betaResults$p.val)),]
    write.table(betaResults, gsub(".tsv$", "betaResults.tsv",output) ,quote=F,row.names=F,col.names=T,sep="\t")
    # Test CpG clusters for differential methylation
    # need to resample individuals and group assignation to creata a null model
    ## Both resampled groups should have the same number of control and treatments
    predictedMethNull <- predictedMeth[, sample(1:ncol(predictedMeth), size = ncol(predictedMeth))]
    colData(predictedMethNull)$group.null <- rep(c(1,2), ncol(predictedMeth)/2)
    ## To shorten the run time, please set mc.cores, if possible!
    print(predictedMethNull)
    if ( regression == "logistic"){
        betaResultsNull  <- logisticRegression(formula = ~group.null, link = "logit", object = predictedMethNull, mc.cores=mc.cores)
    }else{    
        betaResultsNull <- betaRegression(formula = ~group.null,
                                            link = "probit",
                                            object = predictedMethNull,
                                            type="ML", mc.cores=mc.cores)
    }    
    
    vario <- makeVariogram(betaResultsNull)
    ## Auxiliary object to get the pValsList for the test
    ## results of interest:
    vario.aux <- makeVariogram(betaResults, make.variogram=FALSE)
    # Based on the variogram plot we evaluate the sill (usually near 1) of the variogram and smooth the curve:
    vario.sm <- smoothVariogram(vario, sill = 0.9)
    vario.sm$pValsList <- vario.aux$pValsList
    #The correlation of the    Z    scores between two   locations in a cluster can now be estimated
    locCor <- estLocCor(vario.sm)
    #We test each CpG cluster for the presence of at least one differentially methylated site
    clusters.rej <- try(testClusters(locCor, FDR.cluster = FDR.cluster))
    # 
    if (!inherits(clusters.rej, "try-error") )
    {
        print(clusters.rej$clusters.reject)
        #Trim significant CpG clusters
        clusters.trimmed <- trimClusters(clusters.rej,FDR.loc = FDR.loc )
    }else{
        clusters.trimmed <-betaResults 
    }
    DMRs <- findDMRs(clusters.trimmed, max.dist = 100, diff.dir = TRUE)
    write.table(as.data.frame(DMRs), output ,quote=F,row.names=F,col.names=T,sep="\t")
    df <- data.frame(seqnames=seqnames(DMRs),
                     starts=as.integer(unlist(start(DMRs)-1)),
                     ends=as.integer(unlist(end(DMRs))),
                     names=apply(as.data.frame(DMRs)[,-c(1:5)], 1, function (x) paste(paste(names(mcols(DMRs)), x, sep="="), collapse=",")),
                     scores=lapply(mcols(DMRs)["median.meth.diff"], function(x) round(x,3)),
                     strands="+")        
    write.table(df, file=gsub(".tsv$", "_DMRs.bed",output), quote=F, sep="\t", row.names=F, col.names=F)
    
    # Save Bed files
    rrbs<-rawToRel(rrbs)
    track.names <- paste(colnames(rrbs), colData(rrbs)[,group],sep="_")
    writeBED(object = rrbs,
               name = track.names,
               file = paste0(paste(gsub(".tsv$", "",output), track.names, sep="_"), "_raw_data.bed"))
    writeBED(object = predictedMeth,
               name = track.names,
               file = paste0(paste(gsub(".tsv$", "",output), track.names, sep="_"), "_smoothed_data.bed"))
    
    save(predictedMeth,DMRs,rrbs.clust.lim,betaResultsNull , betaResults, file=gsub(".tsv$", ".RData",output))
    # QC PLOTS and methylation plots are done in a separated file    
}

# Read design file
# Create BSSeq object, one per sample or one per input file if 1 input file
pData=read.csv2(design_file, header=T, sep = "\t", na.strings = "0", check.names=F, colClasses =c('character',rep('character',unique(count.fields(design_file))-1)))
rownames(pData)=pData[,1]
print (paste0("NOTICE: Creating a Bsraw object from: ", paste(input_files, collapse=",")))
#BS.all=readMyOutput(input_files,  DataFrame(pData)) 
#print(BS.all)

# Detection of clusters
#rrbs.small <- BS.all[1:1000,]

# Now ready to find differentially methylated regions in samples grouped and contrasted in a design file 
name_sample= as.character(as.vector(pData[,1]))
                         
for (i in 2:ncol(pData)) {
    name_folder = paste(outFile, names(pData[i]), sep="/")
     # Create output directory     
    if (!file.exists(name_folder)) { 
        system(paste("mkdir -p",name_folder,sep=" "))
    }
    current_pData=pData[,i]
    subsampleN=name_sample[!(is.na(current_pData))]
    group=names(pData[i])
    cat("NOTICE: Processing contrast \n")
    cat(paste("Name folder: ",name_folder,"\n",sep=""))
    cat(paste("Design : ",paste(subsampleN, group, sep="=",collapse=" ; "),"\n",spe=""))
    
    BS.all=readMyOutput(input_files[which(sampleNames %in% subsampleN )],  DataFrame(subset(pData,pData[,1] %in% subsampleN)), methylated.counts.cols, coverage.counts.cols, targets.gr) 
    print(BS.all)    
    # Perform Differentially methylated regions
    try(perform_DRM(BS.all, group, cutoff.quantile=0.975,  minCoverage= minCoverage, output=paste(name_folder, "DRM_stats.tsv", sep="/"), FDR.cluster=0.1, FDR.loc=0.05, mc.cores=mc.cores))
    # rrbs=BS.all; group=group; cutoff.quantile=0.975; minCoverage= minCoverage; output=paste(name_folder, "DRM_stats.tsv", sep="/"); plot.top.regions=50 ; FDR.cluster=0.1; FDR.loc=0.05; mc.cores=1
}
