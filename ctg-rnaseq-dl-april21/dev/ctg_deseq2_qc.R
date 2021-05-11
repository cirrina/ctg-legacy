#!/usr/bin/env Rscript

# Main R script for DESeq2 analysis in bulkRNA seq analyses 
# based on nfcore pipeline 3.0 deseq2_qc.r
#

# added functionality 

# Core R functions are defied in source script. ctg_r_functions.R


################################################
################################################
## REQUIREMENTS                               ##
################################################
################################################

# NFCORE 
##  PCA, HEATMAP AND SCATTERPLOTS FOR SAMPLES IN COUNTS FILE
## - SAMPLE NAMES HAVE TO END IN e.g. "_R1" REPRESENTING REPLICATE ID. LAST 3 CHARACTERS OF SAMPLE NAME WILL BE TRIMMED TO OBTAIN GROUP ID FOR DESEQ2 COMPARISONS.
## - PACKAGES BELOW NEED TO BE AVAILABLE TO LOAD WHEN RUNNING R


## ========================================== ## 
## NF CORE vs CTG
## ========================================== ## 
# Diferences bethee

## ========================================== ## 
## LOAD PACKAGES                             ##
## ========================================== ## 
library(optparse)
library(DESeq2)
library(BiocParallel)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)

library(reshape2)
require(dplyr)
require(tidyr)
require(tidyverse)

## ========================================== ## 
## PARSE COMMAND-LINE PARAMETERS              ##
## ========================================== ## 
# sample_sheet added 
# should be add gene/tx info reference file?

option_list <- list(
  make_option(c("-i", "--count_file"    ), type="character", default=NULL    , metavar="path"   , help="Count file matrix where rows are genes and columns are samples."                        ),
  make_option(c("-s", "--sample_sheet"  ), type="character", default=NULL    , metavar="path"   , help="Sample Sheet, colData, where samples are rows and columns are sample annotations."                        ), ## added 
  # make_option(c("-i", "--gene_info_file"     ), type="character", default=NULL    , metavar="path"   , help="Transcript info file where genes/transcripts are rows and columns are sample annotations."        ), ## added 
  
  make_option(c("-f", "--count_col"     ), type="integer"  , default=9       , metavar="integer", help="First column containing sample count data. Detaulf for current featureCounts data is 9"                                             ),
  make_option(c("-d", "--id_col"        ), type="integer"  , default=1       , metavar="integer", help="Column containing identifiers to be used."                                              ),
  ## add filter if only protein coding
  make_option(c("-t", "--protein_coding"           ), type="logical"  , default=TRUE   , metavar="boolean", help="If to run using only protein coding genes. Requires that 'gene_type' column, describing transcipt types, is present in count data file"                                                     ),
  
  # make_option(c("-r", "--sample_suffix" ), type="character", default=''      , metavar="string" , help="Suffix to remove after sample name in columns e.g. '.rmDup.bam' if 'DRUG_R1.rmDup.bam'."),
  make_option(c("-o", "--outdir"        ), type="character", default='./deseq2_qc'    , metavar="path"   , help="Output directory."                                                                      ),
  make_option(c("-p", "--outprefix"     ), type="character", default='deseq2', metavar="string" , help="Output prefix."                                                                         ),
  make_option(c("-v", "--vst"           ), type="logical"  , default=FALSE   , metavar="boolean", help="Run vst transform instead of rlog."                                                     ),
  make_option(c("-c", "--cores"         ), type="integer"  , default=1       , metavar="integer", help="Number of cores."                                                                       )
)

opt_parser <- OptionParser(option_list=option_list)
opt        <- parse_args(opt_parser)

if (is.null(opt$count_file)){
  print_help(opt_parser)
  stop("Please provide a counts file.", call.=FALSE)
}

if (is.null(opt$sample_sheet)){
  print_help(opt_parser)
  stop("Please provide a sample_sheet file.", call.=FALSE)
}


## convert featureCounts output to nfore style count matrix file

## Mach up counts_file with 

## Temp files to work on 
# /Users/david/CTG_projects/bulkRNA/CTG_2020_177/sample_sheet.nf.csv
# opt$count_file <-  "/Users/david/CTG_projects/bulkRNA/CTG_2020_177/Quant/2020_177_geneid.featureCounts.txt"
# opt$sample_sheet <-  "/Users/david/CTG_projects/bulkRNA/CTG_2020_177/samplesheet_2020_177_nfcore_testdeseq.csv"

## ========================================== ## 
## PROCESS COUNTS FILE
## ========================================== ## 
# Counts table is read and hashes removed - 
## comment linesare includedwhen featureCounts produce output - these could be used to double check output


count.table           <- read.delim(file=opt$count_file,header=TRUE,comment="#")
# check if any duplicated ids
## current featureCounts will give --id_col as 1 and --count_col 9
if(any(duplicated(count.table[,opt$id_col]))) stop("Duplicated gene identifiers exist. Please check gene id column number --id_col")
rownames(count.table) <- count.table[,opt$id_col]
rowdata               <- count.table[,1:(opt$count_col-1)]
rownames(rowdata)     <- count.table[,opt$id_col] 
# str(rowdata)
count.table           <- count.table[,opt$count_col:ncol(count.table),drop=FALSE]
# chek if count columns have started in gene info or not
if(!all(unlist(unique(lapply(count.table, class)) %in% c("numeric","integer")))){
  stop("count.table contains non-integer columns. Check --count_col input var")
}
# colnames(count.table) <- gsub(opt$sample_suffix,"",colnames(count.table))
# colnames(count.table) <- gsub(pattern='\\.$', replacement='', colnames(count.table))


## ========================================== ## 
##  Process Sample Sheet (coldata) and match with data
## ========================================== ## 
# currenty, column namings, i.e. corresponding to sample names, should be supplied in sample sheet nf core style under header 'counts_matrix_id'
cat("\n ... fetching samples from 'counts_matrix_id column' in sample sheet.")
coldata <- read.csv(opt$sample_sheet, as.is=T)
## Check sample sheet !! 
## remember that  R will substitute hyphens (and spaces) in columns names to dots
## here get the colnames to match from counts_matrix_id column 
my_columns <- coldata$counts_matrix_id

## doucheck if duplicated samples
if(any(duplicated(my_columns))) stop(cat("There are duplicates in sanple sheet count matrix id column", 
                                         coldataFile, " \n ",
                                         paste(my_columns[duplicated(my_columns)], collapse = "; ")))


## check that all count matrix ids from sample sheet are columns in count matrix file
if(!all(my_columns %in% colnames(count.table))) stop("Mimatch between column names in count matrix file: \n", opt$count_file, "\nand count martrix ids in sample sheet" )

## seelect and match the sample sheet ids with count matrix. only the ids present in sample sheet will be analyzed
count.table <- count.table[,my_columns]
stopifnot(identical(colnames(count.table), coldata$counts_matrix_id))

## now set sample names and column names to original Sample_Names as presented in Sample Sheet
## This for simplicity and easy-use
if(!any(duplicated(coldata$Sample_Name))){
  rownames(coldata) <- coldata$Sample_Name
  colnames(count.table) <- coldata$Sample_Name
}else{
  rownames(coldata) <- coldata$counts_matrix_id
  }


## ===============================================  ## 
##   Barplot aligned Transcript Types - if avilable from GTV/gene info
## ===================================================================== ## 

## Boxplot/Barplot % of different transcript types
## each sample, % of total reads beloning to each ttype

# ttypes
# my_melt <- cbind(gene_type=factor(rowdata[,"gene_type"]), count.table)
# my_melt <- reshape2::melt(my_melt, varnames = c("gene_type"), value.name="counts", variable.name="sample") 
# str(my_melt)
# 
# my_melt <- my_melt %>%
#   group_by(gene_type, sample)  %>%
#   summarise(counts = sum(counts))
# my_melt <- as.data.frame(my_melt) %>% mutate(sample_id=paste(sample, gene_type, sep = "_")) %>% select(sample_id, everything())
# 
# plotViolinScatterRankBarGeneric(
#   x = data.frame(
#     sample_id = my_melt$sample_id,
#     dim1 = my_melt[, "counts"]
#   ),
#   pdata = my_melt,
#   pdata.column = "gene_type",
#   # color.key=color_key, 
#   drop.na.annotation = FALSE,
#   y.title = "", x.title = "",
#   plot.title = "AJCC stage", my.alpha = 1
# )

## ===============================================  ## 
##   PRE FILTER TRANSCRIPT TYPESS - Protein Coding  
## ================================================ ## 


if(opt$protein_coding){
  cat("\n ... Filtering. Keeping only protein coding transcripts")
  if("gene_type" %in% colnames(rowdata)){
    u <- rowdata$gene_type %in% "protein_coding"
    if(any(is.na(u))) stop("Error subsetting protein coding genes")
    count.table <- count.table[u,]
    rowdata <- rowdata[u,]
    opt$outprefix <- paste0(opt$outprefix, ".protein.coding")
  }else{
    stop("gene_type column not present in count matrix. cannot filter on protein coding genes. consider set --protein_coding to FALSE")
  }
}



## ========================================== ## 
##   RUN DESEQ2                               ##
## ========================================== ## 

#  opt$outdir <- "/Users/david/CTG_projects/bulkRNA/CTG_2020_177/r_test"
if (file.exists(opt$outdir) == FALSE) {
  dir.create(opt$outdir,recursive=TRUE)
}
setwd(opt$outdir)

# samples.vec <- sort(colnames(count.table))
# groups      <- sub("_[^_]+$", "", samples.vec)
samples.vec <- colnames(count.table)
groups <- as.factor(coldata$group) 
coldata$condition <- groups
# if (length(unique(groups)) == 1 || length(unique(groups)) == length(samples.vec)) {
#   quit(save = "no", status = 0, runLast = FALSE)
# }


DDSFile <- paste(opt$outprefix,".dds.RData",sep="")
if (file.exists(DDSFile) == FALSE) {
  counts  <- count.table[,samples.vec,drop=FALSE]
  # coldata <- data.frame(row.names=colnames(counts), condition=groups)
   # coldata <- data.frame(row.names=colnames(counts), condition=)
  dds     <- DESeqDataSetFromMatrix(countData=round(counts), colData=coldata, design=~ condition)
  dds     <- DESeq(dds, parallel=TRUE, BPPARAM=MulticoreParam(opt$cores))
  if (!opt$vst) {
    vst_name <- "rlog"
    rld      <- rlog(dds)
  } else {
    vst_name <- "vst"
    rld      <- varianceStabilizingTransformation(dds)
  }
  assay(dds, vst_name) <- assay(rld)
  save(dds,file=DDSFile)
} else {
  load(DDSFile)
  vst_name <- intersect(assayNames(dds), c("vst", "rlog"))
  if (length(vst_name)==0) { # legacy might mean vst was saved as a separate object called rld
    vst_name <- "loaded_rld"
    assay(dds, vst_name) <- assay(rld)
  } else {
    vst_name==vst_name[1]
  }
}


################################################
################################################
## FUNCTIONS                                  ##
################################################
###############################################
##' PCA pre-processeor
##' NF-core v3.0
##' Generate all the necessary information to plot PCA from a DESeq2 object
##' in which an assay containing a variance-stabilised matrix of counts is
##' stored. Copied from DESeq2::plotPCA, but with additional ability to
##' say which assay to run the PCA on, and adds an assessment of how well
##' each PC explains the experimental grouping of the data.
##' 
##' @param object The DESeq2DataSet object.
##' @param intgroup interesting groups: a character vector of names in 'colData(x)' to use for grouping.
##' @param ntop number of top genes to use for principla components, selected by highest row variance.
##' @param assay the name or index of the assay that stores the variance-stabilised data.
##' @return A data.frame containing the projected data alongside the grouping columns.
##' A 'percentVar' attribute is set which includes the percentage of variation each PC explains,
##' and additionally how much the variation within that PC is explained by the grouping variable.
##' @author Gavin Kelly
plotPCA_vst <- function (object, intgroup = "condition", ntop = 500, assay=length(assays(object))) {
  rv         <- rowVars(assay(object, assay))
  select     <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca        <- prcomp(t(assay(object, assay)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = ":"))
  }  else {
    colData(object)[[intgroup]]
  }
  d <- cbind(pca$x, group = group, intgroup.df, name = colnames(object))
  percentFrame <- data.frame(PC=seq(along=percentVar), percentVar=100*percentVar, groupR=0.0)
  for (ipc in seq(along=percentVar)) {
    fit1 <- lm(pca$x[,ipc]  ~ group)
    percentFrame$groupR[ipc] <- 100*summary(fit1)$r.squared
  }
  attr(d, "percentVar") <- percentFrame
  return(d)
}

PlotFile <- paste(opt$outprefix,".plots.pdf",sep="")
if (file.exists(PlotFile) == FALSE) {
  pdf(file=PlotFile,onefile=TRUE,width=7,height=7)
  
  ## PCA
  ntop <- c(500, Inf)
  for (n_top_var in ntop) {
    pca.data      <- plotPCA_vst(dds, assay=vst_name,intgroup=c("condition"),ntop=n_top_var)
    percentVar    <- round(attr(pca.data, "percentVar")$percentVar)
    plot_subtitle <- ifelse(n_top_var==Inf, "All genes", paste("Top", n_top_var, "genes"))
    pl <- ggplot(pca.data, aes(PC1, PC2, color=condition)) +
      geom_point(size=3) +
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance")) +
      labs(title = paste0("First PCs on ", vst_name, "-transformed data"), subtitle = plot_subtitle) + 
      theme(legend.position="top",
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=1))
    print(pl)
    
    pl <- ggplot(attr(pca.data, "percentVar"), aes(x=PC, y=percentVar)) +
      geom_line(aes(colour="explained by PC")) +
      geom_line(aes(y=groupR, colour="of PC explained by condition")) +
      scale_x_continuous(breaks=seq(along=percentVar), minor_breaks=NULL)  +
      labs(title="Diagnostics of PCs", subtitle=plot_subtitle, x="Component", y="Percentage explaned", colour="Percentage variation") +
      theme_bw() +
      theme(legend.position="top")
    print(pl)
    
    pc_r <- order(attr(pca.data, "percentVar")$groupR, decreasing=TRUE)
    pl <- ggplot(pca.data, aes_string(paste0("PC", pc_r[1]), paste0("PC", pc_r[2]), color="condition")) +
      geom_point(size=3) +
      xlab(paste0("PC", pc_r[1], ": ",percentVar[pc_r[1]],"% variance")) +
      ylab(paste0("PC", pc_r[2], ": ",percentVar[pc_r[2]],"% variance")) +
      labs(title = paste0("Group-Explanatory PCs of ", vst_name, "-tranformed data"), subtitle = plot_subtitle) + 
      theme(legend.position="top",
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=1))
    print(pl)
  } # at end of loop, we'll be using the user-defined ntop if any, else all genes
  
  ## WRITE PC1 vs PC2 VALUES TO FILE
  pca.vals           <- pca.data[,1:2]
  colnames(pca.vals) <- paste0(colnames(pca.vals), ": ", percentVar[1:2], '% variance')
  pca.vals           <- cbind(sample = rownames(pca.vals), pca.vals)
  write.table(pca.vals,file=paste(opt$outprefix,".pca.vals.txt",sep=""),row.names=FALSE,col.names=TRUE,sep="\t",quote=TRUE)
  
  ## SAMPLE CORRELATION HEATMAP
  sampleDists      <- dist(t(assay(dds, vst_name)))
  sampleDistMatrix <- as.matrix(sampleDists)
  colors           <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  pheatmap(
    sampleDistMatrix,
    clustering_distance_rows=sampleDists,
    clustering_distance_cols=sampleDists,
    col=colors,
    main=paste("Euclidean distance between", vst_name, "of samples")
  )
  
  ## WRITE SAMPLE DISTANCES TO FILE
  write.table(cbind(sample = rownames(sampleDistMatrix), sampleDistMatrix),file=paste(opt$outprefix,".sample.dists.txt",sep=""),row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)
  
  dev.off()
}


################################################
################################################
## PLOT QC                                    ##
################################################
################################################


PlotFile <- paste(opt$outprefix,".plots.pdf",sep="")
if (file.exists(PlotFile) == FALSE) {
  pdf(file=PlotFile,onefile=TRUE,width=7,height=7)
  
  ## PCA
  ntop <- c(500, Inf)
  for (n_top_var in ntop) {
    pca.data      <- plotPCA_vst(dds, assay=vst_name,intgroup=c("condition"),ntop=n_top_var)
    percentVar    <- round(attr(pca.data, "percentVar")$percentVar)
    plot_subtitle <- ifelse(n_top_var==Inf, "All genes", paste("Top", n_top_var, "genes"))
    pl <- ggplot(pca.data, aes(PC1, PC2, color=condition)) +
      geom_point(size=3) +
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance")) +
      labs(title = paste0("First PCs on ", vst_name, "-transformed data"), subtitle = plot_subtitle) + 
      theme(legend.position="top",
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=1))
    print(pl)
    
    pl <- ggplot(attr(pca.data, "percentVar"), aes(x=PC, y=percentVar)) +
      geom_line(aes(colour="explained by PC")) +
      geom_line(aes(y=groupR, colour="of PC explained by condition")) +
      scale_x_continuous(breaks=seq(along=percentVar), minor_breaks=NULL)  +
      labs(title="Diagnostics of PCs", subtitle=plot_subtitle, x="Component", y="Percentage explaned", colour="Percentage variation") +
      theme_bw() +
      theme(legend.position="top")
    print(pl)
    
    pc_r <- order(attr(pca.data, "percentVar")$groupR, decreasing=TRUE)
    pl <- ggplot(pca.data, aes_string(paste0("PC", pc_r[1]), paste0("PC", pc_r[2]), color="condition")) +
      geom_point(size=3) +
      xlab(paste0("PC", pc_r[1], ": ",percentVar[pc_r[1]],"% variance")) +
      ylab(paste0("PC", pc_r[2], ": ",percentVar[pc_r[2]],"% variance")) +
      labs(title = paste0("Group-Explanatory PCs of ", vst_name, "-tranformed data"), subtitle = plot_subtitle) + 
      theme(legend.position="top",
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=1))
    print(pl)
  } # at end of loop, we'll be using the user-defined ntop if any, else all genes
  
  ## WRITE PC1 vs PC2 VALUES TO FILE
  pca.vals           <- pca.data[,1:2]
  colnames(pca.vals) <- paste0(colnames(pca.vals), ": ", percentVar[1:2], '% variance')
  pca.vals           <- cbind(sample = rownames(pca.vals), pca.vals)
  write.table(pca.vals,file=paste(opt$outprefix,".pca.vals.txt",sep=""),row.names=FALSE,col.names=TRUE,sep="\t",quote=TRUE)
  
  ## SAMPLE CORRELATION HEATMAP
  sampleDists      <- dist(t(assay(dds, vst_name)))
  sampleDistMatrix <- as.matrix(sampleDists)
  colors           <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  pheatmap(
    sampleDistMatrix,
    clustering_distance_rows=sampleDists,
    clustering_distance_cols=sampleDists,
    col=colors,
    main=paste("Euclidean distance between", vst_name, "of samples")
  )
  
  ## WRITE SAMPLE DISTANCES TO FILE
  write.table(cbind(sample = rownames(sampleDistMatrix), sampleDistMatrix),file=paste(opt$outprefix,".sample.dists.txt",sep=""),row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)
  
  dev.off()
}