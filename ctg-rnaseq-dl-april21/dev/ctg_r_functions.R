# 
# ## Load packages
# require(Biobase)
# require(dplyr)
# require(tidyr)
# require(tidyverse)
# require(ggplot2)
# require(plotly)
# require(DESeq2)
# require(sva)
# require(SummarizedExperiment)
# require(vsn)
# require(glmpca)
# require(hexbin)
# require(pheatmap)
# require(RColorBrewer)
# require(PoiClaClu)
# require(plyr)
# require(reshape2)
# require(rlang)
# require(DT)
# require(apeglm)
# require(Rtsne)



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



ttypes <- c(
  "transcribed_unprocessed_pseudogene",
  "unprocessed_pseudogene",
  "miRNA",
  "lncRNA",
  "processed_pseudogene",
  "transcribed_processed_pseudogene",
  "snRNA",
  "protein_coding",
  "TEC",
  "transcribed_unitary_pseudogene",
  "misc_RNA",
  "snoRNA",
  "unitary_pseudogene",
  "scaRNA",
  "polymorphic_pseudogene",
  "rRNA_pseudogene",
  "pseudogene",
  "rRNA",
  "scRNA",
  "IG_C_gene",
  "IG_V_gene",
  "IG_V_pseudogene",
  "translated_processed_pseudogene",
  "TR_V_gene",
  "TR_V_pseudogene",
  "translated_unprocessed_pseudogene",
  "TR_C_gene",
  "TR_J_gene",
  "IG_C_pseudogene",
  "ribozyme",
  "IG_J_gene",
  "IG_J_pseudogene",
  "IG_D_gene")

## Scatterpoints/rankplot with annotations in pdata. variant of the scatterPoint_annotationPlot
#' plotViolinScatterRankBarGeneric
#' @description Generic plot for tSNE (list from tsneWrapper analysis), umap (umap object), or PCA (prcomp object) but also for generic data frames and for plotting Ranks only
#' @description If data.frama supplied. Then:
#' @description - Bar plot: Only when x is data.frame and both dim1 non-numeric characrer/factors. Primary annotation supplied in pdata/pdata.column.
#' @description - Violin: only dim1 numeric OR if eset, gene.id (but not gene.id_2) supplied
#' @description - Rank plot: as for Violin but rank.plot set to TRUE,
#' @description - Scatter plot: dim1 and dim2 numeric OR if eset, gene.id and gene.id.2 supplied
#' @family umap
#' @family tsne
#' @family pca
#' @family Violin
#' @family RankPlot
#' @family Barplot
#' @family ScatterPlot
#' @family expressonset
#' @family PlotWrap
#' @param x results object. list (tSNEres) from run.tsne-analysis, or umap (umap object), or PCA (prcomp object). Can also be a data.frame (sample_id, dim1, dim2) or an expressionSet
#' @param gene.id supplied if x is an ExpressionSet and if scatter or rankplot
#' @param gene.id.2 supplied if x is an ExpressionSet and if scatterplot
#' @param gene.id.signature supplied if x is expression set and if to create a score for violin plot over multiple genes
#' @param pdata Optional. Pairs 'sample_id' with the "annotation' used for group-colors. If not provided an 'annotation' column must be present in x. data frame with, minimally, columns 'sample_id' and a column for annotation (defined by 'annotation.column').
#' @param pdata.column should be a column in the pdata object to use for annotation (colors).
#' @param violin.2nd.group.column specify a columnname if to split violin groups into secondary groups
#' @param drop.na.annotation if to drop NA values from annotation column. defaults TRUE. May be set to FALSE for e.g. barplots to get a better representation of the data.
#' @param my.samples if to highligt samples - character vector of sample names
#' @param plot.title Title of plot
#' @param x.title
#' @param y.title
#' @param axis_names_n if to add number of samples to sample groups (x-axis names). Note: not feasible if downstream facet_wrap
#' @param color.key what color key to use for annotation. Detaults to dlfoo2::color_subtypes (but could be set to e.g. dlfoo2::color_annotations)
#' @param do.gradient set to FALSE if to manually override numeric check for gradient or manual fil
#' @param color.key.gradient what gradient.
#' @param gradient.nomalize TRUE or FALSE if to z normalize annotation value (x-mean(x)) sd(x)
#' @param gradient.squish.probs if to cap squish gradient to remove effects on color scale of outliers. vector of tow percent values. Uses quantile-probs. Defaults to c(0.05, 0.95),i.e. cap scale at 5 and 95 of values
#' @param stripchart.plot set to TRUE if to do stripchart instead of scatterplot. Provided as dim1 (values) and dim2 (groupings).
#' @param rank.plot set to TRUE if to do Rank plot instead of Violin. Defaults to False. Possile if dim1 onöy is  supplied (not dim1 and dim2, then scatterplot).
#' @param rank.decreasing set to True if Rank plot and if rank decerasing on x-scale
#' @param barplot.y.freq TRUE if barplot y-axis should be frequency percent, 0-100 (not total N)
#' @param my.size point size
#' @param my.alpha point alpha
#' @param my.stroke stroke width
#' @param my.stroke.color stroke line colors
#' @param do.boxplot set to true if boxplot insead of violin
#' @param showGuide if to show guide
#' @param coord.fixed if coord_fixed ggplot param
#' @param swap.axes if to swap x and y axes
#' @param xlim.cap  CAPs x - values (dim2) if any outsidde
#' @param xlim.rescale passed on to axis. vector of 2. set upper/lower to NA if only rescaling one value
#' @param ylim.cap  CAPs x - values (dim2) if any outsidde
#' @param ylim.rescale passed on to axis. vector of 2. set upper/lower to NA if only rescaling one value
#' @param zero.lines if to higlight h- and -vlines through origo
#' @param shape.key if to plot annotations in different shapes (connected to annotation)
#' @param zoom.coord list of x and y coordinates to specify/zoom plot coordinates
#' @return a ggplot object
#' @export
plotViolinScatterRankBarGeneric <- function(
  x,
  gene.id=NULL, gene.id.2=NULL, gene.id.signature = NULL,
  pdata=NULL, pdata.column=NULL,
  plot.title=NULL, x.title=NULL, y.title=NULL,
  axis_names_n=T,
  rank.plot = F,
  stripchart.plot =F,
  violin.2nd.group.column = NULL,
  drop.na.annotation = T,
  barplot.y.freq = T,
  rank.decreasing = F,
  my.size=2, my.alpha=0.7, my.stroke=0.05, my.stroke.color="gray10",
  color.key=NULL,
  do.gradient=T, pseudo.gradient=T, gradient.sat = c(-1,1),
  palette.gradient=rev(dlfoo2::palette_gradients[["red2white2blue"]]), gadient.na.col=NA, gradient.normalize=F, gradient.squish.probs=c(0.05, 0.95),
  my.samples=NULL,
  do.boxplot = F,
  swap.axes =F,
  showGuide=F, coord.fixed=F, xlim.cap=NULL, ylim.cap=NULL, xlim.rescale=NULL, ylim.rescale=NULL, zero.lines=F,
  shape.key=NULL, zoom.coord=NULL
){
  
  
  require(umap)
  require(Rtsne)
  require(dplyr)
  require(tidyverse)
  require(tibble)
  require(ggrepel)
  # outputPlotly=T
  #  if plotly text, add  text =  paste("Point: ", annotation, "\n", sample_id, "\n", annotation2)
  
  plot_type <- NULL
  
  ## if umap/tsne/pca then ScatterPlot
  ## -------------
  if(class(x)=="umap"){
    message("Found UMAP results. Begin plotting ... ")
    x_df <- umap2df(x)
    plot_type <- "scatterPlot"
  }
  
  if(class(x)=="prcomp"){
    message("Found PCA results. Begin plotting ... ")
    plot_type <- "scatterPlot"
    x_df <- pca2df(x)
  }
  
  if(class(x)=="list"){
    if(!is.null(x$Y)){
      message("Found list with Y matrix slot  ... ")
      message("... assuming x is a tsne reults list. Begin plotting")
      plot_type <- "scatterPlot"
      x_df <- tsne2df(x)
    }}
  
  
  
  ## Expression set
  ## -------------
  if(class(x)=="ExpressionSet"){
    message("expressionset provided")
    
    x_samples <- intersect(unique(pdata$sample_id), sampleNames(x))
    
    
    gene_i = match(gene.id, featureNames(x))
    
    if(!is.null(gene.id.signature)){ # if multiple genes - gene signature
      gene.id.signature <- unique(gene.id.signature[!is.na(gene.id.signature)])
      gene_i = match(gene.id.signature, featureNames(x))
      
    }
    if(all(is.na(gene_i))){
      message(gene.id, " not found in feature data for pan_gdc")
      return(NULL)
    }
    
    if(!is.null(gene.id.2)){
      # if gene_2 is supplied both dim1 and dim2 has to be numeric and set plot to scatterPlot
      plot_type <- "scatterPlot"
      gene_ii = match(gene.id.2, featureNames(x))
      if(is.na(gene_i)){
        message(gene.id.2, " not found in feature data for pan_gdc")
        return(NULL)
      }
      
      x_df <- data.frame(sample_id=sampleNames(x),
                         dim1 = as.numeric(Biobase::exprs(x[gene_i, x_samples])),
                         dim2 = as.numeric(Biobase::exprs(x[gene_ii, x_samples]))
      )
      x_df <- x_df %>% dplyr::filter(!is.na(dim1)) %>% dplyr::filter(!is.na(dim2))
    }else{
      # if not dim2 is supplied then TWO possible plot types (ViolinPlot or rankPlot) defined by 'rank.plot' param
      ## NOTE: stripchart.plot=T will override rank.plot=T
      if(rank.plot){
        plot_type <- "rankPlot"
        message("... setting plot type to rankPlot")
      } # both dim1 and dim2 has to be neumeric}
      
      if(!rank.plot) {
        plot_type <- "violinPlot"
        message("... setting plot type to violinplot")
      }
      
      ## create pdata_plot with matching sample_id as the eset sample names
      if(length(gene_i)==1){ # if only one gene
        x_df <- data.frame(sample_id=x_samples, dim1 = as.numeric(Biobase::exprs(x[gene_i, x_samples])))
        message("...removing NA-values")
        x_df <- x_df %>% dplyr::filter(!is.na(dim1))
      }
      if(length(gene_i)>1){
        message("\n ... multiple genes supplied - creting gene score based on mean value")
        es_sig <- esetCenterGenes(x[gene_i,])
        x_df <- data.frame(sample_id=x_samples, dim1 = as.numeric(apply(Biobase::exprs(es_sig[,x_samples]), 2, mean, na.rm=T)))
        message("...removing NA-values")
        x_df <- x_df %>% dplyr::filter(!is.na(dim1))
      }
      
    } # end if else
    
    
  } # end if x is expression set
  
  
  
  
  ## Data frame
  ## -------------
  if(class(x) == "data.frame"){
    message("Found data frame ... looking for sample_id and dim1 as columns")
    if(!all(c("sample_id","dim1") %in% colnames(x))) stop()
    
    message("...removing NA-values")
    x_df <- x %>% dplyr::filter(!is.na(dim1)) %>% dplyr::filter(!is.na(sample_id))
    x_samples <- unique(as.character(x_df$sample_id))
    x_df <- x_df[match(x_samples, x_df$sample_id), ]
    
    if("dim2" %in% colnames(x_df)){
      message(" ... found dim2")
      # if dim2 is supplied both dim1 and dim2 has to be numeric and set plot to scatterPlot
      # OR if stripchart.plot=T then allow di2 as not numeric
      if(stripchart.plot){
        plot_type <- "stripChart"
        x_df$dim2 <- as.character(x_df$dim2)
      }else{
        if(is.numeric(x$dim2)) plot_type <- "scatterPlot"
        if(!is.numeric(x$dim2)) stop("dim2 must be numeric. If barPlot then supply binning annotation (y-axis) through the pdata and pdata.colums parameters while dim1 should be the primary grouping annotation (x-axis) ")
      }
      message(" ... setting plot type to:  ",plot_type)
      x_df <- x_df %>% dplyr::filter(!is.na(dim2))
      
    }else{
      # if  not dim2 is supplied then 3 possible plot types
      if(is.numeric(x_df$dim1) && rank.plot) plot_type <- "rankPlot"
      if(is.numeric(x_df$dim1) && !rank.plot) plot_type <- "violinPlot"
      if(!is.numeric(x_df$dim1) && !is.numeric(pdata[,pdata.column])){
        plot_type <- "barPlot"
        if(!is.factor(x_df$dim1)) x_df$dim1 <- factor(x_df$dim1)
        x_df <- droplevels(x_df)
        if(is.null(pdata)) stop("pdata must be provided for barPlot")
      }
      message(plot_type)
    } # end if else
  } # end if dataframe
  
  
  
  # if rank plot then create dim2 (x-axis) rank value
  if(rank.plot && plot_type=="rankPlot"){
    x_df$dim2 <- rank(x_df$dim1, ties.method = "random")
    if(rank.decreasing) x_df$dim2 <- rank((-1)*x_df$dim1, ties.method = "random")
  }
  
  # stop if no x
  if(is.null(plot_type)) stop("x must be an umap, tsne, expressionSet or a data frame with sample_id and dim1 as columns ")
  
  
  
  ## If no pdata provided - then a annotation column must be found in x
  ## ------------
  if(is.null(pdata)){
    message("... no pdata provided\n ... ... checking annotation column in x")
    if(!"annotation" %in% colnames(x_df)) stop("x must contain 'annotation' if no pdata is provided !")
    my_df <- x_df
    ## Check annotation
    if(!is.factor(my_df$annotation) && !is.numeric(my_df$annotation)){
      my_df$annotation <- factor(my_df$annotation)
    }
  }
  
  ## check pdata - add annotation (group - colors) to my_df
  ## ------------
  ## Check x and pdata
  if(!is.null(pdata)){
    message("...  pdata provided\n ... ... checking annotations in x")
    message("... ... integrating pdata into x")
    if(!"sample_id" %in% colnames(pdata)) stop("pdata object must contain 'sample_id' column")
    if(!(pdata.column %in% colnames(pdata))) stop("pdata.column not in selected pdata object")
    
    x_samples <-  unique(intersect(pdata$sample_id, x_df$sample_id))
    my_df <- as.data.frame(
      pdata[match(x_samples, pdata$sample_id),] %>%
        # mutate(annotation = rlang::UQ(rlang::sym(pdata.column))) %>%
        dplyr::rename(annotation=pdata.column) %>% 
        dplyr::select(c(sample_id, annotation, everything())))
    
    ## Check annotation
    if(!is.factor(my_df$annotation) && !is.numeric(my_df$annotation)){
      my_df$annotation <- factor(my_df$annotation)
    }
    suppressWarnings(my_df <- x_df %>% dplyr::filter(sample_id %in% x_samples) %>% dplyr::left_join(my_df, by="sample_id"))
    # Drop NA from annotation colums (default)
    if(drop.na.annotation) my_df <- my_df %>% dplyr::filter(!is.na(annotation))
    if(!drop.na.annotation){
      temp_levels <- levels(my_df$annotation)
      my_df <- my_df %>%
        mutate(annotation=as.character(annotation)) %>%
        mutate(annotation = if_else(is.na(annotation), "_NA_", annotation)) %>%
        mutate(annotation = factor(annotation, levels = c(temp_levels, "_NA_")))
    }
    my_df <- my_df %>% droplevels() %>% arrange(annotation)
    # str(my_df)
    
  } # end integrate pdata
  
  
  ## Check for Inf values
  if(any(my_df$dim1 %in% c("Inf","-Inf"))){
    message("!! WARNING: 'Inf' values found -possibly after log - replacing thse with max/min of non Inf values")
    my_df <- my_df %>% dplyr::mutate(dim1 = if_else(dim1=="Inf", max(my_df$dim1[!my_df$dim1=="Inf"]), dim1))
    my_df <- my_df %>% dplyr::mutate(dim1 = if_else(dim1=="-Inf", min(my_df$dim1[!my_df$dim1=="-Inf"]), dim1))
  }
  
  if("dim2" %in% colnames(my_df)){
    if( any(my_df$dim2 %in% c("Inf","-Inf"))){
      message("Inf values found - replacing thse w max/min")
      my_df <- my_df %>% dplyr::mutate(dim2 = if_else(dim2=="Inf", max(my_df$dim2[!my_df$dim2=="Inf"]), dim2))
      my_df <- my_df %>% dplyr::mutate(dim2 = if_else(dim2=="-Inf", min(my_df$dim2[!my_df$dim2=="-Inf"]), dim2))
    }
  }
  
  
  
  ## Add manual samples
  ## ----------
  # if manually added my.samples Then label these
  ## ??? working?
  if(!is.null(my.samples)){
    if(length(which(my.samples %in% my_df$sample_id)) > 0) {
      my_df$plot_text <- NA
      my_df$plot_text[my_df$sample_id %in% my.samples] <- as.character(my_df$sample_id[my_df$sample_id %in% my.samples] )
    }
    if(length(which(my.samples %in% my_df$sample_id)) < 1) {
      message("cant find any of the supplied sample_ids")
      my.samples <- NULL}
  }
  
  # if fill manual or fill gradient
  if(!is.numeric(my_df$annotation)){
    message("... ... annotation not numeric. setting fill to discrete colors defined by 'color.key'")
    do.fill = "manual"
  }
  
  ## if to do gradient
  ## -------------
  if(is.numeric(my_df$annotation) && do.gradient && !pseudo.gradient){
    message("... ... annotation is numeric. setting fill to gradient as defined by 'color.key.gradient'")
    do.fill = "gradient"
    if(gradient.normalize) my_df$annotation <- my_df$annotation - mean(my_df$annotation, na.rm=T) / sd(my_df$annotation, na.rm=T)
    #if(cap.outliers)
    #color.key.gradient = palette_gradientRamp(x = my_df$annotation, my.pal = palette.gradient, palette.sat = gradient.sat, na.col = gadient.na.col)
    color.key.gradient = palette.gradient
    my.gradient.limits <- quantile(my_df$annotation, probs=gradient.squish.probs, na.rm=T)
  }
  
  
  ## Color.key - Use Pseudo gradient
  ## ----------------------
  # fix for problems with gradient_fill in plotly - define a pseudigradient instaed
  if(pseudo.gradient==T && do.fill == "gradient"){
    message(" ... setting upp pseudo-gradient - manual_fill")
    if(gradient.normalize) my_df$annotation <- my_df$annotation - mean(my_df$annotation, na.rm=T) / sd(my_df$annotation, na.rm=T)
    ## gradient colors is defined here by creating a manual palette for each row
    my_df$annotation_bin <- paste0("annotation",1:nrow(my_df))
    # color_key_pseudo_gradient
    my_df <- my_df %>% dplyr::rename(annotation_num = annotation, annotation = annotation_bin)
    color.key <- dlfoo2::heatcol_creator2(my_df$annotation_num, my.pal = rev(palette.gradient), palette.sat = gradient.sat)$col
    names(color.key) <- my_df$annotation
    do.fill = "manual"
  } # end pseudo.gradient
  
  
  ## if do.shape
  do.shape = F
  if(!is.null(shape.key)){
    if(any(levels(my_df$annotation) %in% names(shape.key))){
      message("... ... setting do.shape to TRUE")
      do.shape = T
      shape_key <- shapeKeyFix(shape.key,  as.character(levels(my_df$annotation)))
      shape_key <- shape_key[order(match(names(shape_key), levels(my_df$annotation)))]
    }}
  
  # zoom coordinates (?? outdated)
  do.zoom = FALSE
  if(!is.null(zoom.coord)){
    message("... ... 'zoom.chord' found. setting do.zoom to T")
    do.zoom=T
    if(class(zoom.coord)!="list") stop("zoom.coord must be list with named character slots x and y character slots of lengtht 2")
    if(length(zoom.coord)!=2) stop("zoom.coord must be list with named character slots x and y character slots of lengtht 2")
    if(!identical(sort(names(zoom.coord)), c("x","y"))) stop("zoom.coord must be list with named character slots x and y character slots of lengtht 2")
  }
  
  ## if capping x/y values
  ## ------------------
  #!!! note that x is dim2 and y dim1 for historic reasons
  if(!is.null(xlim.cap)){
    if(!is.numeric(xlim.cap) | length(xlim.cap)!=2) stop("xlim.cap is not nueric vector of 2 values")
    my_df$dim2[my_df$dim2 < min(xlim.cap)] <- min(xlim.cap)
    my_df$dim2[my_df$dim2 > max(xlim.cap)] <- max(xlim.cap)
    message("x-axis (dim2) is capped to: ", xlim.cap[1], " and ", xlim.cap[2])
  }
  if(!is.null(ylim.cap)){
    if(!is.numeric(ylim.cap) | length(ylim.cap)!=2) stop("xlim.cap is not nueric vector of 2 values")
    my_df$dim1[my_df$dim1 < min(ylim.cap)] <- min(ylim.cap)
    my_df$dim1[my_df$dim1 > max(ylim.cap)] <- max(ylim.cap)
    message("y-axis (dim1) is capped to: ", ylim.cap[1], " and ", ylim.cap[2])
  }
  
  
  
  ## color.key & shape keys (named vectors w colors)
  ## ---------
  if(!drop.na.annotation) color.key['_NA_'] <- NA
  color_key_x <-  colorKeyFix(color.key, as.character(levels(factor(my_df$annotation))))
  
  if(!is.null(violin.2nd.group.column) & plot_type=="violinPlot"){
    if(!c(violin.2nd.group.column %in% colnames(my_df))) stop(paste(violin.2nd.group.column, "not among pdata columns"))
    color_key_x <-  colorKeyFix(color.key, as.character(levels(factor(my_df[,violin.2nd.group.column]))))
  }
  
  ## if swap axes
  ## ------------
  if(swap.axes){
    my_df <- my_df %>%
      dplyr::mutate(dim2_temp=dim1) %>%
      dplyr::mutate(dim1_temp=dim2) %>%
      dplyr::select(-dim1, -dim2) %>%
      dplyr::rename(dim1=dim1_temp, dim2=dim2_temp) %>%
      dplyr::select(sample_id, dim1, dim2, everything())
    
    x.title_temp <- y.title
    y.title <- x.title
    x.title <- x.title_temp
  }
  
  
  ##    Plot ggplot
  ## :::::::::::::::::::
  # if(is.null(plot.title)) plot.title = pdata.column
  
  ## scatterPlot or rankPlot
  ##
  message("... plotting")
  if(plot_type %in% c("scatterPlot","rankPlot")){
    
    g <- ggplot(my_df) +
      aes(x=dim2, y=dim1)
    # text =  paste("Point: ", annotation, "\n", sample_id, "\n", annotation2) :: if plotly
    
    if(do.shape & do.fill == "manual"){
      g <- g +
        geom_point(aes(fill=annotation, shape=annotation), colour=my.stroke.color, stroke=my.stroke, size=my.size, alpha=my.alpha) +
        scale_fill_manual(values=color_key_x, na.value=NA) +
        scale_shape_manual(values = as.numeric(shape_key), na.value = 4) +
        guides(color = guide_legend(override.aes = list(size=7)))
    }
    
    if(!do.shape & do.fill == "manual"){
      g <- g +
        geom_point(aes(fill=annotation), shape=21, stroke=my.stroke, colour=my.stroke.color, size=my.size, alpha=my.alpha) +
        scale_fill_manual(values=color_key_x, na.value=NA)+
        guides(fill=guide_legend(title=pdata.column))
    }
    if(!do.shape & do.fill == "gradient"){
      g <- g +
        geom_point(aes(fill=annotation), shape=21, stroke=my.stroke, colour=my.stroke.color, size=my.size, alpha=my.alpha) +
        # scale_fill_gradientn(colours=color.key.gradient, na.value=NA) +
        scale_fill_gradientn(colours=color.key.gradient, na.value=NA, limits=my.gradient.limits, oob=scales::squish) +
        guides(fill=guide_legend(title=pdata.column))
    }
    
    ## add dim1 and dim2 lines for selected sample name(s)
    if(!is_null(my.samples)){
      # g <- g + geom_label_repel(aes(label = plot_text),  box.padding   = 0.35, point.padding = 0.5, segment.color = 'grey50')
      g <- g +
        geom_hline(yintercept = data.frame(my_df %>% dplyr::filter(!is.na(plot_text)))$dim1) +
        geom_vline(xintercept = data.frame(my_df %>% dplyr::filter(!is.na(plot_text)))$dim2)
    }
    
  } # end if scatter/rankplot
  
  ## Stripchart
  if(plot_type %in% c("stripChart")){
    
    my_sizes <- my_df %>% group_by(dim2) %>% dplyr::summarize(num=n())
    if(axis_names_n){ axis_names <- paste0(my_sizes$dim2, ", ", "n=", my_sizes$num)
    }else{
      axis_names <- paste0(my_sizes$dim2)
    }
    # if(stripchart.box){
    #    g <- ggplot(my_df) +
    #       message(" ... stripchart.box=T, adding box")
    #       aes(x=dim2, y=dim1) +
    #       geom_boxplot() +
    #       theme(legend.position='none') +
    #       geom_jitter(aes(fill=annotation), position=position_jitter(0.2), shape=21, stroke=my.stroke, colour=my.stroke.color, size=my.size, alpha=my.alpha) +
    #       scale_fill_manual(values = color_key_x) +
    #       scale_x_discrete( labels = axis_names)
    # }
    
    g <- ggplot(my_df) +
      aes(x=dim2, y=dim1, fill=annotation) +
      theme(legend.position='none') +
      geom_jitter(position=position_jitter(0.2), shape=21, stroke=my.stroke, colour=my.stroke.color, size=my.size, alpha=my.alpha) +
      scale_fill_manual(values = color_key_x) +
      scale_x_discrete( labels = axis_names)
    
  }
  
  
  # Vioin
  if(plot_type %in% c("violinPlot")){
    
    my_sizes <- my_df %>% group_by(annotation) %>% dplyr::summarize(num=n())
    if(axis_names_n){ axis_names <- paste0(my_sizes$annotation, ", ", "n=", my_sizes$num)
    }else{
      axis_names <- paste0(my_sizes$annotation)
    }
    
    # axis_num <- my_sizes$num
    # my_df %>% left_join(my_sizes, by="annotation") %>% mutate(myaxis = paste0(dim1, "\n", "n=", x.group.num))
    
    if(is.null(violin.2nd.group.column)){
      if(!do.boxplot){
        g <- ggplot(my_df) +
          aes(x=annotation, y=dim1, fill=annotation) +
          # text = paste("Group: ", annotation)
          theme(legend.position='none') +
          geom_violin(alpha=0.5, width=1, trim = TRUE, scale = "width", adjust = 0.5) +
          geom_boxplot(width=0.2, outlier.colour="grey25", notch = FALSE, notchwidth = .75, alpha = my.alpha, colour = "grey10", outlier.size=0.3) +
          scale_fill_manual(values = color_key_x) +
          scale_x_discrete( labels = axis_names)
      }else{
        g <- ggplot(my_df) +
          aes(x=annotation, y=dim1, fill=annotation) +
          # text = paste("Group: ", annotation)
          theme(legend.position='none') +
          geom_boxplot(outlier.colour="grey25", notch = FALSE, notchwidth = .75, alpha = my.alpha, colour = "grey10", outlier.size=0.3) +
          scale_fill_manual(values = color_key_x) +
          scale_x_discrete( labels = axis_names)
      }
    }
    if(!is.null(violin.2nd.group.column)){
      message(" ... ... viollin () plot with 2nd grouping variable")
      g <- ggplot(my_df) +
        aes(x=annotation, y=dim1, fill = UQ(rlang::sym(violin.2nd.group.column))) +
        theme(legend.position='none') +
        # geom_violin(alpha=0.5, width=1, trim = TRUE, scale = "width", adjust = 0.5) +
        geom_boxplot(outlier.colour="grey25", notch = FALSE, notchwidth = .75, alpha = my.alpha, colour = "grey10", outlier.size=0.3) +
        scale_fill_manual(values = color_key_x) +
        scale_x_discrete(labels = axis_names)
    }
  }
  
  # barplot
  if(plot_type %in% c("barPlot")){
    
    ftab <- my_df %>% group_by(dim1, annotation) %>% dplyr::summarise(n=n()) %>% mutate(f=n/sum(n))
    # ftab$annotation <- factor(ftab$annotation, levels=rev(ftab$annotation))
    # ftab <- ftab %>% arrange(!annotation)
    my_sizes <- my_df %>% group_by(dim1) %>% dplyr::summarize(num=n())
    axis_names <- paste0(my_sizes$dim1, ", ", "n=", my_sizes$num)
    
    if(barplot.y.freq == T ){
      # factor(annotation, levels = c(NA, as.character(levels(my_df$annotation))), exclude = NULL)
      
      g <- ggplot(ftab, aes(fill=forcats::fct_rev(annotation), x=dim1, y=n)) +
        geom_bar(position="fill", stat="identity") +
        scale_fill_manual(values =color_key_x) +
        scale_x_discrete(labels = axis_names)
    }
    
    if(barplot.y.freq == F ){
      g <- ggplot(ftab, aes(fill=forcats::fct_rev(annotation), x=dim1, y=n)) +
        geom_col() +
        scale_fill_manual(values =color_key_x) +
        scale_x_discrete(plabels = axis_names)
    }
    # text = paste("Group: ", dim1, "\n annotation: ", annotation)
    
    
  } # end barplot
  
  ## plot stuff
  g <- g + ggtitle(label = plot.title) +
    labs(x=x.title, y=y.title, colour='') +
    theme(axis.title.y = element_text(size=9)) +
    theme(axis.title.x = element_text(size=9)) +
    theme(legend.text = element_text(size=9)) +
    theme(title = element_text(size=9)) +
    theme(axis.text.x=element_text(angle = 45, hjust = 1))
  
  ## if rescaling of x/y
  if(!is.null(xlim.rescale)){
    if(!is.numeric(xlim.rescale[!is.na(xlim.rescale)]) | length(xlim.rescale)!=2) stop("xlim.rescale is not numeric vector of 2 values or NA adn one value ")
    g <- g + xlim(xlim.rescale)
  }
  if(!is.null(ylim.rescale)){
    if(!is.numeric(ylim.rescale[!is.na(ylim.rescale)]) | length(ylim.rescale)!=2) stop("ylim.rescale is not numeric vector of 2 values or NA adn one value ")
    g <- g + ylim(ylim.rescale)
  }
  
  # Guide
  if(showGuide==F){
    g <- g + guides(fill=FALSE, shape=FALSE, color=FALSE) +
      theme(legend.position='none')
  }
  
  if(coord.fixed) {g = g + coord_fixed(ratio = 1)}
  
  if(do.zoom){g = g + coord_cartesian(xlim = zoom.coord$x, ylim = zoom.coord$y)}
  
  if(zero.lines) g <- g + geom_hline(yintercept=0, linetype="dashed", color="black", size=0.5) + geom_vline(xintercept=0, linetype="dashed", color="black", size=0.5)
  message("... done")
  
  return(g)
  
  
} # end violinScatterRankBarPlotGeneric




#' @description Wrapper for ggplot heatmap using geom_raster
#' @description NOTE! cannot export thse as vetor graphics!! Need geom_rect??
#' @description Input (x) must be a matrix with n sample columns and m feature/gene rows
#' @description Note: since plotly cant cope with plotting pdata (character/factors) in geom_raster and custom colorings - groups have to be defined/marked creatively. here by adding vertical line using 'pdata.column.vert.lines'
#' @description also: plotly.js does not (yet) support horizontal legend items. You can track progress here: https://github.com/plotly/plotly.js/issues/53
#' @family Heatmap
#' @family PlotWrap
#' @param my_mat matrix for heatmap with n samples and n features
#' @param pdata if to add sample information to the ggplot object (used for downstream ggplotly tootips)
#' @param pdata.column.vert.lines A pdata column that define groups. if to draw vertical lines that separate groups. Will only work if x is sorted properly and if cluster.samples==F.
#' @param pdata.column.vert.lines.lwd
#' @param normalize.rows if to z normalize rows (divide by sd)
#' @param palette.gradient palette for gradient coloring
#' @param cluster.features if to cluster rows/features
#' @param cluster.samples if to cluster samples/columns
#' @param strip.background if to leave backround blank
#' @param col.labels if to draw ticks and labes for x axis
#' @param nrow.max.for.labels the max rows for x if to plot y-labels
#' @param gradient.squish.probs how palette is squished
#' @param gradient.limits if to manually set the gradient limits (not by squish values in matrix)
#' @param x.title
#' @param y.title
#' @param plot.title
#' @param showGuide defaults to F
#' @param ... parameters passed on to function violinScatterRankBarPlotGeneric
#' @return a ggplot object. Z vill represent the matix values
#' @export
heatmapGeomRaster <- function(
  my_mat =NULL,
  pdata = NULL,
  pdata.column.vert.lines = NULL,
  pdata.column.vert.lines.lwd = 0.5,
  #color.key.pdata = NULL, # not suppored since ggploty cant understand scale_fill_manual
  #pdata.heatmap = F, # not suppored since ggploty cant understand scale_fill_manual
  normalize.rows = F,
  cluster.features = T,
  cluster.samples = T,
  nrow.max.for.labels = 50,
  col.labels = F,
  strip.background = F,
  plot.title=NULL,
  x.title=NULL, y.title=NULL,
  palette.gradient = rev(c(rev(RColorBrewer::brewer.pal(9,"YlOrRd")[1:9]), rep("white",3), RColorBrewer::brewer.pal(9,"Blues")[1:9])),
  gradient.limits = NULL,
  gradient.squish.probs=c(0.025, 0.975),
  showGuide=F,  ...
){
  
  x <- my_mat
  gadient.na.col=NA
  if(!is.matrix(x)) stop("x must be amatrix")
  
  ## match up columns (sample ids) with pdata
  ## -------------------------------
  # pdata must include sample_id column that match colnames of x
  if(!is.null(pdata)){
    if(!"sample_id" %in% colnames(pdata)) stop("pdata must include 'sample_id' column with values matching colnames of x")
    pdata <- pdata %>% dplyr::filter(sample_id %in% colnames(x))
    pdata <- pdata[match(colnames(x), pdata$sample_id),]
    if(ncol(pdata)>10) message("WARNING: pdata contains >10 columns - are these really needed for downstream plotly?")
    #if(!is.null(pdata.column.plotly.text)){
    #   if(!pdata.column.plotly.text %in% colnames(pdata)) stop(paste0("pdata.column.plotly.text not in colnames of pdata:  ", pdata.column.plotly.text))
    # pdata$plotly_text <- as.character( pdata[,pdata.column.plotly.text])
    #}
  }
  
  ## if normalize.rows
  if(normalize.rows){
    ## remove genes that do not vary (with sd=0). Not compatible with hca clustering.
    ## -------------------------
    u <- apply(x, 1, sd)
    if(length(which(u==0))>0){
      message(" .. WARNING: genes with zero standard deviation found. Removing these")
      x <- x[!u==0,]
    }
    message("... normalizing rows")
    x <- t(apply(x, 1, function(y) y/sd(y)))
  }
  
  ## cluster rows and columns
  ## -----------------------
  if(cluster.features){
    ## remove genes that do not vary (with sd=0). Not compatible with hca clustering.
    ## -------------------------
    u <- apply(x, 1, sd)
    if(length(which(u==0))>0){
      message(" .. WARNING: genes with zero standard deviation found. Removing these")
      x <- x[!u==0,]
    }
    message("... clustering features")
    hca_rows <- dlfoo2::hca(t(x), plot.dendogram = F)
    x <- x[hca_rows$hclust$order, ]
  }
  if(cluster.samples){
    message("... clustering samples")
    hca_cols <- dlfoo2::hca(x, plot.dendogram = F)
    x <- x[, hca_cols$hclust$order]
    pdata <- pdata[hca_cols$hclust$order, ]
  }
  
  
  
  ## Matrix as heatmap
  ## :::::::::::::::::::::::::::::
  x_df <- x %>%
    as_tibble() %>%
    mutate(feature_id = rownames(x)) %>%
    gather(key="sample_id", value="Z", -feature_id) %>%
    mutate(Y=factor(feature_id, levels=unique(rownames(x)))) %>%
    mutate(X=factor(sample_id, levels=unique(colnames(x))))
  x_df$Y <- as.numeric(x_df$Y)
  x_df$X <- as.numeric(x_df$X)
  
  # if to add pdata group to plotly text
  # if(!is.null(pdata.column.plotly.text)){
  #    u <- match(x_df$sample_id, pdata$sample_id)
  #    x_df <- x_df %>% mutate(plotly_text = pdata[,pdata.column.plotly.text][u])
  # }
  
  ## add pdata to x_df
  ## -----------------
  if(!is.null(pdata)){
    x_df <- data.frame(x_df) %>% dplyr::left_join(pdata, by="sample_id")
    x_df <- as_tibble(x_df)
  }
  
  
  ## set up color gradient
  ## -------------
  message("... ... setting fill to gradient as defined by 'palette.gradient'")
  #if(gradient.normalize) my_df$annotation <- my_df$annotation - mean(my_df$annotation, na.rm=T) / sd(my_df$annotation, na.rm=T)
  #if(cap.outliers)
  #color.key.gradient = palette_gradientRamp(x = my_df$annotation, my.pal = palette.gradient, palette.sat = gradient.sat, na.col = gadient.na.col)
  color.key.gradient = palette.gradient
  my.gradient.limits <- quantile(x_df$Z, probs=gradient.squish.probs, na.rm=T)
  if(!is.null(gradient.limits)) my.gradient.limits = gradient.limits
  
  ## Axis labels
  ## -------
  y_labs <- unique(x_df$feature_id)
  y.n <- length(y_labs)
  
  
  ## ggPlot Heatmap
  ## ---------------
  message("... plotting")
  g_mat <-  ggplot(data.frame(x_df), aes(x=X, y=Y, fill= Z)) +
    geom_raster() +
    scale_fill_gradientn(colours=color.key.gradient, na.value=NA, limits=my.gradient.limits, oob=scales::squish)
  
  if(!col.labels){
    g_mat <- g_mat +
      theme(axis.ticks.x = element_blank(),
            axis.text.x = element_blank()
      )
  }
  
  
  if(strip.background){
    g_mat <- g_mat + theme_bw() +
      theme(panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            # axis.line = element_line(colour = "black")
            axis.line = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank()
      )
  }
  
  
  
  g_mat <- g_mat +
    ggtitle(label = plot.title) +
    labs(x=x.title, y=y.title, colour='black') +
    #theme(axis.title.y = element_text(size=9)) +
    #theme(axis.title.x = element_text(size=9)) +
    theme(legend.text = element_text(size=9)) +
    theme(title = element_text(size=9))
  # theme(axis.text.x=element_text(angle = 45, hjust = 1))
  # theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
  #theme(panel.background = element_rect(fill = NA, colour = NA)) +
  #theme(panel.grid.major =  element_line(colour = NA), panel.grid.minor =  element_line(colour = NA)) +
  #theme(strip.text = element_text(size=9)) +
  #labs(x="", y="" )
  
  if(nrow(x) < nrow.max.for.labels) g_mat <- g_mat + scale_y_continuous(limits=c(0.5, y.n+0.5), breaks=c(1:y.n), labels=as.character(y_labs), expand=c(0,0))
  
  ## if pdata.column.vert.lines, then identify groups, plot vert lines and add x-labels
  if(!is.null(pdata.column.vert.lines)){
    if(!pdata.column.vert.lines %in% colnames(pdata)) stop("'pdata.column.vert.lines' among pdata columns:  ", pdata.column.vert.lines)
    if(!is.factor(pdata[,pdata.column.vert.lines])) stop("'pdata.column.vert.lines' must be factor:  ", pdata.column.vert.lines)
    my_groups <- levels(pdata[,pdata.column.vert.lines])
    if(any(is.na(my_groups))) stop("'pdata.column.vert.lines' column must not contain NA values")
    my_x_lines <- unlist(lapply(my_groups, function(x){ which(pdata[,pdata.column.vert.lines]==x)[1]}))
    my_x_lines <- c(my_x_lines-0.5, nrow(pdata)+0.5)
    my_x_mid <- my_x_lines[-1] - c(my_x_lines[-1] - my_x_lines[-length(my_x_lines)])/2
    x.n <- length(my_x_mid)
    g_mat <- g_mat + geom_vline(xintercept = my_x_lines, lwd=pdata.column.vert.lines.lwd, col="black")
    g_mat <- g_mat + scale_x_continuous(breaks =my_x_mid, labels=as.character(my_groups), expand=c(0,0))
    
  }
  
  
  if(showGuide==F) g_mat <- g_mat +
    guides(fill=FALSE, shape=FALSE, color=FALSE) +
    theme(legend.position='none')
  
  #if(is.null(pdata.column.plotly.text)) g_mat <- g_mat + (aes(text=paste0(sample_id, "\n",feature_id,"\n", round(Z,2))))
  #if(!is.null(pdata.column.plotly.text)) g_mat <- g_mat + (aes(text=paste0(sample_id, "\n",feature_id,"\n", round(Z,2), "\n", plotly_text)))
  # ggplotly(g_mat)
  
  
  
  return(g_mat)
} # end function heatmapGeomRaster



#' Match up a color key to given a series of given unique annotation
#' Add colors if annots are missing
#' @family plot wraps
#' @family color
#' @param color.key
#' @param annots
#' @return color_key vector where all annotatoions have a color
#' @export
colorKeyFix <- function(color.key = NULL, annots =NULL, na.color = "light gray"){
  if(is.null(annots)) stop("character vector of all annotation values must be provided")
  if(is.null(color.key)){
    color_key <- rainbow(length(annots))
  }else{
    u <- match(annots, names(color.key))
    u.missing <- annots[is.na(u)]
    
    if(all(!is.na(u))){
      color_key <- color.key[u]
    }else{
      temp_key <- rainbow(length(u.missing))
      names(temp_key) <- u.missing
      color_key <- c(color.key[u[!is.na(u)]], temp_key)
    }
  }
  return(color_key)
  
} # end colorKeyFix


#' Match up a shape key to given a series of given unique annotation
#' Add shapes if annots are missing
#' @family plot wraps
#' @family color
#' @param shape.key
#' @param annots
#' @return shape_key vector where all annotatoions have a color
#' @export
shapeKeyFix <- function(shape.key = NULL, annots = NULL, add.shape = "21"){
  annots = sort(unique(as.character(annots)))
  if(is.null(annots)) stop("character vector of all annotation values must be provided")
  if(is.null(shape.key)){
    shape_key <- rep(add.shape, length(annots))
  }else{
    u <- match(annots, names(shape.key))
    u.missing <- annots[is.na(u)]
    
    if(all(!is.na(u))){
      shape_key <- shape.key[u]
    }else{
      temp_key <- rep(add.shape,length(u.missing))
      names(temp_key) <- u.missing
      shape_key <- c(shape.key[u[!is.na(u)]], temp_key)
    }
  }
  # if(re.sort) shape_key <- shape_key[order(names(shape_key))]
  
  return(shape_key)
} # end shapeKeyFix



#' function tsneGex wrapper for RtSNE
#' @description Version 1.0 20190418
#' @description tSNE Aanalysis using \code{\link[Rtsne]{Rtsne}} function in stats package.
#' @family tsne
#' @family analysis
#' @param eset expressionSet object
#' @param runName name of plot/analyses - subfolder
#' @param runFolder higher level directory where to plot/save objects
#' @param perplex vector of integers. Suggested max is n-1  3
#' @param doPlotTaxonomy if to load rcc phenotype from file & plot
#' @param theta Rtsne theta value
#' @param max_iter Rtsne max_iter value
#' @param floor_value Rtsne floor_value
#' @param z_score if to use z-scores (False default)
#' @param log if to log matrix(False default)
#' @param variance numeric between 0 and 1 describing the ratio of Features to keep (0.8 keeps the 80 percent highest vatying Features)
#' @return write tables and tSNEres to file
#' @export
tsneGex <- function(eset, perplex = seq(5, 50, by = 5), runName,
                    theta=0.5, max_iter=5000, floor_value=NULL, z_score=F,
                    log=F, variance=NULL){
  require(Rtsne)
  
  # data('pdata_utc_nt', package="dlfoo2data")
  
  if(class(eset)!="ExpressionSet") stop("eset not of class 'GenoSet'")
  x <- Biobase::exprs(eset)
  if(max(perplex) > ncol(x)/3-1) stop("max perplex should not exeed n/3-1")
  
  # myDir <- file.path(runFolder, runName)
  #if(dir.exists(myDir)) stop("Run already exists. Set alternative 'runName' or delete/move old run from disk")
  #dir.create(myDir)
  
  if(log){
    x<-log(x,2)
    x[which(x<0)]<-0
  }
  if(!is.null(floor_value)){
    x[which(x<floor_value)]<-floor_value
  }
  if(!is.null(variance)){
    yy<-apply(x,1,function(y)var(y))
    ord<-order(yy, decreasing = T)
    ord<-ord[c(1:(length(ord)*variance))]
    x<-x[ord,]
  }
  if(z_score){
    x<-t(apply(x,1,function(y)((y)-mean(y))/sd(y)))
    x[which(is.na(x))]<-0
  }
  
  x <- t(x)
  number <- dim(x)[1]
  #How many samples are in the matrix
  res <- vector(mode = "list", length=length(perplex))
  names(res) <- paste0("perplex_",perplex)
  
  
  for(i in 1:length(perplex)) {
    set.seed(42)
    #set the seed for every iteration
    res[[i]] <- Rtsne::Rtsne(as.matrix(x), theta=theta, max_iter=max_iter, perplexity=perplex[i], check_duplicates = F)
    #Perform Rtsne function
    rownames(res[[i]]$Y)<-rownames(x)
    #Add the samplename names
    
    #Header and result, cat sets the first col
    # cat("\n#H:perplexity_",perplex[i],"\t", file=outfile, append=TRUE)
    #write.table(res[[i]]$Y, file=outfile, col.names = F, quote=F, sep="\t", append=T)
    
  } # end i all peplexes
  
  return(res)
  # saveRDS(res, file=gsub(".txt",".rds",outfile))
  # write.table(data.frame(FeatureID=rownames(x)), file=gsub(".txt","_FeatureIDs.txt",outfile), row.names = F, quote=F)
  
} # end me tSNE



#' Wrapper to Create a Biobase ExpressionSet
#' Create a Biobase \code{\link[Biobase]{ExpressionSet}} from a data matrix, fdata and pddta objects
#' @family eset
#' @family gene expression
#' @param data, a matrix containg expression values
#' @param pdata, data frame with phenotype data (pdata), i.e. sample information.
#' @param fdata,  data frame with feature data (fdata), i.e. reporter/probe information
#' @param sample.names.col, index (number) of what column from pdata to use as sampleNames
#' @param feature.names.col, index (number) of what column from fdata to use as featureNames
#' @param platform, what platform is used - the annotation slot in the ExpressionSet object.
#' @param data_processing, additional information on data processing.
#' @return an \code{\link[Biobase]{ExpressionSet}} object
#' @export
esetCreate<-function(data, pdata, fdata, sample.names.col=1, feature.names.col=1, platform=c(), data_processing=NULL){
  requireNamespace("Biobase")
  
  rownames(pdata)<-as.character(pdata[,sample.names.col])
  p<-new("AnnotatedDataFrame", data=pdata, varMetadata=data.frame(row.names=names(pdata), labelDescription=rep(NA, length(names(pdata)))))
  
  rownames(fdata)<-as.character(fdata[,feature.names.col])
  
  f<-new("AnnotatedDataFrame", data=fdata, varMetadata=data.frame(row.names=names(fdata), labelDescription=rep(NA, length(names(fdata)))))
  
  rownames(data)<-as.character(fdata[,feature.names.col])
  colnames(data)<-as.character(pdata[,sample.names.col])
  
  dp<-new("MIAME")
  if(!is.null(data_processing)){
    preproc(dp)<-as.list(data_processing)
  }
  data <- as.matrix(data)
  my.eset<-new("ExpressionSet", phenoData=p, featureData=f, exprs=data, experimentData=dp)
  #annotation(my.eset)<-platform
  
  print(show(my.eset))
  return(my.eset)
} ## end create.eset


#' Center features/rows in an eset
#' Each row/feature in an eset will be re-rentered based on mean median or midpoint values. Can be performed on only a subset of samples.
#' @family eset
#' @family gene expression
#' @family transformation
#' @param eset, an object of class ExpressionSet
#' @param measure, can be 'mean','median', or 'max.min' (i.e., the midpoint between min and max values)
#' @param subset, specify sample names if measure should be calculated from only a subset of samples. Only these samples will get mean/median of 0.
#' @return an ExpressionSet where each row now is centered to the mean, median or the max.min
#' @export
esetCenterGenes <- function(eset, measure = c("mean", "median", "max.min"), subset = NULL) {
  requireNamespace("Biobase")
  measure <- match.arg(arg = measure, choices = c("mean", "median", "max.min"))
  cat("\n\n ... centering genes using ", measure, " method")
  if (class(eset) %in% c("LumiBatch", "ExpressionSet")) {
    if (is.null(subset)) {
      subset <- c(1:ncol(Biobase::exprs(eset)))
    }
    if (measure %in% c("mean", "median"))
      apa <- esApply(eset, 1, function(x, m, s) x - eval(parse(text = m))(x[s], na.rm = T), m = measure, s = subset)
    if (measure %in% c("max.min"))
      apa <- esApply(eset, 1, function(x, s) x - c(max(x[s], na.rm = T) + min(x[s], na.rm = T))/2, s = subset)
    Biobase::exprs(eset) <- t(apa)
    return(eset)
  } else {
    stop("object is not an eset!")
  }
  
}  # end center genes



#' Variance filter an eset
#' Returns a matrix with mean values for each sample group and for each row (gene). Colums will be named after the annotation provide.
#' @family eset
#' @family gene expression
#' @family filter
#' @param eset, an object of class ExpressionSet
#' @param method, sd (standard deviation), IQR, or
#' @param cutoff, specify the cutoff (sd, number of genes) to be applied
#' @return a filtered eset where only features that meet the variance filter criteria are kept
#' @export
esetVarianceFilter <- function(eset, method = "sd", cutoff = 1) {
  requireNamespace("Biobase")
  requireNamespace("BiocGenerics")
  requireNamespace("stats")
  
  stopifnot(class(eset) %in% c("LumiBatch", "ExpressionSet"))
  stopifnot(method %in% c("sd", "IQR", "sd.rank"))
  match.arg(method, choices = c("sd", "IQR", "sd.rank"))
  
  if(length(which(is.na(Biobase::exprs(eset))))>=1) message("\n  WARNING: NAs present in data matrix: ", length(which(is.na(Biobase::exprs(eset)))), " out of ", length(Biobase::exprs(eset))," values (",round(length(which(is.na(Biobase::exprs(eset))))/length(Biobase::exprs(eset))*100,2), "%)\n")
  
  if (method == "sd.rank") {
    o <- apply(Biobase::exprs(eset), 1, function(x) sd(x, na.rm = T))
    oo <- order(o, decreasing = T, na.last = F)
    eset <- eset[oo, ][1:cutoff]
    cat("\n ... ... Returning eset with top ", cutoff, " sd genes")
    return(eset)
  }
  
  cat("\n ...", method, " filter eSet using cutoff ", cutoff, sep = "")
  o <- apply(Biobase::exprs(eset), 1, function(x) do.call(what = method, args = list(x, na.rm=T))) > cutoff
  
  if (length(which(o == T)) == nrow(Biobase::exprs(eset))) {
    cat("\n ... ... All features passed threshold")
    return(eset)
  }
  if (length(which(o == F)) == nrow(Biobase::exprs(eset))) {
    cat("\n ... ... No features passed threshold")
    return(eset[o, ])
  }
  cat("\n ... ... ", length(which(o == T)), " features remaining", sep = "")
  return(eset[o, ])
}  # end esetVarianceFilter



#' produce a data frame from a tSNE results list
#' @description Version 1.0 20190213
#' @family tsne
#' @param tSNEres results list (tSNEres) from run.tsne-analysis
#' @return data frame
#' @export
#'
tsne2df <- function(tSNEres){
  x <- tSNEres
  data.frame(sample_id=rownames(x$Y), dim1=x$Y[,1], dim2=x$Y[,2], perplexity=x$perplexity, theta=x$theta, row.names=rownames(x$Y))
}

