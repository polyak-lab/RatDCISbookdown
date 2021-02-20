cnFreq_mod=function (x, CN_low_cutoff = 1.5, CN_high_cutoff = 2.5, plot_title = NULL, 
          CN_Loss_colour = "#002EB8", CN_Gain_colour = "#A30000", x_title_size = 12, 
          y_title_size = 12, facet_lab_size = 10, plotLayer = NULL, 
          plotType = "proportion", genome = "rn6", plotChr = NULL, 
          out = "plot") 
{
  
  multi_cytobandRet <- function(genome)
  {
    # Check if RMySQL is installed
    if(!requireNamespace("RMySQL", quietly=TRUE))
    {
      memo <- paste0("The package RMySQL is required for this functionality, 
                       please run the following: install.packages(\"RMySQL\")!")
      stop(memo)
    }
    # connect to UCSC mySQL genome database
    conn <- DBI::dbConnect(RMySQL::MySQL(), user="genome",
                           host="genome-mysql.cse.ucsc.edu",
                           dbname=genome)
    
    # query for chromosome positions and cytogenetic staining
    result <- DBI::dbGetQuery(conn=conn,
                              statement="SELECT chrom,
                              chromStart, chromEnd, name,
                              gieStain FROM cytoBand")
    DBI::dbDisconnect(conn)
    
    return(result)
  }
  
  multi_chrBound <- function(x)
  {
    # Check that input has size
    if(nrow(x) < 1)
    {
      memo <- paste0("input has 0 rows, it is possible that the UCSC",
                     " MySQL query has failed")
      stop(memo)
    }
    
    # Extract the columns needed
    data <- x[,c('chrom' ,'chromStart' , 'chromEnd')]
    
    # Obtain max for each chromosome
    maxChrom <- stats::aggregate(chromEnd ~ chrom, data=data, max)
    maxChrom <- cbind(maxChrom, maxChrom[,2])
    colnames(maxChrom) <- c('chromosome', 'start', 'end')
    
    # Obtain min for each chromosome
    minChrom <- stats::aggregate(chromStart ~ chrom, data=data, min)
    minChrom <- cbind(minChrom, minChrom[,2])
    colnames(minChrom) <- c('chromosome', 'start', 'end')
    
    # bind all the data together
    data <- rbind(maxChrom, minChrom)
    
    return(data)
  }
  
  cnFreq_qual <- function(x)
  {
    # Check that x is a data frame
    if(!is.data.frame(x)){
      memo <- paste0("Did not detect a data frame in argument supplied",
                     " to x... attempting to coerce")
      warning(memo)
      x <- as.data.frame(x)
      x <- droplevels(x)
    }
    
    # Check that x has at least 1 row
    if(nrow(x) < 1){
      memo <- paste0("x needs at least one row")
      stop(memo)
    }
    
    # remove any NA values in the data
    if(any(is.na(x))){
      na_rows_removed <- nrow(x) - nrow(na.omit(x))
      memo <- paste0("Removing ", na_rows_removed, " rows containing NA values")
      message(memo)
      x <- na.omit(x)
    }
    
    if(all(c('chromosome', 'start','end', 'segmean', 'sample') %in% colnames(x))){
      
      # make sure columns are of the correct type
      x$chromosome <- as.factor(x$chromosome)
      x$start <- as.integer(as.character(x$start))
      x$end <- as.integer(as.character(x$end))
      x$segmean <- as.numeric(as.character(x$segmean))
      x$sample <- as.factor(x$sample)
      
      # make sure windows are consistent if not disjoin them
      tmp <- split(x, x$sample)
      tmp_vec <- tmp[[1]]$end
      if(any(!unlist(sapply(tmp, function(x) x[,"end"] %in% tmp_vec), use.names=F))){
        memo <- paste0("Did not detect identical genomic segments for all samples",
                       " ...Performing disjoin operation")
        message(memo) 
        
        # here we split the DF up in an attempt to avoid complaints that lists are to large
        x <- split(x, f=x$chromosome)
        x <- lapply(x, cnFreq_disjoin)
        
        
        x <- do.call(rbind, x)
      }
      rm(tmp)
      rm(tmp_vec)
    } else {
      memo <- paste0("Did not detect correct columns in argument supplied",
                     " to x!")
      stop(memo)
    }
    
    # Check chromosome column in x
    if(!all(grepl("^chr", x$chromosome))){
      memo <- paste0("Did not detect the prefix \"chr\" in the chromosome",
                     " column of x... adding prefix")
      message(memo)
      x$chromosome <- paste0("chr", x$chromosome)
      x$chromosome <- as.factor(x$chromosome)
    } else if(all(grepl("^chr", x$chromosome))) {
      memo <- paste0("Detected \"chr\" in the chromosome column of x...",
                     " proceeding")
      message(memo)
    } else {
      memo <- paste0("Detected unknown or mixed prefixes in the chromosome ",
                     " column of x, should either have a chr prefix or ",
                     "none at all!")
      stop(memo)
    }
    
    return(x)
  }
  
  cnFreq_disjoin <- function(x){
    
    # create the Granges object for the data
    x <- GenomicRanges::GRanges(seqnames=x$chromosome,
                                ranges=IRanges::IRanges(start=x$start, end=x$end),
                                "sample"=x$sample, "segmean"=x$segmean)
    
    # disjoin with grange, get a mapping of meta columns and expand it
    disJoint_x <- GenomicRanges::disjoin(x, with.revmap=TRUE)
    revmap <- GenomicRanges::mcols(disJoint_x)$revmap
    disJoint_x <- rep(disJoint_x, lengths(revmap))
    
    
    # exract the meta columns and map them back to the disJoint GRanges object
    sample <- unlist(IRanges::extractList(GenomicRanges::mcols(x)$sample, revmap))
    segmean <- unlist(IRanges::extractList(GenomicRanges::mcols(x)$segmean, revmap))
    GenomicRanges::mcols(disJoint_x)$sample <- sample
    GenomicRanges::mcols(disJoint_x)$segmean <- segmean
    
    # convert the GRanges Object back to a data frame
    disJoint_x <- as.data.frame(disJoint_x)[,c("seqnames", "start", "end", "width",
                                               "sample", "segmean")]
    colnames(disJoint_x) <- c("chromosome", "start", "end", "width", "sample", "segmean")
    return(disJoint_x)
  }
  
  cnFreq_buildMain <- function(x, plotType, dummy_data, plot_title=NULL,
                               CN_low_colour='#002EB8', CN_high_colour='#A30000',
                               x_lab_size=12, y_lab_size=12, facet_lab_size=10,
                               plotLayer=NULL)
  {
    # Transform losses to be negative values for plotting purposes
    x$lossFrequency <- -1*x$lossFrequency
    x$lossProportion <- -1*x$lossProportion
    
    # Define parameters of plot
    theme <- theme(strip.text.x=element_text(size=facet_lab_size),
                   axis.text.x=element_blank(),
                   axis.ticks.x=element_blank(),
                   legend.position='right',
                   axis.title.x=element_text(size=x_lab_size, face='bold'),
                   axis.title.y=element_text(size=y_lab_size, face='bold'),
                   panel.grid.major.x=element_blank(),
                   panel.grid.minor.x=element_blank())
    facet <- facet_grid(. ~ chromosome, scales='free', space='free')
    xlabel <- xlab('Chromosomes')
    
    # Choose whether to plot aesthetics for proportion or frequency
    if(grepl("^PROP", plotType, ignore.case=TRUE)){
      ylabel <- ylab("Proportion of Copy Number Gains/Losses")
      ymax <- 1
      x$gain <- x$gainProportion
      x$loss <- x$lossProportion 
    } else if(grepl("^FREQ", plotType, ignore.case=TRUE)){
      ylabel <- ylab("Frequency of Copy Number Gains/Losses")
      ymax <- max(as.numeric(as.character(x$sampleFrequency)), na.rm=TRUE)
      x$gain <- x$gainFrequency
      x$loss <- x$lossFrequency 
    } else {
      memo <- paste0("did not recognize plotType ", plotType,
                     ", please specify one of \"proportion\" or \"frequency\"")
      stop(memo)
    }
    
    # Define the initial plot
    p1 <- ggplot(data=dummy_data,
                 mapping=aes_string(xmin='start',
                                    xmax='end',
                                    ymin=-1*ymax,
                                    ymax=ymax)) + geom_rect(alpha=0) +
      scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0))
    
    # add copy number data 
    p1 <- p1 + geom_rect(data=x, mapping=aes_string(xmin='start',
                                                    xmax='end',
                                                    ymin='loss',
                                                    ymax=0), fill=CN_low_colour)
    p1 <- p1 + geom_rect(data=x, mapping=aes_string(xmin='start',
                                                    xmax='end', 
                                                    ymin=0,
                                                    ymax='gain'), fill=CN_high_colour)
    
    p1 <- p1 + geom_hline(aes(yintercept=0), linetype="dotted")

    
    # build the plot
    p1 <- p1 + ylabel + xlabel + facet + theme_bw() + theme
    
    # if there are other layers, add them
    if(!is.null(plotLayer))
    {
      p1 <- p1 + plotLayer
    }
    
    # if title is supplied plot it
    if(!is.null(plot_title))
    {
      p1 <- p1 + ggtitle(plot_title)
    }
    
    return(p1)
  }
  
  multi_selectOut <- function(data, plot, out="plot", draw="FALSE")
  {
    # Decide what to output
    if(toupper(out) == "DATA")
    {
      return(data)
    } else if(toupper(out) == "PLOT" & isTRUE(draw)) {
      return(grid::grid.draw(plot))
    } else if(toupper(out) == "PLOT" & !isTRUE(draw)) {
      return(plot)
    } else if(toupper(out) == "GROB" & isTRUE(draw)) {
      return(plot)
    } else if(toupper(out) == "GROB" & !isTRUE(draw)) {
      return(ggplot2::ggplotGrob(plot))
    } else {
      warning("Did not recognize input to out...")
      if(isTRUE(draw))
      {
        return(grid::grid.draw(plot))
      } else {
        return(plot)
      }
    }
  }
 
  x <- cnFreq_qual(x)
  samples <- unique(x$sample)
  gainFreq <- function(x) {
    length(x[x >= CN_high_cutoff])
  }
  gainFrequency <- aggregate(segmean ~ chromosome + start + 
                               end, data = x, gainFreq)$segmean
  lossFreq <- function(x) {
    length(x[x <= CN_low_cutoff])
  }
  lossFrequency <- aggregate(segmean ~ chromosome + start + 
                               end, data = x, lossFreq)$segmean
  x <- aggregate(segmean ~ chromosome + start + end, data = x, 
                 length)
  colnames(x)[which(colnames(x) %in% "segmean")] <- "sampleFrequency"
  x$gainFrequency <- gainFrequency
  x$lossFrequency <- lossFrequency
  if (max(x$sampleFrequency) > length(samples)) {
    memo <- paste0("Detected additional sample rows after disjoin operation", 
                   " typically this indicates coordinates are 0-based, please convert", 
                   " coordinates to 1-base for accurate results")
    warning(memo)
  }
  x$gainProportion <- x$gainFrequency/length(samples)
  x$lossProportion <- x$lossFrequency/length(samples)
  preloaded <- c("hg38", "hg19", "mm10", "mm9", "rn5")
  if (any(genome == preloaded)) {
    message("genome specified is preloaded, retrieving data...")
    UCSC_Chr_pos <- GenVisR::cytoGeno[GenVisR::cytoGeno$genome == 
                                        genome, ]
    UCSC_Chr_pos <- multi_chrBound(UCSC_Chr_pos)
  }
  else {
    cyto_data <- genome
    UCSC_Chr_pos <- multi_chrBound(cyto_data)
  }

  if (nrow(UCSC_Chr_pos) < 1) {
    memo <- paste0("did not recognize genome ", genome, ", plotting provided data and ignoring chromosome ", 
                   "boundaries! Output could be decieving!")
    warning(memo)
  }
  dummy_data <- lapply(unique(x$sample), function(sample, chr_pos) cbind(chr_pos, 
                                                                         sample), UCSC_Chr_pos)
  dummy_data <- do.call("rbind", dummy_data)
  chr_order <- gtools::mixedsort(unique(dummy_data$chromosome))
  dummy_data$chromosome <- factor(dummy_data$chromosome, levels = chr_order)
  if (!is.null(plotChr)) {
    if (any(!plotChr %in% dummy_data$chromosome)) {
      missingChr <- plotChr[!plotChr %in% dummy_data$chromosome]
      plotChr <- plotChr[!plotChr %in% missingChr]
      memo <- paste0("The following chromosomes: ", toString(missingChr), 
                     ", could not be found! Valid chromosomes are: ", 
                     toString(unique(dummy_data$chromosome)))
      warning(memo)
    }
    dummy_data <- dummy_data[dummy_data$chromosome %in% plotChr, 
                             ]
    dummy_data$chromosome <- factor(dummy_data$chromosome, 
                                    levels = plotChr)
    x <- x[x$chromosome %in% plotChr, ]
    x$chromosome <- factor(x$chromosome, levels = plotChr)
  }
  p1 <- cnFreq_buildMain(x, plotType, dummy_data = dummy_data, 
                         plot_title = plot_title, CN_low_colour = CN_Loss_colour, 
                         CN_high_colour = CN_Gain_colour, x_lab_size = x_title_size, 
                         y_lab_size = y_title_size, facet_lab_size = facet_lab_size, 
                         plotLayer = plotLayer)
  output <- multi_selectOut(data = list(data = x), plot = p1, 
                            out = out)
  return(list(plot=output, data=x, dummy=dummy_data))
}
    

##cn freq main function


cnFreq_buildMain <- function(x, plotType, dummy_data, plot_title=NULL,
                             CN_low_colour='#002EB8', CN_high_colour='#A30000',
                             x_lab_size=12, y_lab_size=12, facet_lab_size=10,
                             plotLayer=NULL, text)
{
  # Transform losses to be negative values for plotting purposes
  x$lossFrequency <- -1*x$lossFrequency
  x$lossProportion <- -1*x$lossProportion
  
  # Define parameters of plot
  theme <- theme(strip.text.x=element_text(size=facet_lab_size),
                 axis.text.x=element_blank(),
                 axis.ticks.x=element_blank(),
                 legend.position='right',
                 axis.title.x=element_text(size=x_lab_size, face='bold'),
                 axis.title.y=element_text(size=y_lab_size, face='bold'),
                 panel.grid.major.x=element_blank(),
                 panel.grid.minor.x=element_blank())
  facet <- facet_grid(. ~ chromosome, scales='free', space='free')
  xlabel <- xlab('Chromosomes')
  
  # Choose whether to plot aesthetics for proportion or frequency
  if(grepl("^PROP", plotType, ignore.case=TRUE)){
    ylabel <- ylab("Proportion of Copy Number Gains/Losses")
    ymax <- 1
    x$gain <- x$gainProportion
    x$loss <- x$lossProportion 
  } else if(grepl("^FREQ", plotType, ignore.case=TRUE)){
    ylabel <- ylab("Frequency of Copy Number Gains/Losses")
    ymax <- max(as.numeric(as.character(x$sampleFrequency)), na.rm=TRUE)
    x$gain <- x$gainFrequency
    x$loss <- x$lossFrequency 
  } else {
    memo <- paste0("did not recognize plotType ", plotType,
                   ", please specify one of \"proportion\" or \"frequency\"")
    stop(memo)
  }
  
  # Define the initial plot
  p1 <- ggplot(data=dummy_data,
               mapping=aes_string(xmin='start',
                                  xmax='end',
                                  ymin=-1*ymax,
                                  ymax=ymax)) + geom_rect(alpha=0) +
    scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0))
  
  # add copy number data 
  p1 <- p1 + geom_rect(data=x, mapping=aes_string(xmin='start',
                                                  xmax='end',
                                                  ymin='loss',
                                                  ymax=0), fill=CN_low_colour)
  p1 <- p1 + geom_rect(data=x, mapping=aes_string(xmin='start',
                                                  xmax='end', 
                                                  ymin=0,
                                                  ymax='gain'), fill=CN_high_colour)
  
  p1 <- p1 + geom_hline(aes(yintercept=0), linetype="dotted") #+geom_label_repel()
  
  # build the plot
  p1 <- p1 + ylabel + xlabel + facet + theme_bw() + theme
  
  # if there are other layers, add them
  if(!is.null(plotLayer))
  {
    p1 <- p1 + plotLayer
  }
  
  # if title is supplied plot it
  if(!is.null(plot_title))
  {
    p1 <- p1 + ggtitle(plot_title)
  }
  
  return(p1)
}