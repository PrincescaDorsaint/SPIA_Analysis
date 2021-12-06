#!/usr/bin/env Rscript
library("optparse")

option_list = list(
  make_option(c("-s", "--input"), type="character", default=NULL, 
              help="SPIA result csv", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="Output PDF file", metavar="character")
); 
  
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

resultTable <- read.csv(opt$input)
pdf(opt$output)

SPIAPlot <- function(SPIAResult, titleplot = NULL, colors = NULL, xlab=NULL, ...){  
  
  #if title is not defined
  if (is.null(titleplot)){
    titleplot<-paste("Pair-wise comparison of ", SPIAResult$N_samples," samples on ", SPIAResult$N_SNPs, " SNPs",sep="")  
  }
  
  if (is.null(colors)){
    colors <- c(rgb(54,144,192,maxColorValue=255),  #uncertain
                rgb(65,174,118,maxColorValue=255),  #similar
                rgb(215,48,31,maxColorValue=255),  #different
                rgb(107,107,107,maxColorValue=255),  #no computed
                rgb(37,37,37,maxColorValue=255)  #no computed
                
    )    
  }
  
  if(is.null(xlab)){
    xlab="Index of Sample Pairs"
  }
  
  baseCol <- rgb(37,37,37,maxColorValue=255)
  
  plot(NA,
       xlim=c(0,dim(SPIAResult)[1]),
       ylim=c(0,1),
       axes=F,xlab="",ylab="",
       ...)
  
  axis(1,col=baseCol,col.axis=baseCol)
  axis(2,col=baseCol,col.axis=baseCol)
  title(main=titleplot[1],ylab="Distance D (% of discordant genotype calls)",xlab=xlab,col.main=baseCol,col.lab=baseCol)
  
  #check if the SPIAresult contain the probabilistic test
  if (SPIAResult$testDone)
    #compute the color depending on the result of the statistical test
  {      
    for(i in c(1:dim(SPIAResult)[1])){
      pointConf <- switch(SPIAResult[i,4], 
                          Uncertain = c(colors[1],19), 
                          Similar = c(colors[2],19), 
                          Different = c(colors[3],20), 
                          LimitsError= c(colors[4],19), 
                          c(colors[5],19))    
      points(i,SPIAResult[i,3],col=pointConf[1],pch=as.integer(pointConf[2]),cex=0.7);     
    }
    #plot the legend that describe the color code
    legendText <- c("SPIA TEST: uncertain", "SPIA TEST: match", "SPIA TEST: different","SPIA TEST: test not well defined", paste("< ",SPIAResult$PercValidCall[1],"% of available calls",sep=""))
    legend(0,1,legendText, col = colors, pch=c(19,19,20,19,19),cex=0.6)                            
  } else {
    #use the same color for each point
    for(i in c(1:dim(SPIAResult)[1])){      
      points(i,SPIAResult[i,3],col="black",pch=20,cex=0.7);     
    }
  }
  
}

SPIAPlot(resultTable)
dev.off()
  
