
###############################################################################
# read the names and age estimates for each of the fossil analyses
###############################################################################
loadFossilFiles <- function(dir, regxp, fossilTable, minLogRows) {

  require(coda)

  height <- list()
  names <- character(0)
  for (i in 1:nrow(fossilTable)) {

    folder <- dir
    filename <- list.files(path=folder,pattern=paste(fossilTable$File.name[i],regexp,sep=""))

    if (length(filename) == 1) {
   
      name <- as.character(fossilTable$Species.name[i])

      logfile <- read.table(paste(folder,filename[1],sep="/"), sep="\t", header=T)

      if (nrow(logfile) >= minLogRows) {

        logfile$X <- NULL
        lastState <- logfile$Sample[nrow(logfile)]    
        thinState <- lastState - logfile$Sample[nrow(logfile)-1]
   
        logfile <- mcmc(logfile[,2:ncol(logfile)],start=0, end=lastState, thin=thinState)

        varName <- paste("height.",name,".",sep="")
        varName <- gsub(" ","_",varName)

        height[[i]] <- logfile[,(varnames(logfile) == varName ), drop=FALSE]
        names[i] <- name
        varnames(height[[i]]) <- name        
      } else {
        print(paste("Didn't process",filename[1]," because expected at least ", minLogRows, " rows but found ",nrow(logfile), "rows."))
      }
    } else {
      print(paste("Didn't process",fossilTable$File.name[i]," because found ",length(filename), "file matching name."))
    }
  }
  list(heights=height, names=names)
}

###############################################################################
# Construct and return main summary data frame of statistics
###############################################################################
createSummaryTable <- function(names, height, fossilTable, sequences, priorRange=160, hpdProb=0.95) {

  require(coda)

  est <- sapply(height,median)
  geo <- rowMeans(fossilTable[,c("Lower","Upper")])
  err <- abs(est - geo)
  rel_err <- err/geo

  hpds <- sapply(height, function(x) {HPDinterval(x,prob=hpdProb)})
  hpd_lower <- hpds[1,]
  hpd_upper <- hpds[2,]
  ess <- sapply(height, effectiveSize)

  post <- numeric(0)
  for (i in 1:length(height)) {
    post[i] <- length(height[[i]][height[[i]] < fossilTable$Upper[i] & height[[i]] > fossilTable$Lower[i]])/length(height[[i]])    
  }
  
  prior <- numeric(0)
  for (i in 1:length(height)) {
    prior[i] <- (fossilTable$Upper[i] - fossilTable$Lower[i]) / priorRange
  }
  
  probInHPD <- numeric(0)
  for (i in 1:length(height)) {
    if (fossilTable$Upper[i] == fossilTable$Lower[i]) {
    	if (geo[i] > hpd_lower[i] & geo[i] <= hpd_upper[i]) {
    		probInHPD[i] = 1.0;
    	} else {
    		probInHPD[i] = 0.0;
    	}
    } else {
      probInHPD[i] <- max(0,(min(hpd_upper[i], fossilTable$Upper[i]) - max(hpd_lower[i], fossilTable$Lower[i])) / (fossilTable$Upper[i] - fossilTable$Lower[i]))
    }
  }
  
  bf <- (post / prior) / ((1.0-post) / (1.0 - prior))

  df <- data.frame(names, post, probInHPD, prior, bf, err, rel_err, est, hpd_lower, hpd_upper, ess, geo)
  df$variance <- sapply(height,var)
  df$precision <- 1/df$variance

  seqclean <- gsub("\\{[0-9]+\\}","N", seq$Sequence)
  s2 <- gsub("-","",seqclean)
  sequences$charCount <- nchar(s2)
  df$characters <- sequences$charCount

  df
}

writeSummaryTable <- function(filename, caption, label, df) {
  require(xtable)
  
  dfshort <- df[,c("post", "bf", "est", "hpd_lower", "hpd_upper","err", "ess")]
  colnames(dfshort) <- c("post", "BF", "phylo age", "lower", "upper", "error", "ESS")
  
  xtab <- xtable(dfshort, 
    caption=caption, label=label, digits=c(2,2,1,1,1,1,1,0))
  print(file=filename, xtab)
}

###############################################################################
# produce the phylo age versus geo age plot and write to pdf file
###############################################################################
createPhyloAgeVsGeoAgePDF <- function(filename, df, fossilTable, labelled, labelPos, min=0, max=65) {

  require(plotrix)

  pdf(file=filename, width=5, height=5, pointsize=10)

  par(mar = c(4,4,1,1) + 0.1)
  plot(c(min,max*0.925),c(min,max*0.925), type="n", xlab="paleontological age (Myr)", ylab="Bayesian phylogenetic age estimate (Myr)", ylim=c(min,max), xlim=c(min,max))

  draw.ellipse(df$geo, df$est, df$geo-fossilTable$Lower, df$est-df$hpd_lower, angle = 0, segment = c(180,360), border =rgb(0.5,0.5,0.5,0.2), col=rgb(0.5,0.5,0.5,0.1), arc.only=FALSE)

  draw.ellipse(df$geo, df$est, df$geo-fossilTable$Lower, df$hpd_upper-df$est, angle = 0, segment = c(0,180), border =rgb(0.5,0.5,0.5,0.2), col=rgb(0.5,0.5,0.5,0.1), arc.only=FALSE)

  for (i in 1: nrow(df)) {
    lines(c(df$geo[i],df$geo[i]), c(df$hpd_lower[i],df$hpd_upper[i]), col="black")
    lines(c(fossilTable$Lower[i],fossilTable$Upper[i]), c(df$est[i],df$est[i]), col="black")
  }

  lines(c(min,max),c(min,max),col="blue")

  points(df$geo, df$est,pch=21, col="white", bg="black")

  text(df$geo[labelled], df$est[labelled], labels=df$names[labelled],pos=labelPos, offset=0.4, cex=0.7)

  lm <- lm(df$est~df$geo)
  
  r2 <- format(summary(lm)$r.squared,digits=3)
  eq <- bquote(bold(R^2 == .(r2)))
    
  text(max*0.05,max*0.925, labels=eq,pos=4)

  dev.off()
}

###############################################################################
# produce an age versus age plot and write to pdf file
###############################################################################
createPhyloAgeVsPhyloAgePDF <- function(filename, df1, df2, name1, name2, labelled, labelPos, min=0, max=65) {

  require(plotrix)

  if (!is.null(filename)) {
    pdf(file=filename, width=5, height=5, pointsize=10)
  }

  par(mar = c(4,4,1,1) + 0.1)
  plot(c(min,max*0.925),c(min,max*0.925), type="n", xlab=paste(name1, "age estimate (Myr)"), ylab=paste(name2,"age estimate (Myr)"), ylim=c(min,max), xlim=c(min,max))

  draw.ellipse(df1$est, df2$est, df1$hpd_upper-df1$est, df2$hpd_upper-df2$est, angle = 0, segment =  c(0,90), border =rgb(0.5,0.5,0.5,0.0), col=rgb(0.5,0.5,0.5,0.1), arc.only=FALSE)
  draw.ellipse(df1$est, df2$est, df1$est-df1$hpd_lower, df2$hpd_upper-df2$est, angle = 0, segment = c(90,180), border =rgb(0.5,0.5,0.5,0.0), col=rgb(0.5,0.5,0.5,0.1), arc.only=FALSE)
  draw.ellipse(df1$est, df2$est, df1$est-df1$hpd_lower, df2$est-df2$hpd_lower, angle = 0, segment = c(180,270), border =rgb(0.5,0.5,0.5,0.0), col=rgb(0.5,0.5,0.5,0.1), arc.only=FALSE)
  draw.ellipse(df1$est, df2$est, df1$hpd_upper-df1$est, df2$est-df2$hpd_lower, angle = 0, segment = c(270,360), border =rgb(0.5,0.5,0.5,0.0), col=rgb(0.5,0.5,0.5,0.1), arc.only=FALSE)

  for (i in 1: nrow(df1)) {
    lines(c(df1$est[i],df1$est[i]), c(df2$hpd_lower[i],df2$hpd_upper[i]), col="black")
    lines(c(df1$hpd_lower[i],df1$hpd_upper[i]), c(df2$est[i],df2$est[i]), col="black")
  }

  lines(c(min,max),c(min,max),col="blue")

  points(df1$est, df2$est,pch=21, col="white", bg="black")

  if (length(labelled) > 0) {
    text(df1$est[labelled], df2$est[labelled], labels=df1$names[labelled],pos=labelPos, offset=0.4, cex=0.7)
  }

  lm <- lm(df2$est~df1$est)
    
  r2 <- format(summary(lm)$r.squared,digits=3)
  eq <- bquote(bold(R^2 == .(r2)))
    
  text(max*0.05,max*0.925, labels=eq,pos=4)
  
  if (!is.null(filename)) {
    dev.off()
  }
}


####################################################################
# produce the density plots and write to pdf files
####################################################################
createDensityPlotFigures <- function (filePrefix, names, height, fossilTable, combinedPlots) {
  xmax <- 65

  if (combinedPlots) {
  
    pdf(file=paste(filePrefix, "_younger.pdf", sep=""), width=8, height=11)
    par(mfrow = c(6,3),
        oma = c(4,4,0,0) + 0.1,
        mar = c(1,1,2,1) + 0.1,
        mgp=c(2.5,0.65, 0))
        plot1 <- dev.cur()
    pdf(file=paste(filePrefix, "_older.pdf", sep=""), width=8, height=11)
     par(mfrow = c(6,3),
        oma = c(4,4,0,0) + 0.1,
        mar = c(1,1,2,1) + 0.1,
        mgp=c(2.5,0.65, 0))
        plot2 <- dev.cur()
  }

  for (i in 1:length(height)) {

    if (!combinedPlots) {
    pdf(file=paste(fossilTable$File.name[i],".pdf",sep=""))
    }
 
    lower <- fossilTable$Lower[i] 
    upper <- fossilTable$Upper[i]
    density <- 1.0/(upper-lower)

    x <- c(lower, lower, upper, upper)
    y <- c(0.0, density, density, 0.0)

    if (combinedPlots) {
    if (upper < 32) {
      dev.set(plot1)
      xmax <- 45
    } else {
      dev.set(plot2)
      xmax <- 65
    }
    }
    densplot(height[[i]], xlim=c(0,xmax),xlab="",main="")
    title(main = list(names[i], cex = 1.0))
    lines(x,y, col="red")

    if (!combinedPlots) {
    dev.off()
    }
  }            
   
  if (combinedPlots) {
    dev.set(plot1)
    dev.off()
    dev.set(plot2)
    dev.off()
  }
}

####################################################################
# produce the hist plots and write to pdf files
####################################################################
createHistPlotFigures <- function (filePrefix, names, height, fossilTable, combinedPlots, xmin=0, xmax=66,binsPerMY=1) {

  if (combinedPlots) {
  
    pdf(file=paste(filePrefix, ".pdf", sep=""), width=8, height=11.5)
    par(mfrow = c(7,3),
        oma = c(4,4,0,0) + 0.1,
        mar = c(1,1,2,1) + 0.1,
        mgp=c(2.5,0.65, 0),lwd=0.5)
        plot1 <- dev.cur()
  }

  for (i in 1:length(height)) {

    if (!combinedPlots) {
      pdf(file=paste(fossilTable$File.name[i],".pdf",sep=""))
    }
    
    hist <- hist(height[[i]], xlim=c(xmin,xmax),breaks=c((0:(xmax*binsPerMY)/binsPerMY),10*xmax),freq=FALSE, xlab="",main="",col="dark gray", border="white")
    title(main = list(names[i], cex = 1.0))
    #rect(lower, 0.0, upper, density, col=rgb(1.0,0.0,0.0,0.5))
    
    lower <- fossilTable$Lower[i] 
    upper <- fossilTable$Upper[i]
    if (lower == upper) {
      density = max(hist$density)
    } else {
      density <- 1.0/(upper-lower)
    }
    x <- c(lower, lower, upper, upper)
    y <- c(0.0, density, density, 0.0)
    
    
    lines(x,y, col="red",lwd=1)

    if (!combinedPlots) {
      dev.off()
    }
  }            
   
  if (combinedPlots) {
    dev.set(plot1)
    dev.off()
  }
}

####################################################################
# produce the precision versus number of known characters plot
####################################################################
createPrecisionVsKnownCharactersFigure <- function(filename, df, labelled, labelPos) {

  pdf(file=filename, width=5, height=5, pointsize=10)
  plot(df$characters,df$precision, pch=16, xlab="Number of known characters", ylab="Precision of phylogenetic age estimate (1/variance)", xlim=c(0,max(df$characters)*1.05))
  
  text(df$characters[labelled], df$precision[labelled], labels=df$names[labelled],pos=labelPos,offset=0.45, cex=0.7)

  lm <- lm(df$precision ~ df$characters)
  abline(lm, col="red")
      
  r2 <- format(summary(lm)$r.squared,digits=3)
  eq <- bquote(bold(R^2 == .(r2)))
    
  text(max(df$characters)*0.05,max(df$precision)*0.925, labels=eq,pos=4)
  
  dev.off()
}
