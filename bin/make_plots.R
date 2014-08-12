#-- make_plots.R ------------------------------------------------------------------#

run_task <- function(task, ...){
    switch(task,
           scr=scr(...),
           kfr=kfr(...),
           unknown_task(task)
           )
}	 


unknown_task <- function(task){
    write(paste("Unknown task:", task), stderr()); 
    quit(status = 1);
}

#-- task scr --#

scr <-function(reads.hash, seeds.hash, do.plot=TRUE){
    
    # from seed reads
    write("Reading seed kmers", stderr());
    seeds <- kmerPeaks(read.table(pipe(paste('jellyfish histo', seeds.hash, sep=" "))), freq.min=0);

    # from entire input
    write("Reading total kmers", stderr());
    raw <- kmerPeaks(read.table(pipe(paste('jellyfish histo', reads.hash, sep=" "))), freq.min=0);

    # TODO: warn/die if no peaks
    
    # get raw peak closest to seed peak
    delta <- abs(raw$peaks$covs - seeds$peaks$covs[1]);
    # get closest peak(s) - minimum distance, two in case of 2 peaks in of equal dist
    p1.i <- which(delta == sort(unique(delta))[1]);
    if(length(p1.i) < 2){ # get second-closest peak(s)
        p2.i <- which(delta == sort(unique(delta))[2])
    }else{
        p2.i <- p1.i[2];
        p1.i <- p1.i[1];
    }

    double.peak <- FALSE;
    # check if we are dealing with double peak
    for (i in p2.i){
        ratio = raw$peaks$covs[p1.i]/raw$peaks$covs[i];
        if(ratio > 1.8 && ratio < 2.2){ # double peak
            double.peak = TRUE;
            p2.i <- p1.i;
            p1.i <- i;
        }else if(ratio > 0.4 && ratio < 0.6){
            double.peak = TRUE;
            p2.i <- i;
        }
    }
        
    p1 <- raw$peaks[p1.i,];
    p2 <- NULL;
    if(double.peak){
        p2 <- raw$peaks[p2.i,];
        pt.size <- round(p1$sizes/p1$covs + p2$sizes/p2$covs)
    }else{
        pt.size <- round(p1$sizes/p1$covs)
    }
    
    # TODO: warn/die if pt.size > 400kb

    write(paste("coverage", round(p1$cov), sep="\t"), stdout());
    write(paste("size", pt.size, sep="\t"), stdout());
    
    if(do.plot){
        
        pdf("scr-seeds.pdf", width=10, height=6);

        write("Plotting seed and total kmers", stderr());

        # get some plotting limits

        xmax <- max(seeds$peaks$covs[1])^1.5
        xmin <- max(seeds$peaks$covs[1])^(1/2)
        ymax <- max(seeds$peaks$freqs[1])*2.5

        
        barplot(
            seeds$data$freqs,
            diff(seeds$data$covs),
            space=0,
            main="kmer-coverage of plastid seed reads",
            xlab="coverage",
            ylab="frequency",
            xlim=c(xmin,xmax),
            ylim=c(0,ymax),
            log="x",
            col=cl[2],
            axes=FALSE
        );

        # scale raw kmer frequency of seed peak surrounding area to fit the plot area
        raw.i <- which.min(abs(raw$data$covs - seeds$peaks$covs[1]))
        raw.sf <- seeds$peaks$freqs[1] / log(max(raw$data$frels[(raw.i-5):(raw.i+5)]))
        raw$peaks.scaled <- raw$peaks
        raw$peaks.scaled$frels <- log(raw$peaks$frels)*raw.sf

        axis(2)
        axis(1, pos=-3);
        lines(raw$data$covs,log(raw$data$frels)*raw.sf, col=cl[4], lwd=3);
        abline(v=seeds$peaks$covs[1], lwd=3, col=cl[5]);
        abline(v=p1$covs[1], lwd=3, col=cl[1]);
        if(! is.null(p2)){
            abline(v=p2$covs[1], lwd=3, lty=2, col=cl[1]);
        }
        peakSizes.add(raw$peaks.scaled, rel=T);        

        legend(
            "topright",
            c("plastid seed kmers", "plastid seed coverage", "total data kmers (30%)", "total plastid coverage"),
            lwd=3,
            lty=c(1,1,1,1),
            seg.len=2, 
            col=cl[c(2,5,4,1)]
        );

        msg <- dev.off()

    }
}


#-- task kfr --#

kfr <- function(coverage=NULL){

    pdf("kfr.pdf", width=10, height=6);
    coverage <- as.integer(coverage)
    # kmer filter

    write("Reading kmers", stderr());
    ## from hist
    # scr <- kmerPeaks(read.table("Dm_GenIl_019-035-PE.trimmed.ada-qual-w10-ml10-histo.tsv", header=F));
    # kfr1 <- kfr2 <- scr;
    scr  <- kmerPeaks(read.table(pipe('jellyfish histo -l 10 scr.jf'), header=F));
    kfr1 <- kmerPeaks(read.table(pipe('jellyfish histo -l 10 kfr1.jf'), header=F));
    kfr2 <- kmerPeaks(read.table(pipe('jellyfish histo -l 10 kfr2.jf'), header=F));
    
    #print(kfr2)
    
    write("Plotting filtered kmers", stderr());
    plot(kfr2$data, 
         type="n", 
         main="kmer-coverage of subsetted and filtered data sets",
         #log="xy",
         xlim=c(0,scr$peak.max.cov*3),
         #xlim=c(0, 3000),
         ylim=c(1,scr$peak.max.freq*1.5),
         #ylim=c(0,1e7),
         xlab="coverage",
         ylab="frequency"
         );
    
    lines(scr$data, col=cl[4], lwd=3);
    #scr.med  <- kmerPeaks(read.table(pipe('jellyfish histo -l 10 scr.jf'), header=F), smooth="medAdj");
    #lines(scr.med$data, col=cl[4], lwd=3);

    lines(kfr1$data, col=cl[3], lwd=3);
#    lines(kfr1.med$data, col=cl[3], lwd=3);

    lines(kfr2$data, col=cl[2], lwd=3);
#    lines(kfr2.med$data, col=cl[2], lwd=3);


    abline(v=scr$peak.max.cov, lty=5, lwd=2, col=cl[1]);
    abline(v=kfr1$peak.max.cov, lty=3, lwd=2, col=cl[1]);
    abline(v=kfr2$peak.max.cov, lty=4, lwd=2, col=cl[1]);

    if(!is.null(coverage)){
        abline(v=coverage, lwd=2, col=cl[5])
    }
    #add_psizes(kfr2.ex);

    legend(
	"topright",
    	c("scr","kfr1","kfr2"),
    	lwd=3,
    	lty=c(1,1,1),
    	seg.len=2, 
    	col=cl[c(4,3,2)]
    );

    dev.off()
  
}


#-- task rrm --#

rrm <- function(){

    # ref map coverage
    rrm<-read.table("rrm-cov.tsv", header=T);
    for (df in split(rrm, rrm$id)){ 
    	id=df[1,1];
	if(id == "genome") next;
    	dt=df[,2:3]
  	dt.ex = get_extrema(dt, peaks=c(50,100,150,200));
    	plot(dt, 
    	     type="n", 
    	     main=paste("per-base coverage of reference (", id, ")", sep=""),
    	     xlab="coverage",
    	     ylab="frequency",
	     xlim=c(1,dt.ex$cov[length(dt.ex$cov)]*3),
	     ylim=c(0,dt[,2][dt.ex$cov]*2)
    	    );
    	lines(dt, col=cl[2], lwd=3);
	add_psizes(dt.ex);
    };

}




#-- shared --#

kmerPeaks <- function(d, cov.dev=0.25, cov.min=15, k=3, freq.min=15, smooth="histGeom", breaks=100, trim.tail=TRUE){
    
    colnames(d) <- c("covs", "freqs");

    peaks.covs <- c()
    peaks.freqs <- c()
    peaks.frels <- c()
    peaks.sizes <- c()
    peaks.err.covs <- c()
    peaks.err.freqs <- c()
    peaks.err.frels <- c()
    peaks.err.freqs <- c()
    peaks.minor.covs <- c()
    peaks.minor.freqs <- c()
    peaks.minor.frels <- c()
    peaks.minor.sizes <- c()

    total.size <- NA

    #pits.covs <- c()
    #pits.freqs <- c()
    #peak.max.cov <- NA
    #peak.max.freq <- NA

    # total.size
    total.size=sum(apply(d, MARGIN=1, FUN=prod))

    # remove last cov - might be a max bin that contains any higher cov freqs as well
    d <- d[-dim(d)[1],]
    
    # trim tail
    if(trim.tail){
        # trim tail
        tail.start <- tail(which(runmed(d$freqs, k=7, endrule="median") > freq.min), n=1);
        if(length(tail.start)){
            d <- d[1:tail.start,]
        }
    }

    # smooth da biatch
    d <- switch(smooth,
                medAdj = smooth.medAdj(d, k=k),
                med = smooth.med(d, k=k),
                hist = smooth.hist(d, breaks=breaks),
                histGeom = smooth.hist(d, breaks=breaks, geom=T),
                none = d,
                smooth.unknown(smooth)
                );

    #print(str(d));
    # get peaks
    maxs.i <- localMaxima(d$freqs)

    for (i in maxs.i){
        freq <- d$freqs[i]
        cov <- d$covs[i]
        d.env <- d[which(d$covs >= cov-cov*cov.dev & d$covs < cov+cov*cov.dev),]

        freq.env=d$freqs[which(d$covs >= cov & d$covs < cov+cov*cov.dev)]

        # eval peaks in da hood
        if( length(freq.env)==1 || freq==max(freq.env) ){
            size <- peakSize(d, cov)
            if(length(peaks.covs) && sum(peaks.covs > cov-cov*cov.dev)){ # the fight is on
                if(peaks.freqs[length(peaks.freqs)] <= freq*1.1){ # demote the looser
                    peaks.minor.covs <- c(peaks.minor.covs, peaks.covs[length(peaks.covs)])
                    peaks.covs[length(peaks.covs)] <- cov

                    peaks.minor.freqs <- c(peaks.minor.freqs, peaks.freqs[length(peaks.freqs)])
                    peaks.freqs[length(peaks.freqs)] <- freq

                    peaks.minor.sizes <- c(peaks.minor.sizes, peaks.sizes[length(peaks.sizes)])
                    peaks.sizes[length(peaks.sizes)] <- size
                }
                # wannabe
                peaks.minor.covs <- c(peaks.minor.covs, cov)
                peaks.minor.freqs <- c(peaks.minor.freqs, freq)
                peaks.minor.sizes <- c(peaks.minor.sizes, size)
            }else{
                peaks.covs <- c(peaks.covs, cov);
                peaks.freqs <- c(peaks.freqs, freq);
                peaks.sizes <- c(peaks.sizes, size);
            }
        }
    }

    # trim min cov
    d <- d[which(d$covs >= cov.min),]
    
    # handle error peak (cov 1 or 2)
    if(length(peaks.covs)){
        peaks.err.i <- peaks.covs < cov.min

        peaks.err.covs <- peaks.covs[peaks.err.i]
        peaks.covs <- peaks.covs[!peaks.err.i]

        peaks.err.freqs <- peaks.freqs[peaks.err.i]
        peaks.freqs <- peaks.freqs[!peaks.err.i]

        peaks.err.sizes <- peaks.sizes[peaks.err.i]
        peaks.sizes <- peaks.sizes[!peaks.err.i]
    }

    # king of da hill
    if(! is.null(peaks.freqs)){
        peaks.frels <- d$frels[d$covs %in% peaks.covs]

        peaks.o <- order(peaks.freqs, decreasing=T)
        peaks.covs <- peaks.covs[peaks.o]    
        peaks.freqs <- peaks.freqs[peaks.o]    
        peaks.sizes <- peaks.sizes[peaks.o]
        peaks.frels <- peaks.frels[peaks.o]
    }
        
    if(! is.null(peaks.minor.freqs)){
        peaks.minor.frels <- d$frels[d$covs %in% peaks.minor.covs]
        
        peaks.minor.o <- order(peaks.minor.freqs, decreasing=T)
        peaks.minor.covs <- peaks.minor.covs[peaks.minor.o]    
        peaks.minor.freqs <- peaks.minor.freqs[peaks.minor.o]    
        peaks.minor.sizes <- peaks.minor.sizes[peaks.minor.o]
        peaks.minor.frels <- peaks.minor.frels[peaks.minor.o] 
        
    }

    if(! is.null(peaks.err.freqs)){
        peaks.err.frels <- d$frels[d$covs %in% peaks.err.covs]
        
        peaks.err.o <- order(peaks.err.freqs, decreasing=T)
        peaks.err.covs <- peaks.err.covs[peaks.err.o]    
        peaks.err.freqs <- peaks.err.freqs[peaks.err.o]    
        peaks.err.sizes <- peaks.err.sizes[peaks.err.o]
        peaks.err.frels <- peaks.err.frels[peaks.err.o] 
    }
    
    return(
        list(
            data = d,
            peaks = data.frame(
                covs=peaks.covs,
                freqs=peaks.freqs,
                sizes=peaks.sizes,
                frels=peaks.frels
            ),
            peaks.minor = data.frame(
                covs=peaks.minor.covs,
                freqs=peaks.minor.freqs,
                sizes=peaks.minor.sizes,
                frels=peaks.minor.frels
            ),
            peaks.error = data.frame(
                covs=peaks.err.covs,
                freqs=peaks.err.freqs,
                sizes=peaks.err.sizes,
                frels=peaks.err.frels
            ),
            total.size=total.size
        )
    );
}



# stolen from http://stackoverflow.com/questions/6836409/finding-local-maxima-and-minima
localMaxima <- function(x) {
  # Use -Inf instead if x is numeric (non-integer)
  y <- diff(c(-Inf, x)) > 0L
  # print(rle(y)$lengths)
  y <- cumsum(rle(y)$lengths)
  y <- y[seq.int(1L, length(y), 2L)]
  if (x[[1]] == x[[2]]) {
    y <- y
  }
  y
}


geomSeries <- function(base, max) {
    base^(0:floor(log(max, base)))
}

plotKmerPeaks <- function(...){

}

peakSize <- function(d, peak.cov){
    peak.i <- which(d$covs > peak.cov*0.5 & d$covs <= peak.cov*1.5)
    size <- sum(apply(d[peak.i,c("covs", "freqs")], MARGIN=1, FUN=prod));
    return(size);
}


smooth.unknown <- function(d, smooth){
    write(paste("unknown smooth algorithm", smooth), stderr());
    quit(status=1);
    return(d);
}

smooth.med <- function(d, k=3){
    d$freqs <- runmed(d$freqs, k=k, endrule="median");
    return(d)
}

smooth.medAdj <- function(d, k=3){
    # smooth frequencies, gradually increase k size
    d$freqs[1:50] <- runmed(d$freqs[1:50], k=3, endrule="median");
    d$freqs[45:100] <- runmed(d$freqs[45:100], k=5, endrule="median");
    d$freqs[90:160] <- runmed(d$freqs[90:160], k=9, endrule="median");
    d$freqs[150:length(d$freqs)] <- runmed(d$freqs[150:length(d$freqs)], k=13, endrule="median");
    return(d)
}

smooth.hist <- function(d, breaks=200, geom=FALSE){

    cov.max <- tail(d$covs, n=1)

    # geometric hist bin size
    if(geom){
        ## # fixed number of bases
        base <- cov.max^(1/breaks);
        breaks <- cumsum(round(c(1, diff(base^(1:breaks)))));
        
        # fixed base
        #breaks <- geomSeries(2^(1/6), cov.max)
        
        # handle too small breaks
        # remove too small bins at the start
        breaks <- breaks[which(c(diff(breaks)>3, TRUE))];
        breaks <- c(0:floor(breaks[1]/3 -1)*3+1, breaks); # add 3-bin

        breaks[length(breaks)] <- cov.max + 0.0001 # make sure, max is contained
    }

    cov.hist <- hist(d$covs, breaks=breaks, plot=F)

    # hist freqs for cov.hist :)
    cumsums <- cumsum(cov.hist$counts)
    freqs <- c();

    for(i in 1:length(cumsums)){
        f <- ifelse(i>1, cumsums[i-1]+1, 0)
        t <- cumsums[i]
        freqs <- c(freqs, sum(d$freqs[f:t]))
    }

    # smooth freqs
    #freqs <- runmed(freqs, k=5, endrule="median")

    # scale
    if(cov.hist$equidist){
        sf <- cov.max/length(cov.hist$breaks)
        frels <- freqs/sf;
    }else{ # scale geom
        break.intervals <- diff(cov.hist$breaks)
        frels <- freqs/break.intervals;                                        
    }

    return(data.frame(covs=cov.hist$mids, freqs=freqs, frels=frels))
}

##---- DEPRECATED ----##

get_extrema <- function(data,peaks){
    
    cmins=c();
    fmins=c();
    for ( i in 1:length(peaks)){
        if(i>1){ 
            p1=peaks[i-1] 
        }else{ 
            p1=0
        }
        p2=peaks[i];
        p2p=data[,1] >= p1 & data[,1] < p2;
        cmins=c(cmins,data[,1][p2p][which.min(data[,2][p2p])])
        fmins=c(fmins,data[,2][p2p][which.min(data[,2][p2p])])
    }
    
                                        # peak max
    cmaxs=c();
    fmaxs=c();
    for ( i in 1:length(cmins) ){
        p1=cmins[i];
        if(i==length(cmins)){ 
            p2=data[,1][length(data[,1])]
        }else{ 
            p2=cmins[i+1]
        }

        p2p=data[,1] >= p1 & data[,1] < p2;
        cmaxs=c(cmaxs,data[,1][p2p][which.max(data[,2][p2p])])
        fmaxs=c(fmaxs,data[,2][p2p][which.max(data[,2][p2p])])
    }
    
    psizes=c();
    for ( i in 1:length(cmins) ){
        p1= cmins[i];
        if(i==length( cmins)){
            p2=data[,1][length(data[,1])]
        }else{ 
            p2= cmins[i+1]
        }

        p2p=data[,1] >= p1 & data[,1] < p2;
        psizes=c(psizes,sum(apply(data[p2p,], MARGIN=1, FUN=prod)))
    }
  
    cov=cmaxs[which.max(fmaxs)]
    freq=max(fmaxs)
    
    tbp=sum(apply(data, MARGIN=1, FUN=prod))
    
    return(list(cov.maxima=cmaxs, freq.maxima=fmaxs, cov.minima=cmins, freq.minima=fmins, peaksizes=psizes, cov=cov, freq=freq, total=tbp));
}

add_psizes <- function(data=extrema){
    cov=data$cov;
    y=data$freq*1.2
    for(i in 1:length(data$cov.maxima)){
        text(x= data$cov.maxima[i], y=y, labels=paste(c(round( data$peaksizes[i] /10^3/cov)),"kbp",sep=" "), vfont=c("sans serif", "bold"))
    }
    text(x=data$cov/10, pos=4, y=y*1.2, labels=paste("total:", round( data$total/10^3/data$cov),"kbp",sep=" "), vfont=c("sans serif", "bold"))
}

peakSizes.add <- function(p, rel=FALSE){

    for(i in 1:length(p$covs)){
        y <- ifelse(rel, p$frels[i], p$freqs[i]);
        points(p$covs[i], y, cex=1.3, pch=17);
        y <- y * 1.2;
        text(x= p$covs[i], y=y, labels=paste(c(round( p$sizes[i] /10^3/p$covs[i])),"kbp",sep=" "), vfont=c("sans serif", "bold"))
    }
#    text(x=max(p$covs)/10, pos=4, y=y*1.2, labels=paste("total:", round( data$total/10^3/data$cov),"kbp",sep=" "), vfont=c("sans serif", "bold"))
}

#-- main ---------------------------------------------------------------------#

# globals
cl=rainbow(5);

params <- commandArgs(trailingOnly=T);

do.call(run_task, as.list(params));

if(! is.null(warnings())) warnings();





    

    
