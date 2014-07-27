#-- make_plots.R ------------------------------------------------------------------#

run_task <- function(task, ...){
    switch(task,
           scr=scr(),
           kfr=kfr(...),
           unknown_task(task)
           )
}	 


unknown_task <- function(task){
    print(paste("Unknown task:", task), quote=F); 
    quit(status = 1);
}

#-- task scr --#

scr <-function(){

    pdf("scr-seeds.pdf", width=10, height=6);

    # from seed reads
    scr <- read.table(pipe('jellyfish dump -c --tab scr-ref.jf | cut -f2'), header=F);
    med <- median(scr[,1]);
    scr <- scr[scr < 4*med]
    scr.hist <- hist(scr, breaks=200, plot=F);
    scr.hist.df <- data.frame(coverage=scr.hist$mids, frequency=scr.hist$counts)
    scr.ex = get_extrema(scr.hist.df, peaks=c(med));
    
    # from entire input
    raw <- read.table(pipe('jellyfish dump -c --tab jf0.jf | cut -f2'), header=F);
    raw <- raw[raw < 4*med]
    raw.hist <- hist(raw, breaks=600, plot=F);
    raw.hist.df <- data.frame(coverage=raw.hist$mids, frequency=raw.hist$counts)
    raw.ex = get_extrema(raw.hist.df, peaks=c(med));
    
    # DEPRECATED: from histo
    #scr.raw <- read.table(pipe('jellyfish histo scr-ref.jf'), header=F);
    #med <- scr.raw[,1][cumsum(scr.raw[,2]) > sum(scr.raw[,2])/2][1];

    # TODO
    # bin the median to 100 for plotting
    # binsize <- med/100
    # cut(scr, breaks=, labels= ...)
    #scr <- scr.raw
    

    plot(scr.hist.df, 
    	      type="n", 
    	      main="kmer-coverage of seed reads",
    	      xlab="coverage",
    	      ylab="frequency",
    	      xlim=c(1,scr.ex$cov*3),
    	      ylim=c(0,raw.ex$freq*1.5)
    );

    lines(scr.hist.df, col=cl[2], lwd=3);
    lines(raw.hist.df, col=cl[4], lwd=3);
    abline(v=med);
    add_psizes(scr.ex);
    add_psizes(raw.ex);

    legend(
    	"topright",
    	c("scr-seeds", "total data (x 1/3)"),
    	lwd=3,
    	lty=c(1,1),
    	seg.len=2, 
    	col=cl[c(2,4)]
    );

    dev.off()
}


#-- task kfr --#

kfr <- function(coverage){

    pdf("kfr.pdf", width=10, height=5);
    coverage <- as.integer(coverage)
    # kmer filter
    scr <- read.table(pipe('jellyfish histo scr.jf'), header=F);
    scr <- scr[scr[,1]<4*coverage,]
    kfr1 <- read.table(pipe('jellyfish histo kfr1.jf'), header=F);
    kfr1 <- kfr1[kfr1[,1]<4*coverage,]
    kfr2 <- read.table(pipe('jellyfish histo kfr2.jf'), header=F);
    kfr2 <- kfr2[kfr2[,1]<4*coverage,]

    kfr2.ex <- get_extrema(kfr2, peaks=c(coverage/2, coverage));

    plot(kfr2, 
    	      type="n", 
    	      main="kmer-coverage of subsetted and filtered data sets",
    	      xlab="coverage",
    	      ylab="frequency",
	      xlim=c(1,coverage*3),
	      ylim=c(0,kfr2.ex$freq*2)
    );

    lines(scr, col=cl[1], lwd=3);
    lines(kfr1, col=cl[4], lwd=3);
    lines(kfr2, col=cl[2], lwd=3);

    abline(v=coverage);
    add_psizes(kfr2.ex);


    legend(
	"topright",
    	c("scr","kfr1","kfr2" ),
    	lwd=3,
    	lty=c(1,1,1),
    	seg.len=2, 
    	col=cl[c(1,4,2)]
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


#-- main ---------------------------------------------------------------------#

# globals
cl=rainbow(5);

params <- commandArgs(trailingOnly=T);
task <- params[1]
params <- params[-1]

run_task(task, params);


















