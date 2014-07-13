#-- methods ------------------------------------------------------------------#

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
    cmins=c(cmins,data[,1][p1:p2][which.min(data[,2][p1:p2])])
    fmins=c(fmins,data[,2][p1:p2][which.min(data[,2][p1:p2])])
  }
  
  # peak max
  cmaxs=c();
  fmaxs=c();
  for ( i in 1:length(cmins) ){
    p1=cmins[i];
    if(i==length(cmins)){ 
      p2=length(data[,2])
    }else{ 
      p2=cmins[i+1]
    }
    cmaxs=c(cmaxs,data[,1][p1:p2][which.max(data[,2][p1:p2])])
    fmaxs=c(fmaxs,data[,2][p1:p2][which.max(data[,2][p1:p2])])
  }
  
  psizes=c();
  for ( i in 1:length(cmins) ){
    p1= cmins[i];
    if(i==length( cmins)){ 
      p2=length(data[,2])
    }else{ 
      p2= cmins[i+1]
    }
    
    psizes=c(psizes,sum(apply(data[p1:p2,], MARGIN=1, FUN=prod)))
  }

  cov=cmaxs[which.max(data[,2][cmaxs])]
  freq=fmaxs[which.max(data[,2][cmaxs])]
  
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

pdf("ptx.pdf", width=10, height=5);

    cl=rainbow(5);

    # cds cluster
    scr<-read.table("scr-cov.tsv", header=F, fill=T);
    plot(
         scr,
	 outline=F,
	 las=2,
	 ylab="coverage",
	 main="CDS cluster coverages"
        );

    # ref map coverage
    rrm<-read.table("rrm-cov.tsv", header=T);
    for (df in split(rrm, rrm$id)){ 
    	id=df[1,1];
	if(id == "genome") next;
    	dt=df[,2:3]
  	dt.ex = get_extrema(dt, peaks=c(200));
    	plot(dt, 
    	     type="n", 
    	     main=paste("per-base coverage of reference (", id, ")", sep=""),
    	     xlab="coverage",
    	     ylab="frequency",
	     xlim=c(1,dt.ex$cov[length(dt.ex$cov)]*3),
	     ylim=c(0,dt[,2][dt.ex$cov]*1.6)
    	    );
    	lines(dt, col=cl[2], lwd=3);
	add_psizes(dt.ex);
    };

    # kmer filter
    scr<-read.table(pipe('jellyfish histo scr.jf'), header=F);
    kfr<-read.table(pipe('jellyfish histo kfr.jf'), header=F);
    kfr2<-read.table(pipe('jellyfish histo kfr2.jf'), header=F);

    kfr2.ex = get_extrema(kfr2, peaks=c(200));
  
    plot(kfr2, 
    	      type="n", 
    	      main="kmer-coverage",
    	      xlab="coverage",
    	      ylab="frequency",
	      xlim=c(1,kfr2.ex$cov[length(kfr2.ex$cov)]*3),
	      ylim=c(0,kfr2[,2][kfr2.ex$cov]*1.6)
    );

    lines(scr, col=cl[1], lwd=3);
    lines(kfr, col=cl[4], lwd=3);
    lines(kfr2, col=cl[2], lwd=3);

    add_psizes(kfr2.ex);


    legend(
	"topright",
    	c("scr","kfr","kfr2" ),
    	lwd=3,
    	lty=c(1,1,1),
    	seg.len=2, 
    	col=cl[c(1,4,2)]
    );
    

dev.off();





