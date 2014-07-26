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

    pdf("scr-seeds.pdf", width=10, height=5);
    
    # cds cluster
    data <- read.table(pipe('jellyfish histo scr-ref.jf'), header=F);

    med = data[,1][cumsum(data[,2]) > sum(data[,2])/2][1];

    data.ex = get_extrema(data, peaks=c(med));

    plot(data, 
    	      type="n", 
    	      main="kmer-coverage of seed reads",
    	      xlab="coverage",
    	      ylab="frequency",
    	      xlim=c(1,data.ex$cov*3),
    	      ylim=c(0,data.ex$freq*1.5)
    );

    lines(data, col=cl[2], lwd=3);
    abline(v=med);
    add_psizes(data.ex);

    legend(
    	"topright",
    	c("scr-seeds"),
    	lwd=3,
    	lty=c(1),
    	seg.len=2, 
    	col=cl[c(2)]
    );

    dev.off()
}


#-- task kfr --#

kfr <- function(coverage){

    pdf("kfr.pdf", width=10, height=5);
    
    # kmer filter
    scr<-read.table(pipe('jellyfish histo scr.jf'), header=F);
    kfr1<-read.table(pipe('jellyfish histo kfr1.jf'), header=F);
    kfr2<-read.table(pipe('jellyfish histo kfr2.jf'), header=F);

    kfr2.ex = get_extrema(kfr2, peaks=c(coverage));
  
    plot(kfr2, 
    	      type="n", 
    	      main="kmer-coverage of subsetted and filtered data sets",
    	      xlab="coverage",
    	      ylab="frequency",
	      xlim=c(1,kfr2.ex$cov*3),
	      ylim=c(0,kfr2.ex$freq*2)
    );

    lines(scr, col=cl[1], lwd=3);
    lines(kfr1, col=cl[4], lwd=3);
    lines(kfr2, col=cl[2], lwd=3);

    abline(v=coverage);
    add_psizes(kfr2.ex);


    legend(
	"topright",
    	c("scr","kfr2","kfr2" ),
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

# globals
cl=rainbow(5);

params <- commandArgs(trailingOnly=T);
task <- params[1]
params <- params[-1]

run_task(task, params);


















