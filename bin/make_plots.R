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
    for (df in split(rrm, "id")){ 
    	plot(df[,2:3], 
    	     type="n", 
    	     main=paste("per-base coverage of reference (", df[1,1], ")", sep=""),
    	     xlab="coverage",
    	     ylab="frequency"
    	    );
    	lines(df[,2:3], col=cl[2], lwd=3);
    };

dev.off();

