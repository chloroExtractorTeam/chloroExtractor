

num_points <- 10000
xdat <- runif(num_points)*100000
ydat <- rnorm(num_points)*1e7

dat<-read.table("Rg_ChlGen_full_histo.jf");
xdat <- dat$V1;
ydat <- dat$V2;

ready="p";
xlimit<-range(xdat)
ylimit<-range(ydat)
peaklocation<-xlimit
while (ready != "c")
{
    plot(xdat, ydat, xlim=xlimit, ylim=ylimit)
    values<-locator(n=2, type="o")
    cat("New xrange for plot (p) or final location (l) or cancel (c)?")
    choice<-scan(n=1, what="character")
    if (choice[1] == "p")
        {
            xlimit <- range(values$x)
            ylimit <- range(ydat[(xdat>=xlimit[1])&(xdat<=xlimit[2])])
        }
    if (choice[1] == "l")
        {
            peaklocation<-range(values$x)
            dev.off()
            break;
        }
    if (choice[1] == "c")
        {
            dev.off();
            break;
        }
}

peaklocation<-round(peaklocation);

cat(paste("Final location of the peak is: ", peaklocation[1], "-", peaklocation[2], "\n", sep=""));

    
