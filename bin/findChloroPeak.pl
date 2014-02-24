#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use File::Spec;
use FindBin;
use lib "$FindBin::Bin/../lib";
use Data::Dumper;
use Statistics::R;

my %options;

=head1 NAME 

findChloroPeak.pl

=head1 DESCRIPTION

Perl script that calls R with correct configuration to find the chloroplast peak in a kmer histogram.

=head1 USAGE

  $ perl findChloroPeak.pl --histo=<FILE> [options]

=head1 OPTIONS

=over 25

=item --histo=<FILE>

path to the kmer histogram file

=cut

$options{'histo=s'} = \(my $opt_histo);

=item --lower=<INT>

lower bound to search for peak
Default: 0

=cut

$options{'lower=i'} = \(my $opt_lower=0);

=item --upper=<INT>

upper bound to search for peak
Default: No upper bound

=cut

$options{'upper=i'} = \(my $opt_upper);

=item --smooth=<FLOAT>

This passes the f option to the R lowess function
	f: the smoother span. This gives the proportion of points in the
        plot which influence the smooth at each value.  Larger values
        give more smoothness
Default: 0.1

=cut

$options{'smooth=f'} = \(my $opt_smooth = 0.1);


=item [--prefix=<STRING>] 

prefix for the output files. Default is current directory and a prefix
 <query_file>_-_<reference_file>_-_<aligner>.

=cut

$options{'prefix=s'} = \(my $opt_prefix);


=item [--[no]verbose] 

verbose is default.

=cut

$options{'verbose!'} = \(my $opt_verbose = 1);

=item [--help] 

show help

=cut

$options{'help|?'} = \(my $opt_help);

=item [--man] 

show man page

=cut

$options{'man'} = \(my $opt_man);

=back


=head1 CODE

=cut


GetOptions(%options) or pod2usage(1);

pod2usage(1) if($opt_help);
pod2usage(-verbose => 99, -sections => "NAME|DESCRIPTION|USAGE|OPTIONS|AUTHORS") if($opt_man);

print("Checking parameter\n");
pod2usage(-msg => "Missing parameter --histo", -verbose => 0) unless ($opt_histo);

$opt_prefix = get_prefix() unless $opt_prefix;
my ($prefix_name,$prefix_dir) = fileparse($opt_prefix);
my $abs_prefix_dir = File::Spec->rel2abs( $prefix_dir );
my $abs_histo = File::Spec->rel2abs( $opt_histo );

print("Starting R\n");
my $R = Statistics::R->new() ;
$R->startR ;
print("Loading data\n");
$R->send(qq`data<-read.table("$abs_histo")`);
print("Loading functions\n");
$R->send(qq`source("$FindBin::Bin/../find_chloro_kmer_peak/find_chloro_kmer_peak_lowess.R")`);
print("Finding peaks\n");
if($opt_upper){
	$R->send(qq`maxLowess<-findMaxLowess(data, lower=$opt_lower, upper=$opt_upper, f=$opt_smooth)`);
}
else{
	$R->send(qq`maxLowess<-findMaxLowess(data, lower=$opt_lower, f=$opt_smooth)`);
}
$R->send(q`xMax<-maxLowess[[1]][which.max(maxLowess[[2]])[1]]`);
$R->send(q`print(xMax)`);
my $ret = $R->read;
print("Max at:\n");
$ret=~s/\[1\]\s+//;
my $max=$ret;
print("$ret\n");

$R->send(q`minLowess<-lowess(data$V1[data$V1<xMax],data$V2[data$V1<xMax])`);
$R->send(q`xMin<-minLowess[[1]][which.min(minLowess[[2]])[1]]`);
$R->send(q`print(xMin)`);
$ret = $R->read;
$ret=~s/\[1\]\s+//;
my $min=$ret;
print("Min at:\n");
print("$ret\n");

open(OUT, ">$abs_prefix_dir"."/"."$prefix_name"."_minmax.tsv") or die "Can't open file $abs_prefix_dir"."/"."$prefix_name"."_minmax.tsv$!";
print OUT "$min\t$max\n";
close OUT or die "$!";
print("Output min-max written to $abs_prefix_dir"."/"."$prefix_name"."_minmax.tsv\n");

my $pdf_file = "$abs_prefix_dir"."/"."$prefix_name"."_fit.pdf";
$R->send(qq`c(pdf("$pdf_file"),plot(data[,1], data[,2], ylim=c(0,2*max(maxLowess[[2]])), xlim=c(0,3*xMax)), lines(maxLowess, col="red", lwd=3),lines(minLowess, col="green", lwd=3), abline(v=xMin, col="green"), dev.off())`);
$R->stopR() ;
print("Output plot written to $abs_prefix_dir"."/"."$prefix_name"."_fit.pdf\n");
print("findChloroPeak.pl finished\n");



=head2 get_prefix

Returns a default prefix if none is specified by the user. Style: <histo_-_> (without .jf/.histo)

=cut

sub get_prefix{
	my ($histo_name,$histo_path,$histo_suffix) = fileparse($opt_histo, qw(.jf .histo));
	return './'.$histo_name.'_-_';
}

=head1 LIMITATIONS

If you encounter a bug, please drop me a line.

=head1 AUTHORS

=over

=item * Markus Ankenbrand, markus.ankenbrand@stud-mail.uni-wuerzburg.de

=back


