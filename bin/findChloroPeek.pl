#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use FindBin;
use lib "$FindBin::Bin/../lib";
use Data::Dumper;
use Statistics::R;

my %options;

=head1 NAME 

findChloroPeek.pl

=head1 DESCRIPTION

Perl script that calls R with correct configuration to find the chloroplast peek in a kmer histogram.

=head1 USAGE

  $ perl findChloroPeek.pl --histo=<FILE> [options]

=head1 OPTIONS

=over 25

=item --histo=<FILE>

path to the kmer histogram file

=cut

$options{'histo=s'} = \(my $opt_histo);

=item [--prefix=<STRING>] 

prefix for the output files. Default is current directory and a prefix
 <query_file>_-_<reference_file>_-_<aligner>.

=cut

$options{'prefix=s'} = \(my $opt_prefix);

=item [--lower_l=<INT>] 

lower constraint on the "lower" argument (default:1000)

=cut

$options{'lower_l=i'} = \(my $opt_lower_l=1000);

=item [--lower_s=<INT>] 

start value of the "lower" argument (default:1500)

=cut

$options{'lower_s=i'} = \(my $opt_lower_s=1500);

=item [--lower_u=<INT>] 

upper constraint on the "lower" argument (default:2000)

=cut

$options{'lower_u=i'} = \(my $opt_lower_u=2000);

=item [--upper_l=<INT>] 

lower constraint on the "upper" argument (default:8000)

=cut

$options{'upper_l=i'} = \(my $opt_upper_l=8000);

=item [--upper_s=<INT>] 

start value of the "upper" argument (default:9000)

=cut

$options{'upper_s=i'} = \(my $opt_upper_s=9000);

=item [--upper_u=<INT>] 

upper constraint on the "upper" argument (default:9999)

=cut

$options{'upper_u=i'} = \(my $opt_upper_u=9999);

=item [--mean_l=<INT>] 

lower constraint on the "mean" argument (default:1000)

=cut

$options{'mean_l=i'} = \(my $opt_mean_l=1000);

=item [--mean_s=<INT>] 

start value of the "mean" argument (default:5000)

=cut

$options{'mean_s=i'} = \(my $opt_mean_s=5000);

=item [--mean_u=<INT>] 

upper constraint on the "mean" argument (default:9999)

=cut

$options{'mean_u=i'} = \(my $opt_mean_u=9999);


=item [--sd_l=<INT>] 

lower constraint on the "sd" argument (default:500)

=cut

$options{'sd_l=i'} = \(my $opt_sd_l=500);

=item [--sd_s=<INT>] 

start value of the "sd" argument (default:1000)

=cut

$options{'sd_s=i'} = \(my $opt_sd_s=1000);

=item [--sd_u=<INT>] 

upper constraint on the "sd" argument (default:1500)

=cut

$options{'sd_u=i'} = \(my $opt_sd_u=1500);

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

print("Starting R\n");
my $R = Statistics::R->new() ;
$R->startR ;
print("Loading data\n");
$R->send(qq`data<-read.table("$opt_histo")`);
$R->send(q`maxLowess<-lowess(data$V1, data$V2)`);
$R->send(q`xMax<-which(maxLowess[[2]]==max(maxLowess[[2]]))`);
$R->send(q`print(xMax)`);
my $ret = $R->read ;
print("Max at:\n");
print("$ret\n");
$R->send(q`minLowess<-lowess(data$V1[data$V1<xMax],data$V2[data$V1<xMax])`);
$R->send(q`xMin<-which(minLowess[[2]]==min(minLowess[[2]]))`);
$R->send(q`print(xMin)`);
print("Min at:\n");
print("$ret\n");
my $pdf_file = "$prefix_dir"."/"."$prefix_name"."_fits.pdf";
$R->send(qq`c(pdf("$pdf_file"),plot(data[,1], data[,2], ylim=c(0,3*maxLowess[[2]][xMax])),lines(maxLowess, col="red", lwd=3),lines(minLowess, col="green", lwd=3),dev.off())`);
$R->stopR() ;
print("Output plot written\n");

print("findChloroPeek.pl finished\n");



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


