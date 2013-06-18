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
$R->send(qq`data<-read.table("$opt_histo")`) ;
print("Loading functions\n");
$R->send(qq`source("$FindBin::Bin/../find_chloro_kmer_peek_ml/find_chloro_kmer_peek_ml.R")`) ;
$R->send(q`parameters<-cbind(c(1000,1500,2000),c(5500,6000,6500),c(3000,4000,4500),c(300,500,1000))`) ;
$R->send(q`opt<-optim(par=parameters[2,], fn=logLhistoMinus, histo=data, lower=parameters[1,], upper=parameters[3,], method="L-BFGS-B")`) ;
my $pdf_file = "$prefix_dir"."/"."$prefix_name"."histo.pdf";
$R->send(qq`pdf("$pdf_file")`);
$R->send(q`plotResult(opt$par, data)`) ;
$R->send(q`plotResultFull(opt$par, data)`) ;
$R->send(q`dev.off()`) ;
my $ret = $R->read ;
print("$ret\n");
$R->stopR() ;


print('findChloroPeek.pl finished');



=head2 get_prefix

Returns a default prefix if none is specified by the user. Style: <reads_-_> (without .fq/.fastq)

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


