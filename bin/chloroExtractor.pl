#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use FindBin;
use lib "$FindBin::Bin/../lib";
use Verbose;
use Data::Dumper;

my %options;

=head1 NAME 

chloroExtractor.pl

=head1 DESCRIPTION

Wrapper to extract chloroplast reads from genomic samples, assemble and annotate them.

=head1 USAGE

  $ perl chloroExtractor.pl --reads=<fasta> [--mates=<fasta>] [options]

=head1 OPTIONS

=over 25

=item --reads=<FASTQ>

path to the reads in fastq format


=cut

$options{'reads=s'} = \(my $opt_reads);

=item --mates=<FASTA>

path to the paired reads in fastq format


=cut

$options{'mates=s'} = \(my $opt_mates);

=item [--prefix=<STRING>] 

prefix for the output files. Default is current directory and a prefix
 <query_file>_-_<reference_file>_-_<aligner>.

=cut

$options{'prefix=s'} = \(my $opt_prefix);

=item [--jellyfish-kmer-size=<INT>] 

desired kmer size for jellyfish

=cut

$options{'jellyfish-kmer-size=i'} = \(my $opt_jellyfish_kmer_size=23);

=item [--velvet-kmer-size=<INT>] 

desired kmer size for velvet

=cut

$options{'velvet-kmer-size=i'} = \(my $opt_velvet_kmer_size=51);
	  
=item [--jellyfish-bin=<FILE>] 

Path to jellyfish binary file. Default tries if jellyfish is in PATH;

=cut

$options{'jellyfish-bin=s'} = \(my $opt_jellyfish_bin = `which jellyfish 2>/dev/null`);

=item [--allpath-correction-bin=<FILE>] 

Path to Allpath-correction executable (ErrorCorrectReads.pl). Default tries if ErrorCorrectReads.pl is in PATH;

=cut

$options{'allpath-correction-bin=s'} = \(my $opt_allpath_correction_bin = `which ErrorCorrectReads.pl 2>/dev/null`);

=item [--velvet-bin=<FILE>] 

Path to the velveth or velvetg binary file. Default tries if velveth is in PATH
The containing folder is used to find both velveth and velvetg.

=cut

$options{'velvet-bin=s'} = \(my $opt_velvet_bin = `which velveth 2>/dev/null`);


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


chomp($opt_jellyfish_bin,$opt_allpath_correction_bin,$opt_velvet_bin);
my $opt_velvet_path = dirname($opt_velvet_bin);

GetOptions(%options) or pod2usage(1);

my $vwga = Verbose->new(
	report_level => $opt_verbose // 0, #/
	format => "[{TIME_ELAPSED}] {MESSAGE}\n",
	line_width => 70
);

my $vbash = Verbose->new(
	report_level => $opt_verbose // 0, #/
	line_width => 70,
	line_delim => "\\\n",
);

my $vplain = Verbose->new(
	report_level => $opt_verbose // 0, #/
	line_width => 70,
);

pod2usage(1) if($opt_help);
pod2usage(-verbose => 99, -sections => "NAME|DESCRIPTION|USAGE|OPTIONS|AUTHORS") if($opt_man);

$vwga->verbose('Checking parameter');
pod2usage(-msg => "Missing parameter reads", -verbose => 0) unless ($opt_reads);
pod2usage(-msg => 'jellyfish not in $PATH and binary (--jellyfish-bin) not specified', -verbose => 0) unless ($opt_jellyfish_bin);
pod2usage(-msg => 'ErrorCorrectReads.pl not in $PATH and binary (--allpath-correction-bin) not specified', -verbose => 0) unless ($opt_allpath_correction_bin);
pod2usage(-msg => 'velvet not in $PATH and binary (--velvet-bin) not specified', -verbose => 0) unless ($opt_velvet_bin);

$opt_prefix = get_prefix() unless $opt_prefix;
my ($prefix_name,$prefix_dir) = fileparse($opt_prefix);

$vwga->verbose('Counting kmers');
$vwga->hline();
my $jellyfish_count_cmd = jellyfish_count_command();
$vbash->verbose( $jellyfish_count_cmd );
my $jellyfish_count_re = qx($jellyfish_count_cmd); 
$vwga->nline();
$vplain->verbose($jellyfish_count_re) if $jellyfish_count_re;
$vwga->exit('ERROR: Counting kmers failed') if $?>> 8;

$vwga->verbose('Merging kmer counts');
$vwga->hline();
my $jellyfish_merge_cmd = jellyfish_merge_command();
$vbash->verbose( $jellyfish_merge_cmd );
my $jellyfish_merge_re = qx($jellyfish_merge_cmd); 
$vwga->nline();
$vplain->verbose($jellyfish_merge_re) if $jellyfish_merge_re;
$vwga->exit('ERROR: Merging kmer counts') if $?>> 8;

$vwga->verbose('Histogramming kmer counts');
$vwga->hline();
my $jellyfish_histo_cmd = jellyfish_histo_command();
$vbash->verbose( $jellyfish_histo_cmd );
my $jellyfish_histo_re = qx($jellyfish_histo_cmd); 
$vwga->nline();
$vplain->verbose($jellyfish_histo_re) if $jellyfish_histo_re;
$vwga->exit('ERROR: Histogramming kmer counts failed') if $?>> 8;

$vwga->verbose('chloroExtractor finished');


=head2 jellyfish_count_command

Returns the command to call jellyfish for counting.

=cut

sub jellyfish_count_command{
	my $cmd = "$opt_jellyfish_bin count -m $opt_jellyfish_kmer_size ";
	$cmd .= "-o $opt_prefix"."_jf_parts ";
	$cmd .= "-s 100000000 -t 20 --both-strands $opt_reads";
	return $cmd;
}

=head2 jellyfish_merge_command

Returns the command to call jellyfish for merging.

=cut

sub jellyfish_merge_command{
	my $cmd = "$opt_jellyfish_bin merge ";
	$cmd .= "-o $opt_prefix"."_full.jf ";
	$cmd .= "$opt_prefix"."_jf_parts*";
	return $cmd;
}

=head2 jellyfish_histo_command

Returns the command to call jellyfish for histogramming.

=cut

sub jellyfish_histo_command{
	my $cmd = "$opt_jellyfish_bin histo ";
	$cmd .= "-o $opt_prefix"."_full_histo.jf ";
	$cmd .= "--threads 20 --high 100000 ";
	$cmd .= "$opt_prefix"."_full.jf";
	return $cmd;
}



=head2 get_prefix

Returns a default prefix if none is specified by the user. Style: <reads_-_> (without .fq/.fastq)

=cut

sub get_prefix{
	my ($reads_name,$reads_path,$reads_suffix) = fileparse($opt_reads, qw(.fq .fastq));
	return './'.$reads_name.'_-_';
}

=head1 LIMITATIONS

This pipeline is meant to give a quick and easy to use possibility to generate and visualize whole genome alignments.
It is not meant for extensive analyses of whole genomes. For this purpose see the documentation of the alignment programs
and consider using them.  

=head1 CHANGELOG

=head2 08.11.2012

=over

=item * The --unmask option is now set by default, there is an additional --mask flag to enforce masking by lastz

=item * There are two additional option --minidentity <INT> and --minlength <INT> that are passed through to the Circos Parsers

=item * A statistic of total and aligned contigs is returned as part of the log

=back


=head1 AUTHORS

=over

=item * Markus Ankenbrand, markus.ankenbrand@stud-mail.uni-wuerzburg.de

=item * Thomas Hackl, thomas.hackl@uni-wuerzburg.de

=back


