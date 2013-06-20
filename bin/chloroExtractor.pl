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

=item --insertsize=<INT>

Insert size of the paired library as passed to downstream programs.


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

=item [--sickle-bin=<FILE>] 

Path to sickle executable. Default tries if sickle is in PATH;

=cut

$options{'sickle-bin=s'} = \(my $opt_sickle_bin = `which sickle 2>/dev/null`);

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
pod2usage(-msg => 'sickle not in $PATH and binary (--sickle-bin) not specified', -verbose => 0) unless ($opt_sickle_bin);
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

$vwga->verbose('Finding chloroplast peak in kmer histogram');
$vwga->hline();
my $findChloroPeak_cmd = findChloroPeak_command();
$vbash->verbose( $findChloroPeak_cmd );
my $findChloroPeak_re = qx($findChloroPeak_cmd); 
$vwga->nline();
$vplain->verbose($findChloroPeak_re) if $findChloroPeak_re;
open(IN, "<$opt_prefix"."_minmax.tsv") or die "Can't open file $opt_prefix"."_minmax.tsv$!";
my ($min, $max) = split(/\t/,<IN>);
chomp $max;
$vwga->exit('ERROR: Chloroplast peak detection failed') if $?>> 8;

$max *= 3; # Take three times the maximal x value (expect IR at double)

$vwga->verbose('Dumping kmers in count range $min - $max');
$vwga->hline();
my $jellyfish_dump_cmd = jellyfish_dump_command();
$vbash->verbose( $jellyfish_dump_cmd );
my $jellyfish_dump_re = qx($jellyfish_dump_cmd); 
$vwga->nline();
$vplain->verbose($jellyfish_dump_re) if $jellyfish_dump_re;
$vwga->exit('ERROR: Dumping kmers failed') if $?>> 8;

$vwga->verbose('Quality trimming raw reads');
$vwga->hline();
my $quality_trimming_cmd = quality_trimming_command();
$vbash->verbose( $quality_trimming_cmd );
my $quality_trimming_re = qx($quality_trimming_cmd); 
$vwga->nline();
$vplain->verbose($quality_trimming_re) if $quality_trimming_re;
$vwga->exit('ERROR: Quality trimming failed') if $?>> 8;

$vwga->verbose('Dumping reads by kmer coverage');
$vwga->hline();
my $initial_read_dump_cmd = initial_read_dump_command();
$vbash->verbose( $initial_read_dump_cmd );
my $initial_read_dump_re = qx($initial_read_dump_cmd); 
$vwga->nline();
$vplain->verbose($initial_read_dump_re) if $initial_read_dump_re;
$vwga->exit('ERROR: Dumping reads by kmer coverage') if $?>> 8;

########### ALLPATH Correction

########### velvet assembly

########### contig Filtering (simple size filter)

########### Iteration

########### IR-Resolution

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

=head2 jellyfish_dump_command

Returns the command to call jellyfish for dumping.

=cut

sub jellyfish_dump_command{
	my $cmd = "$opt_jellyfish_bin dump ";
	$cmd .= "--column --tab ";
	$cmd .= "-o $opt_prefix"."_dump_$min"."_$max".".jf ";
	$cmd .= "--lower-count=$min --upper-count=$max ";
	$cmd .= "$opt_prefix"."_full.jf";
	return $cmd;
}

=head2 findChloroPeak_command

Returns the command to call findChloroPeak.pl for chloroplast peak detection.

=cut

sub findChloroPeak_command{
	my $cmd = "perl $FindBin::Bin/findChloroPeak.pl ";
	$cmd .= "--histo $opt_prefix"."_full_histo.jf ";
	$cmd .= "--prefix $opt_prefix";
	return $cmd;
}

=head2 quality_trimming_command

Returns the command to call sickle for quality trimming of the raw reads.

=cut

sub quality_trimming_command{
	my $cmd = "$opt_sickle_bin pe ";
	# TODO PHRED offset is fixed to sanger (33) at the moment
	$cmd .= "-f $opt_reads -r $opt_mates -t sanger ";
	$cmd .= "-o $opt_prefix"."_reads_trimmed_1.fq ";
	$cmd .= "-p $opt_prefix"."_reads_trimmed_2.fq ";
	$cmd .= "-s $opt_prefix"."_reads_trimmed_singles.fq ";
	$cmd .= "-l 50";
	return $cmd;
}

=head2 initial_read_dump_command

Returns the command to call Kmer.pl for the initial dumping of reads (by kmer coverage).

=cut

sub initial_read_dump_command{
	my $cmd = "perl $FindBin::Bin/Kmer.pl ";
	$cmd .= "--kmers $opt_prefix"."_dump_$min"."_$max".".jf ";
	$cmd .= "--reads $opt_prefix"."_reads_trimmed_1.fq ";
	$cmd .= "--mates $opt_prefix"."_reads_trimmed_1.fq ";
	$cmd .= "--out $opt_prefix"."_reads_trimmed_dumped ";
	$cmd .= "--histo $opt_prefix"."_trusted_kmers.histo ";
	$cmd .= "--cutoff 50 --maxreads 200000 --notrustall";
	return $cmd;
}

=head2 error_correction_command

Returns the command to call ErrorCorrectReads.pl for ErrorCorrection of the dumped reads.

=cut

sub error_correction_command{
	my $cmd = "$opt_allpath_correction_bin ";
	$cmd .= "PHRED_ENCODING=33 READS_OUT=SRR492316_trimmed_dump_corr ";
	$cmd .= "PAIRED_READS_A_IN=SRR492316_trimmed_dump_1.fq PAIRED_READS_B_IN=SRR492316_trimmed_dump_2.fq ";
	$cmd .= "PAIRED_SEP=300 PAIRED_STDEV=100 PLOIDY=1 THREADS=20"
	return $cmd;
}

=head2 get_prefix

Returns a default prefix if none is specified by the user. Style: <reads_-> (without .fq/.fastq)

=cut

sub get_prefix{
	my ($reads_name,$reads_path,$reads_suffix) = fileparse($opt_reads, qw(.fq .fastq));
	return './'.$reads_name.'_-';
}

=head1 LIMITATIONS

If you encounter a bug, please drop me a line.

=head1 AUTHORS

=over

=item * Markus Ankenbrand, markus.ankenbrand@stud-mail.uni-wuerzburg.de

=back


