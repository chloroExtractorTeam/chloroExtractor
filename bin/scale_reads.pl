#!/usr/bin/env perl

=head1 NAME

scale_reads.pl

=head1 DESCRIPTION

Generate fastq library files with sufficient chloroplast coverage (200X). Coverage is deduced from coverage distributions of read mapping onto cluster of conserved CDS sequences. The mapping will be terminated as soon as enough chloroplast sequences have been found and fastq library files subsets will be created accordingly.

=head1 SYNOPSIS

  $ perl scale_reads.pl -1 lib_1.fq -2 lib_2.fq -ref_cluster cds.fa -t 200

=head1 OPTIONS

=over

=item -1|--reads

Input reads file, first of pair.

=item -2|--mates

Input reads file, second of pair

=item -o|--out

Output prefix.

=item -r|--ref-cluster

Reference sequences for mapping. Fasta or already created bowtie2 index prefix.

=item -t|--target-coverage [200]

Return read files with this estimated coverage.

=item -c|--config

Use user customized config file. Superseeds default config.

=item -V|--version

Display version.

=item -h|--help

Display this help.

=back

=head1 CHANGELOG

see git log.

=head1 TODO

=head1 CODE

=cut

#-----------------------------------------------------------------------------#
# Modules

# core
use strict;
use warnings;
no warnings 'qw';

use Carp;
use Getopt::Long qw(:config no_ignore_case);
use Pod::Usage;
use Log::Log4perl qw(:no_extra_logdie_message);
use Log::Log4perl::Level;

use Data::Dumper;
$Data::Dumper::Sortkeys = 1;

use FindBin qw($RealBin $Script);
use lib "$RealBin/../lib/";

use File::Basename;
use File::Copy;


# additional modules
use Cfg;

use Sam::Parser;
use Sam::Seq;
use Sam::Alignment;

use Jellyfish;

#enable in script mapping
use Bowtie2;


#-----------------------------------------------------------------------------#
# Globals

our $VERSION = 0.01;

our $ID = 'scr';

# get a logger
my $L = Log::Log4perl::get_logger();

my $log_cfg = 'log4perl.rootLogger                     = INFO, Screen
log4perl.appender.Screen                = Log::Log4perl::Appender::Screen
log4perl.appender.Screen.stderr         = 1
log4perl.appender.Screen.layout         = PatternLayout
log4perl.appender.Screen.layout.ConversionPattern = [%d{yy-MM-dd HH:mm:ss}] ['.$ID.'] %m%n
';

Log::Log4perl->init( \$log_cfg );



#-----------------------------------------------------------------------------#
# GetOptions

my %opt = (config => []);

GetOptions( # use %opt (Cfg) as defaults
	    \%opt, qw(
                target_coverage|target-coverage|coverage=i
                reads|1=s
                mates|2=s
		ref_cluster|ref-cluster|r=s
		kmer_hash|kmer-hash=s
		max_reads|max-reads=i
		out|o=s
		bowtie2_bin|bowtie2-bin=s
		jellyfish_bin|jellyfish-bin=s
		threads=i
		version|V!
		debug|D!
		help|h!
		config|c=s{,}
	)
    ) or $L->logcroak('Failed to "GetOptions"');

# help
$opt{help} && pod2usage(1);

# version
if($opt{version}){
    print "$VERSION\n"; 
    exit 0;
}

# debug level
$opt{quiet} && $L->level($WARN);
$opt{debug} && $L->level($DEBUG);


##----------------------------------------------------------------------------##
# Config

my %cfg;

# core
my $core_cfg = "$RealBin/../".basename($Script, qw(.pl)).".cfg";

if(-e $core_cfg){
    $opt{core_config} = File::Spec->rel2abs($core_cfg);
    %cfg = (%cfg, Cfg->Read($opt{core_config}, $ID));
}



# read all configs
if (@{$opt{config}}){
    foreach my $cfg ( @{$opt{config}} ){
	# $L->info("Reading config $cfg");
	$cfg=File::Spec->rel2abs($cfg);
	%cfg = (%cfg, Cfg->Read($cfg, $ID));
    }
}

# create template for user cfg
if(defined $opt{create_config}){
    pod2usage(-msg => 'To many arguments', -exitval=>1) if @ARGV > 1;
    my $user_cfg = Cfg->Copy($core_cfg, $opt{create_config}) or $L->logdie("Creatring config failed: $!");
    $L->info("Created config file: $user_cfg");
    exit 0;
}


# Merge opt and cfg
%opt = (%cfg, %opt);


##----------------------------------------------------------------------------##
# required	
for(qw(reads mates ref_cluster kmer_hash)){
    if(ref $opt{$_} eq 'ARRAY'){
	pod2usage("required: --$_") unless @{$opt{$_}}
    }else{
	pod2usage("required: --$_") unless defined ($opt{$_})
    }
};


# debug level
$L->level($DEBUG) if $opt{debug};
$L->debug('Verbose level set to DEBUG');

$L->debug(Dumper(\%opt));



#-----------------------------------------------------------------------------#
# MAIN

# TODO: working dir
my $opt_o1 = $opt{out}."_1.fq";
my $opt_o2 = $opt{out}."_2.fq";



my $bowtie2 = Bowtie2->new(
    path => $opt{bowtie2_bin}
    );

# TODO: bowtie2 generate db -> prevent issues with different indices
# on different architectures or bowtie2 versions

unless(-e $opt{ref_cluster}.'.1.bt2'){
    $L->info("Building bowtie2 index");
    $bowtie2->build($opt{ref_cluster});
}else{
    $L->info("Using existing bowtie2 index $opt{ref_cluster}.*.bt2");
}

$L->info("Running bowtie2");

$bowtie2->run(
    @{$opt{bowtie2_params}},
    "-1" => $opt{reads},
    "-2" => $opt{mates},
    "-x" => $opt{ref_cluster},
    "-p" => $opt{threads} || 1
    );

$L->debug(Dumper($opt{bowtie2_params}));

# read output on the fly
my $sp = Sam::Parser->new(
    fh => $bowtie2->stdout
    );

my $bam = $opt{out}."-cds.bam";
open(BAM, "| samtools view -Sbu /dev/fd/0 > $bam") or $L->logdie($!);
open(BED, ">", $opt{out}."-cds.bed") or $L->logdie($!);
open(FQ, ">", $opt{out}."-cds.fq") or $L->logdie($!);

my %h; 
my $ss;

my %seqs;

Sam::Seq->Trim(0); # disable trimming of read

# read output on the fly
while(%h = $sp->next_header_line('@SQ')){
    $seqs{$h{'SN'}} = Sam::Seq->new(
	id => $h{'SN'},
	len => $h{'LN'}
	);

    print BAM '@SQ'."\tSN:$h{'SN'}\tLN:$h{'LN'}\n";
    print BED "$h{'SN'}\t$h{'LN'}\n";
}

close BED;

my $current_cov;
my $last_id;
my $closest_ref;
my %refs;

my $c;
my $tlen_sum; 
my $seqwithtlen;


while(my $aln = $sp->next_aln()){
    print BAM "$aln\n";
    print FQ "@",$aln->qname,"\n", $aln->seq, "\n+\n",$aln->qual,"\n";
    
    # closest ref
    my $id = $aln->rname();
    my ($ref_id) = $id =~ /(NC_\d+)/;
    $refs{$ref_id}++; # increment mapping on ref count

    # insert size
    if (abs($aln->tlen) > 0){
	$tlen_sum+= abs($aln->tlen); # compute isize
	$seqwithtlen++;
    }

    $c++;
    last if $c >= $opt{max_reads}
}       


$bowtie2->cancel();
close BAM;
close FQ;


$L->info("Recalculating accurate per base coverages");
bam_coverage($opt{out});


if($c < $opt{max_reads}){
    $L->warn("Could not detect sufficient plastid data in input reads");	
    $L->info("You might need to increase the amount of input data");
    $L->info("Also make sure, your library contains plastid reads at all");

    exit 1;
}

$current_cov = estimate_kmer_coverage($c);

# what if not enough coverage in entire file.
if(! $current_cov || $current_cov < $opt{target_coverage}){
    $L->warn("Could not detect sufficient plastid data in input reads");	
    $L->info("You might need to increase the amount of input data");
    $L->info("Also make sure, your library contains plastid reads at all");

    exit 1;
}



$closest_ref = (sort{$refs{$b} <=> $refs{$a}}keys %refs)[0];

print "coverage\t", $current_cov, "\n";
print "insert_size\t", int($tlen_sum/$seqwithtlen), "\n";
print "closest_ref\t", $closest_ref || 'NA', "\n";



$L->info("Creating libraries");



my $current_reads = qx/wc -l <$opt{reads}/ / 4;
my $target_reads = int($current_reads / $current_cov * $opt{target_coverage});

$L->debug("Got $current_reads reads, extracting $target_reads reads");

my $target_lines = $target_reads * 4;

my $sed_cmd1 = "sed '".$target_lines."q' ".$opt{reads}." >".$opt_o1;

$L->debug($sed_cmd1);
qx($sed_cmd1);

my $sed_cmd2 = "sed '".$target_lines."q' ".$opt{mates}." >".$opt_o2;

$L->debug($sed_cmd2);
qx($sed_cmd2);






#-- SUBS ----------------------------------------------------------------------#




sub median {
  my @array = @{$_[0]};
  my $median = (sort{$a<=>$b}@array)[@array/2];
  return $median;
}


sub estimate_kmer_coverage{
    my $c = shift;
    my ($reads) = $opt{out}."-cds.fq";
    
    my $jf = Jellyfish->new(
	$opt{jellyfish_bin} ? (bin => $opt{jellyfish_bin}."/jellyfish") : ()
	);

    $L->debug("Running jellyfish");

    # hashing collapses ident. kmers
    my %counts = $jf->query(["--sequence", $reads, $opt{kmer_hash}]); 

    my %hist;
    while(my ($k,$v) = each %counts){
	
	if($v){
	    $hist{$v}++;
	}else{
	    delete $counts{$k}
	}
    }

    my $total_count = 0;
    $total_count += $_ for values %hist;

    open(HIST, ">", $opt{out}."-cds-kmer-cov.tsv") or $L->logdie($!);
    
    my $cum_count = 0;
    my $med_cov;
    foreach (sort{$a<=>$b}keys %hist){
	$cum_count += $hist{$_};
	$med_cov = $_ if $cum_count < $total_count/2;
	$L->debug($med_cov," ", $cum_count," ",$total_count);
	print HIST $_,"\t",$hist{$_},"\n";
    }

    close HIST;

    my $cov = $med_cov * 1.1;  # kmer cov underestimates

    $L->debug("Estimated coverage of $cov based on ",scalar keys %counts," kmers");

    return $cov;

}



=head2 bam_coverage

Compute per contig coverage using bedtools and write to file.

=cut

sub bam_coverage{
    my ($pre) = @_;
    my $bam = $pre."-cds.bam";
    my $bed = $pre."-cds.bed";
    my $pre_sorted = $pre."-cds-sorted";
    my $bam_sorted = $pre."-cds-sorted.bam";
    my $tsv = $pre."-cds-bp-cov.tsv";
    my $head = '"id\tposition\tcoverage"';
    my $re;
    $re = qx/samtools sort $bam $pre_sorted/ && $L->warn($re);
    $re = qx/echo $head > $tsv/ && $L->warn($re);
    $re = qx/bedtools genomecov -d -ibam $bam_sorted -g $bed >> $tsv/ && $L->warn($re);
}

#-----------------------------------------------------------------------------#

=head1 AUTHOR

Clemens Weiss S<clemens.weiss@stud-mail.uni-wuerzburg.de>

Thomas Hackl S<thomas.hackl@uni-wuerzburg.de>

=cut



