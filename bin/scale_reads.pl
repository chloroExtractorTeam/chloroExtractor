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

#enable in script mapping
use Bowtie2;


#-----------------------------------------------------------------------------#
# Globals

our $VERSION = 0.01;

our $ID = 'scr';

# get a logger
my $L = Log::Log4perl::get_logger();
Log::Log4perl->init( \(q(
	log4perl.rootLogger                     = INFO, Screen
	log4perl.appender.Screen                = Log::Log4perl::Appender::Screen
	log4perl.appender.Screen.stderr         = 1
	log4perl.appender.Screen.layout         = PatternLayout
	log4perl.appender.Screen.layout.ConversionPattern = [%d{yy-MM-dd HH:mm:ss}] [).$ID.q(] %m%n
)));


#-----------------------------------------------------------------------------#
# GetOptions

my %opt = (config => []);

GetOptions( # use %opt (Cfg) as defaults
	\%opt, qw(
                target_coverage|target-coverage|coverage=i
                reads|1=s
                mates|2=s
		ref_cluster|ref-cluster|r=s
		out|o=s
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

##----------------------------------------------------------------------------##
# Config

my %cfg;

# core
my $core_cfg = "$RealBin/../".basename($Script, qw(.pl)).".cfg";

if( -e $core_cfg){
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
# required stuff  
for(qw(reads mates ref_cluster)){
    pod2usage("required: --$_") unless defined ($opt{$_}) || $opt{$_} eq ""
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
    path => $opt{bowtie2_path},
);

# TODO: bowtie2 generate db -> prevent issues with different indices on different architectures or bowtie2 versions

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
    "-p" => $opt{threads} || 1,
);

$L->debug(Dumper($opt{bowtie2_params}));

# read output on the fly
my $sp = Sam::Parser->new(
  fh => $bowtie2->stdout
);

my $bam = $opt{out}.".bam";
open(BAM, "| samtools view -Sbu /dev/fd/0 > $bam") or $L->logdie($!);
open(BED, ">", $opt{out}.".bed") or $L->logdie($!);


my %h; 
my $ss;

my %seqs;

while(%h = $sp->next_header_line('@SQ')){
  $seqs{$h{'SN'}} = Sam::Seq->new(
    id => $h{'SN'},
    len => $h{'LN'},
  );  

  print BAM '@SQ'."\tSN:$h{'SN'}\tLN:$h{'LN'}\n";
  print BED "$h{'SN'}\t$h{'LN'}\n";

}

close BED;

my $current_cov;
my $last_id;
my $closest_ref;

my $c;
my $tlen_sum; 
my $seqwithtlen;


while(my $aln = $sp->next_aln()){
    print BAM "$aln\n";
my $id = $aln->rname();
if (abs($aln->tlen) > 0){
    $tlen_sum+= abs($aln->tlen); # compute isize
    $seqwithtlen++;
}
$L->logdie($id, ": $aln") unless exists $seqs{$id};
$seqs{$id}->add_aln($aln);

$c++;
unless ($c % $opt{coverage_check_interval}){
    $current_cov = estimate_coverage(\%seqs);
    $last_id = $aln->qname;
    if ($current_cov >= $opt{target_coverage}) {
	
	$bowtie2->cancel();
	close BAM;


	my %refwise_coverage;

	foreach my $prot_id (keys %seqs){
	    my ($ref) = $prot_id =~ /^([^_]+_[^_]+)_/;
	    $refwise_coverage{$ref}+= median([$seqs{$prot_id}->coverage]) || 0;
	    # $L->debug(Dumper(\%seqs, \%refwise_coverage));
	}
	$closest_ref = (sort{$refwise_coverage{$b} <=> $refwise_coverage{$a}}keys %refwise_coverage)[0];
	last;
    }
}
}


print "coverage\t", $current_cov, "\n";
print "insert_size\t", int($tlen_sum/$seqwithtlen), "\n";
print "closest_ref\t", $closest_ref || 'NA', "\n";


#what if not enough coverage in entire file.
if(! $current_cov || $current_cov < $opt{target_coverage}){
    close BAM;

    #$L->info("Recalculating accurate per base coverages");
    #bam_coverage($opt{out});

    $current_cov 
	? $L->info("Ran out of reads at estimated coverage of $current_cov")
	: $L->info("Could not detect chloroplast reads in input data");
    $L->info("You might need to increase the amount of input data");
    $L->info("Also make sure, your library contains chloroplast sequences");

    exit 1;
}


#$L->info("Recalculating accurate per base coverages");
#bam_coverage($opt{out});


$L->info("Creating libraries");

$last_id =~ s?/[12]$??;

my $sed_cmd1 = 'sed "1~4 {/^\@'.$last_id.'\(\/[12]\)*\s/{N;N;N;q}}" '.$opt{reads}." >".$opt_o1;

$L->debug($sed_cmd1);
qx($sed_cmd1);

if ($opt{mates}) {

  my $sed_cmd2 = 'sed "1~4 {/^\@'.$last_id.'\(\/[12]\)*\s/{N;N;N;q}}" '.$opt{mates}." >".$opt_o2;

  $L->debug($sed_cmd2);
  qx($sed_cmd2);

}










#-- SUBS ----------------------------------------------------------------------#



sub pairwise_sum {
  my @array1 = @{$_[0]};
  my @array2 = @{$_[1]};

  my @len;
  push @len, scalar @array1;
  push @len, scalar @array2;

  my @sort = sort {$a <=> $b} @len;

  my @added;

  for (my $i=0;$i<=$sort[-1];$i++) {

    if (exists $array1[$i] && exists $array2[$i]) {
      $added[$i] = $array1[$i] + $array2[$i];
    }
    elsif (exists $array1[$i]) {
      $added[$i] = $array1[$i];
    }
    elsif (exists $array2[$i]) {
      $added[$i] = $array2[$i];
    }
  }
  return @added;
}

sub median {
  my @array = @{$_[0]};
  my $median = (sort{$a<=>$b}@array)[@array/2];
  return $median;
}


sub estimate_coverage {

  # get coverages for complete sam file
  # and build hash with protein name as key
  # and protein wise sum of coverage in
  # array

  open(COV, ">", $opt{out}."-cov.tsv") or $L->logdie($!);

  my %protein_wise_coverage;
  for my $ss(values %{$_[0]}){
    my @covs=$ss->coverage();
    my $protein = (split /_/, $ss->id())[-1];
    @{$protein_wise_coverage{$protein}} = pairwise_sum(\@{$protein_wise_coverage{$protein}}, \@covs);
  }

  # print proteinwise coverage // omitting 
  # zeros and empty proteins // also trim
  # both ends

  my %medians;
  my @median_array;

  for (keys%protein_wise_coverage) {
    @{$protein_wise_coverage{$_}} = grep { $_ != 0 } @{$protein_wise_coverage{$_}};
    if (@{$protein_wise_coverage{$_}}-100>200) {
      @{$protein_wise_coverage{$_}} = @{$protein_wise_coverage{$_}}[100..@{$protein_wise_coverage{$_}}-100];
    }
    else { # ignore to short CDS
      @{$protein_wise_coverage{$_}} = (); 
    }
    if (@{$protein_wise_coverage{$_}}) {
      my $prot_median = median(\@{$protein_wise_coverage{$_}});

      for my $cov (@{$protein_wise_coverage{$_}}){ 
	  print COV $_,"\t",$cov,"\n"; 
      }

      $medians{$_} = $prot_median;
      push @median_array, $prot_median if ($prot_median > $opt{min_single_coverage});
    }
  }

  close COV;

  $L->debug(Dumper(\%medians));
  my $median_cov = @median_array > $opt{min_covered_CDS} 
  	? median(\@median_array)
	    : -1;


  $L->info("Cov: ", $median_cov < 0 ? 'NA' : $median_cov, " Prots: ", scalar @median_array, " [@median_array]");
  
  return $median_cov;

}



=head2 bam_coverage

Compute per contig coverage using bedtools and write to file.

=cut

sub bam_coverage{
    my ($pre) = @_;
    my $bam = $pre.".bam";
    my $bed = $pre.".bed";
    my $pre_sorted = $pre."-sorted";
    my $bam_sorted = $pre."-sorted.bam";
    my $tsv = $pre."-cov.tsv";
    my $head = '"id\tposition\tcoverage"';
    $L->info(qx/samtools sort $bam $pre_sorted/);
    $L->info(qx/echo $head > $tsv/);
    $L->info(qx/bedtools genomecov -d -ibam $bam_sorted -g $bed >> $tsv/);
}

#-----------------------------------------------------------------------------#

=head1 AUTHOR

Clemens Weiss S<clemens.weiss@stud-mail.uni-wuerzburg.de>

Thomas Hackl S<thomas.hackl@uni-wuerzburg.de>

=cut



