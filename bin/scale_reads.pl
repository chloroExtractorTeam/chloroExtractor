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
use Getopt::Long qw(:config no_ignore_case bundling);
use Pod::Usage;
use Log::Log4perl qw(:no_extra_logdie_message);
use Log::Log4perl::Level;

use Data::Dumper;
$Data::Dumper::Sortkeys = 1;

use FindBin qw($RealBin);
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

# get a logger
my $L = Log::Log4perl::get_logger();
Log::Log4perl->init( \q(
	log4perl.rootLogger                     = INFO, Screen
	log4perl.appender.Screen                = Log::Log4perl::Appender::Screen
	log4perl.appender.Screen.stderr         = 1
	log4perl.appender.Screen.layout         = PatternLayout
	log4perl.appender.Screen.layout.ConversionPattern = [%d{yy-MM-dd HH:mm:ss}] [scr] %m%n
));


#-----------------------------------------------------------------------------#
# Config

# core
my $core_cfg = "$RealBin/../chloroExtractor.cfg";
my %cfg = Cfg->Read_Cfg($core_cfg); 

# user defaults and overwrite core
my $user_cfg;
for(my $i=0; $i<@ARGV; $i++){
        if($ARGV[$i] =~ /-c$|--config$/){
                $user_cfg = $ARGV[$i+1];
                last;
        }
}

%cfg = (%cfg, Cfg->Read_Cfg($user_cfg)) if $user_cfg; # simple overwrite

my %opt = %{$cfg{scr}};



#-----------------------------------------------------------------------------#
# GetOptions

GetOptions( # use %opt (Cfg) as defaults
	\%opt, qw(
                target_coverage|coverage=i
                reads|1=s
                mates|2=s
		ref_cluster|r=s
		version|V!
		debug|D!
		help|h!
		config|c=s
	)
) or $L->logcroak('Failed to "GetOptions"');

# help
$opt{help} && pod2usage(1);

# version
if($opt{version}){
	print "$VERSION\n"; 
	exit 0;
}


# required stuff  
for(qw(reads mates ref_cluster)){
        pod2usage("required: --$_") unless defined ($opt{$_})
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
    log => $opt{bowtie2_log},
);

# TODO: bowtie2 generate db -> prevent issues with different indices on different architectures or bowtie2 versions

unless(-e $opt{ref_cluster}.'.1.bt2'){
    $L->info("Building bowtie2 index");
    $bowtie2->bowtie2_build(
	ref => $opt{ref_cluster},
	pre => $opt{ref_cluster},
	);
    my $bis = $bowtie2->stdout;
    my $bie = $bowtie2->stderr;

    $L->debug(<$bis>, <$bie>);

    $bowtie2->finish();
}else{
    $L->info("Using existing bowtie2 index $opt{ref_cluster}.*.bt2");
}

$L->info("Running bowtie2");

$bowtie2->bowtie2(
    %{$opt{bowtie2_params}},
    ref => undef,
    "-1" => $opt{reads},
    "-2" => $opt{mates},
    "-x" => $opt{ref_cluster},
    "-p" => $cfg{threads},
);

## DEPRECATED: mapping paired
# WARNING: Bowtie required "-U" twice for paired; not possible in option hash
# Therefore: if input mates option exists, set singleorpair to 2 to map only reads

# read output on the fly
my $sp = Sam::Parser->new(
  fh => $bowtie2->stdout
);

my %h; 
my $ss;

my %seqs;

while(%h = $sp->next_header_line('@SQ')){
  $seqs{$h{'SN'}} = Sam::Seq->new(
    id => $h{'SN'},
    len => $h{'LN'},
  );  
  #print '@SQ'."\tSN:$h{'SN'}\tLN:$h{'LN'}\n" if ($opt{debug});
}

my $current_cov;
my $last_id;

my $c;
my $tlen_sum; 

while(my $aln = $sp->next_aln){
    my $id = $aln->rname();
    $tlen_sum+= abs($aln->tlen); # compute isize
    $seqs{$id}->add_aln($aln);

    $c++;
    unless ($c % $opt{coverage_check_interval}){
	$current_cov = estimate_coverage(\%seqs);
	$L->debug("Coverage: ",$current_cov);
	$last_id = $aln->qname;
	if ($current_cov >= $opt{target_coverage}) { 
	    last;
	}
    }
}

#what if not enough coverage in entire file.
if(! $current_cov || $current_cov < $opt{target_coverage}){
    $current_cov 
	? $L->info("Ran out of reads at estimated coverage of $current_cov")
	: $L->info("Could not detect chloroplast reads in input data");
    $L->info("You might need to increase the amount of input data");
    $L->info("Also make sure, your library contains chloroplast sequences");
    exit 1;
}

$bowtie2->cancel();

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

print "Estimated coverage: ", $current_cov, "\n";
print "Estimated insert size: ", int($tlen_sum/$c), "\n";








# ----- SUBS ------ #



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
    else {
      @{$protein_wise_coverage{$_}} = ();
    }
    if (@{$protein_wise_coverage{$_}}) {
      my $prot_median = median(\@{$protein_wise_coverage{$_}});
      $medians{$_} = $prot_median;
      push @median_array, $prot_median if ($prot_median > $opt{min_single_coverage});
    }
  }

  $L->debug("Coverage Array ","@median_array"," -- ",scalar @median_array);

  if (@median_array > $opt{min_covered_CDS}) {
    return median(\@median_array);
  }
  else {
    return -1;
  }

}



#-----------------------------------------------------------------------------#

=head1 AUTHOR

Clemens Weiss S<clemens.weiss@stud-mail.uni-wuerzburg.de>

Thomas Hackl S<thomas.hackl@uni-wuerzburg.de>

=cut



