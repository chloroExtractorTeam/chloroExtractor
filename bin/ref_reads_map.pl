#!/usr/bin/env perl

=head1 NAME

scale_reads.pl

=head1 DESCRIPTION


=head1 SYNOPSIS

  $ perl ref_reads_map.pl -1 lib_1.fq -2 lib_2.fq --ref ref-genome.fa -o <out-prefix>

=head1 OPTIONS

=over

=item -1|--reads

Input reads file, first of pair.

=item -2|--mates

Input reads file, second of pair

=item -o|--out

Output prefix.

=item -r|--reference

Reference sequences for mapping. Fasta or already created bowtie2 index prefix.

=item -c|--config

Use user customized config file. Superseeds default config.

=item --bowtie2-path=<DIR> []

Path to bowtie2 binaries.

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

use FindBin qw($RealBin $Script);
use lib "$RealBin/../lib/";

use File::Basename;
use File::Copy;

# additional modules
use Cfg;

use Fasta::Parser;

use Sam::Parser;
use Sam::Seq;
use Sam::Alignment ':flags';

#enable in script mapping
use Bowtie2;


#-----------------------------------------------------------------------------#
# Globals

our $VERSION = 0.02;

our $ID = 'rrm';

# get a logger
my $L = Log::Log4perl::get_logger();
Log::Log4perl->init( \(q(
	log4perl.rootLogger                     = INFO, Screen
	log4perl.appender.Screen                = Log::Log4perl::Appender::Screen
	log4perl.appender.Screen.stderr         = 1
	log4perl.appender.Screen.layout         = PatternLayout
	log4perl.appender.Screen.layout.ConversionPattern = [%d{yy-MM-dd HH:mm:ss}] [).$ID.q(] %m%n
)));

# opt: multi-params need to be initiated with ARRAYREF!
my %opt = (
#    reads => [], # multi-reads not yet supported
#    mates => [],
    config => [],
    );

# Setup defaults
my %def = (
    );


#-----------------------------------------------------------------------------#
# GetOptions

GetOptions( # use %opt (Cfg) as defaults
	\%opt, qw(
                reads|1=s
                mates|2=s
		reference|reference|r=s
		out|o=s
		bowtie2_path|bowtie2-path=s
		threads|t=i
		version|V!
		debug|D!
		quiet|Q!
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


# debug level
$opt{quiet} && $L->level($WARN);
$opt{debug} && $L->level($DEBUG);


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
%opt = (%def, %cfg, %opt);

$L->debug("GetOptions:\n", Dumper(\%opt));



##------------------------------------------------------------------------##	
# required	
for(qw(reads mates reference)){
    if(ref $opt{$_} eq 'ARRAY'){
	pod2usage("required: --$_") unless @{$opt{$_}}
    }else{
	pod2usage("required: --$_") unless defined ($opt{$_})
    }
};


#------------------------------------------------------------------------------#
# MAIN


my $bowtie2 = Bowtie2->new(
    path => $opt{bowtie2_path},
);

# TODO: bowtie2 generate db -> prevent issues with different indices on different architectures or bowtie2 versions

unless(-e $opt{reference}.'.1.bt2'){
    $L->info("Building bowtie2 index");
    $bowtie2->build($opt{reference})
}else{
    $L->info("Using existing bowtie2 index $opt{reference}.*.bt2");
}

$L->info("Running bowtie2");

$bowtie2->run(
    @{$opt{bowtie2_params}},
    "-1" => $opt{reads},
    "-2" => $opt{mates},
    "-x" => $opt{reference},
    "-p" => $opt{threads} || 1,
    "-S" => $opt{out}.'.sam',
    )->finish();

my $fp = Fasta::Parser->new(file => $opt{reference});

open(BED, '>', "rrm.bed") or $L->logdie($!);
while(my $fa=$fp->next_seq){
    print BED $fa->id,"\t",length($fa->seq),"\n";
}
close BED;

$L->info("Creating and sorting bam");
my $sam = $opt{out}.'.sam';
my $bampre = $opt{out};
$bampre =~ s/\.sam$/.sorted/;
qx(samtools view -Sbu $sam | samtools sort /dev/fd/0 $bampre);

$L->info("Computing coverage");
my $bam = $bampre.'.bam';
my $bed = $opt{out}.'.bed';
my $tsv = $opt{out}.'-cov.tsv';
my $head = '"id\tcoverage\tfrequency\tlength\tfraction"';
$L->info(qx/echo $head > $tsv/);
$L->info(qx/genomeCoverageBed -ibam $bam -g $bed >> $tsv/);



#Rplot(); # if ..


#-- SUBS ----------------------------------------------------------------------#



#------------------------------------------------------------------------------#

=head1 AUTHOR

Thomas Hackl S<thomas.hackl@uni-wuerzburg.de>

=cut



