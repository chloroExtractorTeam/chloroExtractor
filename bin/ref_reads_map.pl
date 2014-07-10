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

use FindBin qw($RealBin);
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

our $VERSION = 0.01;

# get a logger
my $L = Log::Log4perl::get_logger();
Log::Log4perl->init( \q(
	log4perl.rootLogger                     = DEBUG, Screen
	log4perl.appender.Screen                = Log::Log4perl::Appender::Screen
	log4perl.appender.Screen.stderr         = 1
	log4perl.appender.Screen.layout         = PatternLayout
	log4perl.appender.Screen.layout.ConversionPattern = [%d{yy-MM-dd HH:mm:ss}] [rrm] %m%n
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

my %opt = %{$cfg{rrm}};



#-----------------------------------------------------------------------------#
# GetOptions

GetOptions( # use %opt (Cfg) as defaults
	\%opt, qw(
                reads|1=s
                mates|2=s
		reference|reference|r=s
		out|o=s
		version|V!
		debug|D!
		help|h!
		config|c=s
		bowtie2_path|bowtie2-path=s
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
for(qw(reads mates reference)){
        pod2usage("required: --$_") unless defined ($opt{$_})
};


# debug level
$L->level($DEBUG) if $opt{debug};
$L->debug('Verbose level set to DEBUG');

$L->debug(Dumper(\%opt));


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
    "-p" => $cfg{threads},
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



