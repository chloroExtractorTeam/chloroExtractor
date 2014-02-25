#!/usr/bin/env perl

=head1 NAME

filter_contigs.pl

=head1 DESCRIPTION

Filter contigs using length cutoff and Blast search against CDS.

=head1 SYNOPSIS

  $ perl filter_contigs.pl <OPTIONS> -i contigs.fa -d cds.fa

=head1 OPTIONS

=over

=item -i|--in <FASTA FILE(S)>

Files containing sequences to be filtered.

=item -o|--out <PATHNAME> [basename "--in".fil.fa]

Output file name.

=item -d|--db|--blast-db <BLASTDB>

Reference blastn database sequences

=item --blast-binary

Path to blastn binary. Required unless exported.

=item --blast-options

Pass-through parameter for blast search.
NOTE: enclose in "".

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
	log4perl.appender.Screen.layout.ConversionPattern = [%d{yy-MM-dd HH:mm:ss}] [fic] %m%n
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

my %opt = %{$cfg{fic}};



#-----------------------------------------------------------------------------#
# GetOptions

GetOptions( # use %opt (Cfg) as defaults
	\%opt, qw(
                in|i=s
                out|o=s
                blast_db|blast-db|db|d=s
                blast_binary|blast-binary=s
		blast_options|blast-options=s@
                seqfilter_options=s@
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
for(qw(in blast_db blast_binary)){
        pod2usage("required: --$_") unless defined ($opt{$_})
};


# debug level
$L->level($DEBUG) if $opt{debug};
$L->debug('Verbose level set to DEBUG');

$L->debug(Dumper(\%opt));



#-----------------------------------------------------------------------------#
# MAIN

$opt{out} ||= basename($opt{in}, qw(.fasta .fas .fa .FASTA .FAS .FA)).".fil.fa";

my @blast_cmd  = ($opt{blast_binary},'-db', $opt{blast_db}, '-query', $opt{in}, @{$opt{blast_options}});
$L->debug("Running: @blast_cmd");
my @ids = qx(@blast_cmd);

$L->debug("Passing filter:\n @ids");

my @seqfilter_cmd = ($RealBin."/SeqFilter", $opt{in}, @{$opt{seqfilter_options}}, qw(--ids - --out), $opt{out});
$L->debug("Running: @seqfilter_cmd");

open(SF, "|-", @seqfilter_cmd);
print SF @ids;
close SF;

#-----------------------------------------------------------------------------#

=head1 AUTHOR

Thomas Hackl S<thomas.hackl@uni-wuerzburg.de>

=cut



