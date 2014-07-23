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


#-----------------------------------------------------------------------------#
# Globals

our $VERSION = 0.02;

our $ID = 'fic';

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

# opt: multi-params need to be initiated with ARRAYREF!
my %opt = (
    config => [],
    );

# Setup defaults
my %def = (
    );

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


# debug level
$opt{quiet} && $L->level($WARN);
$opt{debug} && $L->level($DEBUG);

$L->debug(Dumper(\%opt));



##----------------------------------------------------------------------------##
# required	
for(qw(in blast_db blast_binary)){
    if(ref $opt{$_} eq 'ARRAY'){
	pod2usage("required: --$_") unless @{$opt{$_}}
    }else{
	pod2usage("required: --$_") unless defined ($opt{$_})
    }
};



#-----------------------------------------------------------------------------#
# MAIN

$opt{out} ||= basename($opt{in}, qw(.fasta .fas .fa .FASTA .FAS .FA)).".fil.fa";

my @blast_cmd  = ($opt{blast_binary} ? $opt{blast_binary}."/blastn" : "blastn",'-db', $opt{blast_db}, '-query', $opt{in}, @{$opt{blast_options}});
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



