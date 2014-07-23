#!/usr/bin/env perl

package asr;

=head1 NAME


=head1 DESCRIPTION


=head1 SYNOPSIS


=head1 OPTIONS

=over

=item -c|--config

Use user customized config file. Superseeds default config.

=item --debug

Verbose debug messages.

=item -V|--version

Display version.

=item -h|--help

Display this help.

=item 

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
use Cwd; # as we have to change directories
use Fasta::Parser;
use File::Spec;  # is a filename absolut or relative? Is needed for the output filename

#-----------------------------------------------------------------------------#
# Globals

our $VERSION = 0.02;

our $ID = 'asr';

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
    reads => [],
    mates => [],
    config => [],
    );

# Setup defaults
my %def = (
    );

GetOptions( # use %opt (Cfg) as defaults
	\%opt, qw(
		version|V!
		debug|D!
		help|h!
		config|c=s
                out|o=s
                outkmer=i{,}
		reads|1=s{,}
		mates|2=s{,}
		insert_size|insert-size|isize|s=i
		workingdir|wd|w=s
 		velveth_parameter=s
 		velvetg_parameter=s
 		velvet_out=s
		velvet_path=s
 		extendmode
		append!
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
for(qw(reads mates workingdir insert_size)){
    if(ref $opt{$_} eq 'ARRAY'){
	pod2usage("required: --$_") unless @{$opt{$_}}
    }else{
	pod2usage("required: --$_") unless defined ($opt{$_})
    }
};


#-----------------------------------------------------------------------------#
# MAIN

$L->info('Assemble reads');

my $oldpwd;
BEGIN{
    $oldpwd = cwd();
}
## define rechange of working directory on exit
END{
    if (defined $oldpwd) { chdir($oldpwd) }
}

$L->debug("Changing current directory to $opt{workingdir}");
chdir($opt{workingdir}) || $L->logdie("Unable to change to working directory $opt{workingdir}");

## run velveth command
$L->info("Running velveth");
my @cmd = (
    $opt{velvet_path} ?  $opt{velvet_path}."/velveth" : "velveth", 
    $opt{velvet_out}, 
    $opt{velveth_parameter});

foreach my $setnum (0..@{$opt{reads}}-1)
{
    push(@cmd, ("-fastq -shortPaired", $opt{reads}[$setnum], $opt{mates}[$setnum]));
}
$L->debug("Running velveth using the command '@cmd'");

qx(@cmd);

my $errorcode = $?;

if ($errorcode != 0)
{
    $L->logdie("The velveth run returned with errorcode $errorcode. The command was '@cmd'");
}

## search for all velvet_out directories
$L->info("Running velvetg");
my @dir_list = glob($opt{velvet_out}."*");
$L->debug("List of velvet-hash-directories: ".join(", ", @dir_list));

## run velvetg in each directory
foreach my $current_dir (@dir_list)
{
    @cmd = (
	$opt{velvet_path} ?  $opt{velvet_path}."/velvetg" : "velvetg", 
	$current_dir, 
	$opt{velvetg_parameter}, 
	"-ins_length", $opt{insert_size}
	);
    $L->debug("Running velvetg using the command '@cmd'");

    qx(@cmd);

    $errorcode = $?;

    if ($errorcode != 0)
    {
	$L->logdie("The velvetg run returned with errorcode $errorcode. The command was '@cmd'");
    }
}

## combine output from contigs.fa into a single fasta-file
## contigs have to be renamed according to their kmer size
$L->info("Generating output");
# check if the given output filename is relative or absolute
# in case of a relative filename it should be written into the oldpwd-directory
my $outputfasta=$opt{out};
if (File::Spec->file_name_is_absolute($outputfasta))
{
    $L->debug("Filename for output file is absolute");
} else 
{
    $outputfasta = $oldpwd."/".$outputfasta;
    $L->debug(sprintf("Filename for output file is realtive, assuming to put it in the original directory. New absolut filename is '%s'", $outputfasta));
}

# if append is set, the output file will be extended
my $mode = ($opt{append}) ? ">>" : ">";
my $output = Fasta::Parser->new(
    file => $outputfasta,
    mode => $mode
    );

# generate a list of kmers in output
my %wanted_kmer = map {($_, 0)} @{$opt{outkmer}};

## go through @dir_list and search contigs.fa files for combination
foreach my $current_dir (@dir_list)
{
    my ($kmersize) = $current_dir =~ /(\d+)$/;
    $L->debug(sprintf("Kmersize for directory '%s' was %d", $current_dir, $kmersize));

    if (! exists $wanted_kmer{$kmersize})
    {
	$L->debug(sprintf("Skipping kmersize %d from output", $kmersize));
	next;
    }

    # test if the file exists
    my $inputfilename = $current_dir."/contigs.fa";
    if (! -e $inputfilename)
    {
	$L->debug(sprintf("Skipping kmersize %d from output because the file '%s' does not exist.", $kmersize, $inputfilename));
	next;
    }

    my $input = Fasta::Parser->new(
	file => $inputfilename,
	mode => '<'
	);

    while (my $contig=$input->next_seq())
    {
	# add the kmer size to the contig_name
	$contig->id("Kmer".$kmersize."_".$contig->id());
	$output->append_seq($contig);
    }
}


#-----------------------------------------------------------------------------#

=head1 AUTHOR

Frank Foerster NAME S<frank.foerster@biozentrum.uni-wuerzburg.de>

=cut















