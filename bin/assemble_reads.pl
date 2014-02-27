#!/usr/bin/env perl

package asr;

=head1 NAME


=head1 DESCRIPTION


=head1 SYNOPSIS


=head1 OPTIONS

=over

=item --create-config

Create a config file with default settings for user customization.

=item -c|--config

Use user customized config file. Superseeds default config.

=item --debug

Verbose debug messages.

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
use Cwd; # as we have to change directories
use Fasta::Parser;
use File::Spec;  # is a filename absolut or relative? Is needed for the output filename

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
	log4perl.appender.Screen.layout.ConversionPattern = [%d{yy-MM-dd HH:mm:ss}] [%C] %m%n
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
my %opt = %{$cfg{asr}};

#TODO: custom config


#-----------------------------------------------------------------------------#
# GetOptions

GetOptions( # use %opt (Cfg) as defaults
	\%opt, qw(
		version|V!
		debug|D!
		help|h!
		config|c=s
                out|o=s
                outkmer=i@
		reads|1=s
		mates|2=s
		insert_size|insert-size|isize|s=i
		workingdir|wd|w=s
 		velvetparameter=s
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
for(qw(workingdir insert_size mates reads)){
        pod2usage("required: --$_") unless defined ($opt{$_})
};

# debug level
$L->level($DEBUG) if $opt{debug};
$L->debug('Verbose level set to DEBUG');

$L->debug(Dumper(\%opt));



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
my $cmd = join(" ", ($opt{velvet_path}."/velveth", $opt{velvet_out}, $opt{velveth_parameter}, "-fastq -shortPaired", $opt{reads}, $opt{mates}));
$L->debug("Running velveth using the command '$cmd'");

qx($cmd);

my $errorcode = $?;

if ($errorcode != 0)
{
    $L->logdie("The velveth run returned with errorcode $errorcode");
}

## search for all velvet_out directories
$L->info("Running velvetg");
my @dir_list = glob($opt{velvet_out}."*");
$L->debug("List of velvet-hash-directories: ".join(", ", @dir_list));

## run velvetg in each directory
foreach my $current_dir (@dir_list)
{
    $cmd = join(" ", ($opt{velvet_path}."/velvetg", $current_dir, $opt{velvetg_parameter}, "-ins_length", $opt{insert_size}));
    $L->debug("Running velvetg using the command '$cmd'");

    qx($cmd);

    $errorcode = $?;

    if ($errorcode != 0)
    {
	$L->logdie("The velvetg run returned with errorcode $errorcode. The command was '$cmd'");
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

my $output = Fasta::Parser->new(
    file => $outputfasta,
    mode => '>'                    # overwrite an existing file
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















