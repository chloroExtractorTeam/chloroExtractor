#!/usr/bin/env perl

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

use File::Spec;

use File::Basename;
use File::Copy;

# additional modules
use Cfg;
use Fasta::Parser;

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
my %opt = %{$cfg{mir}};


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
                input|in|i=s
	)
) or $L->logcroak('Failed to "GetOptions"');

# help
$opt{help} && pod2usage(1);

# version
if($opt{version}){
	print "$VERSION\n"; 
	exit 0;
}

#required stuff  
for(qw(input out)){
       pod2usage("required: --$_") unless defined ($opt{$_})
};

# debug level
$L->level($DEBUG) if $opt{debug};
$L->debug('Verbose level set to DEBUG');

$L->debug(Dumper(\%opt));



#-----------------------------------------------------------------------------#
# MAIN

$L->info('Make inverted repeat');

## first read the import file, count number of input sequences and their length for later usage

$L->info("Generating intermediate fasta file");

my $absinput = File::Spec->rel2abs($opt{input});
my $absintermediate = $absinput."_intermediate";
my $absoutput = File::Spec->rel2abs($opt{out});

my $input = Fasta::Parser->new(
    file => $absinput,
    mode => "<"
    );

my $intermediate = Fasta::Parser->new(
    file => $absintermediate,
    mode => ">"
    );

my %seq_names_length = ();
my %new_contig_name2old_name = ();

my $contigcounter = 0;

while (my $contig=$input->next_seq())
{
    # rename the contigs
    $contigcounter++;
    my $newname = sprintf("%s_%04d", "contig", $contigcounter);

    # store the new and old name and the length
    $new_contig_name2old_name{$newname} = $contig->id()." ".$contig->desc();
    $seq_names_length{$newname} = length($contig->seq());

    # rename contig id
    $contig->id($newname);
    $contig->desc("");

    # store the fasta-block in the new file
    $intermediate->append_seq($contig);
}

$L->info("Creating BLAST database");

my @cmd = (
    $opt{blast_path}."/makeblastdb", 
    '-in', $absintermediate, 
    '-dbtype', 'nucl'
    );

$L->debug("Running makeblastdb using the command '@cmd'");

qx(@cmd);

my $errorcode = $?;

if ($errorcode != 0)
{
    $L->logdie("The makeblastdb run returned with errorcode $errorcode. The command was '@cmd'");
}

$L->info("Running BLAST search");

@cmd = (
    $opt{blast_path}."/blastn", 
    '-query', $absintermediate,
    '-task', 'blastn',
    '-db', $absintermediate,
    '-evalue', '1e-30',
    '-outfmt', '"6 sseqid sstart send qseqid qstart qend"',
    '-num_threads', '1'
    );	
		
$L->debug("Running blastn using the command '@cmd'");

open(BLAST, join(" ", @cmd)."|") || $L->logdie("Unable to run BLAST using command '@cmd'");

# for storage of 
my @irRegions = ();

while (my $blastline = <BLAST>)
{
    chomp($blastline);
    my ($sseqid, $sstart, $send, $qseqid, $qstart, $qend) = split(/\t/, $blastline);

    ### used from Markus code
    # only IRs on the same contig/chromosom are detected 
    next unless($sseqid eq $qseqid);
    # only rc matches are considered
    next unless(($sstart < $send) != ($qstart < $qend));
    push(@irRegions, {'id'=>$sseqid,'start1'=>($sstart>$send ? $send : $sstart),'end1'=>($sstart>$send ? $sstart : $send),
		      'start2'=>($qstart>$qend ? $qend : $qstart),'end2'=>($qstart>$qend ? $qstart : $qend)});
    ### end of Markus code
}
close(BLAST) || $L->logdie("Unable to finish BLAST run using command '@cmd'");



# 

































#-----------------------------------------------------------------------------#

=head1 AUTHOR

NAME S<MAIL>

=cut















