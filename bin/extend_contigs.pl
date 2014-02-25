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

use File::Basename;
use File::Copy;

# additional modules
use Cfg;
use Fasta::Parser;
use Fastq::Parser;
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
my %opt = %{$cfg{exc}};

#TODO: custom config


#-----------------------------------------------------------------------------#
# GetOptions

GetOptions( # use %opt (Cfg) as defaults
	\%opt, qw(
		version|V!
		debug|D!
		help|h!
		config|c=s
		reads|1=s
		mates|2=s
		insert_size|insert-size|isize|s=i
		border=i
		in|i=s
		out|o=s
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
for(qw(in out insert_size mates reads)){
        pod2usage("required: --$_") unless defined ($opt{$_})
};

## check if the border was set, otherwise use 2x insert_size as default value
unless (exists $opt{border} && $opt{border} > 0)
{
    $opt{border} = $opt{insert_size}*2;
}

# debug level
$L->level($DEBUG) if $opt{debug};
$L->debug('Verbose level set to DEBUG');

$L->debug(Dumper(\%opt));



#-----------------------------------------------------------------------------#
# MAIN

$L->info('Extend contig script');

$L->info("Generating 5' and 3' ends of contigs");

#### first generate a FASTA file for 5' and 3' ends from input file
my $fasta_in = Fasta::Parser->new(
    file => $opt{in}
    );
my $fasta_contig_ends=$opt{in}.'_contig_ends';
my $fasta_out = Fasta::Parser->new(
    file => $fasta_contig_ends,
    mode => '>'                    # overwrite an existing file
    );
# loop through the FASTA file and generate a set of contig ends
my %contig_ends_seen = ();
my %filehandles = ();
while (my $contig=$fasta_in->next_seq())
{
    # filter for a minimum length
    next unless (length($contig->seq()) > $opt{min_seq_length});

    # store 5' end
    store_sequence_and_create_folder(name => $contig->id()."_5prime", seq => substr($contig->seq(), 0, $opt{border}), file_out => $fasta_out, seen_names => \%contig_ends_seen, filehandle => \%filehandles);
    # store 3' end
    store_sequence_and_create_folder(name => $contig->id()."_3prime", seq => substr($contig->seq(), -$opt{border}), file_out => $fasta_out, seen_names => \%contig_ends_seen, filehandle => \%filehandles);
}

$L->info("Building bowtie2 index");

my $bowtie2_db = $opt{in}.'_bowtie2_db';
my $bowtie2 = Bowtie2->new(
    path => $opt{bowtie2_path},
    log => $opt{bowtie2_log},
    ref => $fasta_contig_ends,
    pre => $bowtie2_db,
);

# TODO: bowtie2 generate db -> prevent issues with different indices on different architectures or bowtie2 versions
$bowtie2->bowtie2_build();
my $bis = $bowtie2->stdout;
my $bie = $bowtie2->stderr;

$L->debug(<$bis>, <$bie>);

$bowtie2->finish();

























sub store_sequence_and_create_folder
{
    my %params = @_;
    my $name = $params{name};
    if (exists $params{seen_names}{$name})
    {
	$L->logdie("The folder '$name' already exists! One possibility are multiple occurence of the same sequence name in the input file");
    }
    $params{seen_names}{$name}++;
    mkdir($name) || $L->logdie("Unable to create folder '$name'");

    $params{filehandle}{$params{name}}{reads}=Fasta::Parser->new(
	file => $name.'/'.'reads.fq',
	mode => '>'                    # overwrite an existing file
	);
    $params{filehandle}{$params{name}}{mates}=Fasta::Parser->new(
	file => $name.'/'.'mates.fq',
	mode => '>'                    # overwrite an existing file
	);

    my $seq_obj = Fasta::Seq->new(
	id => $params{name},
	seq => $params{seq}
	);
    $params{file_out}->append_seq($seq_obj);
    return 1;
}

#-----------------------------------------------------------------------------#

=head1 AUTHOR

Frank Foerster NAME S<frank.foerster@biozentrum.uni-wuerzburg.de>

=cut















