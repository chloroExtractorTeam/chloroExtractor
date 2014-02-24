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

#### first generate a FASTA file for 5' and 3' ends from input file
my $fasta_in = Fasta::Parser->new(
    file => $opt{in}
    );
my $fasta_out = Fasta::Parser->new(
    file => $opt{out},
    mode => '>'                    # overwrite an existing file
    );
# loop through the FASTA file and generate a set of contig ends
while (my $contig=$fasta_in->next_seq())
{
    my $start = Fasta::Seq->new(
	id => $contig->id()."_5prime",
	seq => substr($contig->seq(), 0, $opt{border})
	);
    $fasta_out->append_seq($start);
    my $end = Fasta::Seq->new(
	id => $contig->id()."_3prime",
	seq => substr($contig->seq(), -$opt{border}),
	);
    $fasta_out->append_seq($end);
}






































#-----------------------------------------------------------------------------#

=head1 AUTHOR

NAME S<MAIL>

=cut















