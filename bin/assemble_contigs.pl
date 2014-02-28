#!/usr/bin/env perl

package asc;

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
my %opt = %{$cfg{asc}};

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
		workingdir|wd|w=s
		in|i=s
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
for(qw(workingdir in out)){
        pod2usage("required: --$_") unless defined ($opt{$_})
};

# debug level
$L->level($DEBUG) if $opt{debug};
$L->debug('Verbose level set to DEBUG');

$L->debug(Dumper(\%opt));



#-----------------------------------------------------------------------------#
# MAIN

$L->info('Assemble contigs');

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

$L->info("Renaming fasta sequences");

my $input = Fasta::Parser->new(
    file => $opt{in},
    mode => '<'
    );

my $renamedfile = $opt{in}."_renamed";
my $renamed = Fasta::Parser->new(
    file => $renamedfile,
    mode => '>'
    );

my $counter = 0;
while (my $contig=$input->next_seq())
{
    # substitute all period (.) with underscore (_) and add .a to the end of the line
    my $id = $contig->id();
    $id =~ s/\./_/g;
    # count all contigs to garantee unique names
    $counter++;
    $contig->id(sprintf("%05d_%s.a", $counter, $id));
    $renamed->append_seq($contig);
}

$L->info("Running phrap...");
my $cmd = $opt{phrap_path}."/phrap ".$renamedfile;
$L->debug("Run command: '$cmd'");
qx($cmd);

my $errorcode = $?;

if ($errorcode != 0)
{
    $L->logdie("The phrap run returned with errorcode $errorcode. The command was '$cmd'");
}

## last step is copying the output file
my $outputfasta=$opt{out};
if (File::Spec->file_name_is_absolute($outputfasta))
{
    $L->debug("Filename for output file is absolute");
} else 
{
    $outputfasta = $oldpwd."/".$outputfasta;
    $L->debug(sprintf("Filename for output file is realtive, assuming to put it in the original directory. New absolut filename is '%s'", $outputfasta));
}
copy($renamedfile.".contigs", $outputfasta) || $L->logdie(sprintf("Copying the result ('%s') into the output file ('%s') failed: %s", $renamedfile.".contigs", $outputfasta, $!);

#-----------------------------------------------------------------------------#

=head1 AUTHOR

Frank Foerster NAME S<frank.foerster@biozentrum.uni-wuerzburg.de>

=cut















