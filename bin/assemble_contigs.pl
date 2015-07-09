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

our $ID = 'asc';

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
    in => [],
    );

# Setup defaults
my %def = (
    );

GetOptions( # use %opt (Cfg) as defaults
	\%opt, qw(
		version|V!
		debug|D!
		help|h!
		config|c=s{,}
                out|o=s
		workingdir|wd|w=s
		overwrite!
		phrap_path=s
		in|i=s{,}
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
for(qw(in out workingdir)){
    if(ref $opt{$_} eq 'ARRAY'){
	pod2usage("required: --$_") unless @{$opt{$_}}
    }else{
	pod2usage("required: --$_") unless defined ($opt{$_})
    }
};




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

my $renamedfile = join("_", @{$opt{in}})."_renamed";
### check if the file exists and create a warning in this case!
if (-e $renamedfile)
{
    if (! $opt{overwrite})
    {
	$L->logdie("The file '$renamedfile' already exists! Use --overwrite in case you want to overwrite the file");
    } else {
	$L->warn("The file '$renamedfile' already exists and will be overwritten");
    }
}
my $renamed = Fasta::Parser->new(
    file => $renamedfile,
    mode => '>'
    );

my $counter = 0;

foreach my $inputfile (@{$opt{in}})
{
    my $input = Fasta::Parser->new(
	file => $inputfile,
	mode => '<'
    );

    while (my $contig=$input->next_seq())
    {
	# substitute all period (.) with underscore (_) and add .a to the end of the line
	my $id = $contig->id();
	$id =~ s/\./_/g;
	$id = $id.".a";
	# count all contigs to garantee unique names
	$counter++;
	$contig->id(sprintf("%05d_%s", $counter, $id));
	$renamed->append_seq($contig);
    }
}

$L->info("Running phrap...");
my @cmd = (
    $opt{phrap_path} ? $opt{phrap_path}."/phrap" : "phrap", 
    "-retain_duplicates", 
    $renamedfile
    );
$L->debug("Run command: '@cmd'");
qx(@cmd);

my $errorcode = $?;

if ($errorcode != 0)
{
    $L->logdie("The phrap run returned with errorcode $errorcode. The command was '@cmd'");
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


open(my $FH1, "<", $renamedfile.".contigs") || $L->logdie(sprintf("Opening the result ('%s') failed: %s", $renamedfile.".contigs", $!));
open(my $OUT, ">", $outputfasta) || $L->logdie(sprintf("Writing into the output file ('%s') failed: %s", $outputfasta, $!));
while ( <$FH1> )
{
    print $OUT $_;
}
close $FH1;
open my $FH2, "<", $renamedfile.".singlets" || $L->logdie(sprintf("Opening the result ('%s') failed: %s", $renamedfile.".singlets", $!));
my $n = 1;
while ( <$FH2> )
{
    if ( $_ =~ /^>.+/ )
    {
	print $OUT ">ass/contigs.fasta_renamed.sContig".$n."\n";
	$n++;
    }
    else
    {
	print $OUT $_;
    }
}
close $FH2;
close $OUT;

#copy($renamedfile.".contigs", $outputfasta) || $L->logdie(sprintf("Copying the result ('%s') into the output file ('%s') failed: %s", $renamedfile.".contigs", $outputfasta, $!));

#-----------------------------------------------------------------------------#

=head1 AUTHOR

Frank Foerster NAME S<frank.foerster@biozentrum.uni-wuerzburg.de>

=cut















