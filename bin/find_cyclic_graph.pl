#!/usr/bin/env perl

=head1 NAME

find_cyclic_graph.pl

=head1 SYNOPSIS

    $ perl find_cyclic_graph.pl -i graph.fastg

=head1 OPTIONS

=over

=item -i|--in <FASTG FILE>
    
    Fastg Graph File

=item -o|--out <NAME> [basename "--in".fa]

    Output Filename

=item -V|--version

Display version.

=item -h|--help

Display this help.

=back

=cut
#----------------------#
# Modules

use strict;
use warnings;

use Getopt::Long qw(:config no_ignore_case);
use Pod::Usage;
use Log::Log4perl qw(:no_extra_logdie_message);
use Log::Log4perl::Level;


use Data::Dumper;

use FindBin qw($RealBin $Script);
use lib "$RealBin/../lib/";

use File::Basename;
use File::Copy;




#--------------------#
# Globals

our $VERSION = 0.01;

our $ID = 'fcg';

# Log4perl logger

my $L = Log::Log4perl::get_logger();

my $log_cfg = 'log4perl.rootLogger                     = INFO, Screen
log4perl.appender.Screen                = Log::Log4perl::Appender::Screen
log4perl.appender.Screen.stderr         = 1
log4perl.appender.Screen.layout         = PatternLayout
log4perl.appender.Screen.layout.ConversionPattern = [%d{yy-MM-dd HH:mm:ss}] ['.$ID.'] %m%n
';
Log::Log4perl->init( \$log_cfg );
#------------------#
# Get Options

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
               	version|V!
	       	help|h!
                debug|D!
			)

) or $L->logcroak('Failed to "GetOptions"');


# help
$opt{help} && pod2usage(1);

# version
if($opt{version}){
	print "$VERSION\n"; 
	exit 0;
}

#--------------------------#
# Config

my %cfg;

#core
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

#---------------------------#
# required






#------------------------#
# MAIN

$opt{out} ||= basename($opt{in}, qw(.fastg)).".fa";

# Check infile Format
		       
open(FH, "<", $opt{in}) || die "Unable to load file: $opt{in}\n";
my $first_line = <FH>;

unless ($first_line =~ /^>\S+;/){                               
    $L->logdie("InputFile is not a fastg Graph");
}
else {
    $L->info("Start to process $opt{in}");
}



## DO YOUR THING GRAPH!


close FH;

#--------------------------#
# Methodes


=head1 AUTHOR

Frank Förster S<frank.foerster@uni-wuerzburg.de>
Simon Pfaff S<simon.pfaff@stud-mail.uni-wuerzburg.de>
Aaron Sigmund S<aaron.sigmund@stud-mail.uni-wuerzburg.de>

=cut
