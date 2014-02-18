#!/usr/bin/env perl

=head1 NAME


=head1 DESCRIPTION


=head1 SYNOPSIS


=head1 OPTIONS

=over

=item -s|--sam

SAM file. REQUIRED.

=item --create-config

Create a config file with default settings for user customization.

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

use Sam::Parser;
use Sam::Seq;
use Sam::Alignment;


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
my %opt = Cfg->Read_Cfg($core_cfg); 

# user defaults and overwrite core
my $user_cfg;
for(my $i=0; $i<@ARGV; $i++){
        if($ARGV[$i] =~ /-c$|--config$/){
                $user_cfg = $ARGV[$i+1];
                last;
        }
}

%opt = (%opt, Cfg->Read_Cfg($user_cfg)) if $user_cfg; # simple overwrite


#TODO: custom config


#-----------------------------------------------------------------------------#
# GetOptions

GetOptions( # use %opt (Cfg) as defaults
	\%opt, qw(
		version|V!
		debug|D!
		help|h!
		config|c=s
		create_config|create-config
		sam|s=s
	)
) or $L->logcroak('Failed to "GetOptions"');

# help
$opt{help} && pod2usage(1);

# version
if($opt{version}){
	print "$VERSION\n"; 
	exit 0;
}

# create template for user cfg
if($opt{create_config}){
	pod2usage(-msg => 'To many arguments', -exitval=>1) if @ARGV > 1;
	my $user_cfg = @ARGV ? $ARGV[0] : basename($core_cfg);
	copy($core_cfg, $user_cfg) or $L->logdie("Creatring config failed: $!");
	$L->info("Created config file: $user_cfg");
	exit 0;
}

#required stuff  
for(qw(sam)){
        pod2usage("required: --$_") unless defined ($opt{$_})
};

# debug level
$L->level($DEBUG) if $opt{debug};
$L->debug('Verbose level set to DEBUG');

$L->debug(Dumper(\%opt));



#-----------------------------------------------------------------------------#
# MAIN

my $sp = Sam::Parser->new(file => $opt{sam});
my %h; 
my $ss;

my %seqs;

while(%h = $sp->next_header_line('@SQ')){
        $seqs{$h{'SN'}} = Sam::Seq->new(
                id => $h{'SN'},
                len => $h{'LN'},
        );  
}

my $c; 
while(my $aln = $sp->next_aln){
        #last if ++$c > 100;
        my $id = $aln->rname();
        $seqs{$id}->add_aln($aln);
}

#print Dumper($ss);
for my $ss(values %seqs){
        my @covs=$ss->coverage;
        print $ss->id(),"\t@covs\n";
}





































#-----------------------------------------------------------------------------#

=head1 AUTHOR

Clemens Weiss S<clemens.weiss@uni-wuerzburg.de>

Thomas Hackl S<thomas.hackl@uni-wuerzburg.de>

=cut















