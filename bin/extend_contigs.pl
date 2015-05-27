#!/usr/bin/env perl

package exc;

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

use FindBin qw($RealBin $Script);
use lib "$RealBin/../lib/";

use File::Basename;
use File::Copy;
use File::Path qw(make_path);
use File::Find;

# additional modules
use Cfg;
use Fasta::Parser;
use Fastq::Parser;
use Bowtie2;
use Sam::Parser;
use Cwd;

#-----------------------------------------------------------------------------#
# Globals

our $VERSION = 0.02;

our $ID = 'exc';

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
		coverage=s
		threads=i
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
    my $user_cfg = Cfg->Copy($core_cfg, $opt{create_config}) or $L->logdie("Creating config failed: $!");
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
for(qw(in out insert_size mates reads)){
    if(ref $opt{$_} eq 'ARRAY'){
	pod2usage("required: --$_") unless @{$opt{$_}}
    }else{
	pod2usage("required: --$_") unless defined ($opt{$_})
    }
};


## check if the border was set, otherwise use 2x insert_size as default value
unless (exists $opt{border} && $opt{border} > 0)
{
    $opt{border} = $opt{insert_size}*2;
}


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

# create a sequence store for each border combination
# therefore I need to retrieve a list of all sequence names
my @seq_names = sort keys %filehandles;
# now prepare each combination
for (my $i=0; $i < @seq_names - 1; $i++)
{
    for (my $j=$i+1; $j < @seq_names; $j++)
    {
	my $mergename = join("_", ("merged", $seq_names[$i], $seq_names[$j]));
	
	store_sequence_and_create_folder(name => $mergename, file_out => $fasta_out, seen_names => \%contig_ends_seen, filehandle => \%filehandles);
    }
}

#delete empty directories
my $dir = './merged_ass';
my @dirs;
find(\&get_dir,$dir);

rmdir($_) for(reverse(@dirs));

sub get_dir{
  push(@dirs,$File::Find::name) if(-d $_);
}

my $bowtie2 = Bowtie2->new(
    path => $opt{bowtie2_path},
);

$L->info("Building bowtie2 index");

$bowtie2->build($fasta_contig_ends);

##
$L->debug("FILEHANDLES:", Dumper(\%filehandles));

$L->info("Running bowtie2");

$bowtie2->run(
    @{$opt{bowtie2_params}},
    "-1" => $opt{reads},
    "-2" => $opt{mates},
    "-x" => $fasta_contig_ends,
    "-p" => $opt{threads} || 1,
#    "--all" => '',           # generate all alignments
);

# read output on the fly
my $sp = Sam::Parser->new(
  fh => $bowtie2->stdout,
);

# print sequences of read pairs with both reads mapped
while( my ($aln1, $aln2) = $sp->next_pair() ){
    # if no read was aligned, skip this pair
    next if ($aln1->is_unmapped() && $aln2->is_unmapped());

    # two options:
    # 1) both reads are mapped
    # 2) only one read was mapped
    my @contig_name = ();
    my $mappingtype = undef;

    if (!$aln1->is_unmapped() && !$aln2->is_unmapped())
    {
	# both reads are mapped
	# onto same contig?
	if ($aln1->rname() eq $aln2->rname())
	{
	    # both reads on same contig
	    $contig_name[0] = $aln1->rname();
	    $mappingtype = 'complete';
	} else {
	    # both reads on different contigs
	    $L->debug(sprintf("Found joining pair for contigs %s and %s", $aln1->rname, $aln2->rname));
	    $mappingtype = 'overlapping';
	    @contig_name=($aln1->rname(), $aln2->rname());

	    # add a merged name to fill the corresponding folder
	    push(@contig_name, join("_", ("merged", sort @contig_name)));

	    # store the information in the joined_with hash
	    $filehandles{$aln1->rname()}{joined_with}{$aln2->rname()}++;
	    $filehandles{$aln2->rname()}{joined_with}{$aln1->rname()}++;
	}
    } else
    {
	# what is the contig name?
	$contig_name[0] = ($aln1->is_unmapped) ? $aln2->rname() : $aln1->rname();
	$mappingtype = 'half';
    }

    # reorder the alignments
    # SAM always shows the mapped read first, but this one can be first or second read
    if ($aln1->is_second())
    {
	($aln1, $aln2) = ($aln2, $aln1);
    }
    
    # generate a fastq sequence block for the first read
    my $first_read = Fastq::Seq->new(
	"@".$aln1->qname(),
	$aln1->seq(),
	"+",
	$aln1->qual(),
	);

    # generate a fastq sequence block for the second read
    my $second_read = Fastq::Seq->new(
	"@".$aln2->qname(),
	$aln2->seq(),
	"+",
	$aln2->qual(),
	);


    foreach my $contig (@contig_name)
    {
	$L->debug($contig);
	$filehandles{$contig}{reads}->append_seq($first_read);
	$filehandles{$contig}{mates}->append_seq($second_read);
	$filehandles{$contig}{$mappingtype}++;
    }
}

$bowtie2->finish;

$L->info("Border mapping statistics");

### printing the statistics
# generating a overall stastistic
my ($half, $complete, $overlapping) = (0, 0, 0);
foreach my $contig_border (keys %filehandles) {
    # check if the contig is capable of joining other contigs
    my $joined_contigs = 0;
    my $self_joined = 0;

    my $contig = $contig_border;
    $contig =~ s/_[53]prime$//;
    if (scalar (keys %{$filehandles{$contig_border}{joined_with}}) > 0)
    {
	foreach my $joined_contig (keys %{$filehandles{$contig_border}{joined_with}})
	{
	    $joined_contig =~ s/_[53]prime$//;
	    if ($joined_contig eq $contig)
	    {
		$L->debug(sprintf("Contig %s mapped other contig border", $contig));
		$self_joined++;			  
	    } else {
		$L->debug(sprintf("Contig %s mapped other contig %s", $contig, $joined_contig));
		$joined_contigs++;
	    }
	}
    }
    $L->debug(
	sprintf "%s\tmapped reads (half/complete/overlapping): (%d/%d/%d) joining %d contigs", 
	$contig_border, 
	$filehandles{$contig_border}{half}, 
	$filehandles{$contig_border}{complete}, 
	$filehandles{$contig_border}{overlapping},
	$joined_contigs
	);
    $half+=$filehandles{$contig_border}{half};
    $complete+=$filehandles{$contig_border}{complete};
    $overlapping+=$filehandles{$contig_border}{overlapping};
}
$L->info(sprintf("Overall number of mapped reads (half/complete/overlapping): %d/%d/%d", $half, $complete, $overlapping));

## here we need to run a velvet assembly for each folder
$L->info("Running assemble_reads.pl");

my $patchfilename = cwd()."/extended_asc.fa";

my @cmd = (
    $RealBin."/assemble_reads.pl",
    "-c", @{$opt{config}},
    "--workingdir ./",
    "--isize", $opt{insert_size},
    "--extendmode", 
    "--out", $patchfilename,
    "--velvet_out", "velvet_out_extend",
    '--velvetg_parameter', "'-cov_cutoff ".(0.25*$opt{coverage})." -exp_cov ".(0.9*$opt{coverage})."'",
);

foreach my $contig_border (keys %filehandles)
{
    push(@cmd, ("--reads", $contig_border."/reads.fq"));
    push(@cmd, ("--mates", $contig_border."/mates.fq"));
}

$L->debug("Running assemle_reads using the command '@cmd'");
qx(@cmd);
my $errorcode = $?;
if ($errorcode != 0)
{
    $L->logdie("Errorcode for command '@cmd' was not 0!");
}

## filtering of the output sequences
my $final_out = Fasta::Parser->new(
    file => $opt{out},
    mode => '>'                    # overwrite an existing file
    );

## check if output file exists
unless (-e $patchfilename)
{
    $L->debug("Skipping file '$patchfilename' because it does not exist");
}

my $fasta_patches = Fasta::Parser->new(
    file => $patchfilename
    );


while (my $patch=$fasta_patches->next_seq())
{
    # filter for a minimum length (1.5xinsert size)
    unless (length($patch->seq()) > 1.5*$opt{insert_size})
    {
	$L->debug("Patch skipped because it is not long enough");
	next;
    }
    
    $final_out->append_seq($patch);
}

sub store_sequence_and_create_folder
{
    my %params = @_;
    my $name = $params{name};
    if (exists $params{seen_names}{$name})
    {
	$L->logdie("The folder '$name' already exists! One possibility are multiple occurence of the same sequence name in the input file");
    }
    $params{seen_names}{$name}++;
    make_path($name) || $L->logdie("Unable to create folder '$name' $!");

    # generate an empty set for statistics and overlapping mappings
    $params{filehandle}{$name} = {
	overlapping => 0,
	half => 0,
	complete => 0,
	joined_with => {}
    };

    $params{filehandle}{$name}{reads}=Fasta::Parser->new(
	file => $name.'/'.'reads.fq',
	mode => '>'                    # overwrite an existing file
	);
    $params{filehandle}{$name}{mates}=Fasta::Parser->new(
	file => $name.'/'.'mates.fq',
	mode => '>'                    # overwrite an existing file
	);

    # generate a sequence object if name and seq are given
    if ($name ne "" && exists $params{seq} && $params{seq} ne "")
    {
	my $seq_obj = Fasta::Seq->new(
	    id => $params{name},
	    seq => $params{seq}
	    );
	$params{file_out}->append_seq($seq_obj);
    }

    return 1;
}

#-----------------------------------------------------------------------------#

=head1 AUTHOR

Frank Foerster NAME S<frank.foerster@biozentrum.uni-wuerzburg.de>

=cut















