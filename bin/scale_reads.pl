#!/usr/bin/env perl

=head1 NAME

scale_reads.pl

=head1 DESCRIPTION

Determine coverage of a target sequence population within a NGS read
data set. Target sequence typically are host genome, organelle
genomes, bacterial contaminants, plasmids, transposable elements ...

Requirement for the estimation is reference library of representative
sequences. This can be fragments either form the target organism, a
close relative or a collection of homologs. The library should contain
sequences of regular abundance with in the target sequence, abundance
in the reference library does not play a role. Also the sequences
should be unique to the target sequence and at best not occur in any
other sequences within the sample.

=head1 SYNOPSIS

  $ perl ngs_coverage.pl -1 lib_1.fq -2 lib_2.fq -r reference.fa -t 200

=head1 OPTIONS

=over

=item -1|--reads

Input reads file, first of pair.

=item -2|--mates

Input reads file, second of pair

=item -o|--out

Output prefix.

=item -r|--reference

Reference sequences for mapping. Fasta or already created bowtie2
index prefix.

=item -t|--target-coverage [0]

Return subsetted read files with this estimated coverage. 0 disables
the feature.

=item -m|--max-reads [10000]

Use this many mapped reads for coverage computation.

=item --kmer-size [31]

Utilized kmer size.

=item --[no]-ref-cov-hist [TRUE]

Generate a coverage histgram from reads mapped onto
reference.

=item --threads [1]

Number of threads in parallelizable steps.

=item -c|--config

Use user customized config file. Superseeds default config.

=item -V|--version

Display version.

=item --debug

Report debug messages.

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

use Sam::Parser;
use Sam::Seq;
use Sam::Alignment;

use Jellyfish;

#enable in script mapping
use Bowtie2;


#-----------------------------------------------------------------------------#
# Globals

our $VERSION = 0.01;

our $ID = 'scr';

# get a logger
my $L = Log::Log4perl::get_logger();

Log::Log4perl->init( \(q(
	log4perl.rootLogger                     = INFO, Screen
	log4perl.appender.Screen                = Log::Log4perl::Appender::Screen
	log4perl.appender.Screen.stderr         = 1
	log4perl.appender.Screen.layout         = PatternLayout
	log4perl.appender.Screen.layout.ConversionPattern = [%d{yy-MM-dd HH:mm:ss}] [).$ID.q(] %m%n
)));


#-----------------------------------------------------------------------------#
# GetOptions

my %def = (
           threads => 1,
           target_coverage => 0,
           max_reads => 10000,
           ref_cov_hist => 0,
           kmer_size => 31,
);


my %opt = (config => []);

GetOptions( # use %opt (Cfg) as defaults
    \%opt, qw(
                 reads|1=s
                 mates|2=s
                 out|o=s
                 reference|ref-cluster|r=s
                 target_coverage|target-coverage|coverage=i
                 max_reads|max-reads|m=i
                 ref_cov_hist|ref-cov-hist!
                 kmer_size|kmer-size=i
                 threads=i
                 config|c=s{,}
                 version|V!
                 debug|D!
                 help|h!
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
%opt = (%def, %cfg, %opt);


##----------------------------------------------------------------------------##
# required	
for(qw(reads mates reference target_coverage max_reads)){
    if(ref $opt{$_} eq 'ARRAY'){
	pod2usage("required: --$_") unless @{$opt{$_}}
    }else{
	pod2usage("required: --$_") unless defined ($opt{$_})
    }
};

# required binaries
Cfg->Check_binaries(qw(bowtie2 jellyfish samtools bedtools));     


# debug level
$L->level($DEBUG) if $opt{debug};
$L->debug('Verbose level set to DEBUG');

$L->debug(Dumper(\%opt));



#-----------------------------------------------------------------------------#
# MAIN

# TODO: working dir
my $opt_o1 = $opt{out}."_1.fq";
my $opt_o2 = $opt{out}."_2.fq";



my $bowtie2 = Bowtie2->new(
    path => $opt{bowtie2_bin}
    );

# TODO: bowtie2 generate db -> prevent issues with different indices
# on different architectures or bowtie2 versions

unless(-e $opt{reference}.'.1.bt2'){
    $L->info("Building bowtie2 index");
    $bowtie2->build($opt{reference});
}else{
    $L->info("Using existing bowtie2 index $opt{reference}.*.bt2");
}

$L->info("Running bowtie2");

$bowtie2->run(
    @{$opt{bowtie2_params}},
    "-1" => $opt{reads},
    "-2" => $opt{mates},
    "-x" => $opt{reference},
    "-p" => $opt{threads} || 1
    );

$L->debug(Dumper($opt{bowtie2_params}));

# read output on the fly
my $sp = Sam::Parser->new(
    fh => $bowtie2->stdout
    );

my $bam = $opt{out}."-ref.bam";
open(BAM, "| samtools view -Sbu /dev/fd/0 > $bam") or $L->logdie($!);
open(BED, ">", $opt{out}."-ref.bed") or $L->logdie($!);
open(FQ, ">", $opt{out}."-ref.fq") or $L->logdie($!);

my %h; 
my $ss;

my %seqs;

Sam::Seq->Trim(0); # disable trimming of read

# read output on the fly
while(%h = $sp->next_header_line('@SQ')){
    $seqs{$h{'SN'}} = Sam::Seq->new(
	id => $h{'SN'},
	len => $h{'LN'}
	);

    print BAM '@SQ'."\tSN:$h{'SN'}\tLN:$h{'LN'}\n";
    print BED "$h{'SN'}\t$h{'LN'}\n";
}

close BED;

my $current_cov;
my $genome_size;
my $last_id;
my $closest_ref;
my %refs;

my $c;
my $tlen_sum; 
my $seqwithtlen;


while(my $aln = $sp->next_aln()){
    print BAM "$aln";
    print FQ "@",$aln->qname,"\n", $aln->seq, "\n+\n",$aln->qual,"\n";
    
    # closest ref
    my $id = $aln->rname();
    my ($ref_id) = $id =~ /(NC_\d+)/;
    $refs{$ref_id}++; # increment mapping on ref count

    # insert size
    if (abs($aln->tlen) > 0){
	$tlen_sum+= abs($aln->tlen); # compute isize
	$seqwithtlen++;
    }

    $c++;
    last if $c >= $opt{max_reads}
}       

  

$bowtie2->cancel();
close BAM;
close FQ;


if ($opt{ref_cov_hist}){
    $L->info("Recalculating accurate per base coverages");
    bam_coverage($opt{out});
}


if($c < $opt{max_reads}){
    $L->warn("Could not detect sufficient plastid data in input reads");	
    $L->info("You might need to increase the amount of input data");
    $L->info("Also make sure, your library contains plastid reads at all");

    exit 1;
}

($current_cov, $genome_size) = estimate_kmer_coverage($c);

# what if not enough coverage in entire file.
if(! $current_cov || $current_cov < $opt{target_coverage}){
    $L->warn("Could not detect sufficient plastid data in input reads");	
    $L->info("You might need to increase the amount of input data");
    $L->info("Also make sure, your library contains plastid reads at all");

    exit 1;
}

# TODO: check size


$closest_ref = (sort{$refs{$b} <=> $refs{$a}}keys %refs)[0];

print "coverage\t", $current_cov, "\n";
print "genome_size\t", $genome_size, "\n";
print "target_coverage\t", $opt{target_coverage}, "\n";
print "insert_size\t", int($tlen_sum/$seqwithtlen), "\n";
print "closest_ref\t", $closest_ref || 'NA', "\n";


if($opt{target_coverage}){
    $L->info("Creating libraries");

    my $current_reads = qx/wc -l <$opt{reads}/ / 4;
    my $target_reads = int($current_reads / $current_cov * $opt{target_coverage});
    
    $L->info("Got $current_reads reads, extracting $target_reads reads");
    
    my $target_lines = $target_reads * 4;
    
    my $sed_cmd1 = "sed '".$target_lines."q' ".$opt{reads}." >".$opt_o1;
    
    $L->debug($sed_cmd1);
    qx($sed_cmd1);
    
    my $sed_cmd2 = "sed '".$target_lines."q' ".$opt{mates}." >".$opt_o2;
    
    $L->debug($sed_cmd2);
    qx($sed_cmd2);
}
    


#-- SUBS ----------------------------------------------------------------------#




sub median {
  my @array = @{$_[0]};
  my $median = (sort{$a<=>$b}@array)[@array/2];
  return $median;
}


sub estimate_kmer_coverage{
    my $c = shift;
    my $core_reads = $opt{out}."-ref.fq";
    my $jff = $opt{out}."-ref.jf";

    # jellyfish pt counts
    $L->info("Running jellyfish");
    my $jf_count = join(" ",
                        "jellyfish",
			qw(count -s 10M -C -L 20),
                        "-m" => $opt{kmer_size},
			"-t" => $opt{threads},
			"--if" => $core_reads, 
			"--output" => $jff, 
                        $opt{reads},
			$opt{mates}
                       );

    $L->debug($jf_count);
    my $jf_count_re = qx/$jf_count/;
    
    $L->logdie($jf_count_re) if $jf_count_re;

    # dump and R stat counts
    # my $R = q/counts <- read.table(pipe('jellyfish dump -c --tab /.$jff.q/ | cut -f2 '), header=F);/   #   .q/summary(counts[,1])/;

    # my $R = q{data <- read.table(pipe('jellyfish histo -h 1000000 }
    #   .$opt{out}."-ref.jf"
    #     .q{'), header=F);}
    #       .q{med = data[,1][cumsum(data[,2]) > sum(data[,2])/2][1];}
    #         .q{print(med);};

    # $L->debug("Running R: '$R'");
    # my @Rr = qx(echo "$R" | R --vanilla --slave);
    # $L->debug("R:\n", @Rr);

    # $L->logdie(@Rr) unless $Rr[0] =~ /^\[1\]/;
    
    # my ($r, $med) = split(/\s+/, $Rr[0]);

    # $L->info("Estimated coverage\n", @Rr);

    # return $med;

    my $R = "R --vanilla --slave --args scr jf0.jf scr-ref.jf <$RealBin/make_plots.R |";
    $L->debug("Running R: '$R'");
    open(R, $R) or $L->logdie($!);

    my %Rr;
    while(<R>){
	my ($k, @v) = split("\t", $_);
	$Rr{$k} = @v > 1 ? \@v : $v[0];
    }
    $L->debug("R:\n", Dumper(\%Rr));

    return $Rr{coverage}, $Rr{size};
    
}



=head2 bam_coverage

Compute per contig coverage using bedtools and write to file.

=cut

sub bam_coverage{
    my ($pre) = @_;
    my $bam = $pre."-ref.bam";
    my $bed = $pre."-ref.bed";
    my $pre_sorted = $pre."-ref-sorted";
    my $bam_sorted = $pre."-ref-sorted.bam";
    my $tsv = $pre."-ref-bp-cov.tsv";
    my $head = '"id\tposition\tcoverage"';
    my $re;
    $re = qx/samtools sort $bam $pre_sorted/ && $L->warn($re);
    $re = qx/echo $head > $tsv/ && $L->warn($re);
    $re = qx/bedtools genomecov -d -ibam $bam_sorted -g $bed >> $tsv/ && $L->warn($re);
}

#-----------------------------------------------------------------------------#

=head1 AUTHOR

Clemens Weiss S<clemens.weiss@stud-mail.uni-wuerzburg.de>

Thomas Hackl S<thomas.hackl@uni-wuerzburg.de>

=cut



