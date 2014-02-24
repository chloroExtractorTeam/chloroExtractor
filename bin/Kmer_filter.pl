#!/usr/bin/env perl
use warnings;
no warnings 'qw';
use strict;

use Carp;
use Getopt::Long;
use Pod::Usage;
use File::Basename;

use Log::Log4perl qw(:no_extra_logdie_message);
use Log::Log4perl::Level;
Log::Log4perl->init(\<<'CFG');
	log4perl.rootLogger                 = DEBUG, Screen
	log4perl.appender.Screen			= Log::Log4perl::Appender::Screen
	log4perl.appender.Screen.stderr		= 0
	log4perl.appender.Screen.layout		= PatternLayout
	log4perl.appender.Screen.layout.ConversionPattern = [%d{MM-dd HH:mm:ss}] [rkf] %m%n
CFG

use FindBin qw($RealBin);
use lib "$RealBin/../lib/";

use Data::Dumper;
$Data::Dumper::Sortkeys = 1;

use Fastq::Parser;
use Fastq::Seq;
use Jellyfish;
use Kmer;
use Verbose::ProgressBar;


our $VERSION = '0.02';


##------------------------------------------------------------------------##

=head1 NAME 

read_kmer_filter.pl

=cut

=head1 DESCRIPTION

Filter reads based on kmer count information.

=cut


=head1 CHANGELOG

=cut

=head2 0.02

=over

=item [Change] Only retrieve _2 counts if _1 did not fail.

=item [Feature] --kmer-shift - only take every n-th kmer into account

=back

=head2 0.01

=over

=item [Initial]

=back

=cut

=head1 TODO

--help

-q 0 <= q <= 100

#=item --both [OFF]
#
#Only in paired mode. By default, reads pair is discarded if either one read 
# fails the filter. With --both pairs, where both reads fail are discarded.
#
#=item --dump-discarded
#
#Dump discarded reads in extra files.
#
#=item -o|--out <FILE>
#
#Output file name. There will be thread wise temporary output files created
# at --out location and finally merged to --out.
#
#=item [-t|--threads <INT>] [1]
#
#Number of threads to use

=cut

=head1 SYNOPSIS

  read_kmer_filter.pl [<OPTIONS>] -h <HASH> -k <KMERSIZE> -c <NUMBER_OF_GOOD_KMERS> -l <LOWER_BOUND> -u <UPPER_BOUND> -1 <FASTQ> [ <FASTQ> ... ] [-2 <FASTQ> <FASTQ> ... ]
  
=cut

=head1 OPTIONS

=over

=item -k|--kmer-hash <FILE>

Jellyfish kmer count hash to use for count assignment.

=item -1|--reads <FASTQ> [<FASTQ> ...]

Read files in FASTQ format, single or first of pair.

=item -2|-mates <FASTQ> [<FASTQ> ...]

Read files in FASTQ format, second of pair. Use the same order as for --reads.

=item -o|--out

Prefix for the output files

=item -s|--kmer-size

Size of kmers used in Jellyfish hash.

=item -c|--cutoff

Minimum number of kmers in range [--lower .. --upper] for read to be accepted.

=item -u|--upper <INT> [1000]

Upper kmer count cutoff. 0 deactivates cutoff.

=item -l|--lower <INT> [5]

Lower kmer count cutoff. 0 deactivates cutoff.

=item [--[no]-penalize-N] [ON]

"N" containing kmers are set to count "0" by default. Disable --penalize-N to use their acutal counts.

=item [--max-reads <INT>] [0]

Stop dumping after --max-reads have been dumped.

=item [--histogram <FILE>]

Write a histogram of trusted kmers count frequencies.

# TODO:
# currently not working - test assumes only "good" kmer in hash - jellyfish hash however
# can contain bad kmers.
#=item [--[no]-perl-hash] [ON]
#
#By default, kmers are read to a Perl hash structure for fastest access. This hash can become massive on large datasets. To avoid memory clashs you can disable this behaviour. The kmers are then queried directly against the Jellyfish hash database, which is slower,
#but does not require loading all kmers to RAM.

=item [--debug]

Turn on debug messages.

=item [--quiet]

=item [--help]

Show this help screen.

Supress verbose information messages.

=back

=cut

my $L = Log::Log4perl::get_logger();
$L->level($INFO);

my @fq_suffixes=(qw/.fq .fastq .FQ .FASTQ/);
my @opt_reads;
my @opt_mates;
my %opt = (
	reads => \@opt_reads,
	mates => \@opt_mates,
	cutoff => undef,
	lower => 5,
	upper => 1000,
	'penalize-N' => 1,
	'kmer-shift' => 1,
    'perl-hash' => 1,
);

Getopt::Long::Configure("no_ignore_case");
GetOptions(\%opt, qw(
	reads|1=s@{,}
	mates|2=s@{,}
	out|o=s
	kmer-hash|h=s
	kmer-size|k=i
	cutoff|c=i
	lower|l=i
	upper|u=i
	kmer-shift|x=i
	penalize-N!
	perl-hash!
    quiet
	debug
	help
)) or $L->logcroak($!);

pod2usage(1) if $opt{help};

$opt{quiet} && $L->level($WARN);
$opt{debug} && $L->level($DEBUG);

$L->debug("GetOptions:\n", Dumper(\%opt));

my $opt_c = $opt{cutoff};
my $opt_u = $opt{upper};
my $opt_l = $opt{lower};
my $opt_m = $opt{'max-reads'};
##------------------------------------------------------------------------##	
# required	
for(qw(kmer-hash kmer-size cutoff)){
	pod2usage("required: --$_") unless defined ($opt{$_}) 
};

@{$opt{reads}} || pod2usage("required: --reads");

if(@opt_mates){
	$L->logcroak("Number of -1 and -2 files differs!") if @opt_mates != @opt_reads;
}

# check db files -e -s
foreach my $file(@opt_reads, @opt_mates){
	$L->logcroak("Cannot find file: $file ") unless -e $file && -f $file;
}

$L->warn("--lower($opt_l) > --upper($opt_u), are you sure that's what you want?") if $opt_u && $opt_l > $opt_u;

$L->logdie("Multiple input files currently unimplemented") if(@opt_reads>1);
my $out_file1 = $opt{out} ? $opt{out}."_1.fq" : basename($opt_reads[0], @fq_suffixes).".fil.fq";
my $out_file2 = $opt{out} ? $opt{out}."_2.fq" : basename($opt_mates[0], @fq_suffixes).".fil.fq";

##------------------------------------------------------------------------##	

my $jf = Jellyfish->new();
my $km = Kmer->new(
	kmer_size => $opt{'kmer-size'},
	shift_by => $opt{'kmer-shift'}
);

my %H = (); # histogram

##------------------------------------------------------------------------##	

my %K = ();
$L->logdie("--no-perl-hash currently not implemented!") unless $opt{'perl-hash'};
my $khash = $opt{'perl-hash'};

if($khash){
	$L->info("Loading kmer hash");
	my $kfh = $jf->dump([
			qw(-c -t), 						# column, tab separated
			$opt_l ? ("-L", $opt_l) : (), 	# lower cutoff
			$opt_u ? ("-U", $opt_u) : (), 	# upper cutoff
			$opt{'kmer-hash'}				# hash
		]);
	if($opt{'penalize-N'}){
		while(<$kfh>){
			my ($k, $v) = split("\t", $_);
			$K{$k} = $v	unless $k =~ tr/N//;
		}
	}else{
		while(<$kfh>){
			my ($k, $v) = split("\t", $_);
			$K{$k} = $v;
		}
	}
	$L->info(scalar keys %K," distinct kmers loaded");
}


##------------------------------------------------------------------------##	

# --no-perl-hash
#		my @counts = $khash
#				? sort{$a<=>$b} $jf->query(
#						['--both-strands', $opt{'kmer-hash'}], 
#						kmers => [$km->cmerize($fq1->seq)],
#						table=>0 
#					)
#			
#			if(
#				$opt_u && $counts[int($#counts * $opt_quf)] > $opt_u
#				||
#				$counts[int($#counts * $opt_qlf)] < $opt_l
#			){ 
#				# TODO: dumping of discarded ...
#			}else{
#				print FQ1 "$fq1";
#			}



my $FC;
my $rc = 0;
my $rct = 0;

unless(@opt_mates){
	$L->info("Filtering: single end mode");
	
	for($FC=0; $FC < @opt_reads;$FC++){
		
		$L->info("File: $opt_reads[$FC]");
		my $fp1 = Fastq::Parser->new(file => $opt_reads[$FC]);
		open (FQ1, '>', $out_file1) or $L->logcroak("$!");

		my $pgc = 0;
		my $pg = Verbose::ProgressBar->new(
			size => $fp1->fh,
			level => 2,
			report_level => $opt{quiet} ? 0 : 2,
		);

		while(my $fq1 = $fp1->next_seq){
			$rct++;
			$pg->update unless $pgc++%10000;
			
			my $c=0;
			$c+= exists $K{$_} for $km->cmerize($fq1->seq),
			
			$H{$c}++;

			if($c<$opt_c){
				$L->debug("discarded ",$fq1->id());
				# TODO: dumping of discarded ...
			}else{
				print FQ1 "$fq1";
				$rc++;
				if($opt_m && $rc > $opt_m){
					$L->info("--max-reads ",$opt{'max-reads'}," reached");
					last;
				}
			}
			
		}
		$pg->finish;
		close FQ1;
	}	
}else{
	$L->info("Filtering: paired end mode");
	
	for($FC=0; $FC < @opt_reads;$FC++){
		$L->info("Files: $opt_reads[$FC] $opt_mates[$FC]");
		my $fp1 = Fastq::Parser->new(file => $opt_reads[$FC]);
		open (FQ1, '>', $out_file1) or $L->logcroak("$!");
		my $fp2 = Fastq::Parser->new(file => $opt_mates[$FC]);
		open (FQ2, '>', $out_file2) or $L->logcroak("$!");

		my $pgc=0;
		my $pg = Verbose::ProgressBar->new(
			size => $fp1->fh,
			level => 2,
			report_level => $opt{quiet} ? 0 : 2
		);
		
		while(
			(my $fq1 = $fp1->next_seq) &&
			(my $fq2 = $fp2->next_seq)
		){
			$rct++;
			$pg->update unless $pgc++%10000;
			
			my $c1=0;
			my $c2=0;

			$c1+= exists $K{$_} for $km->cmerize($fq1->seq);
			$c2+= exists $K{$_} for $km->cmerize($fq2->seq);
			
			$H{$c1}++;
			$H{$c2}++;

			if($c1 < $opt_c && $c2 < $opt_c){
				$L->debug("discarded ".$fq1->id." / ".$fq2->id);
			}else{
				print FQ1 "$fq1";
				print FQ2 "$fq2";
		
				$rc++;
				if($opt_m && $rc > $opt_m){
					$L->info("--max-reads ",$opt{'max-reads'}," reached");
					last;
				}
			}
		}
		$pg->finish;
		close FQ1;
		close FQ2;
	}
}

$L->info("Kept $rc of $rct (",sprintf("%0.1f%%", $rc/$rct*100),") reads/pairs");

if($opt{histogram}){
	open(HIST, ">$opt{histogram}") or $L->logdie("Can't open file $opt{histogram}!");
	print HIST "$_\t$H{$_}\n" foreach(sort {$a <=> $b} keys %H);
	close HIST or $L->logdie("$!");
}














