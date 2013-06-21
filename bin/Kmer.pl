use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;
use FindBin;
use lib "$FindBin::Bin/../lib";
use Fastq::Parser;
use Fastq::Seq;

use Kmer;

my %options;

=head1 NAME 

Kmer.pl

=head1 DESCRIPTION

This script uses kmer information to draw sequencing reads

=head1 USAGE

  $ perl Kmer.pl --kmers=<file> --reads=<file> [options]

=head1 OPTIONS

=over 25

=item --kmers=<file>

path to the kmer file. 

=cut

$options{'kmers=s'} = \(my $opt_kmers);


=item --reads=<file>

the read file in fastq format

=cut

$options{'reads=s'} = \(my $opt_reads);

=item --mates=<file>

the read file containing the mates of <reads> in fastq format

=cut

$options{'mates=s'} = \(my $opt_mates);

=item --out=<prefix>

prefix for the output reads / read pairs

=cut

$options{'out=s'} = \(my $opt_out);

=item --histo=<file>

desired file to dump the histogram (k-mers per read)

=cut

$options{'histo=s'} = \(my $opt_histo);


=item [--cutoff=<int>]

the cutoff (so much kmers are needed to dump the read)

=cut

$options{'cutoff=i'} = \(my $opt_cutoff=1);

=item [--[no]trustall]

if set, all kmers on this read and its mate are trusted if this read is trusted (default: off)

=cut

$options{'trustall!'} = \(my $opt_trustall=0);


=item [--maxreads=<int>]

the maximum number of reads (pairs) to dump 

=cut

$options{'maxreads=i'} = \(my $opt_maxreads=-1);

=item [--[no]verbose] 

verbose is default.

=cut

$options{'verbose!'} = \(my $opt_verbose = 1);


=item [--help] 

show help

=cut

$options{'help|?'} = \(my $opt_help);


=item [--man] 

show man page

=cut

$options{'man'} = \(my $opt_man);

=back




=head1 CODE

=cut


GetOptions(%options) or pod2usage(1);

pod2usage(1) if($opt_help);
pod2usage(-verbose => 99, -sections => "NAME|DESCRIPTION|USAGE|OPTIONS|LIMITATIONS|AUTHORS") if($opt_man);

pod2usage(-msg => "Missing kmers or reads", -verbose => 0) unless ($opt_kmers && $opt_reads);

my %target_kmer;
open(IN, "<$opt_kmers") or die "Can't open file $opt_kmers$!";
while(<IN>){
	chomp;
	my($kmer, $count) = split(/\t/, $_);
	$target_kmer{$kmer}=1;
}
close IN or die "$!";

my $kmer_size = length((keys %target_kmer)[0]);
print "Detected kmer size: $kmer_size\n";

Kmer->KmerSize($kmer_size);
my $fp =Fastq::Parser->new(file => $opt_reads);
open(OUT1, ">$opt_out"."_1.fq") or die "Can't open output file $opt_out"."_1.fq$!";
my $fp2;
if($opt_mates){
	$fp2 =Fastq::Parser->new(file => $opt_mates);
	open(OUT2, ">$opt_out"."_2.fq") or die "Can't open output file $opt_out"."_2.fq$!";
}

my %histogram;
foreach(0..250){
	$histogram{$_}=0;
}

my $dumped_reads = 0;
while(my $fq =$fp->next_seq()){
	my $fq2;
	$fq2 =$fp2->next_seq() if($opt_mates);
	my $counter = 0;
	foreach my $kmer (Kmer->Kmerize_nr($fq->seq)){
		$counter++ if(exists $target_kmer{$kmer});
	}
	$histogram{$counter}++;
	if($counter>=$opt_cutoff){
		print OUT1 $fq->string();
		print OUT2 $fq2->string() if($opt_mates);
		$dumped_reads++;
		last if($opt_maxreads>0 and $dumped_reads>$opt_maxreads);
		if($opt_trustall){
			$target_kmer{$kmer}=1 foreach my $kmer (Kmer->Kmerize_nr($fq->seq));
			$target_kmer{$kmer2}=1 foreach my $kmer2 (Kmer->Kmerize_nr($fq2->seq));
		}
	}
}

close OUT1 or die "$!";
close OUT2 or die "$!" if($opt_mates);

if($opt_histo){
	open(HIST, ">$opt_histo") or die "Can't open file $opt_histo$!";
	print HIST "$_\t$histogram{$_}\n" foreach(sort {$a <=> $b} keys %histogram);
	close HIST or die "$!";
}

sub reverse_complement(){
	my $seq = shift(@_);
	$seq = reverse($seq);
	$seq =~ tr/ACGTacgt/TGCAtgca/;
	return $seq;
}