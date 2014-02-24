#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use FindBin;
use lib "$FindBin::Bin/../lib";
use Data::Dumper;
use Storable;

use Blast;
use Blast::Parser::TSV;

my %options;

=head1 NAME 

highlight_and_mask_IR.pl

=head1 DESCRIPTION

This script is part of the wgaPipeline and if triggered by the pipeline highlights IR regions 
and/or masks them in a temporary fasta file before the WGA

=head1 USAGE

  $ perl highlight_and_mask_IR.pl --in=<file> [--seqprefix <string> --outprefix=<string> --[no]highlight --[no]mask] [options]

=head1 OPTIONS

=over 25

=item --in=<file>

path to the input fasta file. 

=cut

$options{'in=s'} = \(my $opt_in);


=item [--color=<string>] 

Color of the highlighting (default: black)

=cut

$options{'color=s'} = \(my $opt_color='black');


=item [--outprefix=<string>] 

prefix for the output files

=cut

$options{'outprefix=s'} = \(my $opt_outprefix);

=item [--seqprefix=<string>] 

prefix for the sequences in the highlight file (Q_ for query or R_ for reference)

=cut

$options{'seqprefix=s'} = \(my $opt_seqprefix);

=item [--[no]highlight] 

Trigger highlighting of the IRs

=cut

$options{'highlight!'} = \(my $opt_highlight);


=item [--[no]mask] 

Trigger masking of the IRs

=cut

$options{'mask!'} = \(my $opt_mask);

=item --blast-bin [autodetect]

=cut

$options{'blast-bin=s'} = \(my $blast_bin = `which blastn 2>/dev/null`);
chomp($blast_bin);

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

pod2usage(-msg => "Missing argument --in", -verbose => 0) unless ($opt_in);
$blast_bin || pod2usage(-msg => "Can't locate blast binary, please specify\n", -exitval => 1);

my %SEQS_IN;
my %sequence_only;
my @irRegions;

my ($infile_name,$infile_dir,$infile_suffix) = fileparse($opt_in, ".fa", ".fas", ".fasta", ".txt");

my $db_fa;
my $db_idx;

$db_fa = $opt_in;
$db_idx = $infile_name.'.idx';
my $masked_file = $infile_name.'.masked'.$infile_suffix;

############################################################################
# DB index
############################################################################

my $db;
open (my $db_fh, $db_fa) or die "Can't open FASTA database $db_fa\n";
if (-f $db_idx){
	$db = retrieve($db_idx); 
}else{
	$db = index_db($db_fh, $db_idx);
}


############################################################################
# Input sequences
############################################################################

my $seq_id;

open (IN, $opt_in) or die "can't open initial transcript file $opt_in\n";

while(<IN>){
	if (/^>/){
		$seq_id = fa2id($_);
		$SEQS_IN{$seq_id} = $_;
	}else{
		$SEQS_IN{$seq_id}.= $_;	
		chomp;
		$sequence_only{$seq_id} = "" unless(exists $sequence_only{$seq_id});	
		$sequence_only{$seq_id}.= $_;	
	}
}



############################################################################
# Blast and highlightResults
############################################################################

my ($ir_id, $ir_start, $ir_end);
@irRegions = blast();
open(OUT, '>', $opt_outprefix.".high") or die "Can't open file $!";
foreach my $irRegion(@irRegions){
	if($irRegion->{start1}==1){
#		print OUT $irRegion->{start2}."\n";
		$ir_start = $irRegion->{start2};
		$ir_id = $irRegion->{id};
	}
	elsif($irRegion->{end1}==length($sequence_only{$irRegion->{id}})){
#		print OUT $irRegion->{end2}."\n";
		$ir_end = $irRegion->{end2};
		$ir_id = $irRegion->{id};
	}
}
my $ir_seq = substr($sequence_only{$ir_id}, $ir_start-1, $ir_end-$ir_start+1);
#print OUT ">ir\n$ir_seq\n";
my $ir_seq_rc = $ir_seq;
$ir_seq_rc =~ tr/ACGTacgt/TGCAtgca/;
$ir_seq_rc = scalar(reverse($ir_seq_rc));
print OUT ">ir_rc\n$ir_seq_rc\n";
close OUT or die "$!";


############################################################################
# Subs
############################################################################


sub blast{
	qx(makeblastdb -in $opt_in -dbtype nucl);
	my @irRegions;
	
	# blast
	my $blast = Blast->blast({
		_bin => 			$blast_bin,
		-query => 			$opt_in, 
		-task => 			'blastn',
		-db => 				$opt_in,
		-evalue => 			'1e-30',
		-outfmt => 			'"6 sseqid sstart send qseqid qstart qend"',
		-num_threads => 	1
	});	
		
	my $bp = Blast::Parser::TSV->new(output => $blast->output);		
	
	while (my ($sseqid, $sstart, $send, $qseqid, $qstart, $qend) = $bp->get_next_row()){
		# only IRs on the same contig/chromosom are detected 
		next unless($sseqid eq $qseqid);
		# only rc matches are considered
		next unless(($sstart < $send) != ($qstart < $qend));
		push(@irRegions, {'id'=>$sseqid,'start1'=>($sstart>$send ? $send : $sstart),'end1'=>($sstart>$send ? $sstart : $send),
			'start2'=>($qstart>$qend ? $qend : $qstart),'end2'=>($qstart>$qend ? $qstart : $qend)});
	}
		
	# return new contigs
	return @irRegions;		
};


sub index_db{
	my ($db_fh, $db_idx) = @_;
	my $db;
	while(<$db_fh>){
		if (/^>(\S+)/){
			$db->{$1} = tell()-length($_);
		}
	}
	store($db, $db_idx); 
	seek($db_fh, 0, 0);
	return $db;
}

sub write_fasta{
	my $file = shift;
	my $contigs = shift;
	open(FILE, '>', $file) or die "Can't open $file\n";
	print FILE values %$contigs;
	close FILE;
}


sub fa2id{
	$_[0] =~ m/^>\s*(\S+)/;
	return $1;
}

sub fa2header{
	return ( split /\n/, shift )[0];
}

sub fa2seq{
	(my $seq = shift) =~ s/.*\n//;
	return $seq;
}

sub rc{
	scalar reverse $_[0] =~ tr/ATGC/TACG/r;
}

=head1 LIMITATIONS

If you encounter a bug please drop me a line.

=head1 AUTHORS

=over

=item * Markus Ankenbrand, markus.ankenbrand@stud-mail.uni-wuerzburg.de

=back


