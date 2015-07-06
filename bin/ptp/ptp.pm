package ptp;

use strict;
use warnings;
use Data::Dumper;
use Fasta::Parser;
use Cwd;
use FindBin;
use lib "$FindBin::Bin/../../lib/";
use File::Basename;

my $dir = getcwd();
my $basename = basename($dir);

my @path = qw(Contig2_3prime Contig2_5prime Contig3_3prime Contig3_5prime Contig1_3prime Contig1_5prime Contig3_5prime);
my $contig1;
my $contig2;
my $contigend1;
my $contigend2;
my $contigsdir;
my $overlapdir;
my $seq_contig1;
my $seq_contig2;
my $seq_overlap;
my $output;
my $seqminuslength = 0;
my $consensus;
my $consensusend;
my $bestalign;

$contigsdir = $dir."/asc.fa";

$path[0] =~ /(.+)_[53]prime/;
$contig1 = $1;
$consensus = readseq($contigsdir, $contig1);
$consensus =~ /.+\n(.+)/;
$consensus = $1;

for ( my $i = 0; $i < @path+0; $i++ )
{
    $path[$i] =~ /(.+)_([53]prime)/;
    $contig1 = $1;
    $contigend1 = $2;
    if ( $path[$i+1] )
    {
    $path[$i+1] =~ /(.+)_([53]prime)/;
    $contig2 = $1;
    $contigend2 = $2;
    }

    if ($contig1 eq $contig2)
    {
	next();
    }
    
    $seq_contig1 = readseq($contigsdir, $contig1);
    $seq_contig2 = readseq($contigsdir, $contig2);

    $overlapdir = getoverlapdir($dir, $contig1, $contigend1, $contig2, $contigend2);
    
    $seq_overlap = readseq($overlapdir, "NODE_1");

    #if ( ( $contigend1 eq "3prime" && $contigend2 eq "5prime" ) || ( $contigend1 eq "5prime" && $contigend2 eq "3prime" ) )
    #{
    print STDERR "Aligning $contig1 with Overlap\n";
    $consensusend = substr($consensus, length($consensus)-1000, 1000);
    $output = align( $consensusend, $seq_overlap );
    #$consensusend = substr($consensus, 0, 1000);
    
    $bestalign = find_bestalign();
    print STDERR "choosing $bestalign\n";
    $consensus = substr($consensus, 0, length($consensus)-1000).find_consens($bestalign);
    my $bla = <>;
    print STDERR "Aligning $contig2 with Overlap\n";
    $consensusend = substr($consensus, length($consensus)-1000, 1000);
    $output = align( $consensusend, $seq_contig2 );
    $bestalign = find_bestalign();
    print STDERR "choosing $bestalign\n";
    $consensus = substr($consensus, 0, length($consensus)-1000).find_consens($bestalign);
    print STDERR "DONE\n";
    #}
    #else
    #{
    #     print"Aligning $contig1 with Overlap\n";
    #     $consensusend = substr($consensus, 0, 1000);
    #     $output = align( $consensusend, $seq_overlap );
    #     $consensus = find_consens().substr($consensus, 1000, length($consensus));
    #     print"Aligning $contig2 with Overlap\n";
    #     $consensusend = substr($consensus, 0, 1000);
    #     $output = align( $consensusend, $seq_contig2 );
    #     $consensus = find_consens().substr($consensus, 1000, length($consensus));
    #     print ">Chloroplast_Consensus\n$consensus";
    #     print STDERR "DONE\n";
    #}
    my $bla = <>;
}

open(my $OUT, '>', 'pl_cons.fa') or die "Could not open file 'pl_cons.fa' $!";
print $OUT ">Plastid_Consensus_$basename\n$consensus";
close $OUT;




sub find_bestalign
{
    my $score;
    my $score_rev;

    open(my $FH, '<', 'needle_out_rev.fa') or die "Could not open file 'needle_out_rev.fa' $!";
    while ( <$FH> )
    {
	if ( $_ =~ /Score: (.+)/ )
	{
	    $score_rev = $1;
	}
    }
    close $FH;
    open($FH, '<', 'needle_out.fa') or die "Could not open file 'needle_out.fa' $!";
    while ( <$FH> )
    {
	if ( $_ =~ /Score: (.+)/ )
	{
	    $score = $1;
	}
    }
    close $FH;
    
    print STDERR "needle_out.fa Score: $score\nneedle_out_rev.fa Score: $score_rev\n";

    if ( $score >= $score_rev )
    {
	return "needle_out.fa";
    }
    else
    {
	return "needle_out_rev.fa";
    }

}
    
sub getoverlapdir
{
    my $overlapdir;
    if ( -e $_[0]."/merged_ass/contigs.fasta_renamed.".$_[1]."_".$_[2]."_ass/contigs.fasta_renamed.".$_[3]."_".$_[4]."/extended_asc/scaffolds.fasta" )
    {
	$overlapdir = $_[0]."/merged_ass/contigs.fasta_renamed.".$_[1]."_".$_[2]."_ass/contigs.fasta_renamed.".$_[3]."_".$_[4]."/extended_asc/scaffolds.fasta";
    }
    else
    {
	$overlapdir = $_[0]."/merged_ass/contigs.fasta_renamed.".$_[3]."_".$_[4]."_ass/contigs.fasta_renamed.".$_[1]."_".$_[2]."/extended_asc/scaffolds.fasta";
    }
    return $overlapdir;
}

sub readseq
{
    my $seq_contig;
    my $contig_in = Fasta::Parser->new(
	file => $_[0],
	mode => '<'
	);
    
    while ($seq_contig = $contig_in->next_seq())
    {
	if ( $seq_contig =~ /$_[1]/ )
	{
	    $seq_contig =~ s/\//_/g;
	    return $seq_contig;
	    last();
	}
    }
    
    $contig_in.DESTROY();
}

sub align
{
    open(my $FH, '>', 'forneedle1.fa') or die "Could not open file 'forneedle1.fa' $!";
    print $FH $_[0];
    close $FH;
    open($FH, '>', 'forneedle2.fa') or die "Could not open file 'forneedle2.fa' $!";
    print $FH $_[1];
    close $FH;

    my $output = qx(needle -asequence forneedle1.fa -bsequence forneedle2.fa -aformat3 markx3 -outfile needle_out.fa -gapopen 10 -gapextend 0.5);
    $output = qx(needle -asequence forneedle1.fa -sreverse_bsequence forneedle2.fa -aformat3 markx3 -outfile needle_out_rev.fa -gapopen 10 -gapextend 0.5);
    return $output;
}

sub find_consens
{
    my %clustalseqs;
    my $id;
    open(my $FHCO, '<', $_[0]) or die "Could not open file $_[0] $!";
    while ( <$FHCO> )
    {
	if ( $_ =~ /#/ || $_ =~ /^\n/ )
	{
	    next();
	}
	if ( $_ =~ /^>.+/ )
	{
	    chomp();
	    $id = $_;
	    $clustalseqs{$id} = "";
	}
	else
	{
	    chomp();
	    $clustalseqs{$id} = $clustalseqs{$id}.$_;
	}

    }
    close $FHCO;

    my @keys = keys(%clustalseqs);
    
    my @seq1 = split("", $clustalseqs{$keys[0]});
    my @seq2 = split("", $clustalseqs{$keys[1]});
    
    my $consens;

    if ( @seq1+0 >= @seq2+0 )
    {
	$consens = compare_seqs( \@seq1, \@seq2 );
    }
    else
    {
	$consens = compare_seqs( \@seq2, \@seq1 );
    }

    return $consens;
}

sub compare_seqs
{
    my $cons = "";

    for ( my $i = 0; $i < @{$_[0]}+0; $i++ )
    {
	unless ( ${$_[1]}[$i] )
	{
	    $cons = $cons.${$_[0]}[$i];
	    next();
	}
	if ( ${$_[0]}[$i] eq "-" )
	{
	    $cons = $cons.${$_[1]}[$i];
	    next();
	}
	if ( ${$_[1]}[$i] eq "-" )
	{
	    $cons = $cons.${$_[0]}[$i];
	    next();
	}
	$cons = $cons.${$_[0]}[$i];
    }
    return $cons;
}

1;
