package ptp;

use strict;
use warnings;
use Data::Dumper;
use Fasta::Parser;
use Cwd;
use FindBin;
use lib "$FindBin::Bin/../../lib/";
use lib "$FindBin::Bin";
use File::Basename;
use fir_wo_out;

my $dir = getcwd();
my $basename = basename($dir);

my @path = fir_wo_out::get_path();
@path = @{$path[0]};

print "using path:\n";
print Dumper(@path);

my $contig1;
my $contig2;
my $contigend1;
my $contigend2;
my $contigsdir;
my $overlapdir;
my $seq_contig1;
my $seq_contig2;
my $seq_overlap;
#my $output;
#my $seqminuslength = 0;
my $consensus;

$contigsdir = $dir."/asc.fa";

#init first path element as consensus
$path[0] =~ /(.+)_[53]prime/;
$contig1 = $1;
$consensus = readseq($contigsdir, $contig1);
$consensus =~ /.+\n(.+)/;
$consensus = $1;

#main loop
for ( my $i = 0; $i < @path+0; $i++ )
{
    #get path elements (contigs) which will be connected
    $path[$i] =~ /(.+)_([53]prime)/;
    $contig1 = $1;
    $contigend1 = $2;
    if ( $path[$i+1] )
    {
    $path[$i+1] =~ /(.+)_([53]prime)/;
    $contig2 = $1;
    $contigend2 = $2;
    }
    
    #skip path element to next different contig
    if ($contig1 eq $contig2)
    {
	next();
    }
    
    #get sequence of elements
    $seq_contig1 = readseq($contigsdir, $contig1);
    $seq_contig2 = readseq($contigsdir, $contig2);

    #get the directory of sequence which connects the two path elements
    $overlapdir = getoverlapdir($dir, $contig1, $contigend1, $contig2, $contigend2);
    
    #read its sequence
    $seq_overlap = readseq($overlapdir, "NODE_1");

    #merge consensus with connecting sequence 
    $consensus = do_it($consensus , $seq_overlap);

    #my $bla = <>;
    
    #merge connecting sequnece with next path element (contig)
    $consensus = do_it($consensus , $seq_contig2);
    
    #my $bla = <>;
}

open(my $OUT, '>', 'pl_cons.fa') or die "Could not open file 'pl_cons.fa' $!";
print $OUT ">Plastid_Consensus_$basename\n$consensus";
close $OUT;



sub do_it
{
    my $consensusendup = substr($_[0], length($_[0])-1000, 1000);
    #$output = align( $consensusend, $seq_overlap );
    my $consensusenddown = substr($_[0], 0, 1000);
    my $output = align( $consensusendup, $consensusenddown, $_[1] );
    my $bestalign = find_bestalign();
    print STDERR "choosing $bestalign\n";
    if ( $bestalign =~ /up/ )
    {
	return substr($consensus, 0, length($consensus)-1000).find_consens($bestalign);
    }
    else
    {
	return find_consens($bestalign).substr($consensus, 1000, length($consensus));
    }
}

sub find_bestalign
{
    my %scores;
    my $up_score;
    my $up_score_rev;
    my $down_score;
    my $down_score_rev;

    open(my $FH, '<', 'needle_out_up_rev.fa') or die "Could not open file 'needle_out_up_rev.fa' $!";
    while ( <$FH> )
    {
	if ( $_ =~ /Score: (.+)/ )
	{
	    $scores{"needle_out_up_rev.fa"} = $1;
	}
    }
    close $FH;
    open($FH, '<', 'needle_out_up.fa') or die "Could not open file 'needle_out_up.fa' $!";
    while ( <$FH> )
    {
	if ( $_ =~ /Score: (.+)/ )
	{
	    $scores{"needle_out_up.fa"} = $1;
	}
    }
    close $FH;
    open($FH, '<', 'needle_out_down.fa') or die "Could not open file 'needle_out_down.fa' $!";
    while ( <$FH> )
    {
        if ( $_ =~ /Score: (.+)/ )
        {
            $scores{"needle_out_down.fa"} = $1;
	}
    }
    close $FH;
    open($FH, '<', 'needle_out_down_rev.fa') or die "Could not open file 'needle_out_down_rev.fa' $!";
    while ( <$FH> )
    {
        if ( $_ =~ /Score: (.+)/ )
        {
            $scores{"needle_out_down_rev.fa"} = $1;
	}
    }
    close $FH;


    
    print STDERR Dumper(%scores);

    my @sorted = sort {$scores{$b} <=> $scores{$a}} keys %scores;

    return $sorted[0];

}
    
sub getoverlapdir
{
    my $overlapdir;
    if ( -e $_[0]."/merged_ass/contigs.fasta_renamed.".$_[1]."_".$_[2]."_ass/contigs.fasta_renamed.".$_[3]."_".$_[4]."/extended_asc/scaffolds.fasta" )
    {
	$overlapdir = $_[0]."/merged_ass/contigs.fasta_renamed.".$_[1]."_".$_[2]."_ass/contigs.fasta_renamed.".$_[3]."_".$_[4]."/extended_asc/scaffolds.fasta";
    }
    elsif ( -e $_[0]."/merged_ass/contigs.fasta_renamed.".$_[3]."_".$_[4]."_ass/contigs.fasta_renamed.".$_[1]."_".$_[2]."/extended_asc/scaffolds.fasta" )
    {
	$overlapdir = $_[0]."/merged_ass/contigs.fasta_renamed.".$_[3]."_".$_[4]."_ass/contigs.fasta_renamed.".$_[1]."_".$_[2]."/extended_asc/scaffolds.fasta";
    }
    else
    {
	die "Can not find contig merging $_[1]_$_[2] with $_[3]_$_[4]";
    }

    return $overlapdir;
}

sub readseq
{
    print "Reading $_[1] from $_[0]\n";
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
    open(my $FH, '>', 'forneedle_up.fa') or die "Could not open file 'forneedle1.fa' $!";
    print $FH $_[0];
    close $FH;
    open($FH, '>', 'forneedle_down.fa') or die "Could not open file 'forneedle2.fa' $!";
    print $FH $_[1];
    close $FH;
    open($FH, '>', 'forneedle_overlap.fa') or die "Could not open file 'forneedle2.fa' $!";
    print $FH $_[2];
    close $FH;

    my $output = qx(needle -asequence forneedle_up.fa -bsequence forneedle_overlap.fa -aformat3 markx3 -outfile needle_out_up.fa -gapopen 10 -gapextend 0.5);
    $output = qx(needle -asequence forneedle_up.fa -sreverse_bsequence forneedle_overlap.fa -aformat3 markx3 -outfile needle_out_up_rev.fa -gapopen 10 -gapextend 0.5);

    $output = qx(needle -asequence forneedle_down.fa -bsequence forneedle_overlap.fa -aformat3 markx3 -outfile needle_out_down.fa -gapopen 10 -gapextend 0.5);
    $output = qx(needle -asequence forneedle_down.fa -sreverse_bsequence forneedle_overlap.fa -aformat3 markx3 -outfile needle_out_down_rev.fa -gapopen 10 -gapextend 0.5);
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
