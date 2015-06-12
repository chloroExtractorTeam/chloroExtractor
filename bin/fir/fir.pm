package fir;

use strict;
use warnings;
use File::Find;
use Cwd;
use Data::Dumper;
use Graph::Easy;

my %alldirs = ();
my $dir = cwd();
my $g = Graph::Easy->new();
my $n2inv;
my $n4inv;
my $edge1;
my $edge2;
my $edge3;

#alle ordner in hash damit keine doppelten
find(sub { $alldirs{ $File::Find::dir } = 0; },  $dir."/merged_ass" );

#testen ob reads auf gleiches ende mappen und dann contigpaarung speichern und contig auftreten zählen
foreach (keys(%alldirs))
{
    if ( $_ =~ /.+\.(.+)_([53]prime).+\.(.+)_([53]prime)$/ ) #$1=contig1, $2=ende contig1, §3=contig2, $4 = ende contig2
    {
	print STDERR "$_\n";

	if ( $2 eq "5prime")
	{
	    $n2inv = "3prime"
	}
	else
	{
	    $n2inv = "5prime"
	}
	if ( $4 eq "5prime")
	{
	    $n4inv = "3prime"
	}
	else
	{
	    $n4inv = "5prime"
	}

	$edge1 = 0;
	$edge2 = 0;
	$edge3 = 0;


	$g->add_vertices( $1."_".$2, $3."_".$4, $1."_".$n2inv, $3."_".$n4inv );


	unless ($g->edge($1."_".$2, $3."_".$4) || $g->edge($3."_".$4, $1."_".$2)) {$edge1 = $g->add_edge_once( $1."_".$2, $3."_".$4 )};
	if ($edge1) {$edge1->undirected(1)};
	unless ($g->edge($1."_".$n2inv, $1."_".$2) || $g->edge($1."_".$2, $1."_".$n2inv)) {$edge2 = $g->add_edge_once( $1."_".$n2inv, $1."_".$2 )};
	if ($edge2) {$edge2->undirected(1)};
	unless ($g->edge($3."_".$n4inv, $3."_".$4) || $g->edge($3."_".$4, $3."_".$n4inv)) {$edge3 = $g->add_edge_once( $3."_".$n4inv, $3."_".$4 )};
	if ($edge3) {$edge3->undirected(1)};

    }
}

$g->set_attribute('type','undirected');
$g->set_attribute('title', "$dir");
$g->timeout(60);
print $g->as_svg_file();
print STDERR $g->is_directed() . "\n" . $g->is_undirected() . "\n";
#print Dumper(%irmate);
