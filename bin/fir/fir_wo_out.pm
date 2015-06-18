package fir;

use strict;
use warnings;
use File::Find;
use Cwd;
use Data::Dumper;
use Graph;

my %alldirs = ();
my $dir = cwd();
my $g = Graph::Undirected->new();
my $n2inv;
my $n4inv;
my $edge1;
my $edge2;
my $edge3;

#alle ordner in hash damit keine doppelten
find(sub { $alldirs{ $File::Find::dir } = 0; },  $dir."/merged_ass" );

#GRAPH aufbauen
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

	#Knoten hinzufügen
	$g->add_vertices( $1."_".$2, $3."_".$4, $1."_".$n2inv, $3."_".$n4inv );


	#knoten verbinden und mit anderem contigende verbinden
	unless ($g->has_edge($1."_".$2, $3."_".$4) || $g->has_edge($3."_".$4, $1."_".$2)) {$edge1 = $g->add_edge( $1."_".$2, $3."_".$4 )};
	if ($edge1) {$edge1->undirected(1)};
	unless ($g->has_edge($1."_".$n2inv, $1."_".$2) || $g->has_edge($1."_".$2, $1."_".$n2inv)) {$edge2 = $g->add_edge( $1."_".$n2inv, $1."_".$2 )};
	if ($edge2) {$edge2->undirected(1)};
	unless ($g->has_edge($3."_".$n4inv, $3."_".$4) || $g->has_edge($3."_".$4, $3."_".$n4inv)) {$edge3 = $g->add_edge( $3."_".$n4inv, $3."_".$4 )};
	if ($edge3) {$edge3->undirected(1)};

    }
}


my @knoten = $g->vertices();
my $curknoten = $knoten[0];
my @path = ($curknoten);
my @nachbarn = ();
my $last;
my $curknotenend;
my $nachbarend;
my %pathsgone;
my %ngone;
my $nngone = 0;
my $lastnngone = 1;
my @paths;
my %allpathsnodes;


$allpathsnodes{$curknoten}++;
$ngone{$curknoten}++;

#wie viele knoten müssen durhclaufen werden?
my $nknoten = @knoten+0;

print"==================================\nStarte mit Knoten $curknoten\n==================================\n";
    
while ( $nknoten >= $nngone + 1 )
{
    #wenn wenn schon pfade vorhanden neuen anfangsknoten wählen
    if (@paths && $nngone == 0)
    {
	print"Wähle neuen Knoten\n";
	foreach (@knoten)
	{
	    unless (exists($allpathsnodes{$_}))
	    {
		$curknoten = $_;
		push(@path, $curknoten);
		$allpathsnodes{$curknoten}++;                                                                                                              
		print"Neuer knoten: $curknoten\n";
		print"==================================\nStarte mit Knoten $curknoten\n==================================\n";
		last();
	    }
	}
    }
    

    
    #erstmal zum andren contigende
    if ($curknoten =~ /(.+)5prime$/)
    {
	$last = $curknoten;
	$curknoten = $1."3prime";
	push(@path, $curknoten);
	$allpathsnodes{$curknoten}++;
	$ngone{$curknoten}++;
    }
    elsif ($curknoten =~ /(.+)3prime$/)
    {
	$last = $curknoten;
	$curknoten = $1."5prime";
	push(@path, $curknoten);
	$allpathsnodes{$curknoten}++;
	$ngone{$curknoten}++;
    }
    
    print"Aktueller Knoten: $curknoten\n";
    
    #contigende feststellen um zu passendem ende zu gehen 
    $curknoten =~ /.+([35])prime$/;
    $curknotenend = $1;
    
    @nachbarn = $g->neighbours($curknoten);
    print"Nachbarn: @nachbarn\n";
    
    #$pathsgone{$curknoten};
    
    #nachbarn durchsuchen nach nächstem knoten
    foreach (@nachbarn)
    {
	$_ =~ /.+([35])prime$/;
	$nachbarend = $1;
	
	#knoten auswählen wenn es nicht der letzte knoten ist und das ende passt und der knoten von hier schonmal besucht wurde
	if ( $_ ne $last && $curknotenend != $nachbarend && !exists($pathsgone{$curknoten}{$_}) )
	{
	    $pathsgone{$curknoten}{$_}++;
	    $curknoten = $_;
	    print Dumper(%pathsgone);
	    push(@path, $curknoten);
	    $allpathsnodes{$curknoten}++;
	    $ngone{$curknoten}++;
	    last();
	}
    }
    
    print"Wähle $curknoten\n";
    
    #wenn im kreis dann abbrechen
    $nngone = keys(%allpathsnodes);
    
    if ( $nngone == $lastnngone )
    {
	print"Läuft im Kreis!\n";
	#noch offene knoten
	#$nkoten = $nknoten - $nngone;
	push(@paths, [@path]);
	@path = ();
	$nngone = 0;
	%ngone = ();
	next();
    }
    
    $lastnngone = $nngone;
    
    #print Dumper(%ngone);
    print"Pfad bis jetzt: @path\n";
    print"\nKnoten abgearbeitet: $nngone von $nknoten\n";
    my $pro = <>;
}

#$nknoten = $nknoten - $nngone;
push(@paths, [@path]);
@path = ();


print Dumper(@paths);

#print"\n@path\n";
