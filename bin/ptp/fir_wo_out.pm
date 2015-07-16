package fir_wo_out;

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

sub get_path
{
#alle ordner in hash damit keine doppelten
find(sub { $alldirs{ $File::Find::dir } = 0; },  $dir."/merged_ass" );

#GRAPH aufbauen
#testen ob reads auf gleiches ende mappen und dann contigpaarung speichern und contig auftreten zählen
foreach (keys(%alldirs))
{
    if ( $_ =~ /.+\.(.+)_([53]prime).+\.(.+)_([53]prime)$/ ) #$1=contig1, $2=ende contig1, §3=contig2, $4 = ende contig2
    {
	print "$_\n";

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





#PFAD suchen

my @knoten = $g->vertices();
my @nachbarn = ();
my @nachbarsnachbarn = ();

#isolierte knoten(paare) suchen und löschen

for ( my $i = 0; $i < 1; $i++)
{
    foreach (@knoten)
    {
	@nachbarn = $g->neighbours($_);
	@nachbarsnachbarn = $g->neighbours($nachbarn[0]);
	if ( (@nachbarn+0 == 1 &&  @nachbarsnachbarn+0 == 1) || @nachbarn+0 == 0 )
	{
	    print"lösche Knoten $_. Nur ";
	    print(@nachbarn+0);
	    print" Kanten\n";
	    $g->delete_vertex($_);
	}
	
    }
}

@knoten = $g->vertices();
my $curknoten = $knoten[0];
my @path = ($curknoten);
my $last;
my @last2 = ( 1 , 2 );
my $curknotenend;
my $nachbarend;
my %pathsgone;
my $nngone = 0;
my $lastnngone = 1;
my @paths;
my %allpathsnodes;
my $gleich = 0;

$allpathsnodes{$curknoten}++;

#wie viele knoten müssen durchlaufen werden?
my $nknoten = @knoten+0;

print"\n==================================\nStarte mit Knoten $curknoten\n==================================\n";
    
while ( $nknoten >= $nngone + 1 )
{
    #wenn wenn schon pfade vorhanden neuen anfangsknoten wählen
    if (@paths && $nngone == 0)
    {
	print"Wähle neuen Startknoten\n";
	foreach (@knoten)
	{
	    unless (exists($allpathsnodes{$_}))
	    {
		$curknoten = $_;
		push(@path, $curknoten);
		$allpathsnodes{$curknoten}++;                                                                                                              
		@last2 = ( 1 , 2 );
		print"\n==================================\nStarte mit Knoten $curknoten\n==================================\n";
		last();
	    }
	}
    }
    
    if ( $last2[0] eq $last2[1] )
    {
	$curknoten = $last2[0];
	@path = ();
	push(@path, $curknoten);
	$allpathsnodes{$curknoten}++;
	@last2 = ( 1 , 2 );
	print"\n==================================\nSackgasse: Starte neu mit Knoten $curknoten\n==================================\n";
    }

    
    #erstmal zum andren contigende
    if ( $curknoten =~ /(.+)5prime$/ )
    {
	$pathsgone{$curknoten}{$1."3prime"}++;
	$last = $curknoten;
	pop(@last2);
	unshift(@last2, $curknoten);
	$curknoten = $1."3prime";
	push(@path, $curknoten);
	$allpathsnodes{$curknoten}++;
    }
    elsif ( $curknoten =~ /(.+)3prime$/ )
    {
	$pathsgone{$curknoten}{$1."5prime"}++;
	$last = $curknoten;
	pop(@last2);
	unshift(@last2, $curknoten);
	$curknoten = $1."5prime";
	push(@path, $curknoten);
	$allpathsnodes{$curknoten}++;
    }
    
    print"Aktueller Knoten: $curknoten\n";
    
    #contigende feststellen um zu passendem ende zu gehen 
    $curknoten =~ /.+([35])prime$/;
    $curknotenend = $1;
    
    @nachbarn = $g->neighbours($curknoten);
   
    print"Wähle aus Nachbarn: @nachbarn\n";    
    print"schon gewählte Pfade:\n";
    print Dumper($pathsgone{$curknoten});

    #nachbarn durchsuchen nach nächstem knoten
    foreach (@nachbarn)
    {
	$_ =~ /.+([35])prime$/;
	$nachbarend = $1;
	
	#knoten auswählen wenn es nicht der letzte knoten ist und der knoten nicht von hier schonmal besucht wurde
	if ( $_ ne $last && !exists($pathsgone{$curknoten}{$_}) ) # && $nachbarend ne $curknotenend)
	{
	    $pathsgone{$curknoten}{$_}++;
	    $last = $curknoten;
	    pop(@last2);
	    unshift(@last2, $_);
	    $curknoten = $_;
	    #print Dumper(%pathsgone);
	    push(@path, $curknoten);
	    $allpathsnodes{$curknoten}++;
	    last();
	}
    }
   

    $nngone = keys(%allpathsnodes);

    my @npathsgone = keys %{$pathsgone{$curknoten}};
    print(@npathsgone+0);
    print"\n";
    print(@nachbarn+0);
    print"\n";

    #pfadsuche abbrechen wenns nicht mehr weiter geht
    if ( @nachbarn+0 == @npathsgone+0 && @nachbarn+0 != 0 && @npathsgone+0 != 0 )
    {
	print"DAMN!\n";
	push(@paths, [@path]);
	@path = ();
	$nngone = 0;
    }

    

 
    print"Wähle $curknoten\n";
    
    print"Pfad bis jetzt: @path\n";

    print "\nLetzte beiden: @last2 \n";
    my $bla = <>;

    #wenn im kreis dann abbrechen

    print"\nKnoten abgearbeitet: $nngone von $nknoten\n";

    # if ( $nngone == $lastnngone )
    # {
    # 	$gleich++;
    # }
    
    # if ( $gleich == 10 )
    # {
    # 	print"Läuft im Kreis!\n";
    # 	push(@paths, [@path]);
    # 	@path = ();
    # 	$nngone = 0;
    # 	$gleich = 0;
    # }
    
    # $lastnngone = $nngone;
}

#$nknoten = $nknoten - $nngone;
push(@paths, [@path]);
@path = ();

#pfade nach länge sortieren
my @paths_sort = sort {@{$b}+0 <=> @{$a}+0} @paths;

print Dumper(@paths_sort);

return @paths_sort;

}
#print"\n@path\n";

return 1;
