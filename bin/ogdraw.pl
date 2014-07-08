#!/usr/bin/perl -w
use strict;
use WWW::Mechanize;
use WWW::Mechanize::FormFiller;
use URI::URL;

use File::Spec;
use File::Basename;

use Data::Dumper;

# USAGE perl submit_to_CpGAVAS.pl <chloro.fa> <species>
my $usage = "Usage:\n ogdraw.pl *.gb/*.gff";

my $agent = WWW::Mechanize->new( autocheck => 1 );
my $formfiller = WWW::Mechanize::FormFiller->new();
$agent->env_proxy();

die "File required\n".$usage unless @ARGV;

foreach my $file (@ARGV){
    my $base = basename($gff, qw(.gff .gff3 .gb));

    if($file =~ /\.gff\d?$/){ # convert gff to gb
	my $gff = $file;
	$file = $base.".gb";
	print "Converting gff to genbank ...\n";
	qx( seqret -feature -osformat2 gb -outseq $file $gff );
    }

    

    $file=File::Spec->rel2abs($file);    

    die "File $file not found\n$usage\n" unless(-e $file);

    $agent->get('http://ogdraw.mpimp-golm.mpg.de/cgi-bin/ogdraw.pl');

    $agent->form_number(1) if $agent->forms and scalar @{$agent->forms};
    $formfiller->add_filler( 'seqfile' => Fixed => $file );
    $formfiller->fill_form($agent->current_form);
    print "Submitting $file...\n";
    $agent->submit();

    $agent->form_number(1) if $agent->forms and scalar @{$agent->forms};
    $formfiller->add_filler( 'output_type' => Fixed => "ps" );
    $formfiller->fill_form($agent->current_form);
    print "Submitting options...\n";
    $agent->submit();
    
    my $link=$agent->find_link( text => 'here' )->url_abs;
    print "Retrieving results from $link\n";
    qx(wget $link);
#    print Dumper($agent); 

    
}
