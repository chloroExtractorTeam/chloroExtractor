package Kmer;

# $Id$

use warnings;
use strict;

use Verbose;

our $VERSION = '0.01';
our ($REVISION) = '$Revision$' =~ /(\d+)/;
our ($MODIFIED) = '$Date$' =~ /Date: (\S+\s\S+)/;


##------------------------------------------------------------------------##

=head1 NAME 

Kmer.pm

=head1 DESCRIPTION

Class for handling Kmer stuff.

=head1 SYNOPSIS

=cut

=head1 CHANGELOG

=head2 0.01

=over

=item [Initial]

=back

=cut

=head1 TODO

=over

=back

=cut


##------------------------------------------------------------------------##

=head1 Class Attributes

=cut

our $KmerSize = 19;

our $_Unpack = "(A".$KmerSize."X".($KmerSize-1).")";

##------------------------------------------------------------------------##

=head1 Class METHODS

=cut

=head2 KmerSize

Get/Set the KmerSize the module is currently using

=cut

sub KmerSize{
	my ($class, $ks) = @_;
	if($ks){
		$KmerSize = $ks;
		$_Unpack = "(A".$ks."X".($ks-1).")";
	};
	return $KmerSize;
}


sub Kmerize{
	my ($class,$seq) = @_;
	return unpack($_Unpack.((length $seq) - $KmerSize+1 ), $seq);
}

sub Kmerize_nr{
	my ($class,$seq) = @_;
	map{my $krc = reverse $_; $krc =~ tr/ATGC/TACG/; $_ gt $krc ? $krc : $_}unpack($_Unpack.((length $seq) - $KmerSize+1 ), $seq);
}


##------------------------------------------------------------------------##

=head1 Constructor METHOD

=head2 new

Create a new FASTQ seq object. Either provide "seq_head", "seq", 
 "qual_head", and "qual" as first four parameter or the four lines as one 
 STRING. Additional parameter can be specified in key => value format. 
 C<phred_offset> defaults to undef. 

When used on a Fastq::Seq object, it acts as cloning method, also 
 allowing for additional parameter in key => value format. Does not
 perform deep copy, not required for this object.

  # construct new object
  $fq = Fastq::Seq->new(<4_line_fastq_string>, phred_offset => 33);
  # or
  fq = Fastq::Seq->new(
  	<seq_head_string>,
  	<seq_string>,
  	<qual_head_string>,
    <qual_string>,
  	phred_offset => 33,
  );
  
  # clone object and set a new header at the same time
  $fq_clone = $fq->new(seq_head => 'NEW_HEADER');

=cut

sub new{
	my $proto = shift;
	my $self;
	my $class;
	
	# object method -> clone + overwrite
	if($class = ref $proto){ 
		return bless ({%$proto, @_}, $class);
	}else{ # init empty obj
		$class = $proto;
		$self = {};
	}
	
	# class method -> construct + overwrite
	if(@_){
		if(@_%2){ # create object from string

		}else{ # create object from array

		}
	}
	
	return bless $self, $proto;
}




##------------------------------------------------------------------------##

=head1 Object METHODS

=cut


##------------------------------------------------------------------------##

=head1 Accessor METHODS

=cut

=head1 AUTHOR

Thomas Hackl S<thomas.hackl@uni-wuerzburg.de>

=cut



1;



