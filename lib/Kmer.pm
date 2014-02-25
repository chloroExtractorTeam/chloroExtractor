package Kmer;

use warnings;
use strict;

use Verbose;

our $VERSION = '0.02';


##------------------------------------------------------------------------##

=head1 NAME 

Kmer.pm

=cut

=head1 DESCRIPTION

Extract kmers from sequence strings.

=cut

=head1 SYNOPSIS

  $seq = "ATTATATATCGACTAGCC";
  $kh = Kmer->new(kmer_size => 4);
  
  # get every kmer
  $kh->kmerize($seq);
	ATTA
	TTAT
	TATA
	ATAT
	TATA
	ATAT
	TATC
	ATCG
	TCGA
	CGAC
	GACT
	ACTA
	CTAG
	TAGC
	AGCC

  # get every 4th kmer in canonical representation
  $kh->shift_by(4); 
  $kh->cmerize($seq); 
	'ATTA',	# TAAT
	'TATA',	# TATA
	'TCGA',	# TCGA
	'CTAG',	# CTAG

=cut

=head1 CHANGELOG

=cut

=head2 0.02

=over

=item [Feature] Cache unpack templates for different sequence lengths to 
 boost performance.

=item [Refactoring] Redesigned calculation for "number of kmers" to handle
 'shift_by'

=item [Feature] 'shift_by' bytes between kmers instead of only '1'.

=back

=head2 0.01

=over

=item [BugFix] chomp($seq) in kmerize/cmerize

=item [Refacture] OO handle and OO METHODS instead of Class METHODS

=item [PODfix]

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

##------------------------------------------------------------------------##

=head1 Class METHODS

=cut

##------------------------------------------------------------------------##

=head1 Constructor METHOD

=head2 new

  $kh = Kmer->new(kmer_size => 19);

=cut

sub new{
	my $proto = shift;
	my $self;
	my $class;
	
	# object method -> clone + overwrite
	if($class = ref $proto){ 
		return bless ({%$proto, @_}, $class);
	}

	# class method -> construct + overwrite
	# init empty obj
	$self = {
		kmer_size => undef,
		shift_by => 1,
		@_,
		_u_patt => undef,
		_u_tpl => {}, #  this is a cache for speed
		_u_shift => undef,
		_u_size => undef,
	};
	
	die "kmer_size required" unless $self->{kmer_size};
	die "shift_by needs to be INT >= 1" unless $self->{shift_by} > 0;
	
	
	bless $self, $proto;
	$self->_reset_u;
	
	return $self; 
}



##------------------------------------------------------------------------##

=head1 Object METHODS

=cut

=head2

  $kh->kmer_size()
    # 19
  $kh->kmer_size(4)
    # 4

=cut

sub kmer_size{
	my ($self, $ks) = @_;
	if($ks){
		$self->{kmer_size} = $ks;
		$self->_reset_u;
	};
	return $self->{kmer_size};
}

sub shift_by{
	my ($self, $incr) = @_;
	if (defined ($incr)){
		$self->{shift_by} = $incr;
		$self->_reset_u;
	}
	return $self->{shift_by};
}

=head2 kmerize

Factor a STRING into a LIST of overlapping kmers. Kmers are returned in 
 their literal version.

  $kh->kmerize("ATAGG");
    # ATAG,TAGG

=cut

sub kmerize{
	my ($self,$seq) = @_;
	chomp($seq);
	return unpack($self->{_u_tpl}{length $seq} || $self->_create_u_tpl(length $seq), $seq);
}

=head2 cmerize

Factor a STRING into a LIST of overlapping kmers. Kmers are returned in 
 their canonical representation. The canonical representation is the
 lexically first of the literal and the reverse complement version of any
 kmer.

  $kh->cmerize("ATAGG");
    # ATAG,CCTA

=cut

sub cmerize{
	my ($self,$seq) = @_;
	chomp($seq);
	map{
		my $krc = reverse $_; $krc =~ tr/ATGC/TACG/; $_ gt $krc ? $krc : $_
	}unpack($self->{_u_tpl}{length $seq} || $self->_create_u_tpl(length $seq), $seq);
}

##------------------------------------------------------------------------##
=head1 Private METHODS

=cut

=head2 _u_shift

=cut

sub _u_shift{
	my $self = shift;
	$self->{_u_shift} = $_[0] if defined $_[0];
	return $self->{_u_shift};
}

=head2 _u_patt

=cut

sub _u_patt{
	my $self = shift;
	$self->{_u_patt} = $_[0] if defined $_[0];
	return $self->{_u_patt};
}

=head2 _u_tpl

=cut

sub _u_tpl{
	my $self = shift;
	$self->{_u_tpl} = $_[0] if defined $_[0];
	return $self->{_u_tpl};
}

=head2 _u_size

=cut

sub _u_size{
	my $self = shift;
	$self->{_u_size} = $_[0] if defined $_[0];
	return $self->{_u_size};
}



=head2 _reset_u

=cut

sub _reset_u{
	my $self = shift;
	$self->{_u_shift} = - ($self->kmer_size - $self->shift_by);
	$self->{_u_patt} = $self->_u_shift < 0
		? sprintf("(A%dX%d)%%dA%1\$d", $self->kmer_size, abs($self->_u_shift))
		: sprintf("(A%dx%d)%%dA%1\$d", $self->kmer_size, $self->_u_shift);
	$self->{_u_tpl} = {};
	$self->{_u_size} = $self->kmer_size + $self->_u_shift;
}

=head2 _create_u_tpl

=cut

sub _create_u_tpl{
	my ($self,$len) = @_;
	my $rep = int( ($len + $self->_u_shift) / $self->_u_size ) - 1;
	return $self->{_u_tpl}{$len} = sprintf($self->_u_patt, $rep);
}

=head1 AUTHOR

Thomas Hackl S<thomas.hackl@uni-wuerzburg.de>

=cut



1;



