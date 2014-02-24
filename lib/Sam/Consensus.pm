package Sam::Consensus;

use warnings;
use strict;

# $Id$

# preference libs in same folder over @INC
use lib '../';

use Sam::Parser;
use Sam::Alignment qw(:flags);
use Verbose;

our $VERSION = '0.02';
our ($REVISION) = '$Revision$' =~ /(\d+)/;
our ($MODIFIED) = '$Date$' =~ /Date: (\S+\s\S+)/;



=head1 NAME 

Sam::Consensus.pm

=head1 DESCRIPTION

Module to create consensus fastq sequence from SAM file.

=cut

=head1 SYNOPSIS


=cut

=head1 CHANGELOG

=head2 0.02

=over

=item [Change] Preference libs in same folder over @INC

=item [Change] Added svn:keywords

=back

=over

=item 0.01

Initial Consensus module. Provides Constructor, generic accessor 
 and consensus call method

=back

=cut

=head1 TODO

=over

=item Tests

=back

=cut


############################################################################

# Init a few globals, e.g. class attributes

our $V;



############################################################################


=head1 Constructor METHOD

=head2 new

Initialize a sam consensus object. Takes parameters in key => value 
 format. Calculates the pileup matrix from the sam file.

  $consensus = $Sam::Consensus->new();

=cut

sub new{
	my $class = shift;
	
	# defaults=
	my $self = {
		fh => undef,
		file => undef,
		rname => undef,
		sp => undef,
		@_
	};
	
	unless ($self->{sp}){	
		$self->{sp} = Sam::Parser->new( 
			file => $self->{file},
			fh => $self->{fh},
		);
	};
	
	bless $self, $class;
	
	$self->_state_matrix;
}


############################################################################

=head1 Class ATTRIBUTES

=cut

=head2 $V

Verbose messages are handled using the Verbose.pm module. To 
 customize verbose message behaviour, overwrite the attribute with
 another Verbose object created with the Verbose module.

=cut

$V = Verbose->new();

=head2 $MaxCoverage [100]

=cut

our $MaxCoverage = 100;

=head2 $PhredOffset [33]

=cut

our $PhredOffset = 33;


############################################################################

=head1 Object METHODS

=cut

=head2 consensus

Calculate and the consensus sequence from state matrix of the Sam::Consensus
 object. In SCALAR context, returns the consensus object, wit consensus 
 sequence and quality stored in the object C<< $con->seq, $con->qual >>, in
 LIST context returns C<< rname (id), seq, qual >>, which can be easily 
 concatenated sequence strings of any desired format.

=cut

sub consensus{
	my ($self) = @_;
	
	$self->_consensus;
	return wantarray ? ($self->rname, $self->seq, $self->qual) : $self;
}


=head2 variants

=cut

sub variants{
	my $self = shift;
	my %p = (
		min_prob => 0.1,
		accuracy => 5,
		@_
	);
	
	
	#print Dumper($self->{_state_matrix});
	my @seq;
	my %states_rev = reverse %{$self->{_states}}; # works since values are also unique
	
	foreach my $col (@{$self->{_state_matrix}}){
		# cov
		unless($col){
			push @{$self->{covs}}, 0;
			push @{$self->{vars}}, ['?'];
			push @{$self->{freqs}},[''];
			push @{$self->{probs}},[''];
			next;
		}
		my $cov;
		my %vars;
		# variants
		for(my $i=0; $i<@$col; $i++){
			if (defined(my $v = $col->[$i])){
				$cov+= $v; 
				$vars{$states_rev{$i}} = $v;
			};
		};
		push @{$self->{covs}}, $cov;
		my @vars = sort{$vars{$b} <=> $vars{$a}}keys %vars;
		my @freqs = @vars{@vars};
		my @probs = map{sprintf("%0.".$p{accuracy}."f", $_/$cov)}@freqs;
		my $k = grep{$_>= $p{min_prob}}@probs;
		$k--;
		push @{$self->{vars}}, $k >= 0 ? [@vars[0..$k]] : ['?'];
		push @{$self->{freqs}}, $k >= 0 ? [@freqs[0..$k]] : [''];
		push @{$self->{probs}}, $k >= 0 ? [@probs[0..$k]] : [''];
	}
	# rel freq

	return $self;
}

############################################################################

=head1 Accessor METHODS

=cut

=head2 rname

Get/Set the rname.

=cut

sub rname{
	my ($self, $rname) = @_;
	$self->{rname} = $rname if $rname;
	return $self->{rname};
}

=head2 seq

Get/Set the seq.

=cut

sub seq{
	my ($self, $seq) = @_;
	$self->{seq} = $seq if $seq;
	return $self->{seq};
}

=head2 qual

Get/Set the seq.

=cut

sub qual{
	my ($self, $qual) = @_;
	$self->{qual} = $qual if $qual;
	return $self->{qual};
}

############################################################################

=head1 Private METHODS

=cut

sub _state_matrix{
	my $self = shift;
	
	# state matrix
	my @S;
	# predefined states
	my %states = (
		A => 0,
		T => 1,
		G => 2,
		C => 3,
		'-' => 4,
		N => 5,
		# .. complex states, dynamically added
	);
	
	
	while(my $aln = $self->{sp}->next_aln()){
		# skip all alignments of other references
		if ($self->rname){
			next unless $aln->rname eq $self->rname
		# if no rname is provided, use the first rname found in sam
		}else{
			$self->rname($aln->rname);
		}
		
		# get read seq
		my $seq = $aln->seq;
		
		# get read cigar
		my @cigar = split(/(\d+)/,$aln->cigar);
		shift @cigar;
		
		# reference position
		my $rpos = $aln->pos-1;
		
		my $state; # buffer last match, required if followed by insertion
		for(my $i=0; $i<@cigar;$i+=2){
			if($cigar[$i+1] eq 'M'){
				my @subseq = split(//,substr($seq,0,$cigar[$i],''));
				foreach $_ (@subseq){
					($S[$rpos][$states{$_}])++;  # match states always exist
					$rpos++;
				}
				$state = $subseq[$#subseq];
			}elsif($cigar[$i+1] eq 'D'){
				for(1..$cigar[$i]){
					($S[$rpos][4])++;  # $states{'-'} is always 4 
					$rpos++;
				}
				$state = '-';
			}elsif($cigar[$i+1] eq 'I'){
				#unless ($state){print STDERR $aln->pos," : ",$rpos,"\n"} 
				my $complex_state;
				if($state){
					$complex_state = $state.substr($seq,0,$cigar[$i],'');
					($S[$rpos-1][$states{$state}])--; #
				}else{
					$complex_state = substr($seq,0,$cigar[$i],'');
				}
				# replace by complex state, add state idx to %states if new
				if(exists ($states{$complex_state})){
					#TODO: insertion before first M
					next if ($rpos-1 < 0);
					$S[$rpos-1][$states{$complex_state}]++
				}else{
					next if ($rpos-1 < 0);
					$states{$complex_state} = scalar keys %states;
					$S[$rpos-1][$states{$complex_state}]++;
				}
				#($S[$rpos-1][exists ($states{$complex_state}) ? $states{$complex_state} : $states{$complex_state} = keys %states])++; 
				#$seq[$#seq].= 
			}else{
				$V->exit("Unknown Cigar '".$cigar[$i+1]."'");
			}
		}
	}
	
	# return state matrix
	$self->{_state_matrix} = \@S;
	$self->{_states} = \%states;
	return $self;
	
}

sub _consensus{
	my $self = shift;
	my %states_rev = reverse %{$self->{_states}}; # works since values are also unique
	my $seq;
	my $qual;
	
	foreach my $col (@{$self->{_state_matrix}}){
		# uncovered col
		unless ($col){
			$seq.='n';
			$qual.=chr(0+$PhredOffset);
			next;
		}
		
		my $idx=undef;
		my $v=0;
		my $i;
		# majority vote
		for($i=0; $i<@$col; $i++){
			if (defined($col->[$i]) && $col->[$i] > $v){
				$v = $col->[$i];
				$idx = $i; 
			};
		};
		
		next if $idx == 4; # insertion on reference
		
		# uncovered wouldn't go this far because $col is undef and not array ref
		# my $con = defined ($idx) ? $states_rev{$idx} : '?';
		
		my $con = $states_rev{$idx};
		$seq.= $con;
		$qual.= _phred($col) x length($con);
	}
	$self->{seq} = $seq;
	$self->{qual} = $qual;
	
	return $self;
}


sub _phred{
	my $col = shift;
	my $total = 0;
	my @states = grep{$_}@$col;
	$total += $_ for @states;
	my @probability = map{$_/$total}@states;
	my $Hx;
	$Hx += $_ for map{$_ * log $_}@probability;
	$Hx = -1 if $Hx < -1;
	$Hx = 0 if $Hx > 0;
	# prevent phreds >40 if actual coverage > $MaxCoverage
	# TODO: just a sloppy fix...
	$total = $MaxCoverage if $total > $MaxCoverage;
	return chr(int(((1+$Hx)/($MaxCoverage/$total))*40)+$PhredOffset);
}


=head1 AUTHOR

Thomas Hackl S<thomas.hackl@uni-wuerzburg.de>

=cut



1;
