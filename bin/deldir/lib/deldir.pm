package deldir;

use 5.014002;
use strict;
use warnings;
use File::Find;

require Exporter;

our @ISA = qw(Exporter);

# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.

# This allows declaration	use deldir ':all';
# If you do not need this, moving things directly into @EXPORT or @EXPORT_OK
# will save memory.
our %EXPORT_TAGS = ( 'all' => [ qw(
	
) ] );

our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

our @EXPORT = qw(
	
);

our $VERSION = '0.01';

=pod
=head2 sub add

Expects directory as parameter and returns a list of paths to empty files in dir and subdirs.

=cut


sub getemptyfiledirs
{


    my $dir = $_[0];

    our @allfiles = ();
    
    #find all dirs and files
    find(\&wanted,  $dir);

    sub wanted
    {
	push( @allfiles , "$File::Find::name");
    }


my @emptyfiles = ();
my $dir = $_[0];

#filter out not empty paths and files
foreach my $file (@allfiles)
{
    if ( -s "$file" == 0)
    {
	#remove searchdir from found paths
	$file =~ /^$dir(.+)/;
	push( @emptyfiles , $1 )
    }
	
}

return @emptyfiles;

}

1;
__END__
# Below is stub documentation for your module. You'd better edit it!

=head1 NAME

deldir - Perl extension for deleting empty files and their paths

=head1 SYNOPSIS

  use deldir;
  blah blah blah

=head1 DESCRIPTION

Stub documentation for deldir, created by h2xs. It looks like the
author of the extension was negligent enough to leave the stub
unedited.

Blah blah blah.

=head2 EXPORT

None by default.



=head1 SEE ALSO

Mention other useful documentation such as the documentation of
related modules or operating system documentation (such as man pages
in UNIX), or any relevant external documentation such as RFCs or
standards.

If you have a mailing list set up for your module, mention it here.

If you have a web site set up for your module, mention it here.

=head1 AUTHOR

Maik Guendel, E<lt>s319932@E<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2015 by Maik Guendel

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.14.2 or,
at your option, any later version of Perl 5 you may have available.


=cut
