# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl deldir.t'

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use strict;
use warnings;
use File::Copy::Recursive qw(fcopy rcopy dircopy fmove rmove dirmove);
use File::Path qw(remove_tree);
use Cwd;

#prepare testfiles

my $dir = getcwd;
remove_tree( 't/emptyfiles/' , {verbose => 1}) or die "Cannot remove_tree 't/emptyfiles/' : $!";

warn("$dir\n");
dircopy("t/testdata/emptyfiles","$dir"."/t/emptyfiles") or die "Copy failed: $!";

use Test::More;# tests => 1;
BEGIN { use_ok('deldir') };

#########################



# Insert your test code below, the Test::More module is use()ed here so read
# its man page ( perldoc Test::More ) for help writing this test script.


can_ok('deldir', ('getemptyfilesdir'));




done_testing();