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

dircopy("t/testdata/emptyfiles","$dir"."/t/emptyfiles") or die "Copy failed: $!";

use Test::More;# tests => 1;
BEGIN { use_ok('deldir') };

#########################



# Insert your test code below, the Test::More module is use()ed here so read
# its man page ( perldoc Test::More ) for help writing this test script.

can_ok('deldir', ('getemptyfiledirs'));

my @expectedemptyfiles = ("empty1/file1" , "empty1/file2" , "empty2/file1" , "empty2/file2" , "empty3/file1" , "empty3/file2");

my $testdir = "t/emptyfiles/";

my @emptyfiles = deldir::getemptyfiledirs($testdir);

is_deeply( [sort @emptyfiles], [sort @expectedemptyfiles] , 'Are found empty files correct?');

my $filesize = 0;

foreach my $file (@emptyfiles)
{
	ok( -e "$testdir$file" , "Is $file existing?" );
	$filesize = -s "$file";
	ok( $filesize == 0 , "Is $file empty?" );
}

can_ok('deldir', ('removeemptyfiledirs'));

deldir::removeemptyfiledirs(\@emptyfiles);


done_testing();