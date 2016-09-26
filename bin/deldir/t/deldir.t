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

sub preptest
{
	my $dir = getcwd;
	remove_tree( 't/emptyfiles/' , {verbose => 1}) or die "Cannot remove_tree 't/emptyfiles/' : $!";

	dircopy("t/testdata/emptyfiles","$dir"."/t/emptyfiles") or die "Copy failed: $!";
}

preptest();

use Test::More;# tests => 1;
BEGIN { use_ok('deldir') };

#########################



# Insert your test code below, the Test::More module is use()ed here so read
# its man page ( perldoc Test::More ) for help writing this test script.

can_ok('deldir', ('getemptyfiledirs'));

my @expectedemptyfiles = ("t/emptyfiles/empty1/file1" , "t/emptyfiles/empty1/file2" , "t/emptyfiles/empty2/file1" , "t/emptyfiles/empty2/file2" , "t/emptyfiles/empty3/file1" , "t/emptyfiles/empty3/file2");

my $testdir = "t/emptyfiles/";

my @emptyfiles = deldir::getemptyfiledirs($testdir);

is_deeply( [sort @emptyfiles], [sort @expectedemptyfiles] , 'Are found empty files correct?');

my $filesize = 0;

foreach my $file (@emptyfiles)
{
	is( -e "$file" , 1 , "Is $file existing?" );
	$filesize = -s "$file";
	is( $filesize == 0 , 1 , "Is $file empty?" );
}

can_ok('deldir', ('removeemptyfiles'));

deldir::removeemptyfiles(\@emptyfiles);

foreach my $file (@emptyfiles)
{
        is( -e "$file" , undef , "Is $file existing after removing it?" );
}

can_ok('deldir', ('reremoveemptyfiles'));

my @emptyfolder = ("t/emptyfiles/empty1", "t/emptyfiles/empty2", "t/emptyfiles/empty3");


preptest();

deldir::reremoveemptyfiles("$testdir");

foreach my $file (@emptyfolder)
{
	is( -e "$file", undef , "Are empty dirs removed?" );
}

done_testing();