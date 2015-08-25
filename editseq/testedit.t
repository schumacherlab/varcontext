#
#===============================================================================
#
#         FILE: testedit.t
#
#  DESCRIPTION: 
#
#        FILES: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Arno Velds (), a.velds@nki.nl
#      COMPANY: NKI
#      VERSION: 1.0
#      CREATED: 04/11/2012 11:08:32 AM
#     REVISION: ---
#===============================================================================

use strict;
use warnings;

use Test::More tests => 20;                      # last test to print

BEGIN	{
	use_ok( 'EditSeq' );
}

#simple single edits
my $s = "abcdefghij";
my $es = EditSeq->new($s);
$es->edit_substitute("123",5);
my ($seq, $map) = $es->editedseq;
print join("\n", map {join("=>", $_, $map->[$_])} 1 .. $#$map) ,"\n";
is($seq, "abcd123hij", "substitute");


$s = "abcdefghij";
$es = EditSeq->new($s);
$es->edit_substitute("12345",8);
is($es->editedseq, "abcdefg123", "substitute beyond end");

$s = "abcdefghij";
$es = EditSeq->new($s);
$es->edit_delete("def",4);
is($es->editedseq, "abcghij", "delete");

$s = "abcdefghij";
$es = EditSeq->new($s);
$es->edit_delete("fghijkl",6);
is($es->editedseq, "abcde", "delete beyond end");

$s = "abcdefghij";
$es = EditSeq->new($s);
$es->edit_insert("123",4);
is($es->editedseq, "abcd123efghij", "insert");

$s = "abcdefghij";
$es = EditSeq->new($s);
$es->edit_insert("123",10);
is($es->editedseq, "abcdefghij123", "insert at end");

$s = "abcdefghij";
$es = EditSeq->new($s);
$es->edit_complex("123",4, "def");
is($es->editedseq, "abc123ghij", "even length complex(=substitution)");

$s = "abcdefghij";
$es = EditSeq->new($s);
$es->edit_complex("1234",4, "de");
is($es->editedseq, "abc1234fghij", "uneven complex");




#multiple edits
$s = "abcdefghij";
$es = EditSeq->new($s);
$es->edit_substitute("12",5);
$es->edit_delete("defgh",4);
is($es->editedseq, "abcij", "delete substitution");

$es = EditSeq->new($s);
$es->edit_substitute("12",5);
$es->edit_delete("bcd",2);
is($es->editedseq, "a12ghij", "sub flanking del");

$es = EditSeq->new($s);
$es->edit_substitute("12",5);
$es->edit_delete("h",8);
is($es->editedseq, "abcd12gij", "sub after del");

$es = EditSeq->new($s);
$es->edit_insert("12",5);
$es->edit_delete("h",8);
my ($seq, $map) = $es->editedseq;
print join("\n", map {join("=>", $_, $map->[$_])} 1 .. $#$map) ,"\n";
is($seq, "abcde12fgij", "del after insert");

$es = EditSeq->new($s);
$es->edit_insert("12",5);
$es->edit_delete("f",6);
is($es->editedseq, "abcde12ghij", "del flanking insert");

$es = EditSeq->new($s);
$es->edit_substitute("Q",6);
$es->edit_complex("123",3, "cd");
is($es->editedseq, "ab123eQghij", "complex before substitution");

$es = EditSeq->new($s);
$es->edit_substitute("Q",6);
$es->edit_complex("123",7, "ghi");
is($es->editedseq, "abcdeQ123j", "complex flankin substitution");

$es = EditSeq->new($s);
$es->edit_delete("fgh",6);
$es->edit_complex("123",3, "cd");
is($es->editedseq, "ab123eij", "complex before delete");

$es = EditSeq->new($s);
$es->edit_delete("bc",2);
$es->edit_complex("123",6, "fg");
is($es->editedseq, "ade123hij", "complex after delete");

#try an overlap
$es = EditSeq->new($s);
$es->edit_insert("12",5);
$es->edit_delete("defg",4); #overlap causes this to be discarded
is($es->editedseq, "abcde12fghij", "overlap1");

$es = EditSeq->new($s);
$es->edit_delete("defg",4);
$es->edit_insert("12",5); #overlap causes this to be discarded
is($es->editedseq, "abchij", "overlap2");



