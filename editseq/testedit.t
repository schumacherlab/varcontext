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

use Test::More tests => 36;                      # last test to print

BEGIN	{
	use_ok( 'EditSeq' );
}

#simple single edits
my $s = "abcdefghij";
my $es = EditSeq->new($s);
$es->edit_substitute("123",5);
is($es->apply_edits, "abcd123hij", "substitute");


$s = "abcdefghij";
$es = EditSeq->new($s);
$es->edit_substitute("12345",8);
is($es->apply_edits, "abcdefg123", "substitute beyond end");

$s = "abcdefghij";
$es = EditSeq->new($s);
$es->edit_delete("def",4);
is($es->apply_edits, "abcghij", "delete");

$s = "abcdefghij";
$es = EditSeq->new($s);
$es->edit_delete("fghijkl",6);
is($es->apply_edits, "abcde", "delete beyond end");

$s = "abcdefghij";
$es = EditSeq->new($s);
$es->edit_insert("123",4);
is($es->apply_edits, "abcd123efghij", "insert");

$s = "abcdefghij";
$es = EditSeq->new($s);
$es->edit_insert("123",10);
is($es->apply_edits, "abcdefghij123", "insert at end");

$s = "abcdefghij";
$es = EditSeq->new($s);
$es->edit_complex("123",4, "def");
is($es->apply_edits, "abc123ghij", "even length complex(=substitution)");

$s = "abcdefghij";
$es = EditSeq->new($s);
$es->edit_complex("1234",4, "de");
is($es->apply_edits, "abc1234fghij", "uneven complex");
is($es->substring_ori(3,3), "cde", "substring_ori");
is($es->substring_edited(3,3), "c12", "substring_edited");
is($es->substring_edited(6,3), "fgh", "substring_edited 2");




#multiple edits
$s = "abcdefghij";
$es = EditSeq->new($s);
$es->edit_insert("12",3);
$es->edit_insert("34",6);
is($es->apply_edits, "abc12def34ghij", "two inserts");

$es = EditSeq->new($s);
$es->edit_substitute("12",5);
$es->edit_delete("defgh",4);
is($es->apply_edits, "abcij", "delete substitution");

$es = EditSeq->new($s);
$es->edit_substitute("12",5);
$es->edit_delete("bcd",2);
is($es->apply_edits, "a12ghij", "sub flanking del");

$es = EditSeq->new($s);
$es->edit_substitute("12",5);
$es->edit_delete("h",8);
is($es->apply_edits, "abcd12gij", "sub after del");

$es = EditSeq->new($s);
$es->edit_insert("12",5);
$es->edit_delete("h",8);
is($es->apply_edits, "abcde12fgij", "del after insert");

$es = EditSeq->new($s);
$es->edit_insert("12",5);
$es->edit_delete("f",6);
is($es->apply_edits, "abcde12ghij", "del flanking insert");

$es = EditSeq->new($s);
$es->edit_substitute("Q",6);
$es->edit_complex("123",3, "cd");
is($es->apply_edits, "ab123eQghij", "complex before substitution");

$es = EditSeq->new($s);
$es->edit_substitute("Q",6);
$es->edit_complex("123",7, "ghi");
is($es->apply_edits, "abcdeQ123j", "complex flankin substitution");

$es = EditSeq->new($s);
$es->edit_delete("fgh",6);
$es->edit_complex("123",3, "cd");
is($es->apply_edits, "ab123eij", "complex before delete");

$es = EditSeq->new($s);
$es->edit_delete("bc",2);
$es->edit_complex("123",6, "fg");
is($es->apply_edits, "ade123hij", "complex after delete");

#try an overlap
$es = EditSeq->new($s);
$es->edit_insert("12",5);
$es->edit_delete("defg",4); #overlap causes this to be discarded
is($es->apply_edits, "abcde12fghij", "overlap1");

$es = EditSeq->new($s);
$es->edit_delete("defg",4);
$es->edit_insert("12",5); #overlap causes this to be discarded
is($es->apply_edits, "abchij", "overlap2");

#test the position map
#abcdijklZZZZZmnoqqrstu
#abcdefghijklmnopqrstu
#fix to a-u
$es = EditSeq->new("abcdijklZZZZZmnoqqrstu");
$es->edit_insert("efgh",4);
$es->edit_delete("ZZZZZ",9);
$es->edit_substitute("p",17);
is($es->apply_edits, "abcdefghijklmnopqrstu", "triple edit");

is($es->convert_position_to_original(4),4);
is($es->convert_position_to_original(5),-1);
is($es->convert_position_to_original(11),7);
is($es->convert_position_to_original(21),22);
is($es->convert_position_to_original(6),-1);

is($es->convert_position_to_edit(4),4);
is($es->convert_position_to_edit(5),9);
is($es->convert_position_to_edit(9),12);
is($es->convert_position_to_edit(11),12);
is($es->convert_position_to_edit(13),12);
is($es->convert_position_to_edit(17),16);




