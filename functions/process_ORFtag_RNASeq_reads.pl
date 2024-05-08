#!/usr/bin/perl
use strict;
use warnings;

while (<>) {
    chomp;
    if (/^@(\S+)/) {
        my $read_name = $1;
        my $sequence  = <>;
        my $qual_name = <>;
        my $qual_line = <>;

        my $trimmed_sequence;
        my $trimmed_quality;

        if ($sequence =~ /^[NACGT]{6}CCACGACGGAGACTACAAGGATCATGATATTGATTACAAAGACGATGACGATAAGCAG(.+)$/) {
            $read_name .= "_f0";
            $trimmed_sequence = $1;
            $trimmed_quality  = substr($qual_line, length($sequence) - length($1));
        } elsif ($sequence =~ /^[NACGT]{6}CCACGACGGAGACTACAAGGATCATGATATTGATTACAAAGACGATGACGATAAGGCAG(.+)$/) {
            $read_name .= "_f1";
            $trimmed_sequence = $1;
            $trimmed_quality  = substr($qual_line, length($sequence) - length($1));
        } elsif ($sequence =~ /^[NACGT]{6}CCACGACGGAGACTACAAGGATCATGATATTGATTACAAAGACGATGACGATAAGGCCAG(.+)$/) {
            $read_name .= "_f2";
            $trimmed_sequence = $1;
            $trimmed_quality  = substr($qual_line, length($sequence) - length($1));
        } else {
            $read_name .= "_ambiguous";
            $trimmed_sequence = $sequence;
            $trimmed_quality  = $qual_line;
            chomp $trimmed_sequence;  # Remove the newline character from the original sequence
        }

        print "\@$read_name\n$trimmed_sequence\n$qual_name$trimmed_quality";
    }
}
