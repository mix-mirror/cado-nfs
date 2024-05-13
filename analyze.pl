#!/usr/bin/env perl -w

use strict;
use warnings;

my $producers = {};
my $stats={};

while(defined($_=<>)) {
    /^@(\d+\.\d+) 0x([\da-f]+) (\d) (\w+) (\d+) (\d+) (\d+)$/ && do {
        my $time = $1;
        my $id = $2;
        my $level = $3;
        my $name = $4;
        my $nbuckets = $5;
        my $npools = $6;
        my $N = $7;
        if ($N) {
            my $r = $producers->{$id} or die;
            my $dt  = $time - $r->{'time'};
            my $what = "$level$name updates to $npools*$nbuckets";
            print "$id dispatched $N $what in $dt\n";
            push @{$stats->{$what}}, $dt/$N;
        }
        $producers->{$id} = { 'time'=>$time };
    };
} 


for my $k (keys %$stats) {
    my @c = sort { $a <=> $b } @{$stats->{$k}};

    my $med = int($#c/2);
    print "$k $c[0] .. $c[$med] .. $c[$#c]\n";
}

