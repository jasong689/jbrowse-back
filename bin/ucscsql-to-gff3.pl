#!/usr/bin/perl
use strict;
use warnings;
use DBI;
use Getopt::Long;
use Bio::GFF3::LowLevel qw/ gff3_format_feature /;

#ucsc-to-json.pl
#queries ucsc public MySQL and creates a GFF3 file based on the options below
#options --primaryTable <table name> --secondaryTable <table name> --link <primaryTable column:secondaryTable column> --primaryName <column name>
#link is the condition by which to link the two tables

my $db = "hg19";
my $user = "genome";
my $host = "genome-mysql.cse.ucsc.edu";
my $password = "";

my ($primaryTable,$secondaryTable,$sqlQuery);
my $primaryName = 'name';
my $out = '';
my $link = '';
my $count = 0;
my $verbose = 0;

GetOptions('primaryTable=s' => \$primaryTable,
           'secondaryTable=s' => \$secondaryTable,
           'link=s' => \$link,
           'primaryName=s' => \$primaryName,
           'out=s' => \$out,
           'verbose' => \$verbose);

my @link = split /:/, $link;

#executes MySQL query based on entered options
my $dbh = DBI->connect("DBI:mysql:$db:$host", $user, $password);
if ($secondaryTable) {
    $sqlQuery = "select * from $primaryTable, $secondaryTable where $primaryTable.$link[0] = $secondaryTable.$link[1];";
} else {
    $sqlQuery = "select * from $primaryTable;";
}

my $sth= $dbh->prepare($sqlQuery) or die "Problem preparing: " . DBI->errstr;
my $validquery = $sth->execute or die "Problem executing: " . DBI->errstr;

#fetch all column names
my @fields = @{$sth->{NAME}};

#by default gff3 file is created within the working directory
open my $gff3, ">$out$primaryTable.gff3" or die "Could not create $primaryTable.gff3";

print $gff3 "##gff-version\t3\n";

#fetch rows from database and first writing columns 1-8
#then writes column 9 with leftover fields
while (my $databases = $sth->fetchrow_hashref()){
    my @gff3Required = ("chrom","source","type","score","strand","phase");
    my $hashRef = {
        seq_id => $$databases{chrom},
        source => $$databases{source},
        type => $$databases{type},
        start => $$databases{chromStart} || $$databases{txStart},
        end => $$databases{chromEnd} || $$databases{txEnd},
        score => $$databases{score},
        strand => $$databases{strand},
        phase => $$databases{phase},
    };
    if (exists $$databases{chromStart}) {
        push(@gff3Required,("chromStart","chromEnd"));
    } else {
        push(@gff3Required,("txStart","txEnd"));
    }
    my @leftOver = grep {not $_ ~~ @gff3Required} @fields;
    my %uniqueLeftOver;
    
    #removes duplicate fields returned by SQL query
    $uniqueLeftOver{$_} = $_ foreach @leftOver;
    
    #removes any subfeatures from leftover attributes
    my @subfeatures = grep /.+Start|.+End/, keys %uniqueLeftOver;
    delete $uniqueLeftOver{$_} foreach @subfeatures;
    
    #creates an attributes hashref for attributes and adds to the features hashref
    my %attributes;
    $attributes{$_} = [($$databases{$_})] foreach keys %uniqueLeftOver;
    $attributes{ID} = $$databases{name};
    #switches default name with primaryName
    ($attributes{name},$attributes{$primaryName}) = ($attributes{$primaryName},$attributes{name});

    $hashRef->{attributes} = \%attributes;
    
    print $gff3 gff3_format_feature($hashRef);
    
    #creates the subfeatures
    #sorts subfeatures so that End is always one element before Start
    my @sortedSub = sort @subfeatures;
    for (my $x = 0;$x < scalar @sortedSub;$x++) {
        if ($sortedSub[$x] =~ /(.+)Start/) {
            my @start = split /,/, $$databases{$sortedSub[$x]};
            my @end = split /,/, $$databases{$sortedSub[($x - 1)]};
            for (my $i=0;$i < scalar @start;$i++) {
                my $hashRef = {
                    seq_id => $$databases{chrom},
                    source => $$databases{source},
                    type => $1,
                    start => $start[$i],
                    end => $end[$i],
                    score => undef,
                    strand => $$databases{strand},
                    phase => undef,
                    attributes => {
                        Parent => $$databases{name}
                    }
                };
                print $gff3 gff3_format_feature($hashRef);
            }
        }
    }
    print "Working on row $count\n" if (($count % 1000) == 0) && $verbose;
    $count++;
}


$sth->finish() or die "Problem finishing: " . DBI->errstr;

$dbh->disconnect() or die "Problem disconnecting: " . DBI->errstr;