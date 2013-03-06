use Getopt::Long;
use Bio::GFF3::LowLevel qw/ gff3_format_feature /;
use PerlIO::gzip;
use IO::Handle;

use strict;
use warnings;

my ($primaryTable,$secondaryTable);
my $trackdb = "trackDb";
my $primaryName = 'name';
my $out = '.';
my @link;
my $verbose = '';
my $indir = ".";
my $subfeatures = '';

GetOptions('primaryTable=s' => \$primaryTable,
           'secondaryTable=s' => \$secondaryTable,
           'link=s{2}' => \@link,
           'primaryName=s' => \$primaryName,
           'out=s' => \$out,
           'verbose' => \$verbose,
	   'dir=s' => \$indir,
	   'getSubfeatures' => \$subfeatures);

my %typeMaps =
    (
        "genePred" =>
            ["txStart", "txEnd", "strand", "name", "score", "itemRgb"],
        "bed" =>
            ["chromStart", "chromEnd", "strand", "name", "score", "itemRgb"]
        );

$primaryTable or die "--primaryTable is required to run this script";

#only checks is the primary table is contained in trackDb since secondary table
#does not have to contain track data
my %trackColNames = name2column_map($indir . "/" . $trackdb);
my $tableCol = $trackColNames{tableName};
my $primRow = selectall($indir . "/" . $trackdb, sub {$_[0]->[$tableCol] eq $primaryTable});
if (! $primRow->[0]) {
    die "Table $primaryTable does not exist in $trackdb";
}

my $trackMeta = arrayref2hash($primRow->[0], \%trackColNames);
my @trackTypes = split(/ /,$trackMeta->{type});
my $type = $trackTypes[0];

$typeMaps{$type} or die "This script cannot handle $type tracks";

unless( -f "$indir/$primaryTable.sql" && -f "$indir/$primaryTable.txt.gz" ) {
    die "To format the $primaryTable track, you must have both files $indir/$primaryTable.sql and $indir/$primaryTable.txt.gz\n";
}

#retrieves all data from the primary table
my %primFields = name2column_map($indir . "/" . $primaryTable);
my $primData = selectall($primaryTable, sub {$_});

#retrieves the index number of the linking field in the secondary table
#and checks if sql and txt.gz files exist
my (%secFields,$matchIndex);
if ($secondaryTable) {
    unless( -f "$indir/$secondaryTable.sql" && -f "$indir/$secondaryTable.txt.gz" ) {
	die "To format the $secondaryTable track, you must have both files $indir/$secondaryTable.sql and $indir/$secondaryTable.txt.gz\n";
    }
    %secFields = name2column_map($indir . "/" . $secondaryTable);
    $matchIndex = $secFields{$link[1]} or die "--link is required if --secondaryTable is used";
}

#checks for left over fields to be written as attributes
my @gff3Required = ("chrom","source","type","score","strand","phase");
push (@gff3Required, ($typeMaps{$type}->[0],$typeMaps{$type}->[1]));
my @leftOver = grep {not $_ ~~ @gff3Required} keys %primFields;

open my $gff3, ">$out/$primaryTable.gff3" or die "Could not create $out/$primaryTable.gff3";

foreach my $row (@$primData) {
    my $rowData = arrayref2hash($row, \%primFields);
    my %gffAttr;
    #adds attributes to hash ref from the secondary table if it is specified
    #selects rows which have matching data
    if ($secondaryTable) {
        my $matchData = selectall($indir . "/" . $secondaryTable, sub{$_[0]->[$matchIndex] eq $rowData->{$link[0]}});
        foreach my $linkRow (@$matchData) {
            my $matchRow = arrayref2hash($linkRow, \%secFields);
            $gffAttr{$_} = $matchRow->{$_} foreach keys %$matchRow;
        }
    }
    #adds primary table data to gff hash ref
    my $gffHashRef = {
		    seq_id => $rowData->{chrom},
		    source => $rowData->{source} || ".",
		    type => $rowData->{type} || ".",
		    start => $rowData->{$typeMaps{$type}->[0]},
		    end => $rowData->{$typeMaps{$type}->[1]},
		    score => $rowData->{score},
		    strand => $rowData->{strand},
		    phase => $rowData->{phase} || ".",
    };
    
    #writes the left over fields as attributes
    $gffAttr{$_} = $rowData->{$_} foreach @leftOver;
    $gffHashRef->{attributes} = \%gffAttr;
    print $gff3 gff3_format_feature($gffHashRef);
    
    #retrieves subfeatures eg exons and cds
    if ($subfeatures) {
	my @getSub = ("cds","exon");
	foreach (@getSub) {
	    #some fields have start other have starts
	    #not the best way to go about doing this but it works
	    next unless $rowData->{$_ . "Starts"} || $rowData->{$_ . "Start"};
	    my @start = split /,/, $rowData->{$_ . "Starts"} || $rowData->{$_ . "Start"};
	    my @end = split /,/, $rowData->{$_ . "Ends"} || $rowData->{$_ . "End"};
	    my $subHashRef;
	    for (my $i=0;$i < scalar @start;$i++) {
		$subHashRef = {
		    seq_id => $rowData->{chrom},
		    source => $rowData->{source} || ".",
		    type => $_,
		    start => $start[$i],
		    end => $end[$i],
		    score => ".",
		    strand => $rowData->{strand},
		    phase => $rowData->{phase} || ".",
		    attributes => {
			Parent => $rowData->{name}
		    }
		};
	    print $gff3 gff3_format_feature($subHashRef);
	    }
	}
    }
}
close $gff3 or die "Could not close $out/$primaryTable.gff3";

#Reusing some of Robs subroutines to make parsing the sql and txt.gz
#files easier

# subroutine to crudely parse a .txt.gz table dump, and, for each row,
# apply a given subroutine to a array ref that holds the values for the
# columns of that row
sub for_columns {
    my ($table, $func) = @_;

    # my $gzip = new IO::Uncompress::Gunzip "$table.txt.gz"
    #     or die "gunzip failed: $GunzipError\n";
    my $gzip;
    open $gzip, "<:gzip", "$table.txt.gz"
        or die "failed to open $table.txt.gz: $!\n";

    my $lines = 0;
    my $row = "";
    while (<$gzip>) {
	chomp;
	if (/\\$/) {
            # unescape newline
	    chop;
	    $row .= "$_\n";
	} else {
	    $row .= $_;
            
	    my @data = split /(?<!\\)\t/, $row; # split on unescaped tabs
            map { s/\\\t/\t/g } @data; # unescape escaped tabs
	    &$func (\@data);
	    $row = "";
	}
	if (++$lines % 100000 == 0) { warn "(processed $lines lines)\n" }
    }
    $gzip->close()
        or die "couldn't close $table.txt.gz: $!\n";
}

# processes a table to find all the rows for which $filter returns true.
# returns a list of arrayrefs, where each arrayref represents a row.
sub selectall {
    my ($table, $filter) = @_;
    my @result;
    for_columns($table, sub { push @result, $_[0] if ($filter->($_[0])) });
    return \@result;
}

# subroutine to crudely parse a .sql table description file and return a map from column names to column indices
sub name2column_map {
    my ($table) = @_;
    my $sqlfile = "$table.sql";

    my @cols;
    local *SQL;
    local $_;
    open SQL, "<$sqlfile" or die "$! reading $sqlfile";
    while (<SQL>) { last if /CREATE TABLE/ }
    while (<SQL>) {
	last if /^\)/;
	if (/^\s*\`(\S+)\`/) { push @cols, $1 }
    }
    close SQL;

    return map (($cols[$_] => $_), 0..$#cols);
}

# converts an array ref of values and a hash ref with field name->index mappings
# into a hash of name->value mappings
sub arrayref2hash {
    my ($aref, $fields) = @_;
    my %result;
    foreach my $key (keys %$fields) {
        $result{$key} = $aref->[$fields->{$key}];
    }
    return \%result;
}