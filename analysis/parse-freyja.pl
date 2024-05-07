
open(FILE, $ARGV[0]) or die "cannot open file\n";
print join("\t","file","variant","percentage")."\n";
my $line=<FILE>;
while( $line=<FILE>){
	chomp $line;
	my ($id,$tags)=split("\t",$line);
	#print Dumper($tags);
	$tags =~ s/^\[|\]$//g;
	my @tuples = split(/\), \(/, $tags);
# Process each tuple
foreach my $tuple (@tuples) {
    # Clean up the tuple string
    $tuple =~ s/^\(|\)$//g; # Remove leading and trailing parentheses
    $tuple =~ s/'//g; # Remove single quotes
    # Split the cleaned tuple into label and value
    my ($label, $value) = split(/, /, $tuple, 2);
    print join("\t",$id,$label,$value)."\n";
}
 
}
