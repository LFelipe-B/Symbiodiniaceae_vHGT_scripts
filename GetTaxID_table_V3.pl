#!/usr/bin/perl

### USAGE: perl GetTaxID_table_V2.pl blastout_file [database] (the sequence db (nucleotide or protein) the GI is for)
### Requires perl-XML-Twig
### Requires an internet connection to access the taxonomy information
### Written by Fabien Burki
# -------------------------------
# Customization by L. Felipe Benites 2018-2020 (V2 - 23.04.2020) to include bitscore insted of hsp->freq#, GI/GB, and identity.  print in CSV
# WARNING: Only works in diamond outputs given the GI like this gi|1129193959|gb|OLP98445.1| - need to mod. to use on blast outfmt #
# more info on SearchIO tool https://bioperl.org/howtos/SearchIO_HOWTO.html #


use Bio::SearchIO;
use Bio::DB::Taxonomy;

my $blast_report = new Bio::SearchIO('-format' => 'blasttable', # modify this line according to your blast output format (blasttable if tabular; blast if default, etc)
				     				 '-file' => $ARGV[0]);
				     				 
my $db = Bio::DB::Taxonomy->new (-source => 'entrez', -verbose => $verbose);

open (OUT, ">$ARGV[0]_taxIDs.csv");

#WARNING: set how many hits and hsps you want !

my $nbhits=50;
my $nbhsps=50;

print OUT "Query,Gene_Identifier,Scientific_name,Phylum,Kingdom,Superkingdom,Evalue,Identity,Bitscore\n";

while (my $result = $blast_report->next_result) {
	#output "no hits found" if that's the case
	if ($result -> num_hits == 0) {
		print OUT $result->query_name(), ",";
		print $result->query_name(), ",";
		print OUT "No hits\n";
		print "No hits\n";
   	}    	
   	else {
    	my $counter=0;
    	while (my $hit = $result->next_hit()) {
			if ($counter < $nbhits) {
				while (my $hsp = $hit->next_hsp()) {
					if ($counter < $nbhsps) {
						print OUT $result->query_name(), ",";
                        print $result->query_name(), ","; #text in terminal
						my $hitname = $hit->name();
#						print "$hitname,";
#                       print OUT $hit->accession(), ","; #uncomment to print
						if ($hitname =~ /gi\|(\d+)\|/) {
							my @array;
							my $gi = $1;
#                           print OUT $hit->description, ","; #uncomment to print
#							print "$gi,"; #uncomment to print on terminal
							my $node;
							eval {$node = $db->get_Taxonomy_Node(-gi => $gi, -db => $dbname)};														
							if ($@) {
                                print "$gi Unknown taxid\n";
                                print OUT "$gi Unknown taxid\n";
 	                			next;
    	     				}
							else {
                                print OUT $hit->name(),","; #GI AND GB
                                print $hit->name(),"\n"; #PRINT GI AND GB IN FILE OUT
								print OUT $node->scientific_name;
								print $node->scientific_name, ",";
								my $parent = $db->get_Taxonomy_Node($node->parent_id);
#								print $parent->node_name, ",";
								while (defined $parent && $parent->node_name ne 'root') { 									
									if ($parent->rank eq 'superkingdom') {
										push @array, $parent->node_name;
									}
									elsif ($parent->rank eq 'kingdom') {
										push  @array, $parent->node_name;
									}
									elsif ($parent->rank eq 'phylum') {
										push @array, $parent->node_name;
									}							
									$parent = $db->get_Taxonomy_Node($parent->parent_id);
								}
							}
							my $size = @array;
							if ($size == 1) {
								print OUT ",,,$array[0],";
#								print ",,,$array[0],";
							}
							elsif ($size == 2) {
								print OUT ",,$array[0],$array[1],";
#								print ",,$array[0],$array[1],";
							}
							elsif ($size == 3) {
								print OUT ",$array[0],$array[1],$array[2],";
#								print ",$array[0],$array[1],$array[2],";
							}
                            print OUT $hsp->evalue(),",";
                            print $hsp->evalue(),"\n";
                            print OUT $hsp->percent_identity(),",";
                            print $hsp->percent_identity(),"\n";
                            my $identities = $hsp->bits;
                            my $identities_formatted = sprintf ("%.2f", $identities);
                            print OUT "$identities_formatted\n";
		    			$counter+=1;
		    			}	
					}
				}
	    	}	
		}
	}
 }

print "\n";
exit;
