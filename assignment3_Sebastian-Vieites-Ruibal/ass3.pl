use strict;
use warnings;
use Bio::Perl;
use Bio::SeqFeature::Generic;
use Bio::Seq;
use Bio::SeqIO;
use LWP::Simple;

# This is the third assignment script, where we learn how to use BioPerl.
# Sebastian Vieites Ruibal.
# I just wanted to leave constance that i thought about putting all subroutines inside a package and call them from there,
# but realised that would be easier to correct the code if everything is inside the same file, instead of having to go
# from one to the other, so thats the reason why all subroutines are written before the program core.

########################### SUBROUTINES ###############################

sub openfile{
        chomp (my$filename=<STDIN>);                            # When written the name, you chomp it to delete de \n.
        if (-e $filename){                                      # If the file exists
                open (my$filehandle, "<", "$filename");         # we open it
                my@file_array=<$filehandle>                     # and store it inside an array.
        }
        else{
                die "Wrong filename or $filename doesnt exist.\nThe program will stop execution.\n\n";   # If the file doesnt exist, returns an error message and program ends.
        }
}


sub get_info{
	print "Retrieving gene information from server, please wait...\t";
        my@id_array=@_; 	                 			# The array with our TAIR IDs is an entrance value
	my$id_list=join(",", @id_array);				# that we will use to only have to make ONE accesion to the web, by joining all the IDs.
	my%gene_info_hash;						# We initialize the hash where we will save the information for our objects.
	my$get_URL="http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=ensemblgenomesgene&format=embl&id=".$id_list;
	my$seqIO;
	if (my$info=get($get_URL)){					# If we can retrieve information from the server
		$seqIO=Bio::SeqIO->new(	-string => $info,		# we create objects containing that information
					-format => 'EMBL');
		while (my$current_info=$seqIO->next_seq){		# for all the info downloaded
			my$id=shift(@id_array);				# and all the IDs we have
			$current_info->primary_id($id);			# we get the primary ID
			$gene_info_hash{$id}=$current_info;		# and store that information in a hash.
		}		
	}
	else {
		die "Unable to retrieve information from server.\nThe program will stop execution.\n\n";	# If for any reason the conection fails, returns an error message and program ends.
	}
	print "Information successfully retrieved.\n";
	return \%gene_info_hash;
}


sub pattern_search{
        my$ref=shift;														
        my%info_hash=%$ref;                                                                                                     # The hash where we have the info about our IDs.              
        my@hit_array;                                                                                                           # The array where we will store which IDs have match.     
        print "Searching for requested pattern, please wait...\t";
        foreach my$gene(keys%info_hash){											# For each gene in our information hash
                my@info=$info_hash{$gene}->get_all_SeqFeatures;									# we get all the features
		my$full_sequence=$info_hash{$gene}->seq;									# and store the full sequence
		my$full_sequence_length=length($full_sequence);									# and its length
		my$pattern="CTTCTT";
		my$reverse_pattern="AAGAAG";
		while ($full_sequence =~ /$pattern/g){										# so we can search the pattern.
			my$start=$-[0]+1;											# We have to adjust the start so we really have 6 position for nucleotids.
			my$strand_type=1;											# For the normal pattern, we will have normal strand type.
			foreach my$feature(@info){										# For each feature that we have
				if (($feature->primary_tag eq 'exon')&&($feature->start <= $start)&&($feature->end >= $+[0])){	# if its an exon and the match is inside it
                                        my$pattern_match=Bio::SeqFeature::Generic->new( -seq_id => $gene,			# we make a new feature saving the match
                	                        -source_tag => "Sebastian_Vieites_Ruibal_pattern_match_search_script",
        	                                -primary_tag => "PATTERN_MATCH",
     	            		                -start => $start,
                                	        -end => $+[0],
                        	                -score => ".",
                	                        -strand => $strand_type,
        	                                -phase => ".",
	                                        -display_name => "Match of $pattern inside exon");
                                        $info_hash{$gene}->add_SeqFeature($pattern_match);					# and store it in the hash that contains our gene information.
                                        push(@hit_array, $gene);								# Last, we save the id that had the match inside an array.
				}
			}
		}
                while ($full_sequence =~ /$reverse_pattern/g){									# We repeat the process, but with the reverse pattern
                        my$end=$+[0]-1;												# the adjust will be in the end (because its a reverse chain)
                        my$strand_type=-1;											# and the strand type will be reverse.
                        foreach my$feature(@info){
                                if (($feature->primary_tag eq 'exon')&&($feature->start <= $-[0])&&($feature->end >= $end)){
                                        my$pattern_match=Bio::SeqFeature::Generic->new( -seq_id => $gene,
                                                -source_tag => "Sebastian_Vieites_Ruibal_pattern_match_search_script",
                                                -primary_tag => "PATTERN_MATCH",
                                                -start => $-[0],
                                                -end => $end,
                                                -score => ".",
                                                -strand => $strand_type,
                                                -phase => ".",
                                                -display_name => "Match of $reverse_pattern inside exon");
                                        $info_hash{$gene}->add_SeqFeature($pattern_match);
                                        push(@hit_array, $gene);		
				}
			}
		}
	}
	print "Search finished.\n";
	return (\%info_hash, @hit_array);
}


sub uniq{                       # Subroutine to delete repeated elements on an array.
        my%seen;                # All credits go to Greg Hewgill, as answered in http://stackoverflow.com/questions/7651/how-do-i-remove-duplicate-items-from-an-array-in-perl
        grep !$seen{$_}++, @_;      
}


sub get_report{
	my($ref1,$ref2)=@_;				# We introduce both arrays as references
	my@repeat_hit_array=@$ref1;			# but the one with hits has the hits repeated
	my@hit_array=uniq(@repeat_hit_array);		# so we use a subroutine to delete repeated elements.
	my@full_list=@$ref2;
	my%hit_hash=map{$_=>1} @hit_array;		# We compare both arrays to see what elements are in full array that arent in match array (which means that they have no pattern match)
	my@no_hit=grep { !$hit_hash{$_} } @full_list;	# and store it on a new array. Info retrieved from http://stackoverflow.com/questions/2933347/comparing-two-arrays-using-perl.
	my$outfile='no_match.txt';			# We get the name of the no-match file
	open (FILE, "> $outfile");			# and prepare to write on it.
	foreach my$miss(@no_hit){			# Finally we print a report of the identifiers that have no coincidences
		print "No match coincidences for $miss.\n";
		print FILE "$miss\n";			# and print it in a file.
	}
	close FILE;
}


sub create_GFF3{
	my$ref=shift;
	my%info_hash=%$ref;
	my$outfile_local='GFF_table_local_positions.gff';								# We introduce the name for the file where we will store the local coordinates
	my$outfile_global='GFF_table_global_positions.gff';								# and the global coordinates.
	open (LOCALFILE, "> $outfile_local");										# Prepare opening both files.
	open (GLOBALFILE, "> $outfile_global");
	foreach my$gene(keys%info_hash){                                                               			# For each identifier in our hash
                my@info=$info_hash{$gene}->get_all_SeqFeatures;                                         		# we get all the features that contains
                my$global_info=$info_hash{$gene}->accession_number;                                             	# and retrieve some information using the ID
                my($first, $second, $chromosome, $sequence_start, $sequence_end, $last)=split ":", $global_info;	# like the chromosome where it is, the global start of the sequence and the global end$
                while (@info){												# and while we have accessions
			my$feature=shift(@info);									# we iterate on each one.
			if ($feature->primary_tag eq "PATTERN_MATCH"){							# If the primary tag of our feature is a pattern match
				my$GFF_strand;
				if ($feature->strand eq "1"){								# we transform the format of the strand to GFF
					$GFF_strand="+";
				}
				elsif ($feature->strand eq "-1"){
					$GFF_strand="-";
				}
				else {
					$GFF_strand=".";
				}
				# and write the information inside the local file (NEXT LINE)
				print LOCALFILE $feature->seq_id."\t".$feature->source_tag."\t".$feature->primary_tag."\t".$feature->start."\t".$feature->end."\t".$feature->score."\t".$GFF_strand."\t".$feature->phase."\t".$feature->display_name."\n";
				my$local_start=$feature->start;								# For the global file, we need to make some conversions of the coordinates
				my$global_start=$local_start+$sequence_start;
				my$local_end=$feature->end;
				my$global_end=$local_end+$sequence_start;
				# before we can store them inside the global file (NEXT LINE).
				print GLOBALFILE $chromosome."\t".$feature->source_tag."\t".$feature->primary_tag."\t".$global_start."\t".$global_end."\t".$feature->score."\t".$GFF_strand."\t".$feature->phase."\t".$feature->display_name."\n";
			}
		}
	}
	close LOCALFILE;	# We close both files.
	close GLOBALFILE;
}


########################### PROGRAM CORE ##############################

# 1. Get the identificators inside an array.
print "We are going to start the execution of the third assignment program.\n";
print "Please insert the filename where your identificator list is: ";
my@id_list=&openfile;								# Opening the gene list that we have as an array
chomp@id_list;                                                                	# and chomp each element of the array.
print "File succesfully opened.\n";

# 2. Retrieving information for each identificator we have.
my$genes_information=&get_info(@id_list);					# Retrieving information from the dbfecth server.

# 3. Looking for the pattern asked in the assignment and creating the features in case of coincidence.
my@hit_array;
($genes_information, @hit_array)=&pattern_search($genes_information);		# We update the gene information hash and get a match coincidence array.

# 4. Output the report of the identifiers that have no match.
&get_report(\@hit_array,\@id_list);						# Using the array containing the hits and the identifier list, we can get which identifiers have no match.

# 5. Creating the GFF3 format files.
&create_GFF3($genes_information);						# Using the information that we have for each gene, we make the files asked for in the assignment.
