use strict;
use gene;
use stock;
use POSIX qw(strftime); # To take the date for the subroutine &plant_simulator.

############# Subroutines #################
sub openfile{
	chomp (my$filename=<STDIN>); 				# When written the name, you chomp it to delete de \n.
	if (-e $filename){					# If the file exists
		open (my$filehandle, "<", "$filename");		# we open it
		my@file_array;
		while (my$line=<$filehandle>){				
			chomp $line;				# we chomp the line to remove the \n
			my@line_array=split(/\t/, $line);	# we split the line in an array
			push (@file_array, \@line_array);	# we push that line splitted into the definitive array
		}	
		shift @file_array;				# we remove the first line, that is only the column information
		return @file_array;				# and last we return the array where the file info is stored.
	}
	else{
		die "Wrong filename or $filename doesnt exist.\nThe program will stop execution.\n\n" 	# If the file doesnt exist, returns an error message and program ends.
	}
}	

sub load_gene_data{
	my@own=@_;					# We put the array previously made inside a variable to work with	
	my$size=scalar@own;				# we take the measure of the array size
	my$count=0;					# we initialize a counter
	my%gene_hash;					# and we declare the existance of hour hash.
	while ($count<$size){				# Now we make the condition that our counter is smaller of the array size (to not get errors),
		if ($own[$count][0] =~ m/A[Tt]\d[Gg]\d\d\d\d\d/){ # that the ID requires a certain format,
			my$new_gene=gene->new(		# and we create our instances, storing them on the variable $new_gene.
				id => $own[$count][0],
				name => $own[$count][1],
				phenotype => $own[$count][2],
			);
			$gene_hash{$own[$count][0]}=$new_gene;	# We make that the key for the instance inside our hash is the own id of the instance, and the value the instance that we have created.
			$count ++;				# We make the counter increment in 1
		}
		else{
			die "Some id in your gene data file doesnt match the format required /A[Tt]\d[Gg]\d\d\d\d\d/ so please fix it. Program execution will stop.\n";
		}
	}
	print "Finished creating gene instances\n";
	return \%gene_hash;					# and we return the hash reference. If the expression didnt match, it would have stopped the program with an error.
}

sub load_stock_data{
	my($ref1,$ref2)=@_;						# We have to introduce the array of references in two variables
	my@own=@$ref1;							# so we can make the first reference an array inside the subroutine
	my%gene_hash=%$ref2;						# and the second reference a the hash containing the gene information obtained in the step before.
	my$size=scalar@own;						# We measure the size of the array
	my$count=0;							# start the counter in zero
	my%stock_hash;							# and initialize the hash where our data will be saved.
	while ($count<$size){						# We put the condition that the counter is smaller than the array size to dont get errors
		my$new_stock=stock->new(				# and we start creating our instances.
			seed => $own[$count][0],
			mutant_id => $gene_hash{$own[$count][1]},	# In this case, one of the properties of our instance is other instance.
			last_plantation => $own[$count][2],
			storage => $own[$count][3],
			grams => $own[$count][4],
		);
		$stock_hash{$own[$count][0]}=$new_stock;		# We make the stock ID the key from the hash
		$count ++;						# and increment the counter by 1.
	}
	print "Finished creating stock instances\n";
	return \%stock_hash;						# Finally, we return the reference to our new created hash with our stock info inside.
}

sub plant_simulator{
	my$ref=shift;		# We save the hash referene inside a variable
	my%seed_hash=%$ref;	# and using it we dump the hash inside our subroutine.
	print "\nNow we are gonna simulate planting seeds directly from our stock.\nPlease input a positive integer number as the grams of seeds that you would like to plant: ";
        chomp (my$number=<STDIN>);									# We use an STDIN to ask the user how many grams he wants to plant
	if ($number =~ m/^[+]?\d+$/){									# and we check that is in fact an integer which value is at least 0.
		print "We will plant $number grams.\n";
		my$outfile='updated_seed_stock_data.tsv';						# We introduce the name of the file where we will save the updated data via code,
		open (FILE, "> $outfile") ;								# we open that file
		print FILE "Seed_Stock\tMutant_Gene_ID\tLast_planted\tStorage\tGrams_Remaining\n";	# and we write the first line, which we erased when we stored the info in our array.
		foreach my$key (keys%seed_hash){							# Now, we use a foreach to move through the objects stored in our hash as values
			my$ammount=$seed_hash{$key}->grams;						# we store the grams from the storage inside a variable
			if ($ammount<=$number){								# which we will use to check if the grams asked are bigger than the ones we have in our stock.
				$seed_hash{$key}->grams(0);						# If the ammount asked is bigger than remaining, we store that the new value for grams in that object is 0.
				my$finished=$seed_hash{$key}->seed;					
				print "There is not enough $finished to plant, we could only plant $ammount grams and now the remaining stock is 0.\n"; # And we put a friendly message saying it to the user.
			}
			else{										# If the case is the opposite
				$seed_hash{$key}->grams($ammount-$number);				# we make the substraction and store it as the new value of the property for the object.
			}
			my$new_seed=$seed_hash{$key}->seed;						# We store all the properties' value from each object in variables
			my$new_mutant=$seed_hash{$key}->mutant_id->id;
			my$new_plantation=strftime "%d/%m/%Y", localtime;
			my$new_storage=$seed_hash{$key}->storage;
			my$new_grams=$seed_hash{$key}->grams;
			print FILE "$new_seed\t$new_mutant\t$new_plantation\t$new_storage\t$new_grams\n"; # so we can write them in the file where the info is updated.
		}
	}
	else{
		die "\nThe number isnt a positive integer. Program execution will be stopped.\n";	# This is a check, if the number which we ask for in the begginning isnt an integer which value is at least 0, the program stops and output an error message.
	}
	print "The new state of the stock will be saved in 'updated_seed_stock_data.tsv'\n\n";		# This is just info for the user about where the new info is stored.
}

sub ji_test{
	my$ref1=shift;										# We save the reference to our hash inside a string,
	my%stock_hash=%$ref1;									# and we dereference it in order to use the data.
	my@crossing_data=@_;									# The crossing data is stored inside an array.
        my$size=scalar@crossing_data;								# We measure the size of the array
        my$count=0;										# and initialize a counter variable
	while ($count<$size){									# in order to keep reading the array but without getting errors.
		my$parent1=$crossing_data[$count][0];						# We store the info from the first column of the array (id in the stock of the first parent),
		my$parent2=$crossing_data[$count][1];						# from the second (id in the stock of the second parent),
		my$wildtype_observed=$crossing_data[$count][2];					# the third (number of wildtype descendants),
		my$p1_observed=$crossing_data[$count][3];					# the fourh (number of recessive for a phenotype descendants),
		my$p2_observed=$crossing_data[$count][4];					# the fifth (number of recessive for b phenotype descendants),
		my$p1p2_observed=$crossing_data[$count][5];					# and the sixth (number of double recessive descendants).
		my$total_observed=$wildtype_observed+$p1_observed+$p2_observed+$p1p2_observed;	# We will need the total number of descendants for the test
		my$wildtype_hoped=$total_observed*(9/16);					# so we can calculate the hoped wildtype,
		my$single_recessive_hoped=$total_observed*(3/16);				# the hoped single recessive descendants,
		my$double_recessive_hoped=$total_observed*(1/16);				# and the double recessive_descendants.
		#$test_number holds the value obtained from ji-square test, which we will analyze.
		my$test_number=((($wildtype_observed-$wildtype_hoped)**2)/$wildtype_hoped)+((($p1_observed-$single_recessive_hoped)**2)/$single_recessive_hoped)+((($p2_observed-$single_recessive_hoped)**2)/$single_recessive_hoped)+((($p1p2_observed-$double_recessive_hoped)**2)/$double_recessive_hoped);
		my@linked_to_first;								# We definite two arrays,
		my@linked_to_second;
		if (defined $stock_hash{$parent1}->mutant_id->linked){				# which in case of existing previous info of linked genes
			@linked_to_first=@{$stock_hash{$parent1}->mutant_id->linked};		# will contain that info.
		}
		if (defined $stock_hash{$parent2}->mutant_id->linked){
			@linked_to_second=@{$stock_hash{$parent2}->mutant_id->linked};
		}
		if ($test_number>7.8){								# We use the value 7.8, equivalent to 0'05 significance and 3 liberty-degrees.
			my$first_linked=$stock_hash{$parent1}->mutant_id;			# We take the info of the gene instance of the first parent
			my$second_linked=$stock_hash{$parent2}->mutant_id;			# and we do the same with the second.
			my$first_linked_name=$first_linked->name;				# We save the name of the first parent
			my$second_linked_name=$second_linked->name;				# and from the second one
			print "$first_linked_name and $second_linked_name are ligated with a ji-square value of $test_number\n"; # in order to write a first output for the user.		
			push (@linked_to_first, $second_linked);				# We store the info from the linked genes on the linked array
			push (@linked_to_second, $first_linked);				# for both instances
			my$linked_to_first_reference=\@linked_to_first;				# and we dereference the array 
			$stock_hash{$parent1}->mutant_id->linked($linked_to_first_reference);	# so we can save it as the property 'linked'
			my$linked_to_second_reference=\@linked_to_second;
			$stock_hash{$parent2}->mutant_id->linked($linked_to_second_reference);	# for both instances.
		}		
		$count ++;									# Increase the counter!
	}
}

############## Program core ###############

#1. Getting the files' data in arrays and creating the objects, storing them on hashes and creating the reference to them.
print "\nWe are going to start the execution of the first assignment program.\n";
print "Please insert the filename where your gene data is: ";
my@gene_file=&openfile; # Writting the info of your gene on an array.
my$gene_data=&load_gene_data(@gene_file); # Storing the hash reference where gene data is inside the variable $gene_data.

print "Now please insert the filename where your stock data is: ";
my@stock_file=&openfile; # Writting the info of your stock on an array.
my$stock_data=&load_stock_data(\@stock_file,$gene_data); # Storing the hash reference where stock data is inside the variable $stock_data.

print "Now please insert the filename where your crossing data is: ";
my@cross_data=&openfile; # Writing the info of your crossings on an array.

#2. Making the simulation about planting seeds.
&plant_simulator($stock_data); # We simulate planting seeds taken from the stock that we have.

#3. Computing the search for linked genes and updating the data.
print "Now we will compute the crossing data to search for ligated genes...\n";
&ji_test($stock_data,@cross_data); # We make a ji-square test in order to search for linked genes in our data.

print "\nFinal report of linked genes:\n";						# The last report of linked genes in our data
foreach my$genes(keys%{$gene_data}){							# starts parsing all our 'gene' instances
	if ($gene_data->{$genes}->is_linked){						# and if the search finds some gene with info of linkage
		my$first_gene_name=$gene_data->{$genes}->name;				# gets the name of the gene info which we are parsing
		my$linked_genes=$gene_data->{$genes}->linked;				# and gets the reference of its linked genes stored in the array,
		foreach my$linked(@{$linked_genes}){					# which after dereference and parse through
			my$second_gene_name=$linked->name;				# gets us the name of each linked gene to it
			print "$first_gene_name is linked to $second_gene_name.\n";	# and we output on the screen for the user.
		}
	}
}
