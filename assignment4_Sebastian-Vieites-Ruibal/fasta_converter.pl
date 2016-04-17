use strict;
use warnings;
use Bio::Seq;
use Bio::SeqIO;
# This is a short script to rename the fasta headers because some of them gave problems when running the blast.
# Usage: $fasta_converter.pl infile outfile

################################################################
open (FILEHANDLE, $ARGV[0]);						# We put as argument the name and location of the entrance file
my@updated_list;							# and initialize the array where we will store the sequences with the updated headers.
my$entrance_seqs=Bio::SeqIO->new (					# Here we have the unupdated sequences
	-fh => \*FILEHANDLE,
	-format => 'fasta',);
while (my$seq=$entrance_seqs->next_seq){				# and while we still have sequences without update we will keep compiling names.
	my$seq_name=$seq->primary_id;					
	my$final_name=$seq_name;
	my$what_type=substr($seq_name, 0, 2);
	if ($what_type eq 'SP'){					# This conditional is because to only have the gene ID without information in the case of S.pombe genes we have to use a
		($final_name) = $seq_name =~ /(SP\w*\.\d{1,3}\w?)\|/;	# regular expresion to get it-
	}
	my$final_seq=Bio::Seq->new(					# Once we have the gene names we can update de data
		'-id' => $final_name,
		'-seq' => $seq->seq,);
	push(@updated_list, $final_seq);				# and store it inside the array.
}
my$leaving_seqs=Bio::SeqIO->new (					# Because we have as second argument the filename of the outfile
	-file => '>'.$ARGV[1],
	-format => 'fasta',);
while (@updated_list){
	$leaving_seqs->write_seq(shift(@updated_list));			# we can store there the sequences in fasta format with the updated headers.
}
