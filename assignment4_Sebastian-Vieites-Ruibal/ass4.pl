use strict;
use warnings;
use Bio::SeqIO;
use Bio::Tools::Run::StandAloneBlastPlus;

# This is the fourth assignment script, where we learn how to find orthologues genes using a perl script to blast.
# In the moment of uploading the script there is no report file uploaded because of the long time it takes por the PC to work on it, but it will be up as soon as possible.
################################# PROGRAM CORE ###############################

#1. Creating de databases using system commands.
chdir './data';												# We change the work directory where our information already updates (see fasta_converter.pl) is.
my@is_already_done=`ls`;										# We get the name of all the files in the folder to check if we already have the files needed
my%is_already_done_hash= map {$_ => 1 } @is_already_done;						# to run the blast (see the conditional in the next line).
unless ((exists($is_already_done_hash{"AT_DB.nhr"}))&&(exists($is_already_done_hash{"AT_DB.nin"}))&&(exists($is_already_done_hash{"AT_DB.nsq"}))&&(exists($is_already_done_hash{"SP_DB.phr"}))&&(exists($is_already_done_hash{"SP_DB.pin"}))&&(exists($is_already_done_hash{"SP_DB.psq"}))){
	`makeblastdb -in athal.fas -dbtype 'nucl' -title "Co-Expressed Genes" -out "AT_DB" -input_type fasta`;	# If we dont have the files we create the, 
	`makeblastdb -in spombe.fas -dbtype 'prot' -title "Co-Expressed Genes" -out "SP_DB" -input_type fasta`;
}
chdir '..';												# and return to the main directory.

#2. Preparing the factories and the querys to run the blast.
# The paper where we got the '1e-6' can be found in http://bioinformatics.oxfordjournals.org/content/24/3/319.long at 17/04/2016.
my$athal_factory=Bio::Tools::Run::StandAloneBlastPlus->new(						# We prepare the factory for the Arabidopsis genes
	-db_name => 'AT_DB',
	-db_dir => './data',
	-expect => '1e-6',);
my$athal_query=Bio::SeqIO->new(										# and the query we will use against the S.pombe.
	-file => './data/athal.fas',
	-format => 'fasta',);
my$spombe_factory=Bio::Tools::Run::StandAloneBlastPlus->new(						# Repeat the process for S.pombe.
	-db_name => 'SP_DB',
	-db_dir => './data',
	-expect => '1e-6',);
my$spombe_query=Bio::SeqIO->new(
	-file => './data/spombe.fas',
	-format => 'fasta',);

#3. Core of the blast program.
my%athal_hits;								# We initialize the hashes where we will store the positive matches
my%spombe_hits;
while (my$seq=$athal_query->next_seq){					# and first we do the blast using the Arabidopsis as query and the S.pombe as target.
	my$athalhit=$spombe_factory->blastx( -query => $seq);
	if (my$hit=$athalhit->next_hit){
		$athal_hits{$seq->primary_id}=($hit->name);
		$spombe_hits{$hit->name}="HIT";
	}
	$spombe_factory->cleanup;					# We use the cleanup to remove all the temp files created by the blast.
}
while (my$seq=$spombe_query->next_seq){
	if (exists $spombe_hits{$seq->primary_id}){			# We repeat the process on reverse with the possitive matches.
		my$spombehit=$athal_factory->tblastn( -query => $seq);
		if (my$hit=$spombehit->next_hit){
			$spombe_hits{$seq->primary_id}=($hit->name);
		}
	}
	$athal_factory->cleanup;
}

#4. Writing the positive orthologue results to a file.
my$outfile='orthology_results.tsv';					# We indicate the name for the outfile
open (OUTFILE, "> $outfile");						# open it
print OUTFILE "And this is the final report of the linked genes:\nArabidopsis thaliana genes\tSchizosaccharomyces pombe genes\n";
foreach my$gene(keys%spombe_hits){					# and fill it with the cases where are reciprocal hits
	if ($gene eq $athal_hits{$spombe_hits{$gene}}){
		print OUTFILE "$spombe_hits{$gene}\t$gene\n";
	}
}
close OUTFILE;								# before closing it.
