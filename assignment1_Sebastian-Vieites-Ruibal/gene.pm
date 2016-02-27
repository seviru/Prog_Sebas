package gene;
use Moose;

has 'id' => ( # Id of the mutant gene
	is => 'rw',
	isa => 'Str',
	required => 1,
);

has 'name' => ( # Name of the mutant gene
	is => 'rw',
	isa => 'Str',
	required => 1,
);

has 'phenotype' => ( # Phenotype of the mutant gene
	is => 'rw',
	isa => 'Str',
);

has 'linked' => (
	is => 'rw',
	isa => 'ArrayRef[gene]',
	predicate => 'is_linked',
);

1;
