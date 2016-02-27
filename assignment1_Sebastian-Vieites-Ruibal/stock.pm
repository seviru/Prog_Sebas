package stock;
use Moose;

has 'seed' => ( # Id of the seed
	is => 'rw',
	isa => 'Str',
	required => 1,
);

has 'mutant_id' => ( # Id of the mutant gene in the seed
	is => 'rw',
	isa => 'gene',
	required => 1,
);

has 'last_plantation' => ( # Date of seed's last plantation
	is => 'rw',
	isa => 'Str',
);

has 'storage' => ( # Id of the storage
	is => 'rw',
	isa => 'Str',
);

has 'grams' => ( # Grams remaining of the seed
	is => 'rw',
	isa => 'Int',
	required => 1,
);

1;
