Input data:

Donor plasmid DNA sequence, from which the fragment (the "insert") is to be excised;

Vector plasmid DNA sequence, into which the insert is to be integrated;

A pair of distinct restriction enzymes for the donor plasmid (e.g., EcoRI and XhoI);

A pair of distinct restriction enzymes for the vector plasmid.
Note: The restriction enzymes specified for the donor plasmid and the vector plasmid may be identical.

It is guaranteed that the specified restriction enzyme names are available in standard Python libraries; there is no need to create a custom restriction enzyme database.

Both plasmids are considered circular. The fragment within the donor plasmid is defined as the region between the cleavage sites determined by the specified pair of restriction enzymes. This fragment—the "insert"—must be integrated into the vector plasmid, replacing the region between the cleavage sites determined by the vector's respective pair of restriction enzymes.

Output:
The resulting sequence of the recombinant vector plasmid containing the insert in the same orientation as in the original donor sequence.

The program must validate that:

The restriction sites for the specified enzymes are actually present in their respective plasmids.

Each restriction site occurs exactly once in its respective plasmid.

The ends of the insert generated after restriction digestion are compatible with the corresponding ends of the digested vector plasmid.

The program must support cloning with both blunt ends and all types of complementary sticky ends.
