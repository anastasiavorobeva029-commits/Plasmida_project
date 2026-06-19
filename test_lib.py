from Bio.Restriction import RestrictionBatch
from Bio.Seq import Seq
from Bio import Restriction

# print(dna)
# print(type(dna))
#
# print(dna.complement())
# print(dna.reverse_complement())
# print(dna.reverse_complement_rna())

# dna_1 = Seq("GAATTC" + "ATGC" + "CTCGAG" + "CCCGGG")
#
# # print(dna_1)
#
# eco_enzyme = Restriction.AllEnzymes.get('EcoRI')
# eco_sites = eco_enzyme.search(dna_1)
#
# # print(eco_sites)
#
# xh_enzyme = Restriction.AllEnzymes.get('XhoI')
# xh_sites = xh_enzyme.search(dna_1)
#
# # print(xh_sites)
#
# sm_enzyme = Restriction.AllEnzymes.get('SmaI')
# sm_sites = sm_enzyme.search(dna_1)

# print(sm_sites)

# dna_fragment = Seq('GAATTC')
#
# eco_enzyme = Restriction.AllEnzymes.get('EcoRI')
#
# eco = eco_enzyme.site
#
# print(eco)
# print(eco_enzyme.fst5)
# print(eco_enzyme.fst3)


def find_site_in_circular_dna(sequence, enzyme):

    seq = sequence * 2

    all_sites = enzyme.search(seq)

    original_length = len(sequence)

    valid_sites = [pos for pos in all_sites if pos < original_length]

    return valid_sites


def cut_circle_dna(sequence, enzyme_1, enzyme_2):

    sites_1 = find_site_in_circular_dna(sequence, enzyme_1)
    sites_2 = find_site_in_circular_dna(sequence, enzyme_2)

    pos_1 = sites_1[0]
    pos_2 = sites_2[0]

    start = min(pos_1, pos_2)
    end = max(pos_1, pos_2)

    fragment_1 = sequence[start:end]
    fragment_2 = sequence[end:] + sequence[:start]

    return fragment_1, fragment_2


dna = Seq('AAGGAATTTCCTGTCGAATTCGTCCCCATATATATCTCGAGCGGGGAGAGGA')

eco_enzyme = Restriction.AllEnzymes.get('EcoRI')
xh_enzyme = Restriction.AllEnzymes.get('XhoI')

print(find_site_in_circular_dna(dna, eco_enzyme))

frag_1, frag_2 = cut_circle_dna(dna, eco_enzyme, xh_enzyme)

print(frag_1)
print(frag_2)



