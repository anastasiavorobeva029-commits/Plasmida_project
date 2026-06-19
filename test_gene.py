import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Restriction


def read_fasta_file(file_path):

    sequences = []

    for record in SeqIO.parse(file_path, 'fasta'):
        sequences.append({
            'id': record.id,
            'description': record.description,
            'sequence': str(record.seq),
            'length': len(record.seq)
        })

    df = pd.DataFrame(sequences)

    return df

def find_site_in_circular_dna(sequence_str, enzyme_1, enzyme_2):

    original_length = len(sequence_str)

    seq_obj = Seq(sequence_str * 2)

    eco_sites = enzyme_1.search(seq_obj)
    hind_sites = enzyme_2.search(seq_obj)

    valid_eco_sites = [pos for pos in eco_sites if pos < original_length]
    valid_hind_sites = [pos for pos in hind_sites if pos < original_length]


    return valid_eco_sites, valid_hind_sites


def cut_circle_dna(sequence_str, pos_1, enz_1, pos_2, enz_2):

    cut_1 = pos_1 + enz_1.fst3 - 1
    cut_2 = pos_2 + enz_2.fst3 - 1

    # Определяем внутренний и внешний фрагменты
    if cut_1 < cut_2:
        # Внутренний фрагмент между сайтами
        fragment_inside = sequence_str[cut_1:cut_2]
        # Внешний фрагмент (остов плазмиды)
        fragment_outside = sequence_str[cut_2:] + sequence_str[:cut_1]
        return fragment_inside, fragment_outside
    else:
        # Если сайты идут в обратном порядке
        fragment_inside = sequence_str[cut_1:] + sequence_str[:cut_2]
        fragment_outside = sequence_str[cut_2:cut_1]


        return fragment_inside, fragment_outside


def check_compatibility(enz_donor, enz_vector):
    """
    Проверяет совместимость концов рестрикции для донорского и векторного ферментов.
    """
    # Если это один и тот же фермент - совместимы
    if str(enz_donor) == str(enz_vector):
        return True

    # Если оба тупоконечные - совместимы
    if enz_donor.is_blunt() and enz_vector.is_blunt():
        return True

    # Проверяем оверханги (липкие концы)
    if enz_donor.is_3overhang() == enz_vector.is_3overhang():
        return str(enz_donor.ovhgseq) == str(enz_vector.ovhgseq)

    return False

def main():

    donor_file = r'C:\Users\Nastya\PycharmProjects\plasmida_project\data\pUC19_sequence.fasta'
    vector_file = r'C:\Users\Nastya\PycharmProjects\plasmida_project\data\J01749.1.fasta'

    donor_df = read_fasta_file(donor_file)
    vector_df = read_fasta_file(vector_file)

    donor_seq = donor_df.loc[0, 'sequence']
    vector_seq = vector_df.loc[0, 'sequence']

    print(f"Донор (pUC19):     {len(donor_seq)} п.н. - {donor_df.loc[0, 'id']}")
    print(f"Вектор (pBR322):   {len(vector_seq)} п.н. - {vector_df.loc[0, 'id']}")

    enzyme_1_donor = Restriction.AllEnzymes.get('EcoRI')
    enzyme_2_donor = Restriction.AllEnzymes.get('HindIII')

    enzyme_1_vector = Restriction.AllEnzymes.get('EcoRI')
    enzyme_2_vector = Restriction.AllEnzymes.get('HindIII')

    print(f"Ферменты: {enzyme_1_donor} и {enzyme_2_donor}")

    eco_sites_d, hind_sites_d = find_site_in_circular_dna(
        donor_seq, enzyme_1_donor, enzyme_2_donor
    )

    # Для вектора
    eco_sites_v, hind_sites_v = find_site_in_circular_dna(
        vector_seq, enzyme_1_vector, enzyme_2_vector
    )

    print(f"Донор (pUC19):")
    print(f"  Сайты EcoRI:   {eco_sites_d}")
    print(f"  Сайты HindIII: {hind_sites_d}")

    print(f"Вектор (pBR322):")
    print(f"  Сайты EcoRI:   {eco_sites_v}")
    print(f"  Сайты HindIII: {hind_sites_v}")
    print()

    pos1_d = eco_sites_d[0]
    pos2_d = hind_sites_d[0]
    pos1_v = eco_sites_v[0]
    pos2_v = hind_sites_v[0]

    print(f"Выбранные сайты:")
    print(f"  Донор:  EcoRI={pos1_d}, HindIII={pos2_d}")
    print(f"  Вектор: EcoRI={pos1_v}, HindIII={pos2_v}")
    print()

    if not check_compatibility(enzyme_1_donor, enzyme_1_vector):
        raise ValueError(f"Несовместимы концы: {enzyme_1_donor} и {enzyme_1_vector}")
    if not check_compatibility(enzyme_2_donor, enzyme_2_vector):
        raise ValueError(f"Несовместимы концы: {enzyme_2_donor} и {enzyme_2_vector}")

    insert, _ = cut_circle_dna(donor_seq, pos1_d, enzyme_1_donor, pos2_d, enzyme_2_donor)

    # Из вектора берем остов (внешний фрагмент)
    _, vector_backbone = cut_circle_dna(vector_seq, pos1_v, enzyme_1_vector, pos2_v, enzyme_2_vector)

    print(f"Длина вставки:        {len(insert)} п.н.")
    print(f"Длина остова вектора: {len(vector_backbone)} п.н.")
    print()

    final_plasmid = vector_backbone + insert

    print(f"Длина исходного донора (pUC19):   {len(donor_seq)} п.н.")
    print(f"Длина исходного вектора (pBR322): {len(vector_seq)} п.н.")
    print(f"Длина рекомбинантной плазмиды:    {len(final_plasmid)} п.н.")
    print()

    if len(final_plasmid) == len(vector_seq):
        print("⚠ Длина не изменилась - возможно, вставка не добавилась")
    elif len(final_plasmid) > len(vector_seq):
        print(f"✓ Плазмида увеличилась на {len(final_plasmid) - len(vector_seq)} п.н.")
        print("  Клонирование прошло успешно!")
    else:
        print(f"⚠ Плазмида уменьшилась на {len(vector_seq) - len(final_plasmid)} п.н.")
        print("  Возможно, произошла делеция")

        # Проверяем, не совпадает ли с исходными плазмидами
    if final_plasmid == donor_seq:
        print("⚠ Рекомбинантная плазмида идентична донорской!")
    elif final_plasmid == vector_seq:
        print("⚠ Рекомбинантная плазмида идентична векторной!")
    else:
        print("✓ Рекомбинантная плазмида уникальна")



if __name__ == "__main__":
    main()






