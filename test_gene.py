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

    sites_1 = np.array(enzyme_1.search(seq_obj))
    sites_2 = np.array(enzyme_2.search(seq_obj))

    valid_sites_1 = sites_1[sites_1 < original_length].tolist()
    valid_sites_2 = sites_2[sites_2 < original_length].tolist()

    return valid_sites_1, valid_sites_2

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

def validate_site(sites, enzyme_name, plasmid_name):
    # Проверяем, что сайт встречается не больше одного раза
    if not sites:
        print(f"Сайт {enzyme_name} не найден в {plasmid_name}!")

    if len(sites) > 1:
        print(f"Сайт {enzyme_name} встречается {len(sites)} раз в {plasmid_name}!")

    return sites[0]


def load_plasmids(donor_path, vector_path, donor_name="Донор", vector_name="Вектор"):
    """Загружает плазмиды из файлов с заданными именами."""
    donor_df = read_fasta_file(donor_path)
    vector_df = read_fasta_file(vector_path)

    donor_seq = donor_df.loc[0, 'sequence']
    vector_seq = vector_df.loc[0, 'sequence']

    # Используем ПОНЯТНЫЕ имена, а не ID из файла
    print(f"Донор: {donor_name} - {len(donor_seq)} п.н.")
    print(f"Вектор: {vector_name} - {len(vector_seq)} п.н.\n")

    return donor_seq, vector_seq, donor_name, vector_name


def find_and_validate_sites(seq, enz1, enz2, plasmid_name):
    """Находит и валидирует сайты рестрикции."""
    sites1, sites2 = find_site_in_circular_dna(seq, enz1, enz2)

    pos1 = validate_site(sites1, str(enz1), plasmid_name)
    pos2 = validate_site(sites2, str(enz2), plasmid_name)

    print(f"{plasmid_name}: {enz1}={pos1}, {enz2}={pos2}")
    return pos1, pos2


def perform_cloning(donor_seq, vector_seq, enzyme1_name='EcoRI', enzyme2_name='HindIII'):

    # 1. Получаем ферменты
    enz1 = Restriction.AllEnzymes.get(enzyme1_name)
    enz2 = Restriction.AllEnzymes.get(enzyme2_name)

    if not enz1 or not enz2:
        missing = [name for name, enz in [(enzyme1_name, enz1), (enzyme2_name, enz2)] if not enz]
        print(f"Ферменты не найдены в базе данных: {', '.join(missing)}")

    print(f"Используемые ферменты: {enz1} и {enz2}\n")

    # 2. Проверяем совместимость (так как ферменты одинаковые, это скорее sanity check)
    if not check_compatibility(enz1, enz1) or not check_compatibility(enz2, enz2):
        print("Крипты совместимости не пройдены: ферменты не совместимы сами с собой.")

    print("Концы совместимы\n")


    pos1_d, pos2_d = find_and_validate_sites(donor_seq, enz1, enz2, "Донор")
    pos1_v, pos2_v = find_and_validate_sites(vector_seq, enz1, enz2, "Вектор")
    print()

    # 4. Разрезаем и лигируем
    insert, _ = cut_circle_dna(donor_seq, pos1_d, enz1, pos2_d, enz2)
    _, vector_backbone = cut_circle_dna(vector_seq, pos1_v, enz1, pos2_v, enz2)

    final_plasmid = vector_backbone + insert

    # 5. Красивый вывод результатов
    print(f"Вставка:   {len(insert):>5} п.н.")
    print(f"Остов:     {len(vector_backbone):>5} п.н.")
    print(f"Результат: {len(final_plasmid):>5} п.н.\n")

    return final_plasmid

def print_results(final_plasmid, donor_seq, vector_seq):

    len_final = len(final_plasmid)
    len_vector = len(vector_seq)
    diff = len_final - len_vector

    # 1. Оценка изменения размера
    if diff > 0:
        print(f"Размер увеличился на {diff} п.н.")
    elif diff < 0:
        print(f" Размер уменьшился на {abs(diff)} п.н.")
    else:
        print("Итоговая длина совпадает с исходным вектором.")

    identity_map = {
        donor_seq: "Последовательность идентична донору.",
        vector_seq: "Последовательность идентична вектору"
    }
    status = identity_map.get(final_plasmid, "Создана уникальная рекомбинантная плазмида.")

    print(status)

def save_plasmid(sequence, filepath):
    """Сохраняет плазмиду в FASTA файл."""
    with open(filepath, 'w') as f:
        f.write(f">recombinant_plasmid_{len(sequence)}bp\n")
        for i in range(0, len(sequence), 60):
            f.write(sequence[i:i+60] + "\n")
    print(f"\nСохранено: {filepath}")

def main():

    donor_file = r'C:\Users\Nastya\PycharmProjects\plasmida_project\data\pUC19_sequence.fasta'
    vector_file = r'C:\Users\Nastya\PycharmProjects\plasmida_project\data\pBR322.fasta'

    donor_name = "pUC19"
    vector_name = "pBR322"

    donor_seq, vector_seq, donor_name, vector_name = load_plasmids(
        donor_file, vector_file, donor_name, vector_name
    )

    # Выполняем клонирование
    final_plasmid = perform_cloning(donor_seq, vector_seq)

    # Выводим результаты
    print_results(final_plasmid, donor_seq, vector_seq)

    # Сохраняем результат
    output_file = r'C:\Users\Nastya\PycharmProjects\plasmida_project\data\recombinant_plasmid.fasta'
    save_plasmid(final_plasmid, output_file)



if __name__ == "__main__":
    main()






