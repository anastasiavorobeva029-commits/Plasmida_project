import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Restriction
import itertools

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


def find_site_in_circular_dna(sequence_str, enzyme):

    original_length = len(sequence_str)

    seq_obj = Seq(sequence_str * 2)

    sites = np.array(enzyme.search(seq_obj))

    # Отбираем сайты, начавшиеся в первой копии плазмиды
    valid_sites = sites[sites <= original_length].tolist()

    return valid_sites


def suggest_cloning_enzymes(donor_seq, vector_seq):
    """Автоматически подбирает идеальную пару ферментов для направленного клонирования."""
    common_enzymes = Restriction.CommOnly
    donor_candidates = [str(enz) for enz in common_enzymes if len(find_site_in_circular_dna(donor_seq, enz)) == 1]
    vector_candidates = [str(enz) for enz in common_enzymes if len(find_site_in_circular_dna(vector_seq, enz)) == 1]

    # Находим общие single cutters
    shared_candidates = list(set(donor_candidates) & set(vector_candidates))

    ideal_pairs = []

    # Перебираем все возможные пары ферментов (сочетания по 2)
    for name1, name2 in itertools.combinations(shared_candidates, 2):
        enz1 = Restriction.AllEnzymes.get(name1)
        enz2 = Restriction.AllEnzymes.get(name2)

        # Условие 1: Оба фермента должны быть липкоконечными (это эффективнее)
        if enz1.is_blunt() or enz2.is_blunt():
            continue

        # Условие 2: Ферменты должны иметь разные липкие концы (защита от самолигации вектора)
        if check_compatibility(enz1, enz2):
            continue

        # Если пара идеальна, сохраняем её
        ideal_pairs.append((name1, name2))

    if not ideal_pairs:
       print("Не удалось автоматически найти идеальную пару ферментов для направленного клонирования!")

    best_pair = ideal_pairs[0]

    print(f"Найдено отличных пар для направленного клонирования: {len(ideal_pairs)}")
    print(f"Автоматически выбрана наилучшая пара: {best_pair[0]} + {best_pair[1]}\n")

    return best_pair[0], best_pair[1]


def cut_circle_dna(sequence_str, pos_1, enz_1, pos_2, enz_2):
    cut_1 = pos_1 + enz_1.fst3 - 1
    cut_2 = pos_2 + enz_2.fst3 - 1

    if cut_1 < cut_2:
        fragment_between = sequence_str[cut_1:cut_2]
        fragment_outside = sequence_str[cut_2:] + sequence_str[:cut_1]

    else:
        fragment_between = sequence_str[cut_1:] + sequence_str[:cut_2]
        fragment_outside = sequence_str[cut_2:cut_1]

    return fragment_between, fragment_outside

def check_compatibility(enz_donor, enz_vector):

    if enz_donor.is_blunt() and enz_vector.is_blunt():
        return True

    if enz_donor.is_blunt() != enz_vector.is_blunt():
        return False

    if enz_donor.is_5overhang() != enz_vector.is_5overhang():
        return False

    return str(enz_donor.ovhgseq) == str(enz_vector.ovhgseq)

def validate_site(sites, enzyme_name, plasmid_name):

    if not sites:
        print(f"Критическая ошибка: Сайт {enzyme_name} не найден в {plasmid_name}!")

    if len(sites) > 1:
        print(f"Критическая ошибка: Сайт {enzyme_name} встречается {len(sites)} раз в {plasmid_name}!")

    return sites[0]


def load_plasmids(donor_path, vector_path, donor_name="Донор", vector_name="Вектор"):

    donor_df = read_fasta_file(donor_path)
    vector_df = read_fasta_file(vector_path)

    donor_seq = donor_df.loc[0, 'sequence']
    vector_seq = vector_df.loc[0, 'sequence']

    print(f"Донор: {donor_name} - {len(donor_seq)} п.н.")
    print(f"Вектор: {vector_name} - {len(vector_seq)} п.н.\n")

    return donor_seq, vector_seq, donor_name, vector_name


def find_and_validate_sites(seq, enz1, enz2, plasmid_name):

    sites1, sites2 = find_site_in_circular_dna(seq, enz1)

    pos1 = validate_site(sites1, str(enz1), plasmid_name)
    pos2 = validate_site(sites2, str(enz2), plasmid_name)

    print(f"{plasmid_name}: {enz1}={pos1}, {enz2}={pos2}")
    return pos1, pos2


def perform_cloning(
        donor_seq,
        vector_seq,
        enz_d1_name="EcoRI",
        enz_d2_name="HindIII",
        enz_v1_name="EcoRI",
        enz_v2_name="HindIII",
):
    # 1. Валидация наличия ферментов в базе
    enz_names = np.array([enz_d1_name, enz_d2_name, enz_v1_name, enz_v2_name])
    enz_objects = [Restriction.AllEnzymes.get(name) for name in enz_names]

    if None in enz_objects:
        invalid_mask = np.equal(enz_objects, None)
        print(f"Ферменты не найдены в базе данных: {', '.join(enz_names[invalid_mask])}")

    enz = np.array(enz_objects)

    print(f"Рестриктазы донора:  {enz_d1_name} + {enz_d2_name}")
    print(f"Рестриктазы вектора: {enz_v1_name} + {enz_v2_name}\n")

    # 2. Проверка совместимости концов
    enz_d1, enz_d2, enz_v1, enz_v2 = enz

    comp_left = check_compatibility(enz_d1, enz_v1)
    comp_right = check_compatibility(enz_d2, enz_v2)

    if not (comp_left and comp_right):
        print("Концы несовместимы!")

        if not comp_left:
            print(f"  - Конец после {enz_d1_name} не совместим с концом после {enz_v1_name}")

        if not comp_right:
            print(f"  - Конец после {enz_d2_name} не совместим с концом после {enz_v2_name}")

        return None

    print("Концы совместимы\n")

    # 3. Поиск сайтов (теперь передаем чистые переменные, IDE довольна)
    pos_d1 = validate_site(find_site_in_circular_dna(donor_seq, enz_d1), enz_d1_name, "Донор")
    pos_d2 = validate_site(find_site_in_circular_dna(donor_seq, enz_d2), enz_d2_name, "Донор")
    pos_v1 = validate_site(find_site_in_circular_dna(vector_seq, enz_v1), enz_v1_name, "Вектор")
    pos_v2 = validate_site(find_site_in_circular_dna(vector_seq, enz_v2), enz_v2_name, "Вектор")

    print(f"Донор координаты: {enz_d1_name}={pos_d1}, {enz_d2_name}={pos_d2}")
    print(f"Вектор координаты: {enz_v1_name}={pos_v1}, {enz_v2_name}={pos_v2}\n")

    # 4. Резка плазмид (Движемся строго от Фермента_1 к Ферменту_2)
    insert, _ = cut_circle_dna(donor_seq, pos_d1, enz_d1, pos_d2, enz_d2)
    _, vector_backbone = cut_circle_dna(vector_seq, pos_v1, enz_v1, pos_v2, enz_v2)

    # 5. Лигирование с коррекцией на перекрытие липких концов
    is_blunt_mask = np.array([e.is_blunt() for e in enz])
    ovhg_lengths = np.array([0 if blunt else len(e.ovhgseq) for blunt, e in zip(is_blunt_mask, enz)])

    total_ovhg = ovhg_lengths[0] + ovhg_lengths[1]

    # Склеиваем строки и убираем дублирующийся липкий конец
    final_plasmid = (vector_backbone + insert)[:-total_ovhg] if total_ovhg > 0 else vector_backbone + insert

    # Расчет физических длин фрагментов дуплекса
    real_insert_len = len(insert) - total_ovhg / 2
    real_backbone_len = len(vector_backbone) - total_ovhg / 2

    print(f"Вставка:   {int(real_insert_len):>5} п.н.")
    print(f"Остов:     {int(real_backbone_len):>5} п.н.")
    print(f"Результат: {len(final_plasmid):>5} п.н.\n")

    return final_plasmid

def print_results(final_plasmid, donor_seq, vector_seq):
    if final_plasmid is None:
        return

    diff = len(final_plasmid) - len(vector_seq)

    if diff > 0:
        print(f"Размер увеличился на {diff} п.н.")
    elif diff < 0:
        print(f"Размер уменьшился на {abs(diff)} п.н.")
    else:
        print("Итоговая длина совпадает с исходным вектором.")

    if final_plasmid == donor_seq:
        print("Последовательность идентична донору.")
    elif final_plasmid == vector_seq:
        print("Последовательность идентична вектору.")
    else:
        print("Создана уникальная рекомбинантная плазмида.")

def save_plasmid(sequence, filepath):
    if sequence is None:
        return
    with open(filepath, "w") as f:
        f.write(f">recombinant_plasmid_{len(sequence)}bp\n")
        for i in range(0, len(sequence), 60):
            f.write(sequence[i : i + 60] + "\n")
    print(f"\nСохранено: {filepath}")

def main():
    donor_file = r"C:\Users\Nastya\PycharmProjects\plasmida_project\data\pUC19_sequence.fasta"
    vector_file = r"C:\Users\Nastya\PycharmProjects\plasmida_project\data\pBR322.fasta"
    output_file = r"C:\Users\Nastya\PycharmProjects\plasmida_project\data\recombinant_plasmid.fasta"

    donor_df = read_fasta_file(donor_file)
    vector_df = read_fasta_file(vector_file)
    donor_seq = donor_df.loc[0, "sequence"]
    vector_seq = vector_df.loc[0, "sequence"]

    print(f"Донор: pUC19 - {len(donor_seq)} п.н.")
    print(f"Вектор: pBR322 - {len(vector_seq)} п.н.\n")

    selected_enz1, selected_enz2 = suggest_cloning_enzymes(donor_seq, vector_seq)

    final_plasmid = perform_cloning(
        donor_seq,
        vector_seq,
        enz_d1_name=selected_enz1,
        enz_d2_name=selected_enz2,
        enz_v1_name=selected_enz1,
        enz_v2_name=selected_enz2,
    )

    if final_plasmid:
        print_results(final_plasmid, donor_seq, vector_seq)
        save_plasmid(final_plasmid, output_file)

if __name__ == "__main__":
    main()





