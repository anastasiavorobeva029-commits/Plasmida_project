import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Restriction
import itertools


#Возьмем фаста файл и сделаем из него табличку, с которой будет удобно работать
def read_fasta_file(file_path):

    sequences = []

    for record in SeqIO.parse(file_path, "fasta"):
        sequences.append(
            {
                "id": record.id,
                "description": record.description,
                "sequence": str(record.seq),
                "length": len(record.seq),
            }
        )

    df = pd.DataFrame(sequences)

    return df


# теперь будем искать сайты рестрикции в плазмиде
def find_site_in_circular_dna(sequence_str, enzyme):

    original_length = len(sequence_str)

    seq_obj = Seq(sequence_str * 2) # так как это плазмида, нам лучше удвоить последовательность

    sites = np.array(enzyme.search(seq_obj)) # ищем координаты рестрикции.

    # оставляем только те сайты, которые начались в пределах оригинальной длины. дубликаты брать не будем
    valid_sites = sites[sites <= original_length].tolist()

    return valid_sites


def suggest_cloning_enzymes(donor_seq, vector_seq):

    # попробуем подобрать "хорошую" пару ферментов рестрикции
    # для этого возьмем все фермеенты, что есть в библиотеке
    common_enzymes = Restriction.CommOnly

    # проверяем, точно ли фермент найдет сайт один раз, нам же не нужно, чтобы ДНК развалилось
    donor_candidates = [str(enz) for enz in common_enzymes if len(find_site_in_circular_dna(donor_seq, enz)) == 1]
    vector_candidates = [str(enz) for enz in common_enzymes if len(find_site_in_circular_dna(vector_seq, enz)) == 1]

    # а для этого, найдем ферменты, которые уникальны и для донора, и для вектора одновременно
    shared_candidates = list(set(donor_candidates) & set(vector_candidates))

    ideal_pairs = []

    # проверяем ферменты на их совместимость
    for name1, name2 in itertools.combinations(shared_candidates, 2):

        enz1 = Restriction.AllEnzymes.get(name1)
        enz2 = Restriction.AllEnzymes.get(name2)

        # ставим условие, что нам нужен разрез, который даст липкие концы
        if enz1.is_blunt() or enz2.is_blunt():
            continue

        # ещё одно условие - если концы совместимы, то их нужно пропустить, а то схлопнется сам в себя
        if check_compatibility(enz1, enz2):
            continue

        # Если пара идеальна, сохраняем её
        ideal_pairs.append((name1, name2))

    if not ideal_pairs:
       print("Не удалось найти идеальную пару ферментов")

    best_pair = ideal_pairs[0]

    print(f"Найдено отличных пар для направленного клонирования: {len(ideal_pairs)}")
    print(f"Автоматически выбрана наилучшая пара: {best_pair[0]} + {best_pair[1]}\n")

    return best_pair[0], best_pair[1]

# имитируем гидролиз, пытаясь сохранить топологию кольца
def cut_circle_dna(sequence_str, pos_1, enz_1, pos_2, enz_2):

    # плазмида - это мир-кольцо, так что мы считаем позицию начала сайта узнавания + смещение разреза внутри сайта - 1
    # где 1 это переход к индексации Python
    cut_1 = pos_1 + enz_1.fst3 - 1
    cut_2 = pos_2 + enz_2.fst3 - 1

    if cut_1 < cut_2: # фрагмент между разрезами может быть участком от cut_1 до cut_2
        fragment_between = sequence_str[cut_1:cut_2]

        # для внешнего фрагмента берём хвост от cut_2 до конца и прибавляем начало от 0 до cut_1
        # это как шов у плазмиды
        fragment_outside = sequence_str[cut_2:] + sequence_str[:cut_1]

    else: # фрагмент между разрезами может быть участком от cut_1 до cut_2, но проходящим через начало
        fragment_between = sequence_str[cut_1:] + sequence_str[:cut_2]
        fragment_outside = sequence_str[cut_2:cut_1]

    return fragment_between, fragment_outside

# теперь проверим совместимость концов по гороскопу
def check_compatibility(enz_donor, enz_vector):

    # любые два тупых конца могут быть сшиты лигазой, так как у них нет одноцепочных хвостов
    if enz_donor.is_blunt() and enz_vector.is_blunt():
        return True

    # если этот конец тупой, а второй конец - липкий, у них не может получиться ничего хорошего, может, банан только
    if enz_donor.is_blunt() != enz_vector.is_blunt():
        return False

    # проверяем направление - 5' и 3' направлены в разные стороны и не смогут образовать водородные связи друг с другом
    if enz_donor.is_5overhang() != enz_vector.is_5overhang():
        return False

    # проверяем комплементарность липких концов
    return str(enz_donor.ovhgseq) == str(enz_vector.ovhgseq)

# проверяем нашу криворукость, то есть, безопасность
def validate_site(sites, enzyme_name, plasmid_name):

    if not sites: # если ничего не нашли, пытаемся привлечь к себе внимание
        raise ValueError(f"Сайт {enzyme_name} не найден в {plasmid_name}")

    if len(sites) > 1: # мы не хотим развала плазмиды, поэтому, лучше 7 раз померить, один раз...
        raise ValueError(f"Сайт {enzyme_name} встречается {len(sites)} раз в {plasmid_name}")

    return sites[0]

# жмем на рычаг и бежим клонировать
def perform_cloning(
    donor_seq,
    vector_seq,
    enz_d1_name="EcoRI",
    enz_d2_name="HindIII",
    enz_v1_name="EcoRI",
    enz_v2_name="HindIII",
):
    # собираем все названия ферментов в массив
    enz_names = np.array([enz_d1_name, enz_d2_name, enz_v1_name, enz_v2_name])
    enz_objects = [Restriction.AllEnzymes.get(name) for name in enz_names]

    # проверяем, все ли ферменты найдены в базе
    if None in enz_objects:
        invalid_mask = np.equal(enz_objects, None)
        print(
            f"Ферменты не найдены в базе данных: {', '.join(enz_names[invalid_mask])}"
        )
        return None  # добавили return, чтобы код не падал дальше, если фермент не найден

    enz = np.array(enz_objects)
    enz_d1, enz_d2, enz_v1, enz_v2 = enz

    print(f"Рестриктазы донора:  {enz_d1_name} + {enz_d2_name}")
    print(f"Рестриктазы вектора: {enz_v1_name} + {enz_v2_name}\n")

    # проверка совместимости концов
    comp_left = check_compatibility(enz_d1, enz_v1)
    comp_right = check_compatibility(enz_d2, enz_v2)

    if not (comp_left and comp_right):
        print("Концы несовместимы")
        if not comp_left:
            print(
                f"  - Конец после {enz_d1_name} не совместим с концом после {enz_v1_name}"
            )
        if not comp_right:
            print(
                f"  - Конец после {enz_d2_name} не совместим с концом после {enz_v2_name}"
            )
        return None

    print("Концы совместимы\n")

    # поиск сайтов
    pos_d1 = validate_site(
        find_site_in_circular_dna(donor_seq, enz_d1), enz_d1_name, "Донор"
    )
    pos_d2 = validate_site(
        find_site_in_circular_dna(donor_seq, enz_d2), enz_d2_name, "Донор"
    )
    pos_v1 = validate_site(
        find_site_in_circular_dna(vector_seq, enz_v1), enz_v1_name, "Вектор"
    )
    pos_v2 = validate_site(
        find_site_in_circular_dna(vector_seq, enz_v2), enz_v2_name, "Вектор"
    )

    print(
        f"Донор координаты: {enz_d1_name}={pos_d1}, {enz_d2_name}={pos_d2}"
    )
    print(
        f"Вектор координаты: {enz_v1_name}={pos_v1}, {enz_v2_name}={pos_v2}\n"
    )

    # резка плазмид
    insert, _ = cut_circle_dna(donor_seq, pos_d1, enz_d1, pos_d2, enz_d2)
    _, vector_backbone = cut_circle_dna(
        vector_seq, pos_v1, enz_v1, pos_v2, enz_v2
    )

    # лигирование
    len_d1_ovhg = 0 if enz_d1.is_blunt() else len(enz_d1.ovhgseq)
    len_d2_ovhg = 0 if enz_d2.is_blunt() else len(enz_d2.ovhgseq)
    total_ovhg = len_d1_ovhg + len_d2_ovhg

    # склеиваем строки и убираем дублирующийся липкий конец
    final_plasmid = (
        (vector_backbone + insert)[:-total_ovhg]
        if total_ovhg > 0
        else vector_backbone + insert
    )

    # расчет физических длин фрагментов дуплекса
    real_insert_len = len(insert) - total_ovhg / 2
    real_backbone_len = len(vector_backbone) - total_ovhg / 2

    print(f"Вставка:   {int(real_insert_len):>5} п.н.")
    print(f"Остов:     {int(real_backbone_len):>5} п.н.")
    print(f"Результат: {len(final_plasmid):>5} п.н.\n")

    return final_plasmid

def print_results(final_plasmid, vector_seq):

    if final_plasmid is None:
        return

    diff = len(final_plasmid) - len(vector_seq)

    if diff > 0:
        print(f"Размер увеличился на {diff} п.н.")
    elif diff < 0:
        print(f"Размер уменьшился на {abs(diff)} п.н.")
    else:
        print("Итоговая длина совпадает с исходным вектором.")

    # если длины совпали и последовательности одинаковы - это исходный вектор
    if diff == 0 and final_plasmid == vector_seq:
        print("Последовательность идентична вектору.")
    else:
        print("Создана уникальная рекомбинантная плазмида.")


def save_plasmid(sequence, filepath):

    with open(filepath, "w") as f:
        f.write(f">recombinant_plasmid_{len(sequence)}bp\n")
        for i in range(0, len(sequence), 60):
            f.write(sequence[i: i + 60] + "\n")

    print(f"\nСохранено: {filepath}")


def main():

    donor_file = r"C:\Users\Nastya\PycharmProjects\plasmida_project\data\pUC19_sequence.fasta"
    vector_file = r"C:\Users\Nastya\PycharmProjects\plasmida_project\data\pBR322.fasta"
    output_file = r"C:\Users\Nastya\PycharmProjects\plasmida_project\data\recombinant_plasmid.fasta"

    donor_df = read_fasta_file(donor_file)
    vector_df = read_fasta_file(vector_file)

    donor_seq = donor_df.loc[0, "sequence"]
    vector_seq = vector_df.loc[0, "sequence"]

    # используем уже готовые значения длины из датафрейма вместо повторного вызова len()
    print(f"Донор: pUC19 - {donor_df.loc[0, 'length']} п.н.")
    print(f"Вектор: pBR322 - {vector_df.loc[0, 'length']} п.н.\n")

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
        print_results(final_plasmid, vector_seq)
        save_plasmid(final_plasmid, output_file)


if __name__ == "__main__":
    main()





