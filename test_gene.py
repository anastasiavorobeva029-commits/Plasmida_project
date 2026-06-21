import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Restriction


def read_fasta_file(file_path):
    sequences = []
    for record in SeqIO.parse(file_path, "fasta"):
        sequences.append(
            {
                "id": record.id,
                "description": record.description,
                "sequence": str(record.seq).upper(),
                "length": len(record.seq),
            }
        )
    return pd.DataFrame(sequences)


def find_site_in_circular_dna(sequence_str, enzyme, verbose=False):

    if not sequence_str:
        raise ValueError("Последовательность не может быть пустой")
    if enzyme is None:
        raise ValueError("Фермент не найден в библиотеке")

    original_length = len(sequence_str)
    # Для кольцевой ДНК удваиваем строку, чтобы поймать сайты на стыке краев
    doubled_seq = Seq(sequence_str * 2)

    # Поиск сайтов
    sites = np.array(enzyme.search(doubled_seq))

    # Оставляем только те сайты, которые начинаются в пределах первой (оригинальной) копии
    # Используем < вместо <=, так как сайт не может начинаться на позиции original_length
    valid_sites = sites[sites < original_length].tolist()
    valid_sites = sorted(set(valid_sites))

    if verbose:
        print(f"  Фермент {enzyme}: {len(valid_sites)} сайтов. Позиции: {valid_sites}")

    return valid_sites


def cut_circle_dna(sequence_str, pos, enz_1, enz_2):

    # Переводим 1-based в 0-based
    cut_1 = pos[0] - 1
    cut_2 = pos[1] - 1

    # Сортируем для единообразия
    start = min(cut_1, cut_2)
    end = max(cut_1, cut_2)

    # Вставка - это последовательность между сайтами
    fragment_between = sequence_str[start:end]

    # Остов - всё остальное (кольцевая ДНК)
    fragment_outside = sequence_str[end:] + sequence_str[:start]

    return fragment_between, fragment_outside

# теперь проверим совместимость концов по гороскопу
def check_compatibility(enz_donor, enz_vector):
    # Оба тупые - идеально
    if enz_donor.is_blunt() and enz_vector.is_blunt():
        return True

    # Один тупой, другой липкий - несовместимы
    if enz_donor.is_blunt() != enz_vector.is_blunt():
        return False

    # Проверяем типы оверхангов (5' или 3')
    if np.sign(enz_donor.ovhg) != np.sign(enz_vector.ovhg):
        return False

    # Проверяем длину
    if abs(enz_donor.ovhg) != abs(enz_vector.ovhg):
        return False

    # Для совместимости оверханги должны быть обратно-комплементарны
    donor_overhang = str(enz_donor.ovhgseq)
    vector_overhang = str(enz_vector.ovhgseq)

    # Проверяем, комплементарны ли концы
    # (один из них должен быть обратно-комплементарным другому)
    return donor_overhang == str(Seq(vector_overhang).reverse_complement())


# проверяем нашу криворукость, то есть, безопасность
def validate_site(sites, enzyme_name, plasmid_name):

    if not sites:
        raise ValueError(f"Сайт {enzyme_name} не найден в {plasmid_name}")

    if len(sites) > 1:
        raise ValueError(
            f"Сайт {enzyme_name} встречается {len(sites)} раз в {plasmid_name}. Клонирование невозможно."
        )

    return sites[0]


# жмем на рычаг и бежим клонировать
def perform_cloning(
        donor_seq,
        vector_seq,
        enz_d1_name="EcoRI",
        enz_d2_name="HindIII",
        enz_v1_name="EcoRI",
        enz_v2_name="HindIII",
        insert_orientation="forward"
):
    # Проверка на одинаковые ферменты
    if enz_d1_name == enz_d2_name:
        print("ОШИБКА: Для донора выбраны одинаковые ферменты!")
        return None
    if enz_v1_name == enz_v2_name:
        print("ОШИБКА: Для вектора выбраны одинаковые ферменты!")
        return None

    enz_names = np.array([enz_d1_name, enz_d2_name, enz_v1_name, enz_v2_name])
    enz_objects = [Restriction.AllEnzymes.get(name) for name in enz_names]

    if None in enz_objects:
        invalid_mask = np.equal(enz_objects, None)
        print(f"Ферменты не найдены в БД: {', '.join(enz_names[invalid_mask])}")
        return None

    enz_d1, enz_d2, enz_v1, enz_v2 = enz_objects

    print(f"Рестриктазы донора:  {enz_d1_name} + {enz_d2_name}")
    print(f"Рестриктазы вектора: {enz_v1_name} + {enz_v2_name}\n")

# Валидация сайтов с обработкой ошибок

    pos_d1 = validate_site(
        find_site_in_circular_dna(donor_seq, enz_d1), enz_d1_name, "Донор")

    pos_d2 = validate_site(
        find_site_in_circular_dna(donor_seq, enz_d2), enz_d2_name, "Донор")

    pos_v1 = validate_site(
        find_site_in_circular_dna(vector_seq, enz_v1), enz_v1_name, "Вектор")

    pos_v2 = validate_site(
        find_site_in_circular_dna(vector_seq, enz_v2), enz_v2_name, "Вектор")


    print(f"Донор координаты (1-based): {enz_d1_name}={pos_d1}, {enz_d2_name}={pos_d2}")
    print(f"Вектор координаты (1-based): {enz_v1_name}={pos_v1}, {enz_v2_name}={pos_v2}\n")

    # Проверка совместимости стыков
    comp_left = check_compatibility(enz_d1, enz_v1)
    comp_right = check_compatibility(enz_d2, enz_v2)

    if not (comp_left and comp_right):
        print("Концы несовместимы!")

        if not comp_left:
            print(f"  - Конец после {enz_d1_name} не совместим с {enz_v1_name}")

        if not comp_right:
            print(f"  - Конец после {enz_d2_name} не совместим с {enz_v2_name}")

        return None

    print("Концы совместимы\n")

    # Резка плазмид
    insert, _ = cut_circle_dna(donor_seq, (pos_d1, pos_d2))

    # Изменяем ориентацию вставки если нужно
    if insert_orientation == "reverse":
        insert = str(Seq(insert).reverse_complement())
        print("Вставка реверсирована (обратная ориентация)")

    # Вырезаем остов из вектора
    _, vector_backbone = cut_circle_dna(vector_seq, (pos_v1, pos_v2))

    # Лигазное сшивание
    final_plasmid = vector_backbone + insert

    print(f"Вставка:   {len(insert):>5} п.н.")
    print(f"Остов:     {len(vector_backbone):>5} п.н.")
    print(f"Результат: {len(final_plasmid):>5} п.н.\n")

    return final_plasmid



def print_results(final_plasmid, vector_seq):

    if final_plasmid is None:
        return

    diff = len(final_plasmid) - len(vector_seq)
    if diff > 0:
        print(f"Размер плазмиды увеличился на {diff} п.н.")

    elif diff < 0:
        print(f"Размер плазмиды уменьшился на {abs(diff)} п.н.")

    else:
        print("Итоговая длина совпадает с исходным вектором.")

    if diff == 0 and final_plasmid == vector_seq:
        print("Последовательность идентична исходному вектору.")

    else:
        print("Создана уникальная рекомбинантная плазмида")

def save_plasmid(sequence, filepath):

    with open(filepath, "w") as f:
        f.write(f">recombinant_plasmid_{len(sequence)}bp\n")
        for i in range(0, len(sequence), 60):
            f.write(sequence[i : i + 60] + "\n")
    print(f"Результат успешно сохранен: {filepath}")



def main():
    # Пути к файлам в твоем проекте
    donor_file = r"C:\Users\Nastya\PycharmProjects\plasmida_project\data\pUC19_sequence.fasta"
    vector_file = r"C:\Users\Nastya\PycharmProjects\plasmida_project\data\pBR322.fasta"
    output_file = r"C:\Users\Nastya\PycharmProjects\plasmida_project\data\recombinant_plasmid.fasta"

    # Чтение данных

    donor_df = read_fasta_file(donor_file)
    vector_df = read_fasta_file(vector_file)

    donor_seq = donor_df.loc[0, "sequence"]
    vector_seq = vector_df.loc[0, "sequence"]

    print(f"Донор: {donor_df.loc[0, 'id']} - {donor_df.loc[0, 'length']} п.н.")
    print(f"Вектор: {vector_df.loc[0, 'id']} - {vector_df.loc[0, 'length']} п.н.\n")

    print("Введите ферменты для клонирования:")
    print("Доступные ферменты: EcoRI, HindIII, BamHI, XhoI, SalI, PstI, KpnI, NdeI\n")

    # Считываем ввод и нормализуем регистр названий (делаем первую букву заглавной, остальные как в стандартах,
    # но лучше просто убрать пробелы, а в perform_cloning Biopython сам разберется, если они соответствуют номенклатуре)
    enz_d1 = input("1-й фермент для донора: ").strip()
    enz_d2 = input("2-й фермент для донора: ").strip()
    enz_v1 = input("1-й фермент для вектора: ").strip()
    enz_v2 = input("2-й фермент для вектора: ").strip()

    # Выполнение клонирования
    print("\n" + "=" * 60)
    final_plasmid = perform_cloning(
        donor_seq,
        vector_seq,
        enz_d1,
        enz_d2,
        enz_v1,
        enz_v2
    )

    if final_plasmid:
        print("\n" + "=" * 60)
        print_results(final_plasmid, vector_seq)
        save_plasmid(final_plasmid, output_file)
    else:
        print("\nКлонирование не удалось.")


if __name__ == "__main__":
    main()


