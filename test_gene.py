import os
from tkinter import filedialog, Tk
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Restriction

# чтение FASTA файла в DataFrame
def read_fasta_file(file_path):

    # инициализируем пустой список для хранения данных по каждой последовательности
    sequences = []

    # построчно парсим FASTA-файл с помощью SeqIO из Biopython
    for record in SeqIO.parse(file_path, "fasta"):
        # добавляем словарь с извлеченными данными в наш список
        sequences.append({
                "id": record.id,  # уникальный идентификатор последовательности (после '>')
                "description": record.description,  # полная строка описания
                "sequence": str(record.seq).upper(),  # сама последовательность, приведенная к верхнему регистру
                "length": len(record.seq)}) # длина последовательности (количество нуклеотид)

    return pd.DataFrame(sequences) # преобразуем список словарей в объект pandas DataFrame и возвращаем его

# поиск сайтов рестрикции в кольцевой ДНК с использованием удвоения цепочки
def find_site_in_circular_dna(sequence_str, enzyme, plasmid_type="", verbose=True):
    """Поиск сайтов рестрикции в кольцевой ДНК с использованием удвоения цепочки"""
    if not sequence_str or enzyme is None:
        raise ValueError("Пустая последовательность или фермент не найден")

    original_length = len(sequence_str)
    doubled_seq = Seq(sequence_str * 2)

    sites = np.array(enzyme.search(doubled_seq))
    valid_sites = sorted(list(set(sites[sites <= original_length])))

    if verbose:
        # Получаем имя фермента прямо из объекта
        enzyme_name = enzyme.__name__
        positions_str = ', '.join(str(int(pos)) for pos in valid_sites)
        if plasmid_type:
            print(f"  {plasmid_type}: {enzyme_name} ({enzyme.site}) → {len(valid_sites)} сайт(ов). Позиции: {positions_str}")
        else:
            print(f"  {enzyme_name} ({enzyme.site}): {len(valid_sites)} сайт(ов). Позиции: {positions_str}")

    return valid_sites

#  проверка совместимости липких или тупых концов двух ферментов рестрикции для успешного лигирования
def check_compatibility(enz_donor, enz_vector):

    # оба тупые
    if enz_donor.is_blunt() and enz_vector.is_blunt():
        return True

    # один тупой, другой липкий
    if enz_donor.is_blunt() != enz_vector.is_blunt():
        return False

    # разные типы липких концов
    if np.sign(enz_donor.ovhg) != np.sign(enz_vector.ovhg) or abs(enz_donor.ovhg) != abs(enz_vector.ovhg):
        return False

    # проверка комплементарности оверхангов
    donor_overhang = str(enz_donor.ovhgseq)
    vector_overhang = str(enz_vector.ovhgseq)

    # проверяем, совпадает ли последовательность донора с обратным комплементом вектора,
    # если они комплементарны, то смогут отжечься друг с другом для последующего сшивания лигазой
    return donor_overhang == str(Seq(vector_overhang).reverse_complement())

# валидация уникальности сайта рестрикции
def validate_site(sites, enzyme_name, plasmid_name):

    # проверка на наличие сайта: если список позиций пуст, значит фермент не режет эту ДНК
    if not sites:
        raise ValueError(f"Сайт {enzyme_name} не найден в {plasmid_name}")

    # проверка на уникальность сайта: нужно иметь один сайт рестрикции, если найдено больше - плазмида развалиться
    if len(sites) > 1:
        raise ValueError(f"Сайт {enzyme_name} встречается {len(sites)} раз в {plasmid_name} (нужен уникальный)")

    return sites[0]

# жмем на рычаг и клонируем: in silico рестрикция и лигирование фрагмента ДНК в вектор

def perform_cloning(donor_seq, vector_seq, enz_d1, enz_d2, enz_v1, enz_v2):

    # проверяем, что все переданные объекты ферментов существуют
    if None in [enz_d1, enz_d2, enz_v1, enz_v2]:
        print("Один или несколько ферментов не найдены")
        return None

    # запрещаем использовать одинаковые ферменты на одной плазмиде, чтобы не было одинаковых хвостов и проблем с ориентацией
    if enz_d1.site == enz_d2.site or enz_v1.site == enz_v2.site:
        print("Выбраны одинаковые ферменты для одной плазмиды.")
        return None

    print("Поиск сайтов рестрикции:")

    # поиск и валидация координат сайтов.
    # find_site_in_circular_dna находит координаты в кольцевой молекуле,
    # а validate_site проверяет, что фермент делает ровно один разрез.
    pos_d1 = validate_site(find_site_in_circular_dna(donor_seq, enz_d1, "Донор", verbose=True),
                           enz_d1.site,"Донор")
    pos_d2 = validate_site(find_site_in_circular_dna(donor_seq, enz_d2, "Донор", verbose=True),
                           enz_d2.site,"Донор")
    pos_v1 = validate_site(find_site_in_circular_dna(vector_seq, enz_v1, "Вектор", verbose=True),
                           enz_v1.site,"Вектор")
    pos_v2 = validate_site(find_site_in_circular_dna(vector_seq, enz_v2, "Вектор", verbose=True),
                           enz_v2.site,"Вектор")

    # проверка липких/тупых концов на совместимость.
    # должны совпадать первый фермент донора с первым ферментом вектора, и вторые ферменты между собой.
    if not check_compatibility(enz_d1, enz_v1) or not check_compatibility(
        enz_d2, enz_v2
    ):
        print("Концы несовместимы!")

        return None

    print("Концы совместимы\n")

    # вырезание целевой вставки из донора
    # переходим от 1-индексации Biopython к стандартным 0-based индексам Python
    cut_d1, cut_d2 = pos_d1 - 1, pos_d2 - 1

    # если фрагмент лежит внутри линейной строки:
    if cut_d1 < cut_d2:
        insert = donor_seq[cut_d1:cut_d2]

    # если фрагмент пересекает условную точку начала/конца файла в кольцевой ДНК:
    else:
        insert = donor_seq[cut_d1:] + donor_seq[:cut_d2]

    # вырезание остова вектора
    v1_idx, v2_idx = pos_v1 - 1, pos_v2 - 1

    # направленность сборки учитывает, что вектор открывается от v2 к v1.
    if v1_idx < v2_idx:
        vector_backbone = vector_seq[v2_idx:] + vector_seq[:v1_idx]

    else:
        vector_backbone = vector_seq[v1_idx:] + vector_seq[:v2_idx]

    # лигирование, наконец-то: физически объединяем остов вектора и вставку в новую кольцевую плазмиду
    final_plasmid = vector_backbone + insert

    # вывод результатов расчетов размеров фрагментов
    print(f"Вставка:   {len(insert):>5} п.н.")
    print(f"Остов:     {len(vector_backbone):>5} п.н.")
    print(f"Результат: {len(final_plasmid):>5} п.н.\n")

    # возвращаем словарь со всей метаинформацией о собранной плазмиде
    return {
        "final_plasmid": final_plasmid,
        "backbone_length": len(vector_backbone),
        "insert_length": len(insert),
        "pos_v1": pos_v1,
        "pos_v2": pos_v2,
        "pos_d1": pos_d1,
        "pos_d2": pos_d2,
        "enz_v1": enz_v1,
        "enz_v2": enz_v2,
    }

# вывод результатов анализа стыков
def print_results(final_plasmid, vector_seq, donor_seq,
                  enz_v1_name, enz_v2_name, enz_v1_site, enz_v2_site,
                  pos_v1, pos_v2):

    # если клонирование завершилось неудачей на предыдущем этапе, прерываем выполнение
    if final_plasmid is None:
        return

    # расчет и вывод разницы в размерах между новой конструкцией и исходным вектором
    diff = len(final_plasmid) - len(vector_seq)

    if diff > 0:
        print(f"Размер плазмиды увеличился на {diff} п.н.")

    elif diff < 0:
        print(f"Размер плазмиды уменьшился на {abs(diff)} п.н.")

    else:
        print("Итоговая длина совпадает с исходным вектором.")

    # если вектор просто замкнулся сам на себя без вставки
    print("Последовательность идентична вектору"
          if final_plasmid == vector_seq else "Создана уникальная рекомбинантная плазмида")

    print("\n Анализ стыковки")

    # ANSI-коды для цветового выделения текста в консоли
    GREEN, YELLOW, RESET = "\033[92m", "\033[93m", "\033[0m"

    # переводим координаты Biopython (1-based) в индексы срезов Python (0-based)
    v1_idx, v2_idx = pos_v1 - 1, pos_v2 - 1

    # теперь по отдельности - конец вектора соединяется с началом вставки - левый стык

    # берем 30 нуклеотидов вектора до сайта разреза v1
    # конкатенация срезов защищает от выхода за границы строки, если v1 находится близко к началу кольца
    vector_left_flank = (vector_seq[v1_idx - 30:] + vector_seq[:v1_idx])[-30:].upper()

    # Берем первые 40 нуклеотидов донора
    insert_left_flank = donor_seq[:40].upper()

    print(f"Левый стык (Вектор -> Вставка):")

    # выводим вектор, отрезая последние 6 нуклеотидов ([:-6]), которые красятся желтым как сайт рестрикции,
    # а затем выводим вставку целиком зеленым цветом.
    print(f"{vector_left_flank[:-6]}{YELLOW}{vector_left_flank[-6:]}{RESET} ✖ {GREEN}{insert_left_flank}{RESET}...")

    # печатаем маркер-стрелочку строго под сайтом рестрикции
    print(f"{' ' * 27}⬆ {enz_v1_name} ({enz_v1_site}) сайт")

    # граница, где конец вставки переходит обратно в вектор - правый стык

    # берем последние 40 нуклеотидов вставки
    insert_right_flank = donor_seq[-40:].upper()

    # берем 30 нуклеотидов вектора после сайта разреза v2.
    # Конкатенация срезов защищает от выхода за границы, если v2 находится в самом конце линейного представления кольца.
    vector_right_flank = (vector_seq[v2_idx:] + vector_seq[: v2_idx + 30])[:30].upper()

    print(f"Правый стык (Вставка -> Вектор):")

    # выводим вставку (зеленый), а затем вектор, у которого первые 6 нуклеотидов ([:6]) красятся в желтый (сайт),
    # а остальная часть вектора ([6:]) выводится стандартным цветом.
    print(
        f"...{GREEN}{insert_right_flank}{RESET} ✖ {YELLOW}{vector_right_flank[:6]}{RESET}{vector_right_flank[6:]}"
    )
    # печатаем маркер-стрелочку строго под вторым сайтом рестрикции
    print(f"{' ' * 43}⬆ {enz_v2_name} ({enz_v2_site}) сайт")

# интерактивный выбор пары ферментов через пробел
def select_enzyme_pair(exclude=None):

    # список наиболее часто используемых в лабораториях коммерческих рестриктаз
    enzymes = ["EcoRI", "HindIII", "BamHI", "XhoI", "SalI",
               "PstI", "KpnI", "NdeI", "NotI", "SacI",
               "XbaI", "BglII", "SmaI", "SpeI", "NheI"]

    exclude = set(exclude or [])
    available = [e for e in enzymes if e in Restriction.AllEnzymes and e not in exclude]

    print("\n Доступные ферменты:")
    for i, enzyme in enumerate(available, 1):
        print(f"{i:2}. {enzyme}")

    print("\n Введите номера или названия через пробел")

    while True:
        parts = input("Ваш выбор: ").split()

        if len(parts) != 2:
            print("Нужно ровно два фермента")
            continue

        selected = []

        for part in parts:
            if part.isdigit():
                idx = int(part) - 1
                if idx < 0 or idx >= len(available):
                    print(f"Номер {part} вне диапазона")
                    break

                selected.append(available[idx])
            else:
                if part not in Restriction.AllEnzymes:
                    print(f"Фермент '{part}' не найден")
                    break

                selected.append(part)
        else:
            if selected[0] == selected[1]:
                print("Ферменты должны быть разными")
                continue

            return selected[0], selected[1]


def select_enzymes_for_cloning():
    print(f"\n{'=' * 50}\n Выбор ферментов для клонирования \n{'=' * 50}")

    print("\n Ферменты для донора -")
    donor = select_enzyme_pair()
    if not donor:
        return None

    print("\n--- Ферменты для вектора ---")
    vector = select_enzyme_pair()
    if not vector:
        return None

    print(
        f"\n{'=' * 50}\n Выбраны ферменты :\n  Донор:  {donor[0]} + {donor[1]}\n  Вектор: {vector[0]} + {vector[1]}\n{'=' * 50}\n")

    return donor[0], donor[1], vector[0], vector[1]

def main():

    # окно для работы с диалоговыми окнами файлов
    root = Tk()
    root.withdraw()  # прячем основное его, чтобы оно не мозолило глаза
    root.attributes("-topmost", True)  # поверх всех окон, чтобы диалог выбора файла не прятался на задний план

    print("Выбор файлов плазмид ")
    # запрашиваем у пользователя путь к файлу плазмиды-донора

    donor_file = filedialog.askopenfilename(title="Выберите плазмиду-донор",filetypes=[("FASTA", "*.fasta *.fa *.fna")])

    # запрашиваем у пользователя путь к файлу векторной плазмиды
    vector_file = filedialog.askopenfilename(title="Выберите векторную плазмиду", filetypes=[("FASTA", "*.fasta *.fa *.fna")])

    # если пользователь закрыл одно из окон или нажал "Отмена" - выходим
    if not donor_file or not vector_file:
        print("Выбор файлов отменен.")

        return

    # формируем путь для сохранения результата в ту же папку, где лежит вектор
    output_file = os.path.join(os.path.dirname(vector_file), "recombinant_plasmid.fasta")

    # читаем FASTA-файлы в DataFrame
    donor_df = read_fasta_file(donor_file)
    vector_df = read_fasta_file(vector_file)

    # извлекаем последовательность из первой строки (индекс 0) колонки "sequence"
    donor_seq = donor_df.loc[0, "sequence"]
    vector_seq = vector_df.loc[0, "sequence"]

    print(f"\n Донор: {os.path.basename(donor_file)} - {len(donor_seq)} п.н.")
    print(f"Вектор: {os.path.basename(vector_file)} - {len(vector_seq)} п.н. \n")

    # вызываем интерактивное меню для выбора ферментов
    enzymes = select_enzymes_for_cloning()

    if enzymes is None:
        print("Выбор ферментов отменен.")

        return

    # распаковываем имена выбранных ферментов
    enz_d1_name, enz_d2_name, enz_v1_name, enz_v2_name = enzymes

    # получаем сами объекты ферментов рестрикции из базы данных Biopython Dictionary
    enz_d1 = Restriction.AllEnzymes.get(enz_d1_name)
    enz_d2 = Restriction.AllEnzymes.get(enz_d2_name)
    enz_v1 = Restriction.AllEnzymes.get(enz_v1_name)
    enz_v2 = Restriction.AllEnzymes.get(enz_v2_name)

    # запускаем основной процесс in silico клонирования
    cloning_result = perform_cloning(
        donor_seq, vector_seq, enz_d1, enz_d2, enz_v1, enz_v2
    )
    if cloning_result is None:
        print("Клонирование не удалось.")

        return

    # извлекаем результаты успешной сборки
    final_plasmid = cloning_result["final_plasmid"]
    len_backbone = cloning_result["backbone_length"]

    # выводим в консоль текстовый анализ получившихся стыков
    print_results(
        final_plasmid,
        vector_seq,
        donor_seq,
        enz_v1_name,
        enz_v2_name,
        enz_v1.site,
        enz_v2.site,
        cloning_result["pos_v1"],
        cloning_result["pos_v2"],
    )

    # показываем, что получилось - остов и вставленный фрагмент
    file_ready_seq = (final_plasmid[:len_backbone].upper()+ final_plasmid[len_backbone:].lower())

    # записываем финальную плазмиду в формате FASTA
    with open(output_file, "w") as f:

        # пишем заголовок с указанием итогового размера плазмиды
        f.write(f">recombinant_plasmid_size_{len(final_plasmid)}bp\n")

        # разбиваем длинную строку ДНК на стандартные строки по 60 символов для соответствия стандарту FASTA
        for i in range(0, len(file_ready_seq), 60):
            f.write(file_ready_seq[i : i + 60] + "\n")

    print(f"Результат успешно сохранен: {output_file}")

if __name__ == "__main__":
    main()