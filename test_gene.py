import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Restriction
import os
from tkinter import filedialog, Tk


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


def find_site_in_circular_dna(sequence_str, enzyme, verbose=True):

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


def cut_circle_dna(sequence_str, pos):
    cut_1 = pos[0] - 1  # перевод в 0-based
    cut_2 = pos[1] - 1

    # Получаем два варианта разреза кольца
    if cut_1 < cut_2:
        frag_a = sequence_str[cut_1:cut_2]
        frag_b = sequence_str[cut_2:] + sequence_str[:cut_1]
    else:
        frag_a = sequence_str[cut_1:] + sequence_str[:cut_2]
        frag_b = sequence_str[cut_2:cut_1]

    # Для донора важен строгий порядок от enz_1 к enz_2
    # Но для вектора остов - это всегда большая часть кольца.
    # Поэтому мы вернем оба фрагмента, упорядочив их по длине для удобства.
    if len(frag_a) < len(frag_b):
        short_fragment = frag_a
        long_fragment = frag_b
    else:
        short_fragment = frag_b
        long_fragment = frag_a

    return short_fragment, long_fragment

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
        raise ValueError(f"Сайт {enzyme_name} встречается {len(sites)} раз в {plasmid_name}. Клонирование невозможно.")

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
    # Проверка на одинаковые ферменты
    if enz_d1_name == enz_d2_name:
        print("Для донора выбраны одинаковые ферменты")
        return None

    if enz_v1_name == enz_v2_name:
        print("Для вектора выбраны одинаковые ферменты")
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
    pos_d1 = validate_site(find_site_in_circular_dna(donor_seq, enz_d1), enz_d1_name, "Донор")
    pos_d2 = validate_site(find_site_in_circular_dna(donor_seq, enz_d2), enz_d2_name, "Донор")
    pos_v1 = validate_site(find_site_in_circular_dna(vector_seq, enz_v1), enz_v1_name, "Вектор")
    pos_v2 = validate_site(find_site_in_circular_dna(vector_seq, enz_v2), enz_v2_name, "Вектор")

    print(f"Донор координаты (1-based): {enz_d1_name}={pos_d1}, {enz_d2_name}={pos_d2}")
    print(f"Вектор координаты (1-based): {enz_v1_name}={pos_v1}, {enz_v2_name}={pos_v2}\n")

    # Проверка совместимости стыков
    comp_left = check_compatibility(enz_d1, enz_v1)
    comp_right = check_compatibility(enz_d2, enz_v2)

    if not (comp_left and comp_right):
        print("Концы несовместимы")
        if not comp_left:
            print(f"  - Конец после {enz_d1_name} не совместим с {enz_v1_name}")
        if not comp_right:
            print(f"  - Конец после {enz_d2_name} не совместим с {enz_v2_name}")
        return None

    print("Концы совместимы\n")

    # === 1. ВЫРЕЗАЕМ ВСТАВКУ (от d1 к d2 по часовой стрелке) ===
    # Переводим в 0-based индексы
    cut_d1 = pos_d1 - 1
    cut_d2 = pos_d2 - 1

    if cut_d1 < cut_d2:
        insert = donor_seq[cut_d1:cut_d2]
    else:
        insert = donor_seq[cut_d1:] + donor_seq[:cut_d2]

    # === 2. ВЫРЕЗАЕМ ОСТОВ ВЕКТОРА (ГАРАНТИРОВАННО БОЛЬШОЙ КУСОК) ===
    # Переводим в 0-based индексы
    v1_idx = pos_v1 - 1
    v2_idx = pos_v2 - 1

    # ПРАВИЛЬНАЯ ЛОГИКА ДЛЯ КОЛЬЦЕВОЙ ДНК:
    # Векторный остов - это ВСЁ кольцо, КРОМЕ маленького кусочка между сайтами
    if v1_idx < v2_idx:
        # Случай: v1 раньше v2 (например, EcoRI=400, HindIII=450)
        # Берем всё, ЧТО НЕ попадает в диапазон [v1_idx:v2_idx]
        vector_backbone = vector_seq[v2_idx:] + vector_seq[:v1_idx]
        # Для визуализации: левый стык - v1, правый стык - v2
        left_site_pos = pos_v1
        right_site_pos = pos_v2
    else:
        # Случай: v1 позже v2 (например, HindIII=448, XbaI=424)
        # Берем всё, ЧТО НЕ попадает в диапазон [v2_idx:v1_idx]
        vector_backbone = vector_seq[v1_idx:] + vector_seq[:v2_idx]
        # Для визуализации: левый стык - v1, правый стык - v2
        left_site_pos = pos_v1
        right_site_pos = pos_v2

    # === 3. ЛИГАЗНОЕ СШИВАНИЕ ===
    final_plasmid = vector_backbone + insert

    print(f"Вставка:   {len(insert):>5} п.н.")
    print(f"Остов:     {len(vector_backbone):>5} п.н.")
    print(f"Результат: {len(final_plasmid):>5} п.н.\n")

    # Возвращаем словарь со всеми данными для визуализации
    return {
        'final_plasmid': final_plasmid,
        'backbone_length': len(vector_backbone),
        'insert_length': len(insert),
        'pos_v1': pos_v1,
        'pos_v2': pos_v2,
        'pos_d1': pos_d1,
        'pos_d2': pos_d2,
        'enz_v1_name': enz_v1_name,
        'enz_v2_name': enz_v2_name,
        'enz_d1_name': enz_d1_name,
        'enz_d2_name': enz_d2_name
    }

def print_results(final_plasmid, vector_seq, donor_seq, start_insert,
                  enz_v1_name, enz_v2_name, pos_v1, pos_v2):

    if final_plasmid is None:
        return

    # Основная информация о длине
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

    # --- ИДЕАЛЬНЫЙ АНАЛИЗ ТОЧЕК СТЫКОВКИ ---
    print("\n--- АНАЛИЗ ТОЧЕК СТЫКОВКИ (JUNCTION SITES) ---")

    GREEN = "\033[92m"  # Зеленый для вставки
    YELLOW = "\033[93m"  # Желтый для сайтов рестрикции
    RESET = "\033[0m"  # Сброс цвета

    # Переводим в 0-based индексы
    v1_idx = pos_v1 - 1
    v2_idx = pos_v2 - 1

    print("Просмотр стыков ДНК крупным планом (5' -> 3'):\n")

    # === 1. ЛЕВЫЙ СТЫК (Конец остова вектора → начало вставки) ===
    # Берем фланкирующую последовательность из оригинального вектора вокруг сайта v1
    if v1_idx >= 30:
        vector_left_flank = vector_seq[v1_idx - 30:v1_idx].upper()
    else:
        # Если сайт близко к началу, берем с учетом кольцевости
        vector_left_flank = (vector_seq[v1_idx - 30:] + vector_seq[:v1_idx]).upper()

    # Берем начало вставки из донора (первые 40 нуклеотидов)
    insert_left_flank = donor_seq[:40].upper()

    # Подсвечиваем сайт рестрикции на векторе (последние 6 нуклеотидов перед стыком)
    if len(vector_left_flank) >= 6:
        vector_site = vector_left_flank[-6:]
        vector_prefix = vector_left_flank[:-6]
        print(f"1. Левый стык (Вектор -> Вставка):")
        print(f"   {vector_prefix}{YELLOW}{vector_site}{RESET} ✖ {GREEN}{insert_left_flank}{RESET}...")
        print(f"   {' ' * (len(vector_prefix) + 3)}⬆ {enz_v1_name} сайт")
    else:
        print(f"1. Левый стык (Вектор -> Вставка):")
        print(f"   {vector_left_flank} ✖ {GREEN}{insert_left_flank}{RESET}...")
        print(f"   {' ' * 30} ⬆ {enz_v1_name}\n")

    # === 2. ПРАВЫЙ СТЫК (Конец вставки → начало остова вектора) ===
    # Берем конец вставки из донора (последние 40 нуклеотидов)
    insert_right_flank = donor_seq[-40:].upper()

    # Берем фланкирующую последовательность из оригинального вектора вокруг сайта v2
    if v2_idx + 30 <= len(vector_seq):
        vector_right_flank = vector_seq[v2_idx:v2_idx + 30].upper()
    else:
        # Если сайт близко к концу, берем с учетом кольцевости
        vector_right_flank = (vector_seq[v2_idx:] + vector_seq[:v2_idx + 30 - len(vector_seq)]).upper()

    # Подсвечиваем сайт рестрикции на векторе (первые 6 нуклеотидов после стыка)
    if len(vector_right_flank) >= 6:
        vector_site = vector_right_flank[:6]
        vector_suffix = vector_right_flank[6:]
        print(f"2. Правый стык (Вставка -> Вектор):")
        print(f"...{GREEN}{insert_right_flank}{RESET} ✖ {YELLOW}{vector_site}{RESET}{vector_suffix}")
        print(f"   {' ' * (len(insert_right_flank) + 3)}⬆ {enz_v2_name} сайт")
    else:
        print(f"2. Правый стык (Вставка -> Вектор):")
        print(f"...{GREEN}{insert_right_flank}{RESET} ✖ {vector_right_flank}")
        print(f"   {' ' * 43} ⬆ {enz_v2_name}\n")

    print(f"(Обычным текстом показан остов вектора, {GREEN}Зеленым{RESET} — интегрированный фрагмент гена, "
          f"{YELLOW}Желтым{RESET} — сайты рестрикции, ✖ — место лигирования)")
    print("-----------------------------------------------------\n")


def save_plasmid(sequence, filepath):

    with open(filepath, "w") as f:
        f.write(f">recombinant_plasmid_{len(sequence)}bp\n")
        for i in range(0, len(sequence), 60):
            f.write(sequence[i : i + 60] + "\n")
    print(f"Результат успешно сохранен: {filepath}")


def select_enzyme(prompt):
    """
    Интерактивный выбор фермента из списка популярных или ввод своего

    Args:
        prompt: текст приглашения для пользователя

    Returns:
        название выбранного фермента
    """
    # Список популярных ферментов
    popular_enzymes = [
        "EcoRI", "HindIII", "BamHI", "XhoI", "SalI",
        "PstI", "KpnI", "NdeI", "NotI", "SacI",
        "XbaI", "BglII", "SmaI", "SpeI", "NheI"
    ]

    # Проверяем, какие ферменты доступны в БД
    available = []
    for enz in popular_enzymes:
        if enz in Restriction.AllEnzymes:
            available.append(enz)

    print("\n" + prompt)
    print("Доступные ферменты (введите номер или название):")
    print("-" * 40)

    # Выводим ферменты в две колонки для компактности
    half = len(available) // 2 + len(available) % 2
    for i in range(half):
        left_idx = i
        right_idx = i + half
        left = f"{left_idx + 1:2}. {available[left_idx]:<8}" if left_idx < len(available) else ""
        right = f"{right_idx + 1:2}. {available[right_idx]:<8}" if right_idx < len(available) else ""
        print(f"{left:<15} {right}")

    print("-" * 40)
    print("0. Ввести название фермента вручную")
    print("(или введите полное название фермента)")
    print("-" * 40)

    while True:
        choice = input("Ваш выбор: ").strip()

        if not choice:
            print("Ошибка: введите номер или название фермента.")
            continue

        # Проверяем, является ли ввод числом
        if choice.isdigit():
            num = int(choice)

            if num == 0:
                manual = input("Введите название фермента: ").strip()
                if not manual:
                    print("Ошибка: название не может быть пустым.")
                    continue
                if manual in Restriction.AllEnzymes:
                    return manual
                else:
                    print(f"Ошибка: фермент '{manual}' не найден в базе данных.")
                    continue

            if 1 <= num <= len(available):
                return available[num - 1]
            else:
                print(f"Ошибка: номер {num} вне диапазона (1-{len(available)}).")
                continue

        # Если ввели название фермента
        else:
            if choice in Restriction.AllEnzymes:
                return choice
            else:
                print(f"Ошибка: фермент '{choice}' не найден в базе данных.")
                continue


def select_enzymes_for_cloning():
    """
    Интерактивный выбор всех 4 ферментов для клонирования

    Returns:
        tuple: (enz_d1, enz_d2, enz_v1, enz_v2) или None при отмене
    """
    print("\n" + "=" * 50)
    print("ВЫБОР ФЕРМЕНТОВ ДЛЯ КЛОНИРОВАНИЯ")
    print("=" * 50)

    print("\n--- Ферменты для ДОНОРА ---")
    enz_d1 = select_enzyme("Выберите 1-й фермент для донора:")
    enz_d2 = select_enzyme("Выберите 2-й фермент для донора:")

    print("\n--- Ферменты для ВЕКТОРА ---")
    enz_v1 = select_enzyme("Выберите 1-й фермент для вектора:")
    enz_v2 = select_enzyme("Выберите 2-й фермент для вектора:")

    # Показываем выбранные ферменты
    print("\n" + "=" * 50)
    print("ВЫБРАНЫ ФЕРМЕНТЫ:")
    print(f"  Донор:  {enz_d1} + {enz_d2}")
    print(f"  Вектор: {enz_v1} + {enz_v2}")
    print("=" * 50 + "\n")

    # Подтверждение выбора
    confirm = input("Продолжить с этими ферментами? (y/N): ").strip().lower()
    if confirm != 'y':
        print("Отмена. Выход из программы.")
        return None

    return enz_d1, enz_d2, enz_v1, enz_v2


def main():
    # Создаем скрытое фоновое окно для работы файлового диалога
    root = Tk()
    root.withdraw()
    root.attributes(
        "-topmost", True
    )

    print("=== Выбор файлов плазмид ===")

    # 1. Интерактивный выбор файла донора
    print("Выберите FASTA-файл для плазмиды-донора:")
    donor_file = filedialog.askopenfilename(
        title="Выберите плазмиду-донор",
        filetypes=[("FASTA files", "*.fasta *.fa *.fna"), ("All files", "*.*")],
    )
    if not donor_file:
        print("Выбор файла отменен. Выход из программы.")
        return

    # 2. Интерактивный выбор файла вектора
    print("Выберите FASTA-файл для векторной плазмиды:")
    vector_file = filedialog.askopenfilename(
        title="Выберите векторную плазмиду",
        filetypes=[("FASTA files", "*.fasta *.fa *.fna"), ("All files", "*.*")],
    )
    if not vector_file:
        print("Выбор файла отменен. Выход из программы.")
        return

    # Фиксированное имя файла для автоматического обновления результата.
    vector_dir = os.path.dirname(vector_file)
    output_file = os.path.join(vector_dir, "recombinant_plasmid.fasta")

    # Чтение данных
    donor_df = read_fasta_file(donor_file)
    vector_df = read_fasta_file(vector_file)

    donor_seq = donor_df.loc[0, "sequence"]
    vector_seq = vector_df.loc[0, "sequence"]

    donor_name = os.path.splitext(os.path.basename(donor_file))[0]
    vector_name = os.path.splitext(os.path.basename(vector_file))[0]

    print("\n" + "=" * 40)
    print(f"Донор: {donor_name} - {donor_df.loc[0, 'length']} п.н.")
    print(f"Вектор: {vector_name} - {vector_df.loc[0, 'length']} п.н.")
    print("=" * 40 + "\n")

    # Выбор ферментов
    enzymes = select_enzymes_for_cloning()
    if enzymes is None:
        return

    enz_d1, enz_d2, enz_v1, enz_v2 = enzymes

    # Выполнение клонирования
    print("\n" + "=" * 60)
    cloning_result = perform_cloning(
        donor_seq, vector_seq, enz_d1, enz_d2, enz_v1, enz_v2
    )

    if cloning_result is not None:
        # Распаковываем словарь с результатами
        final_plasmid = cloning_result['final_plasmid']
        len_backbone = cloning_result['backbone_length']
        pos_v1 = cloning_result['pos_v1']
        pos_v2 = cloning_result['pos_v2']

        print("\n" + "=" * 60)
        # Вызываем обновленную функцию print_results с правильными параметрами
        print_results(
            final_plasmid,
            vector_seq,
            donor_seq,
            len_backbone,
            enz_v1,
            enz_v2,
            pos_v1,
            pos_v2
        )

        # Для файла: остов большими, вставка маленькими
        file_ready_seq = (
                final_plasmid[:len_backbone].upper()
                + final_plasmid[len_backbone:].lower()
        )

        with open(output_file, "w") as f:
            f.write(f">recombinant_plasmid_size_{len(final_plasmid)}bp\n")
            for i in range(0, len(file_ready_seq), 60):
                f.write(file_ready_seq[i: i + 60] + "\n")

        print(f"Результат успешно сохранен в файле: {output_file}")
        print("-----------------------------------------------------\n")
    else:
        print("\nКлонирование не удалось.")


if __name__ == "__main__":
    main()


