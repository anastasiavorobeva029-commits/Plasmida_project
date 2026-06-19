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


def get_cut_characteristics(sequence, enzyme, site_position):

    site_length = int(enzyme.size)
    if site_position + site_length > len(sequence):
        site_seq = sequence[site_position:] + sequence[:site_position + site_length - len(sequence)]
    else:
        site_seq = sequence[site_position:site_position + site_length]

    # Проверяем, липкий ли конец
    is_sticky = not enzyme.is_blunt()

    if is_sticky:
        # Получаем смещения разрезов
        # fst5 и fst3 показывают координаты разрезов относительно начала сайта
        c1 = enzyme.fst5
        c2 = enzyme.fst3

        # Нам нужна абсолютная разница между разрезами — это и есть длина липкого конца
        overhang_length = abs(c1 - c2)

        if enzyme.is_5overhang():
            # Для 5'-выступающих концов (как у EcoRI и XhoI)
            # Выступающая часть находится в начале или середине сайта
            start_cut = min(c1, c2) if min(c1, c2) >= 0 else 0
            left_end = site_seq[start_cut: start_cut + overhang_length]
            right_end = left_end  # В палиндромах выступающие концы одинаковы
        else:
            # Для 3'-выступающих концов
            start_cut = max(c1, c2) - overhang_length
            left_end = site_seq[start_cut: start_cut + overhang_length]
            right_end = left_end
    else:
        # Тупой конец
        left_end = Seq("")
        right_end = Seq("")

    return left_end, right_end, is_sticky


def get_enzyme_by_name(name):

    enzyme = Restriction.AllEnzymes.get(name)

    if enzyme is None:
        raise ValueError(f"Фермент '{name}' не найден в библиотеке!")

    return enzyme

def are_ends_compatible(end1, end2, sticky1, sticky2):

    # Если оба тупые - они всегда совместимы
    if not sticky1 and not sticky2:
        return True

    # Если один липкий, а другой тупой - склеить нельзя
    if sticky1 != sticky2:
        return False

    # Если оба липкие, полученные одинаковыми ферментами,
    # их текстовые последовательности хвостиков должны быть равны
    return str(end1) == str(end2)


def simulate_cloning(donor_seq, vector_seq, donor_enzyme1, donor_enzyme2,
                     vector_enzyme1, vector_enzyme2):

    d_enzyme1 = get_enzyme_by_name(donor_enzyme1)
    d_enzyme2 = get_enzyme_by_name(donor_enzyme2)
    v_enzyme1 = get_enzyme_by_name(vector_enzyme1)
    v_enzyme2 = get_enzyme_by_name(vector_enzyme2)

    donor_sites1 = find_site_in_circular_dna(donor_seq, d_enzyme1)
    donor_sites2 = find_site_in_circular_dna(donor_seq, d_enzyme2)

    if not donor_sites1:
        raise ValueError(f"Сайт для {donor_enzyme1} не найден в донорской плазмиде!")
    if not donor_sites2:
        raise ValueError(f"Сайт для {donor_enzyme2} не найден в донорской плазмиде!")
    if len(donor_sites1) > 1:
        raise ValueError(f"Найдено несколько сайтов ({len(donor_sites1)}) для {donor_enzyme1} в доноре!")
    if len(donor_sites2) > 1:
        raise ValueError(f"Найдено несколько сайтов ({len(donor_sites2)}) для {donor_enzyme2} в доноре!")
    if donor_sites1[0] == donor_sites2[0]:
        raise ValueError("Сайты для двух ферментов в доноре совпадают!")

    # 3. Находим сайты в векторе

    vector_sites1 = find_site_in_circular_dna(vector_seq, v_enzyme1)
    vector_sites2 = find_site_in_circular_dna(vector_seq, v_enzyme2)

    if not vector_sites1:
        raise ValueError(f"Сайт для {vector_enzyme1} не найден в векторной плазмиде!")
    if not vector_sites2:
        raise ValueError(f"Сайт для {vector_enzyme2} не найден в векторной плазмиде!")
    if len(vector_sites1) > 1:
        raise ValueError(f"Найдено несколько сайтов ({len(vector_sites1)}) для {vector_enzyme1} в векторе!")
    if len(vector_sites2) > 1:
        raise ValueError(f"Найдено несколько сайтов ({len(vector_sites2)}) для {vector_enzyme2} в векторе!")
    if vector_sites1[0] == vector_sites2[0]:
        raise ValueError("Сайты для двух ферментов в векторе совпадают!")

    # 4. Вырезаем вставку из донора

    insert, donor_remainder = cut_circle_dna(donor_seq, d_enzyme1, d_enzyme2)

    pos_donor1 = donor_sites1[0]
    pos_donor2 = donor_sites2[0]

    # Определяем, какой конец вставки левый, а какой правый
    # В зависимости от порядка сайтов в доноре
    if pos_donor1 < pos_donor2:
        # Первый сайт - левый конец вставки, второй - правый
        left_insert_end, _, left_insert_sticky = get_cut_characteristics(
            donor_seq, d_enzyme1, pos_donor1
        )
        right_insert_end, _, right_insert_sticky = get_cut_characteristics(
            donor_seq, d_enzyme2, pos_donor2
        )
    else:
        # Первый сайт - правый конец вставки, второй - левый
        left_insert_end, _, left_insert_sticky = get_cut_characteristics(
            donor_seq, d_enzyme2, pos_donor2
        )
        right_insert_end, _, right_insert_sticky = get_cut_characteristics(
            donor_seq, d_enzyme1, pos_donor1
        )

    # 7. Получаем характеристики концов для вектора
    pos_vector1 = vector_sites1[0]
    pos_vector2 = vector_sites2[0]

    # Определяем, какой конец вектора левый, а какой правый
    # В зависимости от порядка сайтов в векторе
    if pos_vector1 < pos_vector2:
        # Первый сайт - левый конец вектора, второй - правый
        left_vector_end, _, left_vector_sticky = get_cut_characteristics(
            vector_seq, v_enzyme1, pos_vector1
        )
        right_vector_end, _, right_vector_sticky = get_cut_characteristics(
            vector_seq, v_enzyme2, pos_vector2
        )
    else:
        # Первый сайт - правый конец вектора, второй - левый
        left_vector_end, _, left_vector_sticky = get_cut_characteristics(
            vector_seq, v_enzyme2, pos_vector2
        )
        right_vector_end, _, right_vector_sticky = get_cut_characteristics(
            vector_seq, v_enzyme1, pos_vector1
        )


    # Левый конец вставки должен быть совместим с левым концом вектора
    if not are_ends_compatible(left_insert_end, left_vector_end,
                               left_insert_sticky, left_vector_sticky):
        raise ValueError(f"Левые концы несовместимы: {left_insert_end} vs {left_vector_end}")

    # Правый конец вставки должен быть совместим с правым концом вектора
    if not are_ends_compatible(right_insert_end, right_vector_end,
                               right_insert_sticky, right_vector_sticky):
        raise ValueError(f"Правые концы несовместимы: {right_insert_end} vs {right_vector_end}")


    # Определяем части вектора
    # vector_remainder - это вектор без удаленного участка
    # Нам нужно разделить его на левую и правую части
    # Используем позиции сайтов в векторе

    # Находим позиции разрезов в векторе
    v_pos1 = vector_sites1[0]
    v_pos2 = vector_sites2[0]

    start_vec = min(v_pos1, v_pos2)
    end_vec = max(v_pos1, v_pos2)

    # Левая часть вектора (до удаленного участка)
    vector_left = vector_seq[:start_vec]
    # Правая часть вектора (после удаленного участка)
    vector_right = vector_seq[end_vec:]

    donor_ordered = pos_donor1 < pos_donor2
    vector_ordered = pos_vector1 < pos_vector2

    if donor_ordered != vector_ordered:
        insert = insert.reverse_complement()

    # Собираем новую плазмиду
    new_plasmid = vector_left + insert + vector_right

    return new_plasmid


# Донорная плазмида
donor = Seq('ATGAATTCTCGTCCCCATATATATCTCGAGCGGGGAGAGGACCCGGGAAGGAATTTCCTGTCGA')

# Векторная плазмида (для теста возьмем другую последовательность)
vector = Seq('ATGCGAATTCCCCGGGGGATATATCTCGAGCGGGGAGAGGACCCGGG')

print("Донор:", donor)
print("Вектор:", vector)
print()

result = simulate_cloning(
    donor_seq=donor,
    vector_seq=vector,
    donor_enzyme1='EcoRI',
    donor_enzyme2='XhoI',
    vector_enzyme1='EcoRI',
    vector_enzyme2='XhoI'
)


print(f"Новая плазмида: {result}")
print(f"Длина: {len(result)}")