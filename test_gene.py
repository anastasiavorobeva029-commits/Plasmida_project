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

file_path = r'C:\Users\Nastya\PycharmProjects\plasmida_project\data\ncbi_dataset\data\pUC19_sequence.fasta'
df = read_fasta_file(file_path)

print(df)
print(f"Найдено {len(df)} последовательностей")






