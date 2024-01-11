import glob
import numpy as np
import pandas as pd
from Bio import SeqIO
import concurrent.futures
import matplotlib.pyplot as plt
from tempfile import TemporaryDirectory
from matplotlib.ticker import FuncFormatter
from multiprocessing import cpu_count, Pool
from utilities.apis import EutilsNCBI
from utilities.file_handler import FASTA

def dot_matrix(seq1, seq2, window_size=3, max_workers=cpu_count()):
    '''
    Perform sequence comparison and create a matrix
    '''
    # Define the dot plot matrix
    print(f'{seq1.id} vs {seq2.id}')
    seq1 = str(seq1.seq)
    seq2 = str(seq2.seq)
    matrix = np.zeros((len(seq1), len(seq2)), dtype=int)
    # matrix = np.zeros((max(list(map(len, [seq.seq for seq in haplotype_seqs]))), max(list(map(len, [seq.seq for seq in truth_seqs])))), dtype=int)

    # Define a function to calculate dot plot values for a given range of positions
    def calculate_matrix_values(start_i, end_i):
        for i in range(start_i, end_i):
            for j in range(len(seq2) - window_size + 1):
                match = True
                for k in range(window_size):
                    if seq1[i+k] != seq2[j+k]:
                        match = False
                        break
                if match:
                    for k in range(window_size):
                        matrix[i+k,j+k] = 1

    # Create a thread pool and distribute the work across the threads
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        chunk_size = len(seq1) // executor._max_workers
        futures = []
        for i in range(0, len(seq1), chunk_size):
            futures.append(executor.submit(calculate_matrix_values, i, i+chunk_size))

        # Wait for all threads to complete
        for future in concurrent.futures.as_completed(futures):
            pass
    return matrix

def kilo(x, pos):
    '''
    Convert number to a string with kilo (k) suffix
    '''
    if x < 1000:
        return x
    return '{:.0f}k'.format(x/1000)

def draw_plot(matrices, haplotype_seqs, truth_seqs):
    formatter = FuncFormatter(kilo)

    px = 1/plt.rcParams['figure.dpi']
    fig, ax = plt.subplots(len(haplotype_seqs), len(truth_seqs)+1, figsize=(1920*px, 1080*px), squeeze=False)
    i = 0
    j = 0
    for matrix in matrices:
        this_ax = ax[j][i]
        this_ax.yaxis.set_major_formatter(formatter)
        this_ax.xaxis.set_major_formatter(formatter)
        this_ax.imshow(matrix, cmap='gray', interpolation='nearest')
        this_ax.invert_yaxis()

        this_ax.set_ylabel(haplotype_seqs[j].id, loc='bottom')
        try:
            this_ax.set_xlabel(truth_seqs[i].id, loc='left')
        except IndexError:
            this_ax.set_xlabel(haplotype_seqs[j].id, loc='left')
        i += 1
        if i == len(truth_seqs)+1:
            i = 0
            j += 1
    plt.tight_layout()
    plt.savefig('my_plot.png')
    plt.close()

def main():
    truth_fasta_pool = FASTA.read('/hiv64148/data/blast/db/LosAlamos_db')
    for i, sim in enumerate(glob.glob('/workspace/outputs/test_multi_simulation/simulation_*')):
        matrices = []
        haplotype_fasta = f'{sim}/pipeline_output/filter_by_abun/haplotypes.final.fa'

        # Create truth fasta
        df = pd.read_csv('/workspace/')
        truth_data = df.loc[df['simulation_no'] == f'simulation_{i+1}']['Accession']
        truth_seqs = [ truth_fasta_pool.extract(accession) for accession, _subtype in truth_data ]
        haplotype_seqs = list(SeqIO.parse(haplotype_fasta, 'fasta'))

        window_size = 3
        for haplotype in haplotype_seqs:
            for truth in truth_seqs:
                matrix = dot_matrix(haplotype, truth)
                matrices.append(matrix)
            matrix = dot_matrix(haplotype, haplotype)
            matrices.append(matrix)
