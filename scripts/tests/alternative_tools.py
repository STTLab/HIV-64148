import os
import subprocess
from Bio import SeqIO
from multiprocessing import cpu_count
from utilities.settings import settings
from workflow.components import minimap2

def rvhaplo(input_fastq, reference, output_dir, prefix='rvhaplo', threads=cpu_count()):
    alignment_file = f'{output_dir}/alignment.sam'
    minimap2(input_fastq, reference, alignment_file, fmt='bam')
    os.chdir('/opt/RVhaplo')
    ps = subprocess.run(
        [
            './rvhaplo.sh',
            '--input', alignment_file,
            '--reference', reference,
            '--out', output_dir,
            '--prefix', prefix,
            '-t', str(threads)
        ])
    os.chdir('/hiv64148/scripts')
    return ps.returncode

def reformat_rvhaplo(haplotype_fa, output):
    ori = SeqIO.parse(haplotype_fa, 'fasta')
    reformatted = []
    for seq in ori:
        seq_id_ori = seq.id.split('_')
        seq.id = f'{seq_id_ori[0]}_{seq_id_ori[1]} {int(float(seq_id_ori[11]))}x freq={round(float(seq_id_ori[5]), 3)}'
        reformatted.append(seq)
    SeqIO.write(reformatted, output, 'fasta')
    return
