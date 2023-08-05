import os
import shutil
from tempfile import TemporaryDirectory
import subprocess
from multiprocessing import cpu_count
from Bio import SeqIO
from utilities.settings import settings
from workflow.components import minimap2
from utilities.logger import logger

def rvhaplo(input_fastq, reference, output_dir, prefix='rvhaplo', threads=cpu_count()):
    alignment_file = f'{output_dir}/alignment.sam'
    minimap2(input_fastq, reference, alignment_file)
    os.chdir('/opt/RVhaplo')
    process = subprocess.run(
        [
            './rvhaplo.sh',
            '--input', alignment_file,
            '--reference', reference,
            '--out', output_dir,
            '--prefix', prefix,
            '-t', str(threads)
        ], check=True)
    os.chdir('/hiv64148/scripts')
    return process.returncode

def haplodmf(input_fastq, reference, output_dir, prefix='haplodmf', threads=cpu_count()):
    logger.debug('Begin reconstruction with HaploDMF')
    with TemporaryDirectory() as _tmpdir:
        alignment_file = f'{_tmpdir}/alignment.sam'
        minimap2(input_fastq, reference, alignment_file)
        os.chdir('/opt/HaploDMF')
        process = subprocess.run(
            [
                'micromamba', 'run', '-n', 'haplodmf',
                './haplodmf.sh',
                '--input', alignment_file,
                '-r', reference,
                '--out', f'{_tmpdir}/{prefix}',
                '--prefix', prefix,
                '-t', str(threads)
            ], check=True)
        reformat_rvhaplo(f'{_tmpdir}/{prefix}/{prefix}_consensus.fasta', f'{output_dir}/haplotypes.final.fa')
    os.chdir('/hiv64148/scripts')
    return process.returncode

def reformat_rvhaplo(haplotype_fa, output):
    ori = SeqIO.parse(haplotype_fa, 'fasta')
    reformatted = []
    for seq in ori:
        seq_id_ori = seq.id.split('_')
        seq.id = f'{seq_id_ori[0]}_{seq_id_ori[1]} {int(float(seq_id_ori[11]))}x freq={round(float(seq_id_ori[5]), 3)}'
        reformatted.append(seq)
    SeqIO.write(reformatted, output, 'fasta')
    return

def goldrush(input_fastq, output_dir, threads=cpu_count()):
    base_dir = os.getcwd()
    with TemporaryDirectory() as _tmpdir:
        os.chdir(_tmpdir)
        os.symlink(input_fastq, f'{_tmpdir}/raw_reads.fastq')
        cmd = f'micromamba run -n goldrush goldrush run reads=raw_reads G=9600 m=1000 t={threads} z=6000 P=8'
        process = subprocess.run(cmd.split(), check=True)
        reformat_goldrush('w16_x10_golden_path.goldrush-edit-polished.span2.dist500.tigmint.fa.k40.w250.ntLink-5rounds.fa', f'{output_dir}/haplotypes.final.fa')
    os.chdir(base_dir)
    return process.returncode

def reformat_goldrush(haplotype_fa, output):
    ori = list(SeqIO.parse(haplotype_fa, 'fasta'))
    reformatted = []
    for i,seq in enumerate(ori):
        seq.id = f'contig_{i+1} N/Ax freq=N/A'
        reformatted.append(seq)
    SeqIO.write(reformatted, output, 'fasta')
    return
