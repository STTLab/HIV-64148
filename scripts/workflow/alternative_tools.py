import os
import sys
import re
import shutil
from pathlib import Path
from tempfile import TemporaryDirectory
import subprocess
from multiprocessing import cpu_count
from Bio import SeqIO
from workflow.components import minimap2
from utilities.logger import logger

def rvhaplo(input_fastq, reference, output_dir, prefix='rvhaplo', asm_args=[]):
    with TemporaryDirectory() as _tmpdir:
        alignment_file = f'{_tmpdir}/alignment.sam'
        minimap2(input_fastq, reference, alignment_file)
        os.chdir('/opt/RVHaplo')
        process = subprocess.run(
            [
                'micromamba', 'run', '-n', 'haplodmf',
                'bash', './rvhaplo.sh',
                '--input', alignment_file,
                '-r', reference,
                '--out', f'{_tmpdir}/{prefix}',
                '--prefix', prefix,
                *asm_args
            ], check=True)
        reformat_rvhaplo(
            haplotype_fa=f'{_tmpdir}/{prefix}/{prefix}_consensus.fasta',
            output=f'{output_dir}/haplotypes.final.fa'
        )
    os.chdir(f'{Path(__file__).parent.absolute()}/..')
    return process.returncode

def haplodmf(input_fastq, reference, output_dir, prefix='haplodmf', asm_args=[]):
    logger.debug('Begin reconstruction with HaploDMF')
    with TemporaryDirectory() as _tmpdir:
        alignment_file = f'{_tmpdir}/alignment.sam'
        minimap2(input_fastq, reference, alignment_file)
        os.chdir('/opt/HaploDMF')
        process = subprocess.run(
            [
                'micromamba', 'run', '-n', 'haplodmf',
                'bash', './haplodmf.sh',
                '--input', alignment_file,
                '-r', reference,
                '--out', f'{_tmpdir}/{prefix}',
                '--prefix', prefix,
                *asm_args
            ], check=True)
        reformat_rvhaplo(
            haplotype_fa=f'{_tmpdir}/{prefix}/{prefix}_consensus.fasta',
            output=f'{output_dir}/haplotypes.final.fa'
        )
    os.chdir(f'{Path(__file__).parent.absolute()}/..')
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
    base_dir = '/hiv64148/scripts' # os.getcwd()
    with TemporaryDirectory() as _tmpdir:
        os.chdir(_tmpdir)
        shutil.copy(input_fastq, f'{_tmpdir}/raw_reads.fastq')
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

def flye(input_fastq, output_dir, *args, **kwargs):
    '''
    Run de novo haplotype reconstruction with MetaFlye
    '''
    # kwags processing
    min_overlap = kwargs.get('min_overlap', 1000)
    read_error = kwargs.get('read_error', 0.0)
    keep_haplotypes = kwargs.get('keep_haplotypes', True)
    metagenomic_mode = kwargs.get('metagenomic_mode', True)

    with TemporaryDirectory() as _tmpdir:
        raw_reads = f'{_tmpdir}/raw_reads.fastq'
        os.symlink(input_fastq, raw_reads)
        cmd = [
            'micromamba', 'run', '-n', 'venv',
            'flye', '--nano-hq', raw_reads,
            '--genome-size', '10k',
            '--min-overlap', str(min_overlap),
            '--read-error', str(read_error),
            '--out-dir', _tmpdir,
            '--scaffold'
        ]
        if metagenomic_mode:
            cmd.append('--meta')
        if keep_haplotypes:
            cmd.append('--keep-haplotypes')
        process = subprocess.run(cmd, check=True)

        final_assembly = f'{_tmpdir}/assembly.fasta'
        if os.path.exists(final_assembly):
            reformat_flye_output(final_assembly, f'{output_dir}/haplotypes.final.fa')
        else: sys.exit(1)
        # Move additional final output files
        to_move = (
            'assembly_graph.gfa',
            'assembly_graph.gv',
            'flye.log'
        )
        for file in to_move:
            shutil.move(f'{_tmpdir}/{file}', output_dir)
    return process.returncode

def reformat_flye_output(haplotype_fa, output):
    ori = list(SeqIO.parse(haplotype_fa, 'fasta'))
    reformatted = []
    for i,seq in enumerate(ori):
        seq.id = f'contig_{i+1} N/Ax freq=N/A'
        reformatted.append(seq)
    SeqIO.write(reformatted, output, 'fasta')
    return

def igda(input_fastq, reference, output_dir, threads=cpu_count()):
    context_model = '/hiv64148/data/igda_contextmodel/ont/ont_context_effect_read_qv_8_base_qv_8'
    with TemporaryDirectory() as _tmpdir:
        alignment_file = f'{_tmpdir}/alignment.sam'
        preprocessed_alignment_file = f'{_tmpdir}/cleaned_alignment.sam'
        preprocessed_alignment_bam = f'{_tmpdir}/cleaned_alignment.bam'
        minimap2(input_fastq, reference, alignment_file)
        subprocess.run(
            [
                'micromamba', 'run', '-n', 'igda',
                'igda_align_ont', alignment_file, reference,
                preprocessed_alignment_file, str(threads)
            ],
            check=True
        )
        subprocess.run(
            [
                'micromamba', 'run', '-n', 'igda',
                'sam2bam', preprocessed_alignment_file, str(threads)
            ],
            check=True
        )
        subprocess.run(
            [
                'micromamba', 'run', '-n', 'igda',
                'igda_pipe_detect', '-m', 'ont', '-n', str(threads),
                preprocessed_alignment_bam, reference, context_model,
                f'{_tmpdir}/minor_snvs'
            ],
            check=True
        )
        process = subprocess.run(
            [
                'micromamba', 'run', '-n', 'igda',
                'igda_pipe_phase', '-m', 'ont', '-n', str(threads),
                f'{_tmpdir}/minor_snvs', reference,
                f'{_tmpdir}/phased_snvs'
            ],
            check=True
            )
        reformat_goldrush(f'{_tmpdir}/phased_snvs/contigs.fa', f'{output_dir}/haplotypes.final.fa')
        return process.returncode

def canu(input_fastq, output_dir, prefix='canu', asm_args=[]):
    with TemporaryDirectory() as _tmpdir:
        raw_reads = f'{_tmpdir}/raw_reads.fastq'
        os.symlink(input_fastq, raw_reads)
        cmd = [
            'micromamba', 'run', '-n', 'canu',
            'canu', 'genomeSize=9.8k',
            '-p', prefix, '-d', _tmpdir,
            *asm_args, '-nanopore', raw_reads
        ]
        process = subprocess.run(cmd, check=True)
        final_assembly = f'{_tmpdir}/{prefix}.contigs.fasta'
        if os.path.exists(final_assembly):
            reformat_flye_output(final_assembly, f'{output_dir}/haplotypes.final.fa')
        else:
            raise FileNotFoundError('Assembly cannot be completed.')
        reformat_canu_output(final_assembly, f'{output_dir}/haplotypes.final.fa')
        to_move = (
            f'{prefix}.report',
            f'{prefix}.contigs.layout.readToTig',
            f'{prefix}.contigs.layout.tigInfo',
            f'{prefix}.unassembled.fasta'
        )
        for file in to_move:
            shutil.move(f'{_tmpdir}/{file}', output_dir)
        return process.returncode

def reformat_canu_output(haplotype_fa, output):
    ori = list(SeqIO.parse(haplotype_fa, 'fasta'))
    reformatted = []
    for i,seq in enumerate(ori):
        match = re.search(r'reads=(\d+)', seq.id)
        if match:
            coverage = int(match.group(1))
        else:
            coverage = 'N/A'
        seq.id = f'contig_{i+1} {coverage}x freq=N/A'
        reformatted.append(seq)
    SeqIO.write(reformatted, output, 'fasta')
    return
