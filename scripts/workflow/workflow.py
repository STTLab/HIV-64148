'''
This module controlls the workflow of HIV-64148 assembly pipeline
'''
__all__ = ['Worker',]
__version__ = '0.1'
__author__ = 'Sara Wattanasombat'

import os
import sys
import glob
import time
import shutil
import tracemalloc
from pathlib import Path
from uuid import uuid4, UUID
from datetime import timedelta
from tempfile import TemporaryDirectory
from Bio import SeqIO
from sierrapy.sierraclient import Sequence
from utilities.logger import logger
from utilities.settings import settings
from utilities.apis import hivdb_seq_analysis
from utilities.reporter import report_randerer, context_builder
from utilities.benchmark_utils import log_resource_usage
from utilities.requirements import check_softwares, create_env, clone_github_repository
from .alternative_tools import rvhaplo, goldrush, haplodmf, flye, igda, canu
from .components import BLAST, strainline, nanoplot_qc

class Worker(object):

    def __init__(self, assembler='strainline', asm_args=None) -> None:
        self.assembler = assembler
        self._stat = {
            't_created': time.time(),
            'peak_mem': {}
        }
        self.asm_args = asm_args
        self.make_report = True
        self.output_dir = ''
        self.job_id = uuid4()
        self._input_fastq = ''
        self.reference = None

    def set_reference(self, reference_path):
        self.reference = reference_path

    def get_peak_mem(self) -> dict:
        return self._stat.get('peak_mem', {})

    def get_runtime(self) -> int:
        '''
        Return runtime for this worker in seconds.
        '''
        return timedelta(seconds=time.time() - self._stat.get('t_start', time.time())).seconds

    def assign_job(self, input_fastq, output_dir, overwrite=False) -> UUID | None:
        '''
        Assign a job to a worker with provided parameters.
        '''
        self._input_fastq = input_fastq
        if os.path.exists(output_dir):
            if overwrite:
                shutil.rmtree(output_dir)
                os.makedirs(output_dir)
            else:
                logger.info('Directory existed at %s. Overwrite with --overwrite', output_dir)
                sys.exit(73)
        self.output_dir = output_dir
        return self.job_id

    def check_assembler(self):
        match self.assembler:
            case 'canu':
                if not check_softwares('micromamba run -n canu canu'):
                    create_env('canu', '/hiv64148/environments/canu.yaml')
            case 'rvhaplo':
                if not os.path.exists('/opt/RVHaplo/rvhaplo.sh'):
                    clone_github_repository('https://github.com/dhcai21/RVHaplo.git', '/opt/RVHaplo')
                if not check_softwares('micromamba run -n haplodmf /opt/RVHaplo/rvhaplo.sh'):
                    create_env('haplodmf', '/hiv64148/environments/haplodmf.yaml')
            case 'haplodmf':
                if not os.path.exists('/opt/HaploDMF/haplodmf.sh'):
                    clone_github_repository('https://github.com/dhcai21/HaploDMF.git', '/opt/HaploDMF')
                if not check_softwares('micromamba run -n haplodmf /opt/HaploDMF/haplodmf.sh'):
                    create_env('haplodmf', '/hiv64148/environments/haplodmf.yaml')

    @classmethod
    def run_qc(cls, input_file, output_dir):
        match Path(input_file).suffix:
            case '.fa' | '.fasta':
                input_type = 'fasta'
            case '.fastq':
                input_type = 'fastq'
            case _:
                return
        return nanoplot_qc(input_file, input_type, output_dir)

    @log_resource_usage(interval=0.1, output_file="resource_usage_log.csv")
    def run_workflow(self):
        self.check_assembler()
        tracemalloc.start()
        self._stat['t_start'] = time.time()

        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
        self.run_qc(self._input_fastq, f'{self.output_dir}/qc')

        # Haplotype assembly
        match self.assembler:
            case 'rvhaplo':
                # RVHaplo
                rvhaplo(
                    self._input_fastq,
                    self.reference,
                    self.output_dir,
                    self.asm_args
                )
                logger.info('Finished - RVHaplo')
            case 'haplodmf':
                haplodmf(
                    self._input_fastq,
                    self.reference,
                    self.output_dir,
                    self.asm_args
                )
                logger.info('Finished - HaploDMF')
            case 'goldrush':
                goldrush(
                    self._input_fastq,
                    self.output_dir,
                    self.asm_args
                )
            case 'metaflye':
                flye(
                    self._input_fastq,
                    self.output_dir,
                    self.asm_args
                )
                logger.info('Finished - MetaFlye')
            case 'igda':
                igda(
                    self._input_fastq,
                    self.reference,
                    self.output_dir,
                    self.asm_args
                )
                logger.info('Finished - iGDA')
            case 'canu':
                canu(
                    self._input_fastq,
                    self.output_dir,
                    self.asm_args
                )
                logger.info('Finished - Canu')
            case _:
                # Strainline
                with TemporaryDirectory() as _tmpdir:
                    os.symlink(self._input_fastq, f'{_tmpdir}/raw_reads.fastq')
                    strainline(
                        input_fastq=f'{_tmpdir}/raw_reads.fastq',
                        output_dir=_tmpdir,
                        asm_args = self.asm_args
                    )
                    logger.info('Moving Strainline output to %s', self.output_dir)
                    for file in glob.glob(f'{_tmpdir}/*'):
                        shutil.move(file, self.output_dir)
                    logger.info('Finished - Strainline')

        # Log memory
        _, _peak = tracemalloc.get_traced_memory()
        self._stat['peak_mem'][self.assembler] = round(_peak/(1024^2),3) # Memory Mib
        tracemalloc.reset_peak()

        # BLAST
        blast = BLAST.blast_nucleotide(
            f'{self.output_dir}/haplotypes.final.fa',
            settings['data']['blast']['dbtitle'],
            self.output_dir
        )
        logger.info('Finished - BLASTN')

        # Log memory
        _, _peak = tracemalloc.get_traced_memory()
        self._stat['peak_mem']['blast'] = round(_peak/(1024^2),3) # Memory Mib
        tracemalloc.reset_peak()

        if self.make_report:
            logger.info('Generating report')
            gql_schema_file = f'{Path(__file__).parent.absolute()}/../utilities/gql/sequence_analysis.gql'
            with open(gql_schema_file, 'r', encoding='utf-8') as gql_schema:
                gql = gql_schema.read()
            seqs = [
                        Sequence(header=record.id, sequence=str(record.seq)) \
                        for record in SeqIO.parse(f'{self.output_dir}/haplotypes.final.fa', 'fasta')
                    ]
            try:
                hivdb_result = hivdb_seq_analysis(seqs, gql)
            except ConnectionError:
                logger.warning('Connection error. Drug resistant profile will not be reported.')
                hivdb_result = None

            report_ctx = context_builder.context_builder(
                haplotype_fa=f'{self.output_dir}/haplotypes.final.fa',
                nanoplot_html=f'qc/{Path(self._input_fastq).stem}_NanoPlot-report.html',
                worker_info={
                    'job_id': self.job_id,
                    'run_stats':{
                        'assembler': self.assembler,
                        'runtime': self.get_runtime(),
                        'blast_result': blast['output']['blast_result'],
                        'peak_mem': self.get_peak_mem()
                    }
                },
                hivdb_result=hivdb_result
            )
            report_randerer.render_report(report_ctx, f'{self.output_dir}/hiv-64148_report.html')
        return None
