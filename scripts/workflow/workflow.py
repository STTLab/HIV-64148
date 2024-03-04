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
from tests.alternative_tools import rvhaplo, goldrush, haplodmf, flye, igda, canu
from .components import BLAST, strainline, nanoplot_qc

class Worker(object):

    def __init__(self, handler=None, additionals=None, reconstructor='strainline') -> None:
        self.handler = handler
        self.reconstructor = reconstructor
        self._stat = {
            't_created': time.time(),
            'peak_mem': {}
        }
        self.additionals = additionals
        self.make_report = True
        self.output_dir = ''
        self.job_id = uuid4()
        self._input_fastq = ''
        self.reference = None

    def set_reference(self, reference_path):
        self.reference = reference_path

    def get_peak_mem(self) -> dict:
        return self._stat.get('peak_mem', {})

    def get_runtime(self):
        return timedelta(seconds=time.time() - self._stat.get('t_start', time.time()))

    def assign_job(self, input_fastq, output_dir, overwrite=False, ) -> UUID | None:
        '''
        Assign a job to a worker with provided parameters.
        '''
        self._input_fastq = input_fastq
        if os.path.exists(output_dir):
            if overwrite:
                shutil.rmtree(output_dir)
                os.makedirs(output_dir)
                os.makedirs(f'{output_dir}/logging')
            else:
                logger.info('Directory existed at %s. Please allow overwrite with --overwrite to proceed.', output_dir)
                sys.exit(73)
        self.output_dir = output_dir
        return self.job_id

    def assign_job_cli(self):
        '''
        Command line interface for creating a job.
        '''
        fq = input('Path to FASTQ: ')
        od = input('Output directory: ')
        if os.path.exists(od):
            if input('Output directory existed, overwrite? [y/N]: ').lower() != 'y':
                return
            shutil.rmtree(od)
        return self.assign_job(fq, od, True)

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
        tracemalloc.start()
        self._stat['t_start'] = time.time()

        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
        self.run_qc(self._input_fastq, f'{self.output_dir}/qc')

        # Haplotype reconstruction
        match self.reconstructor:
            case 'rvhaplo':
                # RVHaplo
                rvhaplo(
                    self._input_fastq,
                    self.reference,
                    self.output_dir,
                )
                logger.info('Finished - RVHaplo')
            case 'haplodmf':
                haplodmf(
                    self._input_fastq,
                    self.reference,
                    self.output_dir,
                )
                logger.info('Finished - HaploDMF')
            case 'goldrush':
                goldrush(
                    self._input_fastq,
                    self.output_dir
                )
            case 'metaflye':
                flye(
                    self._input_fastq,
                    self.output_dir
                )
                logger.info('Finished - MetaFlye')
            case 'igda':
                igda(
                    self._input_fastq,
                    self.reference,
                    self.output_dir
                )
                logger.info('Finished - iGDA')
            case 'canu':
                canu(
                    self._input_fastq,
                    self.output_dir
                )
                logger.info('Finished - Canu')
            case _:
                # Strainline
                with TemporaryDirectory() as _tmpdir:
                    os.symlink(self._input_fastq, f'{_tmpdir}/raw_reads.fastq')
                    strainline(
                        input_fastq=f'{_tmpdir}/raw_reads.fastq',
                        output_dir=_tmpdir,
                    )
                    logger.info('Moving Strainline output to %s', self.output_dir)
                    for file in glob.glob(f'{_tmpdir}/*'):
                        shutil.move(file, self.output_dir)
                    logger.info('Finished - Strainline')

        # Log memory
        _, _peak = tracemalloc.get_traced_memory()
        self._stat['peak_mem'][self.reconstructor] = round(_peak/(1024^2),3) # Memory Mib
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
            with open('/hiv64148/scripts/utilities/gql/sequence_analysis.gql', 'r', encoding='utf-8') as gql_schema:
                gql = gql_schema.read()
            seqs = [ Sequence(header=record.id, sequence=str(record.seq)) for record in SeqIO.parse(f'{self.output_dir}/haplotypes.final.fa', 'fasta') ]
            try:
                hivdb_result = hivdb_seq_analysis(seqs, gql)
            except ConnectionError:
                hivdb_result = None

            logger.info('Generating report')
            report_ctx = context_builder.context_builder(
                haplotype_fa=f'{self.output_dir}/haplotypes.final.fa',
                nanoplot_html=f'qc/{Path(self._input_fastq).stem}_NanoPlot-report.html',
                worker_info={
                    'job_id': self.job_id,
                    'run_stats':{
                        'reconstructor': self.reconstructor,
                        'runtime': self.get_runtime(),
                        'blast_result': blast['output']['blast_result'],
                        'peak_mem': self.get_peak_mem()
                    }
                },
                hivdb_result=hivdb_result
            )
            report_randerer.render_report(report_ctx, f'{self.output_dir}/hiv-64148_report.html')
        return None
