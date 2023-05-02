
import os
import glob
import time
import shutil
import tracemalloc
from Bio import SeqIO
from pathlib import Path
from uuid import uuid4, UUID
from datetime import timedelta
from tempfile import TemporaryDirectory
from sierrapy.sierraclient import Sequence
from utilities.logger import logger
from utilities.settings import settings
from utilities.file_handler import FASTA
from .components import BLAST, strainline, snippy, nanoplot_qc
from utilities.apis import hivdb_seq_analysis
from utilities.reporter import report_randerer, context_builder

class Worker(object):

    def __init__(self, handler=None) -> None:
        self.handler = handler
        self._stat = {
            't_created': time.time(),
            'peak_mem': {}
        }

    def get_peak_mem(self) -> dict:
        return self._stat.get('peak_mem', {})

    def get_runtime(self):
        return timedelta(seconds=time.time() - self._stat.get('t_start', time.time()))

    def assign_job(self, input_fastq, output_dir, overwrite=False) -> UUID | None:
        '''
        Assign a job to a worker with provided parameters.
        '''
        self.job_id = uuid4()
        self._input_fastq = input_fastq
        if os.path.exists(output_dir):
            if overwrite:
                shutil.rmtree(output_dir)
                os.makedirs(output_dir)
        self.output_dir = output_dir
        return self.job_id

    def assign_job_cli(self):
        '''
        Command line interface for creating a job.
        '''
        fq = input('Path to FASTQ: ')
        od = input('Output directory: ')
        if os.path.exists(od):
            if input('Output directory existed, overwrite? [y/N]: ').lower() != 'y': return
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

    def run_workflow(self):
        tracemalloc.start()
        self._stat['t_start'] = time.time()

        if not os.path.exists(self.output_dir): os.makedirs(self.output_dir)
        self.run_qc(self._input_fastq, f'{self.output_dir}/qc')

        # Strainline
        with TemporaryDirectory() as _tmpdir:
            os.symlink(self._input_fastq, f'{_tmpdir}/raw_reads.fastq')
            strainline(
                input_fastq=f'{_tmpdir}/raw_reads.fastq',
                output_dir=_tmpdir
            )
            logger.info(f'Moving Strainline output to {self.output_dir}')
            [ shutil.move(file, self.output_dir) for file in glob.glob(f'{_tmpdir}/*') ]
            logger.info('Finished - Strainline')

        # Log memory
        _, _peak = tracemalloc.get_traced_memory()
        self._stat['peak_mem']['strainline'] = round(_peak/(1024^2),3) # Memory Mib
        tracemalloc.reset_peak()

        # BLAST
        blast = BLAST.blast_nucleotide(f'{self.output_dir}/haplotypes.final.fa', settings['data']['blast']['dbtitle'], self.output_dir)
        logger.info('Finished - BLASTN')

        # Log memory
        _, _peak = tracemalloc.get_traced_memory()
        self._stat['peak_mem']['blast'] = round(_peak/(1024^2),3) # Memory Mib
        tracemalloc.reset_peak()

        logger.info('Perform variant calling...')
        # Extract sequence according to BLAST result
        blast_result = BLAST.BLASTResult.read(blast['output']['blast_result']).get_iden()[['qseqid', 'sseqid']].to_numpy()
        haplotype_seqs = FASTA.read(f'{self.output_dir}/haplotypes.final.fa')
        rep_seqs = FASTA.read(settings['data']['variant_calling']['rep_fasta'])
        rep_ids = settings['data']['variant_calling']['rep_ids']
        os.mkdir(f'{self.output_dir}/variants')
        for hap, iden in blast_result:
            with TemporaryDirectory() as _temp:
                seq = haplotype_seqs.extract(hap)
                SeqIO.write(seq, f'{_temp}/{hap}.fasta', 'fasta')
                subtype = BLAST.get_subtypes(settings['data']['blast']['dbtitle'], iden)
                logger.info(f'Haplotype {hap} identified as {iden}, {subtype}')
                rep_seqs.extract(rep_ids[subtype], save_to=f'{_temp}/{subtype}.fasta')
                snippy(
                    input_file=f'{_temp}/{hap}.fasta',
                    input_reference=f'{_temp}/{subtype}.fasta',
                    input_type='contigs',
                    output_dir=f'{self.output_dir}/variants/{hap}_vs_{subtype}',
                    tmp_dir=_temp
                )
        logger.info('Finished - Variant calling')

        # Log memory
        _, _peak = tracemalloc.get_traced_memory()
        self._stat['peak_mem']['snippy'] = round(_peak/(1024^2),3) # Memory Mib
        tracemalloc.stop()

        gql = open('../utilities/gql/sequence_analysis.gql').read()
        seqs = [
            Sequence(header=record.id, sequence=str(record.seq)) for record in SeqIO.parse(f'{self.output_dir}/haplotypes.final.fa', 'fasta')
        ]
        hivdb_result = hivdb_seq_analysis(seqs, gql)

        logger.info('Generating report')
        report_ctx = context_builder.context_builder(
            haplotype_fa=f'{self.output_dir}/haplotypes.final.fa',
            nanoplot_html=f'qc/{Path(self._input_fastq).stem}_NanoPlot-report.html',
            worker_info={
                'job_id': self.job_id,
                'run_stats':{
                    'runtime': self.get_runtime(),
                    'blast_result': blast['output']['blast_result'],
                    'peak_mem': self.get_peak_mem()
                }
            },
            hivdb_result=hivdb_result
        )
        report_randerer.render_report(report_ctx, f'{self.output_dir}/hiv-64148_report.html')
