
import os
import glob
import shutil
from Bio import SeqIO
from uuid import uuid4, UUID
from tempfile import TemporaryDirectory
from utilities.logger import logger
from utilities.settings import settings
from utilities.file_handler import FASTA
from .components import BLAST, strainline, minimaap2

class Worker(object):

    def __init__(self, handler=None) -> None:
        self.handler = handler

    def assign_job(self, input_fastq, output_dir, overwrite=False) -> UUID | None:
        '''
        Assign a job to a worker with provided parameters.
        '''
        self.job_id = uuid4()
        self._input_fastq = input_fastq
        if os.path.exists(output_dir):
            if overwrite:
                os.rmdir(output_dir)
                os.makedirs(output_dir)
                return
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
    def run_qc(cls, input_fastq):
        pass

    def run_workflow(self):
        if not os.path.exists(self.output_dir): os.makedirs(self.output_dir)
        with TemporaryDirectory() as _tmpdir:
            os.symlink(self._input_fastq, f'{_tmpdir}/raw_reads.fastq')
            strainline(
                input_fastq=f'{_tmpdir}/raw_reads.fastq',
                output_dir=_tmpdir
            )
            logger.debug(f'Moving Strainline output to {self.output_dir}')
            [shutil.move(file, self.output_dir) for file in glob.glob(f'{_tmpdir}/*')]
            logger.info('Finished - Strainline')

        # BLAST
        blast = BLAST.blast_nucleotide(f'{self.output_dir}/haplotypes.final.fa', '32hiv1_default_db', self.output_dir)
        logger.info('Finished - BLASTN')

        # Extract sequence according to BLAST result
        blast_result = BLAST.BLASTResult.read(blast['output']['blast_result']).get_iden()[['qseqid', 'sseqid']].to_numpy()
        haplotype_seqs = FASTA.read(f'{self.output_dir}/haplotypes.final.fa')
        rep_seqs = FASTA.read(settings['data']['variant_calling']['rep_fasta'])
        for hap, iden in blast_result:
            with TemporaryDirectory() as _temp:
                seq = haplotype_seqs.extract(hap)
                SeqIO.write(seq, f'{_temp}/{hap}.fasta', 'fasta')
                seq.description

                match subtype:
                    case _:
                        pass
                minimaap2(
                    input_reads=f'{_temp}/{hap}.fasta',
                    reference='',
                    output=f'{_temp}/{hap}.bam'
                )
