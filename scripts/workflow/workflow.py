
import os, sys
import glob
import shutil
import argparse
from uuid import uuid4, UUID
from tempfile import TemporaryDirectory
from utilities.logger import logger 
from .components import BLAST, strainline
from utilities.requirements import setup_workflow

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
        blast_result = blast['output']['blast_result']
        

def main():
    parser = argparse.ArgumentParser(
                prog='HIV-64148 Pipeline',
                description='What the program does',
                epilog='Text at the bottom of help')
    parser.add_argument('function')

    args = parser.parse_args()
    match args.function:
        case 'run_cli':
            worker = Worker()
            job = worker.assign_job_cli()
            logger.info(f'Job created (id:{job})')
            worker.run_workflow()
        case 'setup':
            setup_workflow()
        case _: parser.print_help()


if __name__ == '__main__':
    main()
