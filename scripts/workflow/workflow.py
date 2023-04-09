
import os
import shutil
from uuid import uuid4, UUID
from tempfile import TemporaryDirectory
from . import components

class Worker(object):

    def __init__(self, handler=None) -> None:
        self.handler = handler

    def assign_job(self, input_fastq, output_dir, overwrite=False) -> UUID:
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

    @classmethod
    def assign_job_cli(self):
        '''
        Command line interface for creating a job.
        '''
        fq = input('Path to FASTQ: ')
        od = input('Output directory: ')
        if os.path.exists(od):
            if input('Output directory existed, overwrite? [y/N]: ').lower() != 'y': return
        return self.assign_job(self, fq, od, True)

    @classmethod
    def run_qc(self, input_fastq):
        pass

    def run_workflow(self):
        with TemporaryDirectory() as _tmpdir:
            os.symlink(self._input_fastq, f'{_tmpdir}/raw_reads.fastq')
            haplotype = components.strainline(
                input_fastq=f'{_tmpdir.name}/raw_reads.fastq',
                output_dir=_tmpdir
            )
            while haplotype.poll() == None:
                stdout, stderr = haplotype.communicate()
                print(stdout.decode())

            # if haplotype['return_code' == 0]:
            shutil.move(f'{_tmpdir}/haplotype.final.fa', self.output_dir)
            print('Finished')

def main():
    worker = Worker()
    job = worker.assign_job_cli()
    print(f'Job created (id:{job})')
    worker.run_workflow()

if __name__ == '__main__':
    main()
