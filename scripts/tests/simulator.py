
import subprocess
from tempfile import TemporaryDirectory

class Simulator(object):
    def __init__(self) -> None:
        # Create temporary directory
        self._tmpdir = TemporaryDirectory()

    def simulate_metagenome(self, inputs: list, num_reads: int, model_prefix, perfect=True):
        '''
        Simulate ONT reads for linear metagenome sample.
        '''
        # Check if inputs are valid
        if num_reads < 0: raise ValueError('Reads must be more than 0.')
        _abun_sum = 0
        for item in inputs:
            if len(item) != 3: raise ValueError()
            _abun_sum += item[2]
        if _abun_sum != 100: raise ValueError('Abundance must sum up to 100.')

        # Create metadata
        genome_list = open(f'{self._tmpdir}/ref_list.tsv', 'w')
        abun_list = open(f'{self._tmpdir}/abun_list.tsv', 'w')
        dna_type_list = open(f'{self._tmpdir}/dna_type_list.tsv', 'w')
        abun_list.write(f'Size\t{num_reads}')

        for name, fasta_path, abun in inputs:
            genome_list.write(f'{name}\t{fasta_path}\n')
            abun_list.write(f'{name}\t{abun}\n')
            dna_type_list.write(f'{name}\tlinear\n')
        genome_list.close()
        abun_list.close()

        genome_list = f'{self._tmpdir}/ref_list.tsv'
        abun_list = f'{self._tmpdir}/abun_list.tsv'
        dna_type_list = f'{self._tmpdir}/dna_type_list.tsv'

        cmd = [
            'simulator.py', 'metagenome',
            '--genome_list', genome_list,
            '--abun', abun_list,
            '--dna_type_list', dna_type_list,
            '--model_prefix', model_prefix,
            '--output', f'{self._tmpdir}/simulated',
            '--max_len', '12000',
            '--min_len', '3500',
            '--basecaller', 'guppy'
            '--fastq'
        ]
        if perfect:
            cmd.append('--perfect')
        subprocess.Popen(cmd)

    def _check(self):
        process = subprocess.Popen(['read_analysis.py', '--help'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        if stdout:
            print(stdout.decode())
        process = subprocess.Popen(['simulator.py', '--help'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        if stdout:
            print(stdout.decode())

if __name__ == '__main__':
    sim = Simulator()
    sim._check()
