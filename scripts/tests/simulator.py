
import subprocess
from tempfile import TemporaryDirectory

class simulator(object):
    def __init__(self) -> None:
        pass

    def simulate_metagenome(self, inputs: list, num_reads: int, model_prefix):
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

        # Create temporary directory
        with TemporaryDirectory() as _temp:
            # Create metadata
            genome_list = open(f'{_temp}/ref_list.tsv', 'w')
            abun_list = open(f'{_temp}/abun_list.tsv', 'w')
            dna_type_list = open(f'{_temp}/dna_type_list.tsv', 'w')
            abun_list.write(f'Size\t{num_reads}')

            for name, fasta_path, abun in inputs:
                genome_list.write(f'{name}\t{fasta_path}\n')
                abun_list.write(f'{name}\t{abun}\n')
                dna_type_list.write(f'{name}\tlinear\n')
            genome_list.close()
            abun_list.close()

            genome_list = f'{_temp}/ref_list.tsv'
            abun_list = f'{_temp}/abun_list.tsv'
            dna_type_list = f'{_temp}/dna_type_list.tsv'

            cmd = [
                'simulator.py', 'metagenome',
                '--genome_list', genome_list,
                '--abun', abun_list,
                '--dna_type_list', dna_type_list,
                '--model_prefix', model_prefix,
                '--output', f'{_temp}/simulated'
            ]
            subprocess.Popen(cmd)
    