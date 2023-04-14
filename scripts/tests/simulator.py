
import subprocess
from ..utilities.logger import logger
from ..utilities.settings import settings
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

    @classmethod
    def _check(cls) -> int:
        try: subprocess.run(['read_analysis.py', '--help'], check=True)
        except subprocess.CalledProcessError as e:
            logger.error('`read_analysis.py` Error when called.')
            logger.error(e)
            return 1
        try: subprocess.run(['simulator.py', '--help'], check=True)
        except subprocess.CalledProcessError as e:
            logger.error('`simulator.py` Error when called.')
            logger.error(e)
            return 1
        return 0
    
    @classmethod
    def install_nanosim(cls):
        '''
        Install NanoSim simulator in micromamba environment.
        Requires micromamba to function. If NanoSim read_analysis.py and simulator.py able to run
        and their return code are both 0; do nothing and return `None`.
        '''
        # If able to run; Do nothing...
        if cls._check() == 0: return

        try: subprocess.run(['micromamba', '--help'], check=True)
        except subprocess.CalledProcessError:
            logger.error('Micormamba not installed. Install Micromamba and try again.')
            exit()
        
        subprocess.run(['micromamba', 'create', '-n', 'nanosim', '-c', 'conda-forge', 'python=3.6'])
        logger.info('Created micromamba environment \'nanosim\' with Python 3.6')
        subprocess.run(['micromamba', 'install', '-n', 'nanosim', 'nanosim'])
        logger.info('Install nanosim in micromamba environment \'nanosim\'')

        # Check installation
        if cls._check() != 0:
            logger.error('Installation error. You may want to manually install the tools and specify path to program in the settings')
            exit()

if __name__ == '__main__':
    sim = Simulator()
    sim._check()
