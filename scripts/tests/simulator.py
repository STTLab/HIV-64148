
import os, tarfile
import requests
import subprocess
from ..utilities.logger import logger
from ..utilities.settings import settings
from tempfile import TemporaryDirectory

class Simulator(object):
    pre_trained_dir = './tests/pre-trained_models'
    def __init__(self) -> None:
        # Create temporary directory
        self._tmpdir = TemporaryDirectory()

    @classmethod
    def simulate_metagenome(cls, inputs: list, num_reads: int, model_prefix, perfect=True):
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
                '--output', f'{_temp}/simulated',
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
        try:
            shell = subprocess.run([*settings['softwares']['nanosim']['read_analysis'].split(), '--version'], check=True, stdout=subprocess.PIPE)
            logger.debug(f'Found `read_analysis.py` for {shell.stdout}')
        except subprocess.CalledProcessError as e:
            logger.error('`read_analysis.py` Error when called.')
            logger.error(e)
            return 1
        try:
            shell = subprocess.run([*settings['softwares']['nanosim']['simulator'].split(), '--version'], check=True, stdout=subprocess.PIPE)
            logger.debug(f'Found `simulator.py` for {shell.stdout}')
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
            exit(1)
        
        subprocess.run(['micromamba', 'create', '-n', 'nanosim', '-c', 'conda-forge', 'python=3.7', '-y'])
        logger.info('Created micromamba environment \'nanosim\' with Python 3.7')
        subprocess.run(['micromamba', 'install', '-n', 'nanosim', '-c', 'bioconda', '-c', 'conda-forge', 'nanosim=3.1.0', '-y'])
        logger.info('Installed nanosim in micromamba environment \'nanosim\'')

        with TemporaryDirectory() as _temp:
            logger.debug('Downloading pretrained model...')
            res = requests.get('https://github.com/bcgsc/NanoSim/raw/master/pre-trained_models/human_NA12878_DNA_FAB49712_guppy.tar.gz', allow_redirects=True)
            open(f'{_temp}/human_NA12878_DNA_FAB49712_guppy.tar.gz', 'wb').write(res.content)
            res = requests.get('https://github.com/bcgsc/NanoSim/raw/master/pre-trained_models/human_NA12878_DNA_FAB49712_albacore.tar.gz', allow_redirects=True)
            open(f'{_temp}/human_NA12878_DNA_FAB49712_albacore.tar.gz', 'wb').write(res.content)

            if not os.path.exists(cls.pre_trained_dir): os.makedirs(cls.pre_trained_dir)
            logger.info('Available pretrained model \'human_NA12878_DNA_FAB49712_guppy\'')
            tarfile.open(f'{_temp}/human_NA12878_DNA_FAB49712_guppy.tar.gz').extractall(f'{cls.pre_trained_dir}')
            logger.info('Available pretrained model \'human_NA12878_DNA_FAB49712_albacore\'')
            tarfile.open(f'{_temp}/human_NA12878_DNA_FAB49712_albacore.tar.gz').extractall(f'{cls.pre_trained_dir}')

        # Check installation
        if cls._check() != 0:
            logger.error('Installation error. You may want to manually install the tools and specify path to program in the settings')
            exit(0)

if __name__ == '__main__':
    sim = Simulator()
    sim._check()
