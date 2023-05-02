
import os, tarfile, shutil, glob
import requests
import subprocess
import numpy as np
import pandas as pd
from pathlib import Path
from tempfile import TemporaryDirectory
from utilities.logger import logger
from utilities.settings import settings
from utilities.file_handler import FASTA

class Simulator(object):
    pre_trained_dir = './tests/pre-trained_models'
    def __init__(self) -> None:
        pass

    @classmethod
    def simulate_metagenome(cls, inputs: list, num_reads: int, model_prefix, run_dir, perfect=True):
        '''
        Simulate ONT reads for linear metagenome sample.
        '''
        # Check if inputs are valid
        if num_reads < 0: raise ValueError('Reads must be more than 0.')
        _abun_sum = 0
        for item in inputs:
            if len(item) != 3: raise IndexError()
            _abun_sum += item[2]
        if _abun_sum != 100: raise ValueError('Abundance must sum up to 100.')

        # Create metadata
        genome_list = open(f'{run_dir}/ref_list.tsv', 'w')
        abun_list = open(f'{run_dir}/abun_list.tsv', 'w')
        dna_type_list = open(f'{run_dir}/dna_type_list.tsv', 'w')
        abun_list.write(f'Size\t{num_reads}\n')

        for name, fasta_path, abun in inputs:
            genome_list.write(f'{name}\t{fasta_path}\n')
            abun_list.write(f'{name}\t{abun}\n')
            dna_type_list.write(f'{name}\t{name}\tlinear\n')
        genome_list.close()
        abun_list.close()
        dna_type_list.close()

        genome_list = f'{run_dir}/ref_list.tsv'
        abun_list = f'{run_dir}/abun_list.tsv'
        dna_type_list = f'{run_dir}/dna_type_list.tsv'

        cmd = [
            *settings['softwares']['nanosim']['simulator'].split(), 'metagenome',
            '--genome_list', genome_list,
            '--abun', abun_list,
            '--dna_type_list', dna_type_list,
            '--model_prefix', model_prefix,
            '--output', f'{run_dir}/simulated',
            '--max_len', '12000',
            '--min_len', '3500',
            '-t', '32',
            '--basecaller', 'guppy',
            '--fastq'
        ]
        if perfect:
            cmd.append('--perfect')
        process = subprocess.run(cmd)
        return process.returncode

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

class Simulation(object):
    def __init__(self):
        self.mode = ''
        self.dataset = ''

    @classmethod
    def test_random(cls, path_to_metadata, output_dir, path_to_fasta:str|None=None, prob=None, perfect=True):

        # Random simulation mode
        cls.mode = 'metagenome' # np.random.choice(('genome', 'metagenome'), 1)[0]

        if Path(path_to_metadata).suffix == '.tsv':
            delim = '\t'
        else:
            delim = ','
        df = pd.read_csv(path_to_metadata, sep=delim, )

        # Random number of genomes to include in the simulation
        if cls.mode == 'genome':
            n = 1
        else:
            n = np.random.randint(2,6)

        # Sampling from dataset
        choices = list(pd.unique(df['Subtype']))
        choosen = np.random.choice(choices, n, replace=False, p=prob)
        selected = []
        for item in choosen:
            selected.append(df.loc[df['Subtype'] == item].sample())
        selected = pd.concat(objs=selected).reset_index(drop=True)

        fasta = FASTA.read(path_to_fasta)

        # Start extracting data
        with TemporaryDirectory() as _temp:
            seqs_path = [ fasta.extract(id, save_to=f'{_temp}/{id}') for id in selected['ID'] ]

            while True:
                # Randomizing abundance value
                abun = np.random.dirichlet(np.ones(len(selected)),size=1)[0]
                abundances = [n*100 for n in abun]
                if sum(abundances) == 100:
                    break

            to_simulate = list(zip(selected['ID'].values, seqs_path, abundances))
            # Start simulation
            Simulator.simulate_metagenome(
                inputs=to_simulate,
                num_reads=20000,
                model_prefix=settings['data']['nanosim']['model'],
                run_dir=_temp, perfect=perfect
            )
            # Clean up
            if os.path.exists(output_dir):
                shutil.rmtree(output_dir)
            os.makedirs(output_dir)
            selected.to_csv(f'{output_dir}/pre-simulation-seq.tsv', sep='\t')
            keep = [f'{_temp}/abun_list.tsv', *glob.glob(f'{_temp}/*error_profile'), *glob.glob(f'{_temp}/*.fastq')]
            for file in keep:
                if os.path.exists(f'{output_dir}/{Path(file).name}'):
                    os.remove(f'{output_dir}/{Path(file).name}')
                shutil.move(file, f'{output_dir}')
            with open(f'{output_dir}/simulated_all_reads.fastq', 'w') as all_reads:
                for file in glob.glob(f'{output_dir}/*.fastq'):
                    reads = open(file, 'r').read()
                    all_reads.write(reads)

    def test_single_genome(self):
        pass

    def test_metagenome(self):
        pass
