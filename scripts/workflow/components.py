
import os
import glob
import shutil
import tempfile
import subprocess
import pandas as pd
from Bio import SeqIO

from utilities.apis import EutilsNCBI
from multiprocessing import cpu_count
from utilities.settings import settings
from utilities.logger import logger

THREADS = cpu_count()

def strainline(
        input_fastq,
        output_dir,
        platform:str = 'ont',
        mintrimlen:int = 1000,
        topk:int = 50,
        minoverlap:int = 1000,
        miniden:float = 0.99,
        minseedlen:int = 3000,
        maxoh:int = 30,
        iter:int = 2,
        maxgd:float = 0.01,
        maxld:float = 0.001,
        maxco:int = 5,
        min_abun:float = 0.02,
        rm_mis_asm: bool = False,
        err_cor:bool = True,
        threads:int = THREADS
    ) -> int:
    strainline_exe = settings['softwares']['strainline']
    cmd = [
        strainline_exe,
        '-i', input_fastq,
        '-o', output_dir,
        '-p', platform,
        '--minTrimmedLen', str(mintrimlen),
        '--topk', str(topk),
        '--minOvlpLen', str(minoverlap),
        '--minIdentity', str(miniden),
        '--minSeedLen', str(minseedlen),
        '--maxOH', str(maxoh),
        '--iter', str(iter),
        '--maxGD', str(maxgd),
        '--maxLD', str(maxld),
        '--maxCO', str(maxco),
        '--minAbun', str(min_abun),
        '--rmMisassembly', str(rm_mis_asm),
        '--correctErr', str(err_cor),
        '--threads', str(threads)
    ]
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return process.returncode

class BLAST:
    db_path = settings['data']['blast']['db_dir']
    def __init__(self, db_path=settings['data']['blast']['db_dir']) -> None:
        if not os.path.exists(db_path): os.makedirs(db_path)
        self.db_path = db_path

    @classmethod
    def create_db(cls, dbtitle, input_file, dbtype) -> int:
        '''
        Create Blast database in the same path with the input_file.

        Parameters:
            - dbtitle (str): Database title
            - input_file (str): Path to FASTA file to be use as reference for databse building
            - dbtype (str): Database type either "nulc" for nucleotide or "prot" for protein

        Returns:
            (int): Return code of a `makeblsetdb` subprocess
        '''
        if dbtype not in ('nucl', 'prot'):
            raise ValueError('Only "nucl" or "prot" is allowed. Plese refers to BLAST documentation.')
        if not os.path.exists(cls.db_path): os.makedirs(cls.db_path)
        cmd = [
            *settings['softwares']['blast']['makedb'].split(),
            '-title', dbtitle,
            '-dbtype', dbtype,
            '-in', input_file,
            '-parse_seqids',
            '-blastdb_version', '5',
        ]
        logger.debug(f'Creating database {dbtitle}')
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        if stdout:
            logger.debug(stdout.decode())
        if stderr:
            logger.warn(stderr.decode())
        return process.returncode
    
    @classmethod
    def accession_to_db(cls, dbtitle: str, accession_list: list[str] | tuple[str], dbtype: str) -> None:
        '''
        Create Blast database from a list of accession number(s).
        Fasta sequences will be fetched from NCBI database (requires internet connection).

        Parameters:
            - dbtitle (str): Database title
            - accession_list (list | tuple): List or Tuple of accession number in string
            - dbtype (str): Database type either "nulc" for nucleotide or "prot" for protein

        Returns:
            None
        '''
        if dbtype not in ('nucl', 'prot'):
            raise ValueError('Only "nucl" or "prot" is allowed. Plese refers to BLAST documentation.')

        with tempfile.TemporaryDirectory() as _temp:
            seqs = EutilsNCBI.fetch_fasta_parallel(accession_list)
            SeqIO.write(seqs, f'{_temp}/{dbtitle}', 'fasta')
            if not os.path.exists(cls.db_path): os.makedirs(cls.db_path)
            shutil.move(f'{_temp}/{dbtitle}', f'{cls.db_path}/{dbtitle}')
        cls.create_db(dbtitle, f'{cls.db_path}/{dbtitle}', dbtype)
        return

    @classmethod
    def blast_nucleotide(cls, query_fasta, dbtitle, output_dir, threads: int=THREADS) -> dict:
        '''
        Performs Blast analysis on provided sequences based on user-specified database.
        '''
        if not os.path.exists(cls.db_path):
            raise Exception(f'Database {dbtitle} not found in {cls.db_path}.')
        output_file = f'{output_dir}/haplotype.blast.csv'
        cmd =  settings['softwares']['blast'] + f'-db {dbtitle} -query {query_fasta} -out {output_file} -outfmt 18', '-num_threads', str(threads)
        cmd =  settings['softwares']['blast']['blastn'] + f' -db {cls.db_path}/{dbtitle} -query {query_fasta} -out {output_file} -outfmt 6 -num_threads {str(threads)}'
        ps1 = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = ps1.communicate()
        if stdout:
            logger.debug(stdout.decode())
        if stderr:
            logger.warn(stderr.decode())
            raise Exception('Blast error.')
        return {
            'return_code': ps1.returncode,
            'output': {
                'output_dir': output_dir,
                'blast_result': output_file
            }
        }
    
    class BLASTResult(object):
        def __init__(self, data: pd.DataFrame) -> None:
            self.data = data
        
        @classmethod
        def read(cls, file: str):
            '''
            Read file produced by BLASTN in tabular format (tab-seperated).
            Input must be produced with parameter `-outfmt 6` with default output fields.
            Plese refers to BLAST documentation.
            '''
            delim = '\t'
            col_names = [
                'qseqid', 'sseqid', 'piden','length',
                'mismatch', 'gapopen', 'qstart', 'qend',
                'sstart', 'send', 'evalue', 'bitscore'
            ]
            data = pd.read_csv(file, delimiter=delim, header=None, names=col_names)
            return BLAST.BLASTResult(data)
        
        def get_iden(self, qseqid: str|None=None) -> pd.DataFrame:
            '''
            Retrieve data from BLAST result.
            '''
            df = self.data
            res = df.loc[df.reset_index().groupby(['qseqid'])['bitscore'].idxmax()]
            if qseqid: return res.loc[res['qseqid'] == qseqid].reset_index(drop=True)
            return res.reset_index(drop=True)

def minimap2(input_reads, reference, output, platform: str = 'map-ont', fmt='bam') -> int:
    fmt = fmt.lower()
    match fmt:
        case 'bam':
            cmd = [settings['softwares']['minimap2'], '-ax', platform, reference, input_reads]
            alignment = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            process = subprocess.Popen(" ".join(["samtools", "sort", "-o", output]), stdin=alignment.stdout, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
            stdout, stderr = process.communicate()
        case 'sam':
            with open(output, 'w') as f:
                cmd = [settings['softwares']['minimap2'], '-ax', platform, reference, input_reads]
                process = subprocess.Popen(cmd, stdout=f, stderr=subprocess.PIPE)
        case _: raise ValueError(f'{fmt} is not a correct output format. (Only sam or bam is allowed)')
    return process.returncode

def snippy(
        input_file, input_reference, input_type,
        output_dir=os.getcwd(),
        snippy_params:dict = {},
        threads:int = THREADS,
        max_ram:int = -1,
        tmp_dir:str = settings.get('tmp_dir', './tmp')
    ) -> dict:
    '''
    Finds SNPs between a haploid reference genome and your NGS sequence reads using Snippy.

    Parameters:
        - input_bam (str): Alignment file in bam format.
        - input_reference (str): Reference genome. Supports FASTA, GenBank, EMBL (not GFF)
        - output_dir (str): Output directory. (default=os.getcwd())
        - snippy_params (dict): Snippy parameter in Python dictionary, defaults will be use if not specified.
        - threads (int): Maximum number of CPUs to use. (default=multiprocessing.THREADS)
        - max_ram (int): Maximum RAM in Gb. (default=-1: AUTO)
        - tmp_dir (str): Directory to store temporary files. (default='./tmp')

    Returns:
        None
    '''
    if input_type not in ('bam', 'contigs'):
        raise Exception('Incorrect input type. Only "bam" or "contig" is allowed')
    cmd = [
        *settings['softwares']['snippy'].split(),
        '--outdir', output_dir,
        '--ref', input_reference,
        '--mapqual', snippy_params.get('mapqual', '60'),
        '--basequal', snippy_params.get('basequal', '13'),
        '--mincov', snippy_params.get('mincov', '10'),
        '--minfrac', snippy_params.get('minfrac', '0'),
        '--minqual', snippy_params.get('minqual', '100'),
        '--maxsoft', snippy_params.get('maxsoft','10' ),
        '--cpus', str(threads),
        '--tmpdir', tmp_dir,
        '--quiet'
    ]
    if max_ram > 0:
        cmd.append('--ram')
        cmd.append(str(max_ram))

    if input_type == 'bam':
        cmd.extend(['--bam', input_file])
    elif input_type == 'contigs':
        cmd.extend(['--ctgs', input_file])

    try:
        ps = subprocess.run(cmd)
        output = {
            'return_code': ps.returncode,
            'output': {
                'output_dir': output_dir,
                'snps': f'{output_dir}/snps.vcf',
                'consensus': f'{output_dir}/snps.consensus.fa'
            }
        }
        return output
    except subprocess.CalledProcessError as e:
        return e.output
