import re
import os
import glob
import shutil
import tempfile
import subprocess
from pathlib import Path
from functools import lru_cache
from multiprocessing import cpu_count
import pandas as pd
from Bio import SeqIO

from utilities.apis import EutilsNCBI
from utilities.settings import settings
from utilities.logger import logger
from utilities.file_handler import FASTA
from utilities.benchmark_utils import run_command_with_logging

THREADS = str(cpu_count()-1)

def nanoplot_qc(input_file, input_type, output_dir, **filtering_options) -> int:
    '''
    Run NanoPlot, a plotting tool for long read sequencing data and alignments.

    Args:
        input_file (str): Path to the input file.
        input_type (str): Type of the input file. Supported types are:
        'fastq', 'fasta', 'fastq_rich', 'fastq_minimal',
        'summary', 'bam', 'ubam', 'cram'.
        output_dir (str): Directory to save the output files.
        **filtering_options: Additional filtering options for NanoPlot.
        Each option should be passed as a keyword argument.

    Returns:
        int: The return code of the NanoPlot process.

    Raises:
        RuntimeWarning: If the input_type is not one of the supported types.

    Note:
        This function requires NanoPlot to be installed and accessible via the command line.
        Please follow the instruction on how to install NanoPlot 
        from https://github.com/wdecoster/NanoPlot.

    Example:
        `return_code = nanoplot_qc('input.fastq', 'fastq', 'output_dir', min_length=1000, max_length=5000)`
    '''
    if input_type not in [
        'fastq', 'fasta', 'fastq_rich',
        'fastq_minimal', 'summary', 'bam',
        'ubam', 'cram'
    ]:
        raise RuntimeWarning('Format invalid input file\'s suffix may result in an error.')
    cmd = [
        *settings["softwares"]["NanoPlot"].split(),
        '--no_static', '--N50',
        '--prefix', f'{Path(input_file).stem}_',
        '--outdir', output_dir,
        '--threads', THREADS
    ]
    # If file is larger than 1Gib, declare as a huge input.
    if os.stat(input_file).st_size > (1024^3):
        cmd.append('--huge')
    if filtering_options:
        for key, val in zip(filtering_options.keys(), filtering_options.values()):
            cmd.append(f'--{key}')
            cmd.append(val)
    cmd.append(f'--{input_type}')
    cmd.append(input_file)

    process = subprocess.run(cmd, check=False)
    return process.returncode

def strainline(
        input_fastq,
        output_dir,
        clean_up = True,
        platform:str = 'ont',
        mintrimlen:int = 1000,
        topk:int = 50,
        minoverlap:int = 1000,
        miniden:float = 0.99,
        minseedlen:int = 3000,
        maxoh:int = 30,
        iterations:int = 2,
        maxgd:float = 0.01,
        maxld:float = 0.001,
        maxco:int = 5,
        min_abun:float = 0.02,
        rm_mis_asm: bool = False,
        err_cor:bool = True,
        threads:str|int = THREADS
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
        '--iter', str(iterations),
        '--maxGD', str(maxgd),
        '--maxLD', str(maxld),
        '--maxCO', str(maxco),
        '--minAbun', str(min_abun),
        '--rmMisassembly', str(rm_mis_asm),
        '--correctErr', str(err_cor),
        '--threads', str(threads)
    ]
    process, _ = run_command_with_logging(cmd, save_to=f'{output_dir}/logging/strainline_usage.csv')
    if clean_up:
        logger.info('Removing intermediate files.')
        to_remove = [
            'iter*',
            '*reads*',
            'corrected*.fa',
            'tmp',
            'haplotypes.fa',
            'haplotypes.final.fa',
            'haplotypes.rm_misassembly.fa',
            'contig_list.txt'
        ]
        for things in to_remove:
            for item in glob.glob(f'{output_dir}/{things}'):
                if os.path.isdir(item):
                    shutil.rmtree(item)
                else:
                    os.remove(item)
                logger.debug('%s...Removed', item)
        for file in ['haplotypes.final.fa', 'haps.bam', 'haps.depth']:
            shutil.move(f'{output_dir}/filter_by_abun/{file}', output_dir)
        shutil.rmtree(f'{output_dir}/filter_by_abun')
    return process.returncode

class BLAST:
    db_path = settings['data']['blast']['db_dir']
    def __init__(self, db_path=settings['data']['blast']['db_dir']) -> None:
        if not os.path.exists(db_path):
            os.makedirs(db_path)
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
            raise ValueError('Only "nucl" or "prot" is allowed. \
                            Plese refers to BLAST documentation.')
        if not os.path.exists(cls.db_path):
            os.makedirs(cls.db_path)
        cmd = [
            *settings['softwares']['blast']['makedb'].split(),
            '-title', dbtitle,
            '-dbtype', dbtype,
            '-in', input_file,
            '-parse_seqids',
            '-blastdb_version', '5',
        ]
        logger.debug('Creating database %s', dbtitle)
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        if stdout:
            logger.debug(stdout.decode())
        if stderr:
            logger.warning(stderr.decode())
        return process.returncode

    @classmethod
    def accession_to_db(cls, dbtitle: str, accession_list: list[str] | tuple[str], dbtype: str):
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
            raise ValueError('Only "nucl" or "prot" is allowed. \
                            Plese refers to BLAST documentation.')

        with tempfile.TemporaryDirectory() as _temp:
            seqs = EutilsNCBI.fetch_fasta_parallel(accession_list)
            SeqIO.write(seqs, f'{_temp}/{dbtitle}', 'fasta')
            if not os.path.exists(cls.db_path):
                os.makedirs(cls.db_path)
            shutil.move(f'{_temp}/{dbtitle}', f'{cls.db_path}/{dbtitle}')
        cls.create_db(dbtitle, f'{cls.db_path}/{dbtitle}', dbtype)

    @classmethod
    def blast_nucleotide(cls, query_fasta, dbtitle, output_dir, threads: str=THREADS) -> dict:
        '''
        Performs Blast analysis on provided sequences based on user-specified database.
        '''
        if not os.path.exists(cls.db_path):
            raise FileNotFoundError(f'Database {dbtitle} not found in {cls.db_path}.')
        output_file = f'{output_dir}/haplotype.blast.csv'
        cmd = [
            *(settings['softwares']['blast']['blastn'].split()),
            '-db', f'{cls.db_path}/{dbtitle}',
            '-query', query_fasta,
            '-out', output_file,
            '-outfmt', '6',
            '-num_threads', threads
        ]
        process, _ = run_command_with_logging(cmd, save_to=f'{output_dir}/logging/blast_usage.csv')
        return {
            'return_code': process.returncode,
            'output': {
                'output_dir': output_dir,
                'blast_result': output_file
            }
        }
    @classmethod
    @lru_cache(maxsize=100)
    def get_subtypes(cls, dbtitle, accession):
        iden_seq = FASTA.read_and_extract(f'{cls.db_path}/{dbtitle}', accession)
        subtype_regex = re.compile(
            r'CRF[0-9]{2}_[A-Z]{2}|CRF[0-9]{2}_[A-Z][0-9][A-Z]|subtype_[A,B,C,D,F,F1,F2,G]'
        )
        subtype = subtype_regex.findall(iden_seq.description)[0] # pyright: ignore reportGeneralTypeIssues=false
        return subtype

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

        def get_top(self, num):
            '''
            Get the top 'num' records for each group based on the 'bitscore' column.

            Args:
                num (int): Number of top records to retrieve for each group.

            Returns:
                pandas.DataFrame: Sorted DataFrame with the top 'n' records for each group.
                                The DataFrame is sorted based on the original index.

            Note:
                This method assumes that `self.data` is a DataFrame containing the data.

            Example:
                top_records = get_top(5)
            '''
            df_copy = self.data.copy()
            # sort by bitscore, grouped by `qseqid` and get top `n` records from each group.
            sorted_df = df_copy.sort_values('bitscore', ascending=False).groupby('qseqid').head(num)
            # Sort again by the original index in order to group the same group togeather.
            sorted_df = sorted_df.reset_index().sort_values('index')
            # Drop the original index.
            sorted_df = sorted_df.reset_index(drop=True).drop(['index'], axis=1)
            return sorted_df

        def get_iden(self, qseqid: str|None=None) -> pd.DataFrame:
            '''
            Retrieve data from the BLAST result.

            Args:
                qseqid (str | None): Optional. The 'qseqid' to filter the result.
                                    If provided, only rows with the specified 'qseqid' 
                                    will be included. Default is None.

            Returns:
                pd.DataFrame: DataFrame containing the rows with the highest 'bitscore' 
                value for each unique 'qseqid'. If 'qseqid' is provided, the result is
                filtered to include only rows with that 'qseqid'.

            Example:
                # Get all rows with the highest 'bitscore' value for each 'qseqid'
                
                `top_records = get_iden()`

                # Get rows with the highest 'bitscore' value for a specific 'qseqid'
                
                `specific_records = get_iden('my_qseqid')`
            '''
            df_copy = self.data.copy()
            res = df_copy.loc[df_copy.reset_index().groupby(['qseqid'])['bitscore'].idxmax()]
            if qseqid:
                return res.loc[res['qseqid'] == qseqid].reset_index(drop=True)
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
            with open(output, 'w', encoding='utf-8') as file:
                cmd = [settings['softwares']['minimap2'], '-ax', platform, reference, input_reads]
                process = subprocess.Popen(cmd, stdout=f, stderr=subprocess.PIPE)
        case _: raise ValueError(f'{fmt} is not a correct output format. (Only sam or bam is allowed)')
    return process.returncode

def snippy(
        input_file, input_reference, input_type,
        output_dir=os.getcwd(),
        snippy_params = None,
        threads:int|str = THREADS,
        max_ram:int = -1,
        tmp_dir:str = settings.get('tmp_dir', './tmp')
    ) -> dict:
    '''
    Finds SNPs between a haploid reference genome and your NGS sequence reads using Snippy.

    Parameters:
        - input_bam (str): Alignment file in bam format.
        - input_reference (str): Reference genome. Supports FASTA, GenBank, EMBL (not GFF)
        - output_dir (str): Output directory. (default=os.getcwd())
        - snippy_params (dict): Snippy parameter in Python dictionary, 
                                defaults will be use if not specified.
        - threads (int): Maximum number of CPUs to use. (default=multiprocessing.THREADS)
        - max_ram (int): Maximum RAM in Gb. (default=-1: AUTO)
        - tmp_dir (str): Directory to store temporary files. (default='./tmp')

    Returns:
        None
    '''
    if input_type not in ('bam', 'contigs'):
        raise RuntimeError('Incorrect input type. Only "bam" or "contig" is allowed')
    if snippy_params is None:
        snippy_params = {}
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
    }
    return output
