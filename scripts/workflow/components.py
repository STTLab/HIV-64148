
import os
import shutil
import tempfile
import subprocess
from multiprocessing import cpu_count
from utilities.settings import settings

THREADS = cpu_count()

def reads_alignment(input_reads, reference, output, platform: str = 'map-ont', fmt='bam') -> dict:
    fmt = fmt.lower()
    if fmt == 'bam':
        cmd = ['~/Tools/3-assembler/minimap2-2.24_x64-linux/minimap2', '-ax', platform, reference, input_reads]
        alignment = subprocess.Popen(" ".join(cmd), shell=True, stdout=subprocess.PIPE)
        ps = subprocess.run(" ".join(["samtools", "sort", "-o", output]), stdin=alignment.stdout, shell=True, check=True)
    elif fmt == 'sam':
        with open(output, 'w') as f:
            cmd = ['~/Tools/3-assembler/minimap2-2.24_x64-linux/minimap2', '-ax', platform, reference, input_reads]
            ps = subprocess.run(" ".join(cmd), shell=True, check=True, stdout=f)
    else:
        raise ValueError(f'{fmt} is not a correct output format. (Only sam or bam is allowed)')
    return {
        'retrn_code': ps.returncode,
        'output': output
    }

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
    ) -> dict:
    strainline_exe = '/opt/Strainline/src/strainline.sh'
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
    ps = subprocess.run(cmd, shell=True, check=True)
    return {
        'return_code': ps.returncode,
        'output': {
            'corrected_reads': f'{output_dir}/corrected.{iter}.fa',
            'haplotype': f'{output_dir}/haplotype.final.fa'
        }
    }

def rvhaplo(input_reads, reference, prefix, output):
    rvhaplo_script = settings["softwares"]["rvhaplo"]
    os.chdir(settings["softwares"]["rvhaplo"])
    temp1, temp_name1 = tempfile.mkstemp(suffix='.sam')
    reads_alignment(input_reads, reference, temp_name1, fmt='sam')
    with tempfile.TemporaryDirectory() as tmp_dir:
        cmd = f'{rvhaplo_script} -i {temp_name1} -r {reference} -o {tmp_dir} -p {prefix} -t {THREADS}'
        print(cmd)
        ps = subprocess.run(cmd, shell=True)
        shutil.move(f'{tmp_dir}/{prefix}_haplotypes.fasta', output)
    os.close(temp1)
    return {
        'return_code': ps.returncode,
        'output': {
            'haplotype': f'{output}'
        }
    }


class BLAST:

    def __init__(self) -> None:
        pass

    def blast_db_exists(self) -> bool:
        pass

    def create_blast_db(self, input_file):
        '''
        Create Blast database if not existed.
        '''
        if self.blast_db_exists():
            return
        return subprocess.run(['makeblastdb', '-dbtype', 'nucl', '-in', input_file, '-parse_seqids'])

    def blast_sequences(self, query_fasta, database, output_dir, threads: int=THREADS) -> dict:
        '''
        Performs Blast analysis on provided sequences based on user-specified database.
        '''
        if self.blast_db_exists:
            raise Exception(f'Daatabase {database} not found.')
        output_file = f'{output_dir}/haplotype.blast.csv'
        cmd =  settings['softwares']['blast'] + f'-db {database} -query {query_fasta} -out {output_file} -outfmt 18', '-num_threads', str(threads)
        ps1 = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        if os.path.exists(output_file):
            ps2 = subprocess.run(
                ['sed', '-i',
                "'1 i\qseqid,sseqid,pident,length,mismatch,gapopen,\
                    qstart,qend,sstart,send,evalue,bitscore'", output_file
                ],shell=True, check=True
            )
            return {
                'return_code': [ps1.returncode, ps2.returncode],
                'output': {
                    'output_dir': output_dir,
                    'blast_result': output_file
                }
            }
        raise Exception('Blast error.')

def snippy(
        input_file, input_reference, input_type,
        output_dir=os.getcwd(),
        snippy_params:dict = {},
        threads:int = THREADS,
        max_ram:int = -1,
        tmp_dir:str = settings.get('tmp_dir', './tmp')
    ) -> None:
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
        settings['softwares']['snippy'],
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
        ps = subprocess.Popen(cmd, stdout=subprocess.PIPE)
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
