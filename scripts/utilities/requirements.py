
import glob
import os, sys
import requests
import subprocess
from pathlib import Path
from workflow.components import BLAST
from .settings import settings
from .logger import logger
from .apis import EutilsNCBI

def return_code(code: int):
    '''
    Convert return code to its text representation.
    '''
    if sys.version_info[1] < 10:
        return code

    match code:
        case 0: return 'OK'
        case 1: return 'ERROR'
        case 127: return 'MISSING'
        case _: return code

def check_softwares():
    result = {}
    softwares = settings.get('softwares', {})
    for prog in softwares.keys():
        if isinstance(softwares[prog], dict):
            for command in softwares[prog].keys():
                cmd = softwares[prog][command] + ' -h'
                shell = subprocess.Popen(
                    cmd.split(), stdout=subprocess.DEVNULL, stderr=subprocess.PIPE
                )
                output, error = shell.communicate()
                if error:
                    logger.warning(error.decode())
                else: logger.debug(f'{prog} {command} - {return_code(shell.returncode)}')
                result[f'{prog}_{command}'] = shell.returncode
            continue
        else:
            cmd = softwares[prog] + ' -h'
            shell = subprocess.Popen(
                cmd.split(), stdout=subprocess.DEVNULL, stderr=subprocess.PIPE
            )
            output, error = shell.communicate()
            if error:
                logger.warning(error.decode())
            else: logger.debug(f'{prog} - {return_code(shell.returncode)}')
            result[prog] = shell.returncode
    return result

def check_init_files():
    # Check BLAST DB
    db_path = settings['data']['blast']['db_dir']
    if os.path.exists(db_path): logger.info(f'Found BLAST database directory at {db_path}')
    else: logger.warning(f'BLAST database not found at {db_path}')

    # Check representative sequences
    rep_fasta = settings['data']['variant_calling']['rep_fasta']
    if os.path.exists(rep_fasta): logger.info(f'Found HIV-1 representative sequences at {rep_fasta}')
    else: logger.warning(f'HIV-1 representative sequences not found at {rep_fasta}')

    return (os.path.exists(db_path), os.path.exists(rep_fasta))

def check_requirements(repair: bool=False):
    '''
    Check for return code 0 of the required softwares by calling their help command.
    '''
    check_softwares()
    if not all(check_init_files()):
        logger.warning('Not all required files are present. Attempting to create...')
        if repair:
            setup_workflow()

def setup_workflow():
    for file in glob.glob('/opt/Strainline/src/*'):
        try:
            os.symlink(file, '/usr/local/bin', target_is_directory=True)
        except FileExistsError:
            pass
    rep_ids: dict = settings['data']['variant_calling']['rep_ids']
    rep_fasta = settings['data']['variant_calling']['rep_fasta']
    if not os.path.exists(rep_fasta):
        logger.info('Retrieving representative sequences for each HIV-1 subtype.')
        try: os.makedirs(Path(rep_fasta).parent)
        except FileExistsError: pass
    # for subtype, rep_ids in rep_ids.items():
    EutilsNCBI.fetch_fasta_parallel(rep_ids.values(), save_to=rep_fasta)
    logger.info('Downloading LosAlamos HIV-1 sequences.')
    response = requests.get('https://storage.googleapis.com/open-to-plublic/hiv_db_LosAlamos_FullGenome_utf8.fasta', allow_redirects=True)
    os.makedirs(settings['data']['blast']['db_dir'], exist_ok=True)
    with open(f"{settings['data']['blast']['db_dir']}/{settings['data']['blast']['dbtitle']}", 'wb') as file:
        file.write(response.content)
    logger.info('Creating database...')
    BLAST.create_db(settings['data']['blast']['dbtitle'], f"{settings['data']['blast']['db_dir']}/{settings['data']['blast']['dbtitle']}", 'nucl')
