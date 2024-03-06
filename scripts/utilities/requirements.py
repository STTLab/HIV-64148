'''
'''
import os
import sys
import subprocess
from .settings import settings
from .logger import logger


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

def check_softwares(executable) -> bool:
    cmd = executable + ' -h'
    with subprocess.Popen(
        cmd.split(), stdout=subprocess.DEVNULL, stderr=subprocess.PIPE
    ) as shell:
        _, error = shell.communicate()
    if error:
        logger.warning(error.decode())
        return False
    logger.debug('%s (exit: %s)', executable, shell.returncode)
    return True

def create_env(env_name, env_yaml):
    logger.debug('Create new environment %s', env_name)
    cmd = f'micromamba create -n {env_name} -y -f {env_yaml}'
    with subprocess.Popen(
        cmd.split(), stdout=subprocess.DEVNULL, stderr=subprocess.PIPE
    ) as shell:
        _, error = shell.communicate()
    if error:
        logger.warning(error.decode())

def clone_github_repository(repo_url, destination):
    '''
    This function takes the GitHub repository URL and the destination folder
    where you want to clone the repository as input arguments. 
    '''
    try:
        subprocess.run(['git', 'clone', repo_url, destination], check=True)
        logger.debug('Repository cloned successfully to %s', destination)
    except subprocess.CalledProcessError as e:
        logger.debug('Error cloning repository: %s', e)

def check_init_files():
    # Check BLAST DB
    db_path = settings['data']['blast']['db_dir']
    if os.path.exists(db_path): logger.info('Found BLAST database directory at %s', db_path)
    else: logger.warning('BLAST database not found at %s', db_path)

    # Check representative sequences
    rep_fasta = settings['data']['variant_calling']['rep_fasta']
    if os.path.exists(rep_fasta): logger.info('Found HIV-1 representative sequences at %s', rep_fasta)
    else: logger.warning('HIV-1 representative sequences not found at %s', rep_fasta)

    return (os.path.exists(db_path), os.path.exists(rep_fasta))
