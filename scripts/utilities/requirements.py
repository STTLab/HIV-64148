
import sys
import subprocess
from .settings import settings
from .logger import logger
from ..workflow.components import BLAST

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

def check_requirements() -> dict:
    '''
    Check for return code 0 of the required softwares by calling their help command.
    '''
    result = {}
    softwares = settings.get('softwares', {})
    for prog in softwares.keys():
        if isinstance(softwares[prog], dict):
            for command in softwares[prog].keys():
                cmd = softwares[prog][command] + ' -h'
                shell = subprocess.Popen(
                    cmd.split(), stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, shell=True
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
                cmd.split(), stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, shell=True
            )
            output, error = shell.communicate()
            if error:
                logger.warning(error.decode())
            else: logger.debug(f'{prog} - {return_code(shell.returncode)}')
            result[prog] = shell.returncode
    return result

def setup_workflow():
    check_requirements()
    BLAST.accession_to_db(
        '32hiv1_default_db',
        [
            'AF164485.1', 'AF259954.1', 'AF259955.1', 'KU168309.1', 'U51189.1',
            'AB049811.1', 'AB231895.1', 'KC492737', 'KC492738', 'KC503852',
            'KC503853', 'KC503854', 'KC503855', 'KF234628', 'AF385934',
            'AF385935', 'AF385936', 'AB703607.1', 'EF158040.1', 'AF324493.2',
            'MH705157.1', 'AY835759.1', 'K03455.1', 'X01762.1', 'AY772699.1',
            'AB485648.1', 'AJ320484.1', 'KY392769.1', 'AB231893.1', 'AB485662.1',
            'AF084936.1', 'MH705134.1'
        ],
        'nucl'
    )

def _test():
    '''
    Check module funcionality.
    '''
    check_requirements()


if __name__ == '__main__':
    _test()
