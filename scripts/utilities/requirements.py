
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

def check_requirements() -> dict:
    '''
    Check for return code 0 of the required softwares by calling their help command.
    '''
    result = {}
    softwares = settings.get('softwares', {})
    for prog in softwares.keys():
        cmd = softwares[prog] + ' -h'
        shell = subprocess.Popen(
            cmd.split(), stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, shell=True
        )
        output, error = shell.communicate()
        if error:
            logger.warning(error.decode())
        else: logger.debug(f'{prog} - {return_code(shell.returncode)}')
            
        result[prog] = return_code(shell.returncode)
    return result

def _test():
    '''
    Check module funcionality.
    '''
    check_requirements()


if __name__ == '__main__':
    _test()
