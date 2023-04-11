
import sys
from utilities import apis
from workflow import workflow

PYTHON_VERSION = sys.version_info
VERSION = "3.10.6"
PROGRAM = "HIV64148"
AUTHOR = "Sara Wattanasombat (Faculty of Medicine, Chiang Mai University, Thailand)"
CONTACT = "sara_watt@cmu.ac.th"


def main():
    # workflow.main()
    apis._test()
    
if __name__ == '__main__':
    main()
