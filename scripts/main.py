
import sys
from workflow import workflow

PYTHON_VERSION = sys.version_info
VERSION = "3.10.6"
PROGRAM = "HIV64148"
AUTHOR = "Sara Wattanasombat (Faculty of Medicine, Chiang Mai University, Thailand)"
CONTACT = "sara_watt@cmu.ac.th"


def main():
    workflow.main()
    
if __name__ == '__main__':
    main()
