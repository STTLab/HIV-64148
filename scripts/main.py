
import sys
import utilities.requirements

PYTHON_VERSION = sys.version_info
VERSION = "3.10.6"
PROGRAM = "NanoSim"
AUTHOR = "Sara Wattanasombat (Faculty of Medicine, Chiang Mai University, Thailand)"
CONTACT = "sara_watt@cmu.ac.th"


def main():
    utilities.requirements.check_requirements()

if __name__ == '__main__':
    main()
