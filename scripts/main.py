#!/bin/python3
import sys
import argparse
from workflow.workflow import Worker
from utilities.logger import logger 
from utilities.requirements import setup_workflow

PYTHON_VERSION = sys.version_info
VERSION = "3.10.6"
PROGRAM = "HIV64148"
AUTHOR = "Sara Wattanasombat (Faculty of Medicine, Chiang Mai University, Thailand)"
CONTACT = "sara_watt@cmu.ac.th"


def main():
    parser = argparse.ArgumentParser(
                prog='HIV-64148 Pipeline',
                description='What the program does',
                epilog='Text at the bottom of help')
    parser.add_argument('function')
    parser.add_argument('-i', '--input', type=str, required=False)
    parser.add_argument('-o', '--output_dir', type=str, required=False)
    
    args = parser.parse_args()
    match args.function:
        case 'run':
            worker = Worker()
            job = worker.assign_job(args.input, args.output_dir, True)
            logger.info(f'Job created (id:{job})')
            worker.run_workflow()
        case 'run_cli':
            worker = Worker()
            job = worker.assign_job_cli()
            logger.info(f'Job created (id:{job})')
            worker.run_workflow()
        case 'setup':
            setup_workflow()
        case _: parser.print_help()
    
if __name__ == '__main__':
    main()
