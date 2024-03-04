#!/bin/python3
'''
This is the main entry of HIV-64148 pipeline
'''

__all__ = ['main',]
__version__ = '0.1'
__author__ = 'Sara Wattanasombat'

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
    '''
    This is the main function for handling command line argument and pipeline execution.
    '''
    parser = argparse.ArgumentParser(
                prog='HIV-64148 Pipeline',
                description='What the program does',
                epilog='Text at the bottom of help')
    parser.add_argument(
        'function',
        choices=('run', 'report')
    )
    parser.add_argument(
        '-i', '--input',
        type=str,
        required=True,
        help='Path to input file in FASTQ format.'
    )
    parser.add_argument(
        '-o', '--output_dir',
        type=str,
        required=True,
        help='Path to Output directory.'
    )
    parser.add_argument(
        '-r', '--reference',
        type=str,
        required=False,
        help='Path to reference genome, required for reference-based assemblers.'
    )
    parser.add_argument(
        '-rcon', '--reconstructor',
        type=str,
        required=False,
        choices=(
            'canu', 'strainline', 'goldrush',
            'metaflye', 'rvhaplo', 'haplodmf', 'igda'
        ),
        default='strainline',
        help='Assembler selection',
    )
    parser.add_argument(
        '--no-report',
        default=False,
        action=argparse.BooleanOptionalAction,
        help='Don\'t generate HTML report.'
    )
    parser.add_argument(
        '--overwrite',
        required=False,
        default=False,
        action=argparse.BooleanOptionalAction,
        help='Force overwrite the output dirtectory. (CANNOT BE RECOVERED)',
    )
    args = parser.parse_args()
    match args.function:
        case 'run':
            worker = Worker(reconstructor=args.reconstructor)
            job = worker.assign_job(args.input, args.output_dir, args.overwrite)
            logger.info('Job created (id: %s)', job)
            worker.run_workflow()
        case 'setup':
            setup_workflow()
        case _: parser.print_help()

if __name__ == '__main__':
    main()
