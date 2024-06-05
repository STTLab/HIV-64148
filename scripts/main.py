#!/bin/python3
'''
    HIV-64148, a pipeline for analysis of HIV-1 genomic data based on
    long-read Oxford Nanopore Sequencing.
    Copyright (C) 2024  Sara Wattanasombat

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License or any
    later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

    This module is the main entry of HIV-64148 pipeline.
'''

__all__ = ['main',]
__version__ = '0.1'
__author__ = 'Sara Wattanasombat'

import sys
import argparse
from workflow.workflow import Worker
from utilities.logger import logger

PYTHON_VERSION = sys.version_info
PYTHON_DEV_VERSION = "3.10.6"
PROGRAM = "HIV64148"
AUTHOR = "Sara Wattanasombat (Faculty of Medicine, Chiang Mai University, Thailand)"
CONTACT = "sara_watt@cmu.ac.th"

def main():
    '''
    This is the main function for handling command line argument and pipeline execution.
    '''
    parser = argparse.ArgumentParser(
                prog='HIV-64148 Pipeline',
                description='''About HIV-64148, an integration of multiple long-read genome assemblers
                               with a pipeline for analysis of HIV-1 genomic data from Oxford Nanopore
                               Sequencing Technology or PacBio Real-Time (SMRT) Sequencing technology.''',
                epilog='Citing our pipeline use https://doi.org/10.12688/f1000research.149577.1')
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
        default=None,
        help='Path to reference genome, required for reference-based assemblers.'
    )
    parser.add_argument(
        '-a', '--assembler',
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
        '-ag', '--assembler-args',
        dest='assembler_args',
        type=str,
        required=False,
        help='A quoted string of custom parameters for the selected assembler\n\
            Requires equal sign (=) after the argument.\n\
            Example: -ag="--minTrimmedLen 500 --minOvlpLen 1000 -t 8"',
    )
    parser.add_argument(
        '--no-report',
        dest='no_report',
        default=False,
        action=argparse.BooleanOptionalAction,
        help='Don\'t generate HTML report.',
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
            worker = Worker(assembler=args.assembler, asm_args=args.assembler_args.split())
            worker.check_assembler()
            worker.set_reference(args.reference)
            job = worker.assign_job(args.input, args.output_dir, args.overwrite)
            logger.info('Job created (id: %s)', job)
            worker.run_workflow()
        case _: parser.print_help()

if __name__ == '__main__':
    main()
