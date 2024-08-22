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
__version__ = '0.2'
__author__ = 'Sara Wattanasombat'

import sys
import argparse
from utilities.logger import logger
try:
    from workflow.workflow import Worker
except ModuleNotFoundError as err:
    logger.fatal(err)

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
                prog='hiv64148',
                description='About HIV-64148, an integration of multiple long-read genome \
                            assemblers with a pipeline for analysis of HIV-1 genomic data \
                            from Oxford Nanopore Sequencing Technology or PacBio Real-Time \
                            (SMRT) Sequencing technology.',
                epilog='Citing our pipeline use https://doi.org/10.12688/f1000research.149577.1\n'
                'along with an appropriate citataion of the selected assembler.',
                formatter_class=argparse.RawTextHelpFormatter
    )
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
        '-tech', '--tech',
        dest='tech',
        type=str,
        choices=('nanopore', 'pacbio'),
        default='nanopore',
        help='Technology used to generate the reads. (default: nanopore)'
    )
    parser.add_argument(
        '-qual', '--quality',
        dest='quality',
        type=str,
        choices=('raw', 'corrected', 'hifi', 'hq'),
        default='raw',
        help='Expected error rate (default: raw)\n'
        '- raw: <20%% error\n'
        '- corrected: <3%% error\n'
        '- hifi: <1%% error (PacBio only)\n'
        '- hq: <5%% berror (ONT only, Guppy5+ SUP or Q20)'
    )
    parser.add_argument(
        '-a', '--assembler',
        type=str,
        required=False,
        choices=(
            'canu', 'strainline', 'goldrush',
            'flye', 'rvhaplo', 'haplodmf', 'igda'
        ),
        default='strainline',
        help='Assembler selection (default: strainline)',
    )
    parser.add_argument(
        '-r', '--reference',
        type=str,
        required=False,
        default=None,
        help='Path to reference genome, required for reference-based assemblers.'
    )
    parser.add_argument(
        '-g', '--genome-size',
        dest='genome_size',
        type=str,
        required=False,
        default=None,
        help='An estimate of the size of the genome. Requried for running Canu. (default: not set)'
    )
    parser.add_argument(
        '-ag', '--assembler-args',
        dest='assembler_args',
        type=str,
        default='',
        required=False,
        help='Path to a yaml file for assembler settings.',
    )
    parser.add_argument(
        '--non-hiv',
        dest='non_hiv',
        default=False,
        action=argparse.BooleanOptionalAction,
        help='This is a non-HIV sequence. HTML report related to HIV will not be generated.',
    )
    parser.add_argument(
        '--overwrite',
        required=False,
        default=False,
        action=argparse.BooleanOptionalAction,
        help='Force overwrite the output dirtectory. (CANNOT BE RECOVERED)',
    )
    parser.add_argument(
        '-db', '--blast_list',
        dest='db',
        type=str,
        required=False,
        default='LosAlamos_db',
        help='Path to a FASTA file containing sequences you wish to use for BLAST'
        '\nOR `-db remote` if you want to use a remote NCBI BLAST server. (default: LosAlamos_db)'
    )
    args = parser.parse_args()
    match args.function:
        case 'run':
            worker = Worker(
                assembler=args.assembler,
                tech=args.tech,
                quality=args.quality,
                asm_args=args.assembler_args,
                kwargs={'non_hiv':args.non_hiv}
            )
            worker.check_assembler()
            if args.reference:
                worker.set_reference(args.reference)
            if args.db:
                if args.non_hiv and args.db == 'LosAlamos_db':
                    logger.warning('Database not provided for non-HIV sequences, '
                                    'switch to remote NCBI BLAST server.')
                    worker.set_blast_db('remote')
                elif args.db != 'LosAlamos_db':
                    worker.set_blast_db(args.db)
            if args.genome_size:
                worker.set_genome_size(args.genome_size)
            else:
                if args.assembler == 'canu':
                    if args.non_hiv:
                        logger.critical('Canu requires genome size input.'
                                        'Please set via -g/--genome-size argument. Example: -g 9.8k'
                                        )
                        sys.exit(1)
                    else:
                        worker.set_genome_size('9.8k')

            job = worker.assign_job(args.input, args.output_dir, args.overwrite)
            logger.info('Job created (id: %s)', job)
            worker.run_workflow()
        case _:
            parser.print_help()

if __name__ == '__main__':
    main()
