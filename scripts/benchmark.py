import argparse
from utilities.logger import logger
from workflow.workflow import Worker
from tests.simulator import Simulator, Simulation

def main():
    parser = argparse.ArgumentParser(
                prog='HIV-64148 Pipeline',
                description='What the program does',
                epilog='Text at the bottom of help')
    parser.add_argument('function')
    parser.add_argument(
        '-meta', '--metadata', type=str, required=True,
        help='A CSV or TSV file with headers, containing metadata for sequences in Multi-fasta file.'
    )
    parser.add_argument('-refs', '--references', type=str, required=True, help='Multi-fasta file')
    parser.add_argument('-o', '--output_dir', type=str, required=True, help='output directory')
    parser.add_argument('--perfect', action=argparse.BooleanOptionalAction)

    args = parser.parse_args()
    match args.function:
        case 'setup':
            Simulator.install_nanosim()
        case 'simulate_random':
            Simulation.test_random(
                path_to_metadata=args.metadata,
                output_dir=f'{args.output_dir}/data',
                path_to_fasta=args.references,
                prob=[0.0593007943587783, 0.0019115585574104752, 0.9356866743129009, 0.0013168514506605496, 0.0002548744743213967, 0.00038231171148209506, 0.0, 0.0, 4.247907905356612e-05, 0.0, 0.0, 0.001104456055392719],
                perfect=args.perfect
            )
            worker = Worker()
            job = worker.assign_job(f'{args.output_dir}/data/simulated_all_reads.fastq', f'{args.output_dir}/pipeline_output', True)
            logger.info(f'Job created (id:{job})')
            worker.run_workflow()

        case 'simulate_file':
            Simulation.test_from_csv(
                path_to_input=args.metadata,
                output_dir=f'{args.output_dir}/data',
                path_to_fasta=args.references,
                perfect=args.perfect
            )


if __name__ == '__main__':
    main()