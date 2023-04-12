
from sierrapy import SierraClient
from sierrapy.sierraclient import Sequence
import requests
from Bio import SeqIO
from io import StringIO
from workflow.gen_report import read_haplotype_fa
from utilities.logger import logger
from utilities.settings import secrets

class EutilsNCBI:

    API_KEY = secrets.get('API_KEYS', {}).get('NCBI', '____YOUR_KEY____')
    END_POINT = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils'
    def __init__(self) -> None:
        pass

    @classmethod
    def fetch_fasta(cls, accession, save_to=None):
        params = {
            'db': 'nucleotide',
            'id': accession,
            'retmode': 'fasta',
            'rettype': 'fasta'
        }
        if cls.API_KEY != '____YOUR_KEY____': params['api_key'] = cls.API_KEY
        logger.debug(f'Fetching sequence {accession}')
        res = requests.get(f'{cls.END_POINT}/efetch.fcgi', params=params)
        try:
            seq = SeqIO.read(StringIO(res.content.decode()), 'fasta')
            summary: dict = cls.fetch_summery(accession)
            seq.description = f'{summary.get("subname","_|_").split("|")[0]}; {summary.get("title", "")}'
        except Exception as e:
            raise e
        if save_to:
            SeqIO.write(seq, save_to, 'fasta')
            return
        return seq
    
    @classmethod
    def fetch_summery(cls, accession):
        params = {
            'db': 'nucleotide',
            'id': accession,
            'retmode': 'json'
        }
        if cls.API_KEY != '____YOUR_KEY____': params['api_key'] = cls.API_KEY
        else: logger.info('NCBI key not setup. Please input your key to increase API limit.')
        logger.debug(f'Fetching summary for accession number: {accession}')
        res = requests.get(f'{cls.END_POINT}/esummary.fcgi', params=params)
        try:
            data = res.json()
            data = data['result'][data['result']['uids'][0]]
        except Exception as e:
            raise e
        return data

def hivdb_seq_analysis(haplotypes):
    sequences = []
    query = open('./scripts/utilities/gql/hiv1_sequence_analysis_default.gql', 'r').read()
    for record in haplotypes:
        sequences.append(Sequence(header=record.id, sequence=str(record.seq)))
    
    client: SierraClient = SierraClient()
    return client.sequence_analysis(sequences, query)
        

def _test():
    print(EutilsNCBI.fetch_fasta('KC492737'))
    # print(hivdb_seq_analysis(read_haplotype_fa('scripts\\tests\\mock\\HIV1_Strainline_result.fa')))

if __name__ == '__main__':
    _test()
