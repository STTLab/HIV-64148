
from sierrapy import SierraClient
from sierrapy.sierraclient import Sequence
import requests
from Bio import SeqIO
from io import StringIO
from workflow.gen_report import read_haplotype_fa

class EutilsNCBI():
    def __init__(self) -> None:
        pass

    @classmethod
    def fetch_fasta(self, accession, save_to=None):
        url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils'
        params = {
        'db': 'nucleotide',
        'id': accession,
        'retmode': 'fasta',
        'rettype': 'fasta'
        }
        res = requests.get(f'{url}/efetch.fcgi', params=params)
        try:
            seq = SeqIO.read(StringIO(res.content.decode()), 'fasta')
        except Exception as e:
            raise e
        if save_to:
            SeqIO.write(seq, save_to, 'fasta')
            return
        return seq


def hivdb_seq_analysis(haplotypes):
    sequences = []
    query = open('./scripts/utilities/gql/hiv1_sequence_analysis_default.gql', 'r').read()
    for record in haplotypes:
        sequences.append(Sequence(header=record.id, sequence=str(record.seq)))
    
    client: SierraClient = SierraClient()
    return client.sequence_analysis(sequences, query)
        

def _test():
    print(EutilsNCBI.fetch_fasta('AF164485.1'))
    # print(hivdb_seq_analysis(read_haplotype_fa('scripts\\tests\\mock\\HIV1_Strainline_result.fa')))

if __name__ == '__main__':
    _test()
