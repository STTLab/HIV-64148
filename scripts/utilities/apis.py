
import requests
import pandas as pd
from Bio import SeqIO
from io import StringIO
import multiprocessing
from sierrapy import SierraClient
from sierrapy.sierraclient import Sequence
from utilities.logger import logger
from utilities.settings import secrets

class EutilsNCBI:

    API_KEY = secrets.get('API_KEYS', {}).get('NCBI', '____YOUR_KEY____')
    END_POINT = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils'
    RATE_LIMIT = 3

    if API_KEY != '____YOUR_KEY____':
        RATE_LIMIT = 8

    def __init__(self) -> None:
        pass

    @classmethod
    @RateLimiter(max_calls=RATE_LIMIT, period=1.5)
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
            info = dict(zip(summary.get("subtype","").split('|'), summary.get("subname","").split("|")))
            if 'subtype' in info.keys():
                subtype = 'subtype_' + info['subtype']
            else: subtype = 'subtype_' + hivdb_seq_analysis(
                                Sequence(header=seq.id, sequence=str(seq.seq)),
                                'inputSequence { header, SHA512 }, strain { name }, \
                                bestMatchingSubtype { displayWithoutDistance, distance }'
                            )[0]['bestMatchingSubtype']['displayWithoutDistance']

            seq.description = f'{subtype}; {summary.get("title", "")}'
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
        except Exception as err:
            raise err
        return data

    @classmethod
    def fetch_fasta_parallel(cls, *accession, save_to=None):
        if cls.API_KEY != '____YOUR_KEY____': nproc = 10
        else:
            nproc = 3
            logger.info('NCBI key not setup. Please input your key to increase API limit.')
        with multiprocessing.Pool(processes=nproc) as pool:
            seqs = (pool.map(cls.fetch_fasta, accession[0]))
        if save_to:
            SeqIO.write(seqs, save_to, 'fasta')
            return
        if len(seqs) == 1: return seqs[0]
        return seqs

def hivdb_seq_analysis(sequences, gql_query):
    client: SierraClient = SierraClient()
    conv_to_result = [ SequenceAnalysisResult(data=result) for result in client.sequence_analysis([*sequences], gql_query)]
    return conv_to_result

class SequenceAnalysisResult(object):
    def __init__(self, data:dict, **kwargs) -> None:
        self.inputSequence = data['inputSequence']
        self.availableGenes = [''.join(x.values()) for x in data['availableGenes']]
        self.bestMatchingSubtype = data['bestMatchingSubtype']
        self.drugResistance = self.DrugResistantProfile(data['drugResistance'])
        self.mutations = self.MutationProfile(data['mutations'])
        self.alignGeneSequences = [ self.AlignSequences(aln['gene']['name'], aln['matchPcnt'], aln['prettyPairwise']) for aln in data['alignedGeneSequences']]
        self.validationResults = data['validationResults']
        return

    class DrugResistantProfile(object):
        def __init__(self, data: dict) -> None:
            results = []
            for record in data:
                for item in record['drugScores']:
                    drug_data = item['drug']
                    data_dict = {
                        'gene': record['gene']['name'],
                        'drug': str.capitalize(drug_data['fullName']),
                        'drug_class': drug_data['drugClass']['fullName'],
                        'score': item['score'],
                        'level': item['level'],
                        'text': item['text'],
                        'SIR': item['SIR']
                    }
                    results.append(data_dict)
            self.results = results

        def to_dataframe(self) -> pd.DataFrame:
            return pd.DataFrame(self.results)

    class MutationProfile(object):
        def __init__(self, data: dict) -> None:
            results = []
            for record in data:
                data_dict = {
                    'ID': record['text'],
                    'REF': record['reference'],
                    'POS': record['position'],
                    'ALT': record['AAs'],
                    'Type': record['primaryType'],
                    'gene': record['gene']['name'],
                    'is_Apobec': record['isApobecMutation'],
                    'is_DRM': record['isDRM'],
                    'is_ApobecDRM': record['isApobecDRM'],
                    'is_Unusual': record['isUnusual'],
                    'has_Stop': record['hasStop'],
                    'comments': record['comments'][0]['text'] if len(record['comments'])>0 else ''
                }
                results.append(data_dict)
            self.results = results

        def to_dataframe(self):
            return pd.DataFrame(self.results)

    class AlignSequences(object):
        def __init__(self, gene, match_pcnt, align_data) -> None:
            self.gene = gene
            self.match_pcnt = match_pcnt
            self.align = self.PrettyPairwise(align_data)
        class PrettyPairwise(object):
            def __init__(self, prettyPairwise) -> None:
                self.prettyPairwise = prettyPairwise
                self.text = str(self)

            def __str__(self):
                items = tuple(zip(
                            self.prettyPairwise['positionLine'],
                            self.prettyPairwise['refAALine'],
                            self.prettyPairwise['alignedNAsLine'],
                            self.prettyPairwise['mutationLine']
                            )
                        )
                increment = 20
                f = 0
                t = f+increment
                pretty = ''
                while True:
                    for row in [1,2,3]:
                        if row == 2:
                            pretty += str(f+1).center(4)
                        else:
                            pretty += ' '*4
                        for item in items[f:t]:
                            pretty += item[row].center(3)
                            pretty += ' '
                        pretty += '\n'
                    pretty += '\n'
                    f=t
                    t+=increment
                    if t > len(items): t = len(items)
                    if f>=len(items): break
                return pretty


def _test():
    accession_list = [
            'AF164485.1', 'AF259954.1', 'AF259955.1', 'KU168309.1', 'U51189.1',
            'AB049811.1', 'AB231895.1', 'KC492737', 'KC492738', 'KC503852',
            'KC503853', 'KC503854', 'KC503855', 'KF234628', 'AF385934',
            'AF385935', 'AF385936', 'AB703607.1', 'EF158040.1', 'AF324493.2',
            'MH705157.1', 'AY835759.1', 'K03455.1', 'X01762.1', 'AY772699.1',
            'AB485648.1', 'AJ320484.1', 'KY392769.1', 'AB231893.1', 'AB485662.1',
            'AF084936.1', 'MH705134.1'
        ]
    # print(EutilsNCBI.fetch_fasta_parallel(accession_list[:5]))
    # print(EutilsNCBI.fetch_fasta_parallel('AF164485.1'))
    # seqs = [EutilsNCBI.fetch_fasta(accession) for accession in accession_list]
    # SeqIO.write(seqs, f'32_HIV1_for_BLAST.fasta', 'fasta')
    # print(hivdb_seq_analysis(SeqIO.parse('scripts\\tests\\mock\\HIV1_Strainline_result.fa')))

if __name__ == '__main__':
    _test()
