
from io import FileIO
from Bio import SeqIO
# from .logger import logger

class FASTA(object):
    def __init__(self) -> None:
        pass

    @classmethod
    def read(cls, file) -> dict:
        if isinstance(file, FileIO): filename = file.name
        else: filename = file
        # logger.debug(f'Reading sequences from {filename}')
        cls.sequences: dict[str, SeqIO.SeqRecord] = SeqIO.to_dict(SeqIO.parse(file, 'fasta'))
        return cls.sequences
    
    @classmethod
    def read_and_extract(cls, file, id):
        return cls.read(file)[id]

    def extract(self, id: str) -> SeqIO.SeqRecord:
        return self.sequences[id]

x = FASTA()
y = x.read('C:\\Users\\sarat\\Downloads\\32_HIV1_for_BLAST.fasta')
print(len(y))
print(x.sequences.keys())
