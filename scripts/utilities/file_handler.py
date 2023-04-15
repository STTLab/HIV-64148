
import pandas as pd
from io import FileIO
from Bio import SeqIO
# from .logger import logger

class FASTA(object):
    def __init__(self) -> None:
        pass

    @classmethod
    def read(cls, file) -> dict[str, SeqIO.SeqRecord]:
        if isinstance(file, FileIO): filename = file.name
        else: filename = file
        # logger.debug(f'Reading sequences from {filename}')
        cls.sequences: dict[str, SeqIO.SeqRecord] = SeqIO.to_dict(SeqIO.parse(file, 'fasta'))
        return cls.sequences
    
    @classmethod
    def read_and_extract(cls, file, id, save_to: str|None = None) -> SeqIO.SeqRecord:
        seq = cls.read(file)[id]
        if save_to: SeqIO.write(seq, save_to, 'fasta')
        return seq

    def extract(self, id: str) -> SeqIO.SeqRecord:
        return self.sequences[id]
