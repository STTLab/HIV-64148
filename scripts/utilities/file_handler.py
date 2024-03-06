'''
HIV-64148  Copyright (C) 2024  Sara Wattanasombat
This program comes with ABSOLUTELY NO WARRANTY.
This is free software, and you are welcome to redistribute it.
'''

__version__ = '0.1'
__author__ = 'Sara Wattanasombat'

import pandas as pd
from io import FileIO
from Bio import SeqIO
# from .logger import logger

class FASTA(object):
    def __init__(self, sequences) -> None:
        self.sequences = sequences

    @classmethod
    def read(cls, file):
        if isinstance(file, FileIO): filename = file.name
        else: filename = file
        # logger.debug(f'Reading sequences from {filename}')
        cls.sequences: dict[str, SeqIO.SeqRecord] = SeqIO.to_dict(SeqIO.parse(file, 'fasta'))
        return FASTA(cls.sequences)
    
    @classmethod
    def read_and_extract(cls, file, id, save_to: str|None = None) -> str | SeqIO.SeqRecord:
        seq = cls.read(file).extract(id, save_to)
        return seq

    def extract(self, id: str, save_to: str|None = None) -> SeqIO.SeqRecord | str:
        seq = self.sequences[id]
        if save_to:
            SeqIO.write(seq, save_to, 'fasta')
            return save_to
        return seq
    
    def list_ids(self):
        return tuple(self.sequences.keys())

    @classmethod
    def write(cls, sequences, output):
        SeqIO.write(sequences, output, 'fasta')
