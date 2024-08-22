'''
HIV-64148  Copyright (C) 2024  Sara Wattanasombat
This program comes with ABSOLUTELY NO WARRANTY.
This is free software, and you are welcome to redistribute it.
'''

__version__ = '0.2'
__author__ = 'Sara Wattanasombat'

import os
import yaml
from io import FileIO
from Bio import SeqIO
# from .logger import logger

def cmd_args_builder(params: dict, arg_prefix:str='--', value_assign_sign:str=None):
    args = []
    for key, value in params.items():
        if value:
            if value_assign_sign:
                args.append(f'{arg_prefix}{key}{value_assign_sign}{str(value)}')
            else:
                args.extend([f'{arg_prefix}{key}', str(value)])
    return args

def args_loader(asm_settings: str|os.PathLike, arg_prefix:str='--', value_assign_sign:str=None):
    with open(asm_settings, 'r', encoding='utf8') as p_file:
        return cmd_args_builder(yaml.safe_load(p_file), arg_prefix, value_assign_sign)
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
