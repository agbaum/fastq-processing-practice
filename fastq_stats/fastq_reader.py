'''Tools for working with FASTQ files.

Contains data structures for paired and single fastq records,
and iterators for reading through paired and unparied fastq files.

'''

import gzip as gz
import io
from itertools import zip_longest
from typing import List

class FASTQError(Exception): pass

def qual_str2num(qual_str: str) -> List[int]:
    '''Converts a base quality string into a list of ints.
    
    Currently just takes single digits, 
    which isn't how actual fastq quality strings work.
    This is simpler as practice implimentation for now.
    
    '''
    return [int(c) for c in qual_str]

def qual_num2str(qual_list: List[int]) -> str:
    '''Converts a base quality string into a list of ints.'''
    return ''.join([str(i) for i in qual_list])


class FASTQRecord:

    '''Data structure for one FASTQ read
    
    Attributes:
        name, sequence, quality: from the three meaningful lines in 
            a fastq file.
    
    '''

    def __init__(self, name: str, sequence: str, quality: str):
        self.name = name
        self.sequence = sequence
        self.quality = qual_str2num(quality)

    def __str__(self):
        s = (self.name + "\n"
             + self.sequence + "\n" 
             + "+\n" 
             + qual_num2str(self.quality) + "\n"
        )
        return (s)

class FASTQReader:

    '''Iterates through a fastq file, returning records.'''

    def __init__(self, file: io.TextIOBase):
        self._file = file

    def __iter__(self):
        return self

    def __next__(self) -> FASTQRecord:
        try:
            name = next(self._file)[:-1]
        except StopIteration:
            raise StopIteration

        try:
            sequence = next(self._file)[:-1]
            plus = next(self._file)[:-1]
            quality = next(self._file)[:-1]
        except StopIteration:
            raise FASTQError("Incomplete FASTQ file")

        return(FASTQRecord(name, sequence, quality))


class PairedFASTQRecord:

    '''Data structure for one paired-end FASTQ read
    
    Attributes: 
        read1, read2: FASTQ records

    '''

    def __init__(self, read1: FASTQRecord, read2: FASTQRecord):
        if read1.name.split(' ')[0] != read2.name.split(' ')[0]:
            raise FASTQError("Mismatched paired-end records")

        self.read1 = read1
        self.read2 = read2


class PairedFASTQReader:

    '''
    Iterates through paired FASTQ files and returns 
    paired records.
    '''

    def __init__(self, read1_file: io.TextIOBase, read2_file: io.TextIOBase):
        self._reader1 = FASTQReader(read1_file)
        self._reader2 = FASTQReader(read2_file)
        self._zipped_reader = zip_longest(self._reader1, self._reader2)


    def __iter__(self):
        return self

    def __next__(self) -> PairedFASTQRecord:
        read1_record, read2_record = next(self._zipped_reader)
        
        if read1_record is None or read2_record is None:
            raise FASTQError("Different length FASTQ files")

        return PairedFASTQRecord(read1_record, read2_record)


    
    
