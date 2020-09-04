import gzip as gz
from itertools import zip_longest

class FASTQError(Exception): pass

class FASTQRecord:
    def __init__(self, name: str, sequence: str, quality: str):
        self.name = name
        self.sequence = sequence
        self.quality = quality

    def to_string(self):
        # TODO
        pass

class FASTQReader:
    def __init__(self, file):
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
    def __init__(self, read1: FASTQRecord, read2: FASTQRecord):
        if read1.name.split(' ')[0] != read2.name.split(' ')[0]:
            raise FASTQError("Mismatched paired-end records")

        self.read1 = read1
        self.read2 = read2


class PairedFASTQReader:

    def __init__(self, read1_file, read2_file):
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


    
    
