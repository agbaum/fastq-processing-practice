'''Tools for calculating statistics on FASTQ files.

Contains classes that calculate statistics on paired FASTQ records,
and an aggregator for collecting the stats over a whole file.

'''

from typing import List, Dict, Iterable
from abc import ABC, abstractmethod
import re
from . import fastq_reader
import pandas as pd
import numpy as np


class PairedStatsCalculator(ABC):

    '''Abstract Class for calculators of statistics on paired FASTQ records.'''

    def __init__(self):
        pass

    def calc_paired_stats(self, 
                          record: fastq_reader.PairedFASTQRecord) -> pd.Series:
        stats1 = self._calc_single_stats(record.read1)
        stats2 = self._calc_single_stats(record.read2)

        return(pd.concat([stats1, stats2], keys = [1, 2]))
    
    @ abstractmethod
    def _calc_single_stats(self, record: fastq_reader.FASTQRecord) -> pd.Series:
        pass

class PairedStatsAggregator():

    '''Aggregator for multiple stats collectors over paired FASTQ files'''

    def __init__(self, calculators = List[PairedStatsCalculator]):
        if not calculators:
            raise ValueError("Empty calculators list.")
        self._calculators = calculators

    def _single_read_stats(self, 
                           record: fastq_reader.PairedFASTQRecord) -> pd.Series:
        stats = [calc.calc_paired_stats(record) for calc in self._calculators]
        return pd.concat(stats)

    def aggregate_stats(self, 
                        reader: fastq_reader.PairedFASTQReader) -> pd.DataFrame:
        stats = []
        for record in reader:
            stats.append(self._single_read_stats(record))

        return(pd.DataFrame(stats)) 


class SequenceMatcher(PairedStatsCalculator):

    '''Matches read sequences to patterns'''

    def __init__(self, patterns: List[str]):
        if not patterns:
            raise ValueError("Empty patterns list.")
        self._patterns = {p: re.compile(p) for p in patterns}

    def _query(self, s: str) -> pd.Series:
        matches = pd.Series(np.nan, index = self._patterns)
        for pattern, compiled in self._patterns.items():
            match = compiled.search(s)
            if match:
                matches[pattern] = match.start()
        return matches

    def _calc_single_stats(self, record: fastq_reader.FASTQRecord) -> pd.Series:
        return self._query(record.sequence)
            
class QualityAverager(PairedStatsCalculator):

    '''Averages the per-base quality'''

    def _calc_single_stats(self, record: fastq_reader.FASTQRecord) -> pd.Series:
        return pd.Series(np.mean(record.quality), index = ["avg_qual"])