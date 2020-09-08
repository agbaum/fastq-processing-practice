'''Tools for ways of outputting the results of FASTQ stats collection.'''

import pandas as pd
from abc import ABC, abstractmethod

class OutputWriter(ABC):

    '''Abstract class for stats output writers.'''

    def __init__(self, output_file: str):
        self._output_file = output_file

    @abstractmethod
    def output_stats(self, stats: pd.DataFrame):
        pass

class CsvWriter(OutputWriter):

    '''Writes stats to a csv file'''

    def output_stats(self, stats: pd.DataFrame):
        stats.to_csv(self._output_file)
