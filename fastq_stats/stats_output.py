'''Tools for ways of outputting the results of FASTQ stats collection.'''

from abc import ABC, abstractmethod
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import axes

class OutputWriter(ABC):

    '''Abstract class for stats output writers.'''

    def __init__(self, output_file: str):
        self._output_file = output_file

    @abstractmethod
    def output_stats(self):
        pass

class CsvWriter(OutputWriter):

    '''Writes stats to a csv file'''

    def output_stats(self, stats):
        stats.to_csv(self._output_file)

class HistogramPlotter(OutputWriter):

    '''Makes Histograms of the stats, and saves them to a file.'''

    def output_stats(self, stats):
        nstats = stats.shape[1]
        fig = plt.figure(figsize=(4, 2*nstats))
        fig.subplots_adjust(hspace=0.4, wspace=0.4)
        for i, stat in enumerate(stats.columns):
            ax = fig.add_subplot(nstats, 1, i + 1)
            filtered_stat = [s for s in stats[stat] if np.isfinite(s)]
            if len(filtered_stat) != 0:
                sns.distplot(filtered_stat, label = stat, ax=ax)
                ax.set_title(stat)

        fig.savefig(self._output_file)