'''Tools for running full analysis'''

from typing import List
import argparse
from . import fastq_reader
from . import stats_collectors
from . import stats_output


class AnalysisRunner:

    '''Given a set of parameters, runs a particular analysis.
    
    Parameters:
        * Input files
        * Stats calculators
        * Outputs
    '''

    def __init__(self, read1_path: str, read2_path: str,
                 stats: List[stats_collectors.PairedStatsCalculator],
                 outputters: List[stats_output.OutputWriter]):
        self._read1_path = read1_path
        self._read2_path = read2_path
        self._stats_aggregator = stats_collectors.PairedStatsAggregator(stats)
        self._outputters = outputters

    def run_analysis(self):
        with open(self._read1_path, mode='r') as read1, \
             open(self._read2_path, mode='r') as read2:
            reader = fastq_reader.PairedFASTQReader(read1, read2)

            stats = self._stats_aggregator.aggregate_stats(reader)

            for output in self._outputters:
                output.output_stats(stats)


class CLI:

    '''Command Line Interface for running analyses'''
    
    def __init__(self):
        self.parser = argparse.ArgumentParser(
            description = "Collect stats on FASTQ files. " 
                          + "At least one STATISTIC and one OUTPUT required.")
        self.parser.add_argument("input_r1",
                                 help = "FASTQ file for read 1.")
        self.parser.add_argument("input_r2",
                                 help = "FASTQ file for read 1.")
        self.parser.add_argument("--seq", nargs = "+",
                                 help = "STATISTIC: Sequences to match against.")
        self.parser.add_argument(
            "--quality", 
            action="store_true",
            help = "STATISTIC: Average per-base quality scores.")
        self.parser.add_argument('--csv',
                                 help = "OUTPUT: Csv file path.")

    def _setup_runner(self):
        stats = []
        if self.args.seq:
            stats.append(stats_collectors.SequenceMatcher(self.args.seq))
        if self.args.quality:
            stats.append(stats_collectors.QualityAverager())
        if not stats:
            raise ValueError("No statistics specified.")
        
        outputs = []
        if self.args.csv:
            outputs.append(stats_output.CsvWriter(self.args.csv))
        if not outputs:
            raise ValueError("No outputs specified.")

        self.runner = AnalysisRunner(self.args.input_r1,
                                     self.args.input_r2,
                                     stats, outputs)

    def main(self, args=None):
        if args:
            self.args = self.parser.parse_args(args)
        else:
            self.args = self.parser.parse_args()

        self._setup_runner()
        self.runner.run_analysis()
