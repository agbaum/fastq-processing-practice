from typing import List
import fastq_reader
import stats_collectors
import stats_output


class AnalysisRunner:

    def __init__(self, input_file_paths: List[str], 
                 stats_aggregator: stats_collectors.PairedStatsAgregator,
                 outputters: List[stats_output.OutputWriter]):
        self._read1_path = input_file_paths[0]
        self._read2_path = input_file_paths[1]
        self._stats_aggregator = stats_aggregator
        self._outputters = outputters

    def run_analysis(self):
        with open(self._read1_path, mode='r') as read1, \
             open(self._read2_path, mode='r') as read2:
            reader = fastq_reader.PairedFASTQReader(read1, read2)

            stats = self._stats_aggregator.aggregate_stats(reader)

            for output in self._outputters:
                output.output_stats(stats)