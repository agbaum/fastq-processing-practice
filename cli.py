import argparse
import stats_collectors
import stats_output
import run_analysis


class CLI:
    
    def __init__(self):
        parser = argparse.ArgumentParser(description = "Collect stats on FASTQ files.")
        parser.add_argument("input_r1",
                            help = "FASTQ file for read 1.")
        parser.add_argument("input_r2",
                            help = "FASTQ file for read 1.")
        parser.add_argument("--seq", nargs = "+",
                            help = "Sequences to match against.")
        parser.add_argument("--quality", action="store_false",
                            help = "Option to report average per-base quality scores.")
        parser.add_argument('--csv',
                            help = "Output csv file path.")
        self.args = parser.parse_args()

    def main(self):
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

        runner = run_analysis.AnalysisRunner(self.args.read1_path,
                                             self.args.read2_path,
                                             stats, outputs)


if __name__ == "__main__":
    cli = CLI()
    cli.main()

        

        