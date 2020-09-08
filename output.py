import pandas as pd

class csv_writer:

    def __init__(self, output_file: str):
        self._output_file = output_file

    def output_stats(self, stats: pd.DataFrame):
        stats.to_csv(self._output_file)
