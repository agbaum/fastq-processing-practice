import unittest
import stats_collectors
import numpy as np
import fastq_reader
import io

record1_R1 = fastq_reader.FASTQRecord("test1 R1", "abcd", "0678")
record1_R2 = fastq_reader.FASTQRecord("test1 R2", "bxdxyz", "98543")
record2_R1 = fastq_reader.FASTQRecord("test2 R1", "axyz", "5678")
record2_R2 = fastq_reader.FASTQRecord("test2 R2", "bada", "0968")
paired_record1 = fastq_reader.PairedFASTQRecord(record1_R1, record1_R2)
paired_record2 = fastq_reader.PairedFASTQRecord(record2_R1, record2_R2)


class TestSequenceMatcher(unittest.TestCase):

    def setUp(self):
        self.patterns = [r"abc", r"xyz", r"b.d"]
        self.seq_matcher = stats_collectors.SequenceMatcher(self.patterns)

    def test_no_blank(self):
        with self.assertRaises(ValueError):
            stats_collectors.SequenceMatcher([])

    def test_query(self):
        matches = self.seq_matcher._query("abcd")
        for i in range(len(self.patterns)):
            self.assertEqual(self.patterns[i], matches.index[i])

        self.assertEqual(matches[0], 0)
        self.assertTrue(np.isnan(matches[1]))
        self.assertEqual(matches[2], 1)
        
    def test_single(self):
        matches = self.seq_matcher._calc_single_stats(record1_R1)
        
        self.assertEqual(matches[0], 0)
        self.assertTrue(np.isnan(matches[1]))
        self.assertEqual(matches[2], 1)

    def test_paired(self):
        stats = self.seq_matcher.calc_paired_stats(paired_record1)
        self.assertEqual(stats[1, r'abc'], 0)
        self.assertTrue(np.isnan(stats[1, r'xyz']))
        self.assertEqual(stats[1, r'b.d'], 1)
        self.assertTrue(np.isnan(stats[2, r'abc']))
        self.assertEqual(stats[2, r'xyz'], 3)
        self.assertEqual(stats[2, r'b.d'], 0)
        
class TestQualityAverager(unittest.TestCase):

    def setUp(self):
        self.qual_averager = stats_collectors.QualityAverager()

    def test_paired(self):
        qual = self.qual_averager.calc_paired_stats(paired_record1)
        self.assertEqual(qual[1, "avg_qual"], np.mean(record1_R1.quality))
        self.assertEqual(qual[2, "avg_qual"], np.mean(record1_R2.quality))


class TestAggregator(unittest.TestCase):

    def setUp(self):
        patterns = [r"abc", r"xyz", r"b.d"]
        seq_matcher = stats_collectors.SequenceMatcher(patterns)
        qual_averager = stats_collectors.QualityAverager()
        aggregator = stats_collectors.PairedStatsAgregator([seq_matcher,
                                                                 qual_averager])
        sequence1_R1 = str(record1_R1)
        sequence1_R2 = str(record1_R2)
        sequence2_R1 = str(record2_R1)
        sequence2_R2 = str(record2_R2)
        read1 = sequence1_R1 + sequence2_R1
        read2 = sequence1_R2 + sequence2_R2
        
        with io.StringIO(read1) as file1, io.StringIO(read2) as file2:
            paired_reader = fastq_reader.PairedFASTQReader(file1, file2)
            self.stats = aggregator.aggregate_stats(paired_reader)    
            
    def test_quality(self):
        self.assertEqual(self.stats.loc[0, (1, "avg_qual")],
                         np.mean(record1_R1.quality))
        self.assertEqual(self.stats.loc[0, (2, "avg_qual")],
                         np.mean(record1_R2.quality))
        self.assertEqual(self.stats.loc[1, (1, "avg_qual")],
                         np.mean(record2_R1.quality))
        self.assertEqual(self.stats.loc[1, (2, "avg_qual")],
                         np.mean(record2_R2.quality))
        
    def test_abc(self):
        self.assertEqual(self.stats.loc[0, (1, r'abc')], 0)
        self.assertTrue(np.isnan(self.stats.loc[0, (2, r'abc')]))
        self.assertTrue(np.isnan(self.stats.loc[1, (1, r'abc')]))
        self.assertTrue(np.isnan(self.stats.loc[1, (2, r'abc')]))

    def test_xyz(self):
        self.assertTrue(np.isnan(self.stats.loc[0, (1, r'xyz')]))
        self.assertEqual(self.stats.loc[0, (2, r'xyz')], 3)
        self.assertEqual(self.stats.loc[1, (1, r'xyz')], 1)
        self.assertTrue(np.isnan(self.stats.loc[1, (2, r'xyz')]))

    def test_b_d(self):
        self.assertEqual(self.stats.loc[0, (1, r'b.d')], 1)
        self.assertEqual(self.stats.loc[0, (2, r'b.d')], 0)
        self.assertTrue(np.isnan(self.stats.loc[1, (1, r'b.d')]))
        self.assertEqual(self.stats.loc[1, (2, r'b.d')], 0)

        
        