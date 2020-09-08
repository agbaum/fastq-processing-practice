import unittest
import stats_collectors
import numpy as np
import fastq_reader

record1 = fastq_reader.FASTQRecord("test1", "abcd", "0678")
record2 = fastq_reader.FASTQRecord("test1", "bxdxyz", "98543")
paired_record = fastq_reader.PairedFASTQRecord(record1, record2)


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
        matches = self.seq_matcher._calc_single_stats(record1)
        
        self.assertEqual(matches[0], 0)
        self.assertTrue(np.isnan(matches[1]))
        self.assertEqual(matches[2], 1)

    def test_paired(self):
        stats = self.seq_matcher.calc_paired_stats(paired_record)
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
        qual = self.qual_averager.calc_paired_stats(paired_record)
        self.assertEqual(qual[1, "avg_qual"], np.mean(record1.quality))
        self.assertEqual(qual[2, "avg_qual"], np.mean(record2.quality))



