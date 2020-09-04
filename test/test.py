import unittest
import demultiplex_fastq

class TestFASTQReader(unittest.TestCase):
    def test_read_fastq(self):
        with demultiplex_fastq.FASTQReader("test_data/A_R1.fastq") as reader:
            first_record = next(reader)
            self.assertEqual(first_record.name, "TESTA:1 R1")
            self.assertEqual(first_record.sequence, "AGCT")
            self.assertEqual(first_record.quality, "ABCD")

            second_record = next(reader)
            self.assertEqual(second_record.name, "TESTA:2 R1")
            self.assertEqual(second_record.sequence, "TCGA")
            self.assertEqual(second_record.quality, "EFGH")

            with self.assertRaises(StopIteration):
                next(reader)


