import unittest
import fastq_reader
import io

read1_R1 = "TESTA:1 R1\nAGCT\n+\nABCD\n"
read2_R1 = "TESTA:2 R1\nTCGA\n+\nEFGH\n"

class TestFASTQReader(unittest.TestCase):

    def test_read_fastq(self):

        with io.StringIO(read1_R1 + read2_R1) as input_file:
            reader = fastq_reader.FASTQReader(input_file)
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

    def test_incomplete_record(self):
        l1 = "TESTA:1 R1\n"
        l2 = "AGCT\n"
        l3 = "+\n"

        with io.StringIO(l1) as input_file:
            reader = fastq_reader.FASTQReader(input_file)
            with self.assertRaises(fastq_reader.FASTQError):
                next(reader)

        with io.StringIO(l1 + l2) as input_file:
            reader = fastq_reader.FASTQReader(input_file)
            with self.assertRaises(fastq_reader.FASTQError):
                next(reader)

        with io.StringIO(l1 + l2 + l3) as input_file:
            reader = fastq_reader.FASTQReader(input_file)
            with self.assertRaises(fastq_reader.FASTQError):
                next(reader)

