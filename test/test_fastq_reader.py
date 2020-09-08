import unittest
import fastq_reader
import io

sequence1_R1 = "TESTA:1 R1\nAGCT\n+\n1234\n"
sequence2_R1 = "TESTA:2 R1\nTCGA\n+\n5678\n"

sequence1_R2 = "TESTA:1 R2\nCCGT\n+\n1483\n"
sequence2_R2 = "TESTA:2 R2\nGTCT\n+\n0968\n"

class TestFASTQRecord(unittest.TestCase):

    def test_to_string(self):
        record = fastq_reader.FASTQRecord("TESTA:1 R1", "AGCT", "1234")
        self.assertEqual(str(record), sequence1_R1)

class TestFASTQReader(unittest.TestCase):

    def test_read_fastq(self):
        with io.StringIO(sequence1_R1 + sequence2_R1) as input_file:
            reader = fastq_reader.FASTQReader(input_file)
            first_record = next(reader)
            self.assertEqual(str(first_record), sequence1_R1)

            second_record = next(reader)
            self.assertEqual(str(second_record), sequence2_R1)

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

class TestPairedFASTQ(unittest.TestCase):
    
    def test_correct_pair(self):
        with io.StringIO(sequence1_R1) as file1, io.StringIO(sequence1_R2) as file2:
            reader1 = fastq_reader.FASTQReader(file1)
            reader2 = fastq_reader.FASTQReader(file2)

            record1 = next(reader1)
            record2 = next(reader2)
            paired = fastq_reader.PairedFASTQRecord(record1, record2)

            self.assertEqual(paired.read1, record1)
            self.assertEqual(paired.read2, record2)

    def test_bad_pair(self):
        with io.StringIO(sequence1_R1) as file1, io.StringIO(sequence2_R1) as file2:
            reader1 = fastq_reader.FASTQReader(file1)
            reader2 = fastq_reader.FASTQReader(file2)

            record1 = next(reader1)
            record2 = next(reader2)

            with self.assertRaises(fastq_reader.FASTQError):
                paired = fastq_reader.PairedFASTQRecord(record1, record2)

class TestPairedFASTQReader(unittest.TestCase):

    def test_read_pair(self):
        read1 = sequence1_R1 + sequence2_R1
        read2 = sequence1_R2 + sequence2_R2
        with io.StringIO(read1) as file1, io.StringIO(read2) as file2:
            paired_reader = fastq_reader.PairedFASTQReader(file1, file2)
            
            paired_record1 = next(paired_reader)
            self.assertEqual(str(paired_record1.read1), sequence1_R1)
            self.assertEqual(str(paired_record1.read2), sequence1_R2)

            paired_record2 = next(paired_reader)
            self.assertEqual(str(paired_record2.read1), sequence2_R1)
            self.assertEqual(str(paired_record2.read2), sequence2_R2)

            with self.assertRaises(StopIteration):
                next(paired_reader)

    def test_mismatched(self):
        read1 = sequence1_R1 + sequence1_R2
        read2 = sequence1_R2
        
        with io.StringIO(read1) as file1, io.StringIO(read2) as file2:
            paired_reader = fastq_reader.PairedFASTQReader(file1, file2)

            paired_record = next(paired_reader)
            with self.assertRaises(fastq_reader.FASTQError):
                next(paired_reader)
            

