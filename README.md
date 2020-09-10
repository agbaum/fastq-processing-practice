# FASTQ Pracitce
A project for practicing OOP in python.

A tool for reading through paired fastq files, compiling various statistics on them, and outputting those statstics in various ways (plots, tables, etc.).

## Module Structure
* fastq_reader.py - tools for reading fastq files
* stats_collectors.py - different classes for generating stats from fastq records, and tools for aggregating those stats over a whole file.
* stats_output.py - different classes for outputting statistics results in different ways.
* run_analysis.py - A class for running a whole analysis (read -> stats -> output) and a command line interface for using it.

## Example
test_out.csv can be generated with `python -m fastq_stats --seq AG TCC --quality --csv .\test\test_data\test_out.csv .\test\test_data\test_R1.fastq .\test\test_data\test_R2.fastq`

# TODOs
* Write plotting output
* Figure out how/if to do test for outputs
* Figure out integration/application level testing