#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 24 13:34:34 2024

@author: siobhan
"""

import gzip

def sort_sequences_by_prefix(file1, file2, prefix1="GTTC", prefix2="CCTT"):
    # Initialize dictionaries to hold paired sequences
    paired_sequences = {}

    # Function to score the similarity of a sequence to a prefix
    def score_similarity(sequence, prefix):
        return sum(1 for a, b in zip(sequence[:4], prefix) if a == b)

    # Function to process each file and store pairs
    def process_file(filename, read_type):
        with gzip.open(filename, 'rt', encoding='utf-8') as file:
            while True:
                header = file.readline().strip()
                if not header:
                    break
                sequence = file.readline().strip()
                plus = file.readline().strip()
                quality = file.readline().strip()

                read_id = header.split(' ')[0]
                if read_id not in paired_sequences:
                    paired_sequences[read_id] = {}

                paired_sequences[read_id][read_type] = (header, sequence, plus, quality)

    # Process both files
    process_file(file1, 'R1')
    process_file(file2, 'R2')

    # Initialize lists to hold sorted sequences
    prefix1_sequences = []
    prefix2_sequences = []

    # Sort pairs into appropriate output files
    for read_id, reads in paired_sequences.items():
        if 'R1' in reads and 'R2' in reads:
            seq1 = reads['R1'][1]
            seq2 = reads['R2'][1]

            if seq1.startswith(prefix1) or seq2.startswith(prefix2):
                prefix1_sequences.append(reads['R1'])
                prefix2_sequences.append(reads['R2'])
            elif seq1.startswith(prefix2) or seq2.startswith(prefix1):
                prefix1_sequences.append(reads['R2'])
                prefix2_sequences.append(reads['R1'])
            else:
                score1 = max(score_similarity(seq1, prefix1), score_similarity(seq2, prefix2))
                score2 = max(score_similarity(seq1, prefix2), score_similarity(seq2, prefix1))
                
                if score1 >= score2:
                    prefix1_sequences.append(reads['R1'])
                    prefix2_sequences.append(reads['R2'])
                else:
                    prefix1_sequences.append(reads['R2'])
                    prefix2_sequences.append(reads['R1'])

    # Write sorted sequences to new FASTQ files
    def write_sequences_to_file(sequences, output_file):
        with gzip.open(output_file, 'wt', encoding='utf-8') as file:
            for header, sequence, plus, quality in sequences:
                file.write(f"{header}\n{sequence}\n{plus}\n{quality}\n")

    write_sequences_to_file(prefix1_sequences, f"sequences_starting_with_{prefix1}.fastq.gz")
    write_sequences_to_file(prefix2_sequences, f"sequences_starting_with_{prefix2}.fastq.gz")

# Example usage
file1 = "/Users/siobhan/Downloads/SPK_YL011_2_L001_R2_001.fastq.gz"
file2 = "/Users/siobhan/Downloads/SPK_YL011_2_L001_R1_001.fastq.gz"
sort_sequences_by_prefix(file1, file2)




