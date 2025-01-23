# Credit to Template Code for parse_ref_file: https://github.com/eeskin/CM122_starter_code/blob/master/HP1/basic_aligner.py

import sys
import argparse
import numpy as np
import time
import zipfile
    
def parse_fasta(file_path):
    """
    :param file_path: the file containing all of the reads
    :return: outputs a list of all paired-end reads
    """
    reads = []
    with open(file_path, 'r') as file:
        current_read = ""
        for line in file:
            line = line.strip()
            # If the line starts with '>', it's the header of a new read
            if line.startswith(">"):
                if current_read:  # Add the previous read to the list
                    reads.append(current_read)
                current_read = ""  # Reset the current read
            else:
                current_read += line  # Add sequence to the current read
        if current_read:  # Append the last read
            reads.append(current_read)
    return reads

def parse_ref_file(ref_fn):
    """
    :param ref_fn: the file containing the reference genome
    :return: a string containing the reference genome
    """
    try:
        with open(ref_fn, 'r') as gFile:
            # print("Parsing Ref")
            first_line = True
            ref_genome = ''
            for line in gFile:
                if first_line:
                    first_line = False
                    continue
                ref_genome += line.strip()
        return ref_genome
    except IOError:
        print("Could not read file: ", ref_fn)
        return None