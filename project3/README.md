### Roland Yang, UID: 506053914, C122 Winter 2025

# Project 3 C122

This project is a motif identifier (de novo motif discovery). In a Chromatin Immunoprecipitation seqeuencing (ChIP-seq) experiement, this project identifies the peak summit (which is where the transcription factor is most likely binding DNA aka motif) in each sequence presented.

## Usage

```bash
python project3a.py -s [pwm-sequences-filepath] -i [input-filepath] -o [output-filepath]
```

- -s, --sequences: The input sequences to build the PWM on
- -i, --input: The input data to identify peaks on using the PWM
- -o, --output: The output filepath, should be a .txt file

## Project 3a

Given a set of sequences where the motif is likely centered towards the middle, each 201bp, identifies the peak summit in another set of sequences with offsets such that the motif is not always in the middle.

### Score

Overall: 5730/10227

### Implementation

Generates a PWM using the EM (Expectation-Maximum) Algorithm for de novo motif discovery. Found acceptable results using a kmer size of 20 and running the algorithm for 100 iterations.

Used the PWM generated to calculate the most probable kmer in each sequence in the second set of sequences, then found the index of where that kmer first appears in each of the sequences.