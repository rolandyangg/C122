
### Roland Yang, UID: 506053914, C122 Winter 2025

# Project 1a & 1b C122

This project takes in a reference genome sequence along with reads from a donor genome sequence, and identifies variants including: substitutions, insertions, and deletions.

## Usage

```bash
python project1X.py -g [reference-genome-filepath] -r [reads-filepath] -o [output-filepath]
```

- -g, --genome: The reference genome filepath
- -r, --reads: The reads filepath
- -o, --output: The output filepath, should be a .txt file

## Project 1a

Identifies strictly subsitutions.

### Score

Overall: 0.728
Substitutions: 0.8444

### Implementation

Utilizes a Sliding Window approach that aligns if the amount of mismatches between a read and the genome is under a certain threshold. From there the mismatches/substitutions in the aligned read are stored in a dictionary. The dictionary containing the mutations is filtered by throwing away all mutations that do not meet the threshold of occurences necessary.

## Project 1b

Identifies insertions and deletions in addition to substitutions.

### Score

Overall: 0.919
Substitutions: 0.9457
Insertions: 0.783
Deletions: 0.9412

### Implementation

For alignment utilizes a kmer dictionary, and splits the reads into 3 kmers of size 15. At least one of them should match closely. Using the newly identified candidate positions that matched, load the appropriate section from the reference genome and run the Dynamic Programming Edit Distance alignment algorithm on the aligned read and reference genome. Use backtracking during the edit distance function to identify the exact mutations taken to optimize the edit distance. Out of all the candidate windows, take the one with the lowest edit distance and put those mutations in a dictionary. The mutations dictionary is later filtered, throwing away all mutations that do not meet the threshold of occurences necessary.

