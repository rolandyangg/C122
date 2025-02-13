### Roland Yang, UID: 506053914, C122 Winter 2025

# Project 2 C122

This project is a genome assembler. From a set of reads or a spectrum of kmers, it will reconstruct the genome.

## Usage

```bash
python project2X.py -i [input-filepath] -o [output-filepath]
```

- -i, --input: The input data path
- -o, --output: The output filepath, should be a .txt file

## Project 2a

Generates/reconstructs a sequence of reads from a spectrum. 

### Score

Overall: 0.8376

### Implementation

Create Bruijn graph out of a spectrum. Find a Eulerian path in the Bruijn graph. Reconstruct the genome from the Eulerian path.

Use a modified version of project1b.py solution to map the reads back to the reconstructed genome (dynamic programming edit distance method). The reads are then sorted by their identified potential position in the genome.