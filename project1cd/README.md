### Roland Yang, UID: 506053914, C122 Winter 2025

# Project 1c C122

This project performs a metagenomic analysis on a sample. There are many genome sequences which may in the sample and many reads from the sample as well. Our project predicts which reads are from which genome in the sample.

## Usage

```bash
python project1X.py -i [input-folder-filepath] -o [output-filepath] -N [num-genomes]
```

- -i, --input: The input data folder path (*must be the prefix of the file names!*)
- -o, --output: The output filepath, should be a .txt file
- -N, --num: Number of genomes in the folder

Note: Folder name should be the prefix of all the files in the folder in order for these implementations to work.

## Project 1c

100 Genomes, 400,000 Reads

### Score

Overall: 0.9738

### Implementation

Builds a kmer dictionary for every genome in the folder. Then for each read splits it into 3 kmers of size 15. Identify candidate positions in each genome. For each read identifies the potential candidate positions in each genome. Store the genomes with the most candidate positions in another array for each read. Afterwards, filter through this extra array, by going through each genomes potential genomes and choosing the most frequent genome that occured. This makes sense because it is more statistically likely that out of the reads that match multiple genomes, that it matches the one that is most common.

## Project 1d

5,000 Genomes, 4,000,000 Reads

### Score

Overall: 0.9747

### Implementation

Builds a bloom filter for every genome in the folder. Splits each read into 3 kmers of size 15. Identifies the amount of those split kmers found in the bloom filter for each genome. Store the genomes with the most appearances in another array for each read. Afterwards, similarly as last implementation, filter through extra array going through each genomes potential genomes and choosing most frequent genome that occured.

The main innovation to solving this problem was to identify the actual set of genomes that was present in the sample. Running this implementation on a small set of reads (70k), I was able to narrow it down to 49 genomes that were actually present. From there, I ran this implementation with all 4,000,000 reads and those 49 genomes in approximately 10 minutes.