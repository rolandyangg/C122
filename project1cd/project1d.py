from parse import *
from collections import Counter
from tqdm import tqdm
from pybloom_live import BloomFilter
import argparse
import os

# Process reads and divide into kmers
def divide_into_kmers(read, k=15):
    return [read[i:i+k] for i in range(0, len(read) - k + 1, k)]

def build_bloom_filter(genome, k=15, bloom_capacity=1000000, error_rate=0.01):
    bloom_filter = BloomFilter(capacity=bloom_capacity, error_rate=error_rate)
    for i in range(len(genome) - k + 1):
        kmer = genome[i:i+k]
        bloom_filter.add(kmer)
    return bloom_filter

def check_read_in_bloom_filter(read, bloom_filter, k=15):
    kmers = divide_into_kmers(read, k)
    matches = 0
    for kmer in kmers:
        if kmer in bloom_filter:
            matches += 1
    return matches

def format_output(predictions):
    output = ""
    for i in range(len(predictions)):
        output += f">read_{i}\tGenome_Number_{predictions[i]}\n"
    return output

def main(folder, output_file, N):
    bloom_filters = []
    index_to_genome = [] # Have to map the indices back to the correct genome number when we use the filtered genomes

    print("Building Bloom Filters...")
    i = 0
    for filename in tqdm(os.listdir(folder)):
        prefix = f"{folder}_genome_"
        if filename.startswith(prefix):
            genome_number = int(filename[len(prefix):-len(".fasta")])
            genome = parse_reference_genome(f"{folder}/{folder}_genome_{genome_number}.fasta")
            index_to_genome.append(genome_number)
            bloom_filter = build_bloom_filter(genome)
            bloom_filters.append(bloom_filter)
            i += 1

    # Map each one
    reads = parse_reads(f"{folder}/{folder}_reads.fasta")
    predictions = []
    freq = []

    print("Mapping reads...")
    for i in tqdm(range(len(reads))):
        read = reads[i]
        genome_mappings = []
        for j in range(N):
            bloom_filter = bloom_filters[j]
            genome_mappings.append(check_read_in_bloom_filter(read, bloom_filter))
        
        most_mappings = max(genome_mappings)
        potential_genomes = []
        for j in range(N):
            if genome_mappings[j] == most_mappings:
                potential_genomes.append(j)
                freq.append(j)

        predictions.append(potential_genomes)

        # print(f"Read {i}: {potential_read_mappings}")

    freq = Counter(freq)
    print(freq)
    temp_freq = []

    # Handle reads that have multiple potential genomes
    # Select the most frequent genome since statistically more probable it came from that one
    for i in range(len(reads)):
        # There are more than one option
        if len(predictions[i]) > 1:
            most_freq_genome = max(predictions[i], key=lambda x: freq[x])
            predictions[i] = index_to_genome[most_freq_genome] # Map back to the correct genome when listing answers
        else:
            predictions[i] = index_to_genome[predictions[i][0]]
        temp_freq.append(predictions[i])
    
    print(Counter(temp_freq))

    output = format_output(predictions)
    with open(output_file, "w") as file:
        file.write(output)

    # with open("counters2.txt", "w") as file:
    #     file.write(str(Counter(temp_freq)))
    
    # with open("counters.txt", "w") as file:
    #     file.write(str(freq))
    
    # print(output)
    print(f"Successfully outputted to {output_file}")
    
if __name__ == "__main__":
    # Handle CLI
    parser = argparse.ArgumentParser(description="Project 1D Solution")
    parser.add_argument("-i", "--input", required=True, dest="input", help="Input data folder path (should include genomes and a reads file)")
    parser.add_argument("-o", "--output", required=True, dest="output", help="Output file path (should be a txt)")
    parser.add_argument("-N", "--num", required=True, dest="num", help="Number of genomes in folder")

    # NOTE: This Implementation assumes your genomes and reads file use the folder name as the prefix
    # To run it with the input data you must change the folder name of the data to "project1c"

    args = parser.parse_args()
    input_folder = args.input
    output_file = args.output
    N = int(args.num)

    if not os.path.exists(input_folder):
        print("Input folder does not exist")
        exit(1)

    main(input_folder, output_file, N)