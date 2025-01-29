from parse import *
from collections import defaultdict, Counter
from tqdm import tqdm
import argparse
import os

# Sliding Window to build kmer dictionary
def build_kmer_dict(genome, k=15):
    kmer_dict = defaultdict(list)
    for i in range(len(genome) - k + 1):
        kmer = genome[i:i+k]
        kmer_dict[kmer].append(i)
    return kmer_dict

# Process reads and divide into kmers
def divide_into_kmers(read, k=15):
    return [read[i:i+k] for i in range(0, len(read) - k + 1, k)]

# Get candidate positions for reads
def get_candidate_positions(read, kmer_dict, k=15):
    kmers = divide_into_kmers(read, k)
    candidate_positions = []
    for i in range(len(kmers)):
        if kmers[i] in kmer_dict:
            candidate_positions.append((kmer_dict[kmers[i]][0], i))
    return candidate_positions

def format_output(predictions):
    output = ""
    for i in range(len(predictions)):
        output += f">read_{i}\tGenome_Number_{predictions[i]}\n"
    return output

def main(folder, output_file, N):
    dictionaries = []

    # Build the kmer dictionaries
    print("Building Dictionaries...")
    for i in tqdm(range(N)):
        genome = parse_reference_genome(f"{folder}/{folder}_genome_{i}.fasta")
        kmer_dict = build_kmer_dict(genome)
        dictionaries.append(kmer_dict)

    # Map each one
    reads = parse_reads(f"{folder}/{folder}_reads.fasta")
    predictions = []
    freq = []

    print("Mapping reads...")
    for i in tqdm(range(len(reads))):
        read = reads[i]
        genome_mappings = []
        for j in range(N):
            kmer_dict = dictionaries[j]
            candidate_positions = get_candidate_positions(read, kmer_dict)
            genome_mappings.append(candidate_positions)
        
        most_mappings = max(len(_) for _ in genome_mappings)
        potential_genomes = []
        for j in range(N):
            if len(genome_mappings[j]) == most_mappings:
                potential_genomes.append(j)
                freq.append(j)

        predictions.append(potential_genomes)

        # print(f"Read {i}: {potential_read_mappings}")

    freq = Counter(freq)

    # Handle reads that have multiple potential genomes
    # Select the most frequent genome since statistically more probable it came from that one
    for i in range(len(reads)):
        # There are more than one option
        if len(predictions[i]) > 1:
            most_freq_genome = max(predictions[i], key=lambda x: freq[x])
            predictions[i] = most_freq_genome        
        else:
            predictions[i] = predictions[i][0]

    output = format_output(predictions)
    with open(output_file, "w") as file:
        file.write(output)
    
    print(output)
    print(f"Successfully outputted to {output_file}")
    
if __name__ == "__main__":
    # Handle CLI
    parser = argparse.ArgumentParser(description="Project 1C Solution")
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