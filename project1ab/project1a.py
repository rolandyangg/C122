from parse import *
from tqdm import tqdm
import argparse
import os

MISMATCH_THRESHOLD = 3 # Amount of mismatches allowed in a read
MUTATION_THRESHOLD = 3 # If the mutation occured less than this many times, it is thrown out

def sliding_window_solution(reference_genome, reads, max_mismatches=3, mutation_threshold=3):
    mutations = {}

    for read in tqdm(reads):
        # Use sliding window to find the aligned window for each read
        for i in range(len(reference_genome) - len(read)):
            window = reference_genome[i:i+len(read)]
            mismatches = 0
            mismatch_values = []

            for j in range(len(window)):
                if read[j] != window[j]:
                    mismatches += 1
                    mismatch_values.append((i + j, window[j], read[j]))
            if mismatches <= max_mismatches:
                # Update mutations dictionary
                for mismatch in mismatch_values:
                    if mismatch not in mutations:
                        mutations[mismatch] = 0
                    mutations[mismatch] += 1
                break

    # Filter out mutations that occur less than the mutation threshold
    predictions = []

    for mutation in mutations:
        if mutations[mutation] >= mutation_threshold:
            predictions.append(mutation)
    
    return predictions # Outputs a list of tuples --> [(index, reference, mutation), ...]

def format_output(predictions):
    output = ""
    for prediction in predictions:
        output += f">S{prediction[0]} {prediction[1]} {prediction[2]}\n"
    return output

def main(reference_genome, reads, output_file):
    print("Running sliding window and identfying mutations...")
    predictions = sliding_window_solution(reference_genome, reads, MISMATCH_THRESHOLD, MUTATION_THRESHOLD)

    # Sort to be in order of index
    predictions = sorted(predictions, key=lambda x: x[0])
    output = format_output(predictions)

    # Dump the output to a file
    with open(output_file, "w") as file:
        file.write(output)
    print(output)
    print(f"Successfully outputted to {output_file}")

if __name__ == "__main__":
    # Handle CLI
    parser = argparse.ArgumentParser(description="Project 1A Solution")
    parser.add_argument("-g", "--genome", required=True, dest="genome", help="Reference genome file path")
    parser.add_argument("-r", "--reads", required=True, dest="reads", help="Reads file path")
    parser.add_argument("-o", "--output", required=True, dest="output", help="Output file path (should be a txt)")

    args = parser.parse_args()
    reference_genome_file = args.genome
    reads_file = args.reads
    output_file = args.output

    if not os.path.exists(reference_genome_file):
        print("Reference genome file does not exist")
        exit(1)
    
    if not os.path.exists(reads_file):
        print("Reads file does not exist")
        exit(1)

    reference_genome = parse_reference_genome(reference_genome_file)
    reads = parse_reads(reads_file)

    main(reference_genome, reads, output_file)