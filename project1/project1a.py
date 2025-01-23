from parse import *
from collections import Counter

def solution(reference_genome, reads, max_mismatches, mutation_threshold):
    mutations = {}

    for read in reads:
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

def main():
    MISMATCH_THRESHOLD = 3 # Amount of mismatches allowed in a read
    MUTATION_THRESHOLD = 3 # If the mutation occured less than this many times, it is thrown out

    reference_genome = parse_ref_file("data/project1a_reference_genome.fasta")
    reads = parse_fasta("data/project1a_with_error_paired_reads.fasta")

    predictions = solution(reference_genome, reads, MISMATCH_THRESHOLD, MUTATION_THRESHOLD)

    # Sort to be in order of index
    predictions = sorted(predictions, key=lambda x: x[0])

    output = format_output(predictions)
    print(output)

    # Dump the output to a file
    with open("output/1a/predictions.txt", "w") as file:
        file.write(output)

if __name__ == "__main__":
    main()
