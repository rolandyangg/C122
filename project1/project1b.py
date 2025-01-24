from parse import *
from collections import defaultdict
from tqdm import tqdm
import argparse
import os

# Dynamic programming method
def edit_distance(a, b, offset):
    m, n = len(a), len(b)
    # Create the DP table
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    operations = [[None] * (n + 1) for _ in range(m + 1)]
    
    # Fill the base cases
    for i in range(m + 1):
        dp[i][0] = i
        operations[i][0] = "I"  # Insertion
    for j in range(n + 1):
        dp[0][j] = j
        operations[0][j] = "D"  # Deletion
    
    # Fill the DP table and track operations
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if a[i - 1] == b[j - 1]:  # Match
                dp[i][j] = dp[i - 1][j - 1]
                operations[i][j] = "M"
            else:
                # Choose the best operation (substitution, insertion, deletion)
                options = [
                    (dp[i - 1][j] + 1, "I"),  # Insertion
                    (dp[i][j - 1] + 1, "D"),  # Deletion
                    (dp[i - 1][j - 1] + 1, "S")  # Substitution
                ]
                dp[i][j], operations[i][j] = min(options)
    
    # Backtrace to identify operations
    mutations = []
    i, j = m, n
    while i > 0 or j > 0:
        if operations[i][j] == "M":  # Match
            i -= 1
            j -= 1
        elif operations[i][j] == "S":  # Substitution
            mutations.append(("S", i-1 + offset, f"{b[j-1]} {a[i-1]}"))
            i -= 1
            j -= 1
        elif operations[i][j] == "I":  # Insertion
            mutations.append(("I", i-1 + offset, f"{a[i-1]}"))
            i -= 1
        elif operations[i][j] == "D":  # Deletion
            mutations.append(("D", j-1 + offset, f"{b[j-1]}"))
            j -= 1
    
    return dp[m][n], mutations[::-1]  # Return edit distance and mutations

# Sliding Window to build kmer dictionary
def build_kmer_dict(genome, k=15):
    kmer_dict = defaultdict(list)
    for i in tqdm(range(len(genome) - k + 1)):
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

def main(reference_genome, reads, output_file):
    print("Building kmer dictionary...")
    kmer_dict = build_kmer_dict(reference_genome)
    total_mutations = {}
    print("Identifying mutations...")
    for read in tqdm(reads):
        # Divide kmer of size 50 into 3 pieces of size 15
        spots = get_candidate_positions(read, kmer_dict)
        temp = []
        for i in range(len(spots)):
            start = spots[i][0] - spots[i][1] * 15
            reference_window = reference_genome[start:start + 50] # Get the aligned reference window
            distance, mutations = edit_distance(read, reference_window, start) # Dyanmic Programming on aligned read
            temp.append((distance, mutations))
        
        # Choose the best scoring window
        if len(temp) > 0:
            best_mutations = min(temp, key=lambda x: x[0])[1]

            for mutation in best_mutations:
                if mutation not in total_mutations:
                    total_mutations[mutation] = 0
                total_mutations[mutation] += 1
    
    predictions = []

    # Filter
    for mutation in total_mutations:
        if total_mutations[mutation] >= 6:
            predictions.append(mutation)
    
    # Format
    priority = {"S": 0, "I": 1, "D": 2}
    predictions = sorted(predictions, key=lambda x: (priority[x[0]], x[1]))

    output = ""
    for prediction in predictions:
        output += f">{prediction[0]}{prediction[1]} {prediction[2]}\n"

    # Output
    with open(output_file, "w") as file:
        file.write(output)
    print(output)
    print(f"Successfully outputted to {output_file}")

if __name__ == "__main__":
    # Handle CLI
    parser = argparse.ArgumentParser(description="Project 1B Solution")
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