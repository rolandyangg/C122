from parse import *
from collections import defaultdict, deque, Counter
from tqdm import tqdm
import argparse
import os

def build_de_bruijn(patterns):
    graph = defaultdict(list)
    i = 0
    for pattern in tqdm(patterns):
        prefix = pattern[:-1]
        suffix = pattern[1:]
        graph[prefix].append(suffix)

        i += 1
    return graph

def find_eulerian(graph):
    indegree = defaultdict(int)
    outdegree = defaultdict(int)
    
    for u in graph:
        outdegree[u] = len(graph[u])
        for v in graph[u]:
            indegree[v] += 1
    
    start = next(iter(graph))
    for node in set(indegree.keys()).union(set(outdegree.keys())):
        if outdegree[node] == indegree[node] + 1:
            start = node
            break
    
    stack = [start]
    path = deque()
    
    while stack:
        while graph[stack[-1]]:
            stack.append(graph[stack[-1]].pop())
        path.appendleft(stack.pop())
    
    return path

# Dynamic programming method
def edit_distance(a, b):
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
    
    return dp[m][n] # Return edit distance

# Sliding Window to build kmer dictionary
def build_kmer_dict(genome, k=12):
    kmer_dict = defaultdict(list)
    for i in tqdm(range(len(genome) - k + 1)):
        kmer = genome[i:i+k]
        kmer_dict[kmer].append(i)
    return kmer_dict

# Process reads and divide into kmers
def divide_into_kmers(read, k=12):
    return [read[i:i+k] for i in range(0, len(read) - k + 1, k)]

# Get candidate positions for reads
def get_candidate_positions(read, kmer_dict, k=12):
    kmers = divide_into_kmers(read, k)
    candidate_positions = []
    for i in range(len(kmers)):
        if kmers[i] in kmer_dict:
            candidate_positions.append((kmer_dict[kmers[i]][0], i))
    # print(candidate_positions)
    return candidate_positions

######################################################

def generate_spectrum(reads, k=12, threshold=3):
    kmer_freq = Counter()

    for read in reads:
        for i in range(len(read) - k + 1):
            kmer = read[i:i + k]
            kmer_freq[kmer] += 1
    
    spectrum = []
    for kmer, count in kmer_freq.items():
        if count > threshold:
            spectrum.append(kmer)

    return spectrum

######################################################

def main(reads, output_file):
    k = 12
    threshold = 3

    spectrum = generate_spectrum(reads, k, threshold)

    print("Building De Bruijn Graph...")
    graph = build_de_bruijn(spectrum)

    path = find_eulerian(graph)
    # print(path)

    genome = path[0] + "".join(node[-1] for node in path[1:])

    # Map the reads spectrum back to the genome
    kmer_dict = build_kmer_dict(genome, k)
    read_to_positions = {}
    print("Mapping reads back to genome...")
    read_index = 0
    for read in tqdm(reads):
        # Divide kmer of size 50 into 3 pieces of size 12
        spots = get_candidate_positions(read, kmer_dict, k)
        temp = []
        for i in range(len(spots)):
            start = spots[i][0] - spots[i][1] * k
            reference_window = genome[start:start + 50] # Get the aligned reference window
            distance = edit_distance(read, reference_window) # Dyanamic Programming on aligned read
            temp.append((distance, start))

        # Choose the best scoring window
        if len(temp) > 0:
            best_position = min(temp, key=lambda x: x[0])[1]
            read_to_positions[read_index] = best_position
        else:
            read_to_positions[read_index] = 0 # Default, shouldn't occur but just worst case...

        read_index += 1

    # Sort
    res = sorted(read_to_positions.keys(), key=lambda x: read_to_positions[x])

    # Format
    output = ""
    for r in res:
        output += f">read_{r}\n"

    # Output
    with open(output_file, "w") as file:
        file.write(output)
    print(genome)
    print(f"Successfully outputted to {output_file}")

if __name__ == "__main__":
    # Handle CLI
    parser = argparse.ArgumentParser(description="Project 2b Solution")
    parser.add_argument("-i", "--input", required=True, dest="input", help="Reads file path")
    parser.add_argument("-o", "--output", required=True, dest="output", help="Output file path (should be a txt)")

    args = parser.parse_args()
    reads_file = args.input
    output_file = args.output

    if not os.path.exists(reads_file):
        print("Spectrum file does not exist")
        exit(1)

    reads = parse_reads(reads_file)

    main(reads, output_file)