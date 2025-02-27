from parse import *
from tqdm import tqdm
import argparse
import os
import json

# SETUP

def motifs_to_profile(motifs, k, pseudocount=1):
    profile = [{"A": 0, "C": 0, "G": 0, "T": 0} for _ in range(k)]
    
    # Get counts
    for motif in motifs:
        for i in range(len(motif)):
            nucleotide = motif[i]
            profile[i][nucleotide] += 1
    
    # Normalize
    t = len(motifs)
    for i in range(len(profile)):
        for nucleotide in profile[i]:
            profile[i][nucleotide] += pseudocount
            # print(profile[i][nucleotide],  t + pseudocount * 4) # MAYBE DEBUG HERE
            profile[i][nucleotide] /= t + pseudocount * 4
    
    return profile

def initialize_motifs(sequences, k=20):
    # Motif is first kmer in every sequence
    res = []
    for sequence in sequences:
        res.append(sequence[:k])
    return res

######################################################################

# EM ALGORITHM

def e_step(sequences, pwm, k=20):
    # The size of each row in the array is the amount of kmers in each sequence
    # sequences x (amount of kmers in each sequence)
    probabilities = [[0] * (len(sequences[0]) - k + 1) for _ in range(len(sequences))]

    for i in range(len(sequences)):
        sequence = sequences[i]
        denominator = 0
        temp_probs = [] # Store each probability calculated to avoid repeats

        # Compute the denominator
        for j in range(len(sequence) - k + 1):
            kmer = sequence[j:j + k]
            kmer_prob = 1
            for z in range(len(kmer)):
                prob = pwm[z][kmer[z]]
                kmer_prob *= prob
            denominator += kmer_prob
            temp_probs.append(kmer_prob)

        # Go back and compute the probability for each kmer
        for j in range(len(temp_probs)):
            # print(len(pwm))
            # print(i, j)
            probabilities[i][j] = temp_probs[j] / denominator

    return probabilities

def m_step(sequences, pwm, probabilities, k=20, pseudocount=1):
    for i in range(len(sequences)):
        sequence = sequences[i]
        for j in range(len(sequence) - k + 1):
            kmer = sequence[j:j + k]
            for z in range(len(kmer)):
                nucleotide = kmer[z]
                pwm[z][nucleotide] += probabilities[i][j]

    # Normalize
    t = len(sequences)
    for i in range(len(pwm)):
        for nucleotide in pwm[i]:
            pwm[i][nucleotide] += pseudocount
            pwm[i][nucleotide] /= t + pseudocount * 4

    return pwm

######################################################################

# MAPPING

def profile_most_probable_kmer(text: str, k: int,
                               profile: list[dict[str, float]]) -> str:
    res = ""
    max_score = -1
    for i in range(len(text) - k + 1):
        temp_score = 1
        for j in range(k):
            curr_nucleotide = text[i + j]
            temp_score *= profile[j][curr_nucleotide]
        if temp_score > max_score:
            max_score = temp_score
            res = text[i:i + k]
    return res

######################################################################

def main(pwm_sequences, sequences, output_file):
    # Generate the PWM using EM Algorithm
    K = 20
    iterations = 100

    motifs = initialize_motifs(pwm_sequences, K) # Set to the first kmer in every sequence
    profile = motifs_to_profile(motifs, K)

    for i in tqdm(range(iterations)):
        probabilities = e_step(pwm_sequences, profile, K)
        profile = m_step(pwm_sequences, profile, probabilities, K)
        print(profile)

    # Use the PWM to find peaks in the sequences
    predictions = []
    for seq in sequences:
        temp = profile_most_probable_kmer(seq, K, profile)
        predictions.append(seq.find(temp))

    print(predictions)

    with open(output_file, 'w') as file:
        for i in range(len(sequences)):
            file.write(f"seq{i + 1}\t{predictions[i]}\n")
    
    print("Success!")

if __name__ == "__main__":
    # Handle CLI
    parser = argparse.ArgumentParser(description="Project 3a Solution")
    parser.add_argument("-s", "--sequences", required=True, dest="sequences", help="Sequences file path to generate PWM")
    parser.add_argument("-i", "--input", required=True, dest="input", help="Sequences file path to find peaks on using PWM")
    parser.add_argument("-o", "--output", required=True, dest="output", help="Output file path (should be a txt)")

    args = parser.parse_args()
    pwm_sequences_file = args.sequences
    input_file = args.input
    output_file = args.output

    if not os.path.exists(pwm_sequences_file):
        print("PWM Sequence file does not exist")
        exit(1)

    if not os.path.exists(input_file):
        print("Sequence file to identify peaks does not exist")
        exit(1)

    pwm_sequences = parse_reads(pwm_sequences_file)
    sequences = parse_reads(input_file)
    main(pwm_sequences, sequences, "predictions.txt")