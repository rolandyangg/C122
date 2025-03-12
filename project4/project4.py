from parse import *
from tqdm import tqdm
import argparse
import os
import json
import numpy as np

######################################################################

def forward_alg(x, alphabet, states, T, E):
    n = len(x)
    num_states = len(states)
    symbol = {s: i for i, s in enumerate(alphabet)}

    forward = np.zeros((num_states, n))
    scaling_factors = np.zeros(n)

    for k in range(num_states):
        forward[k, 0] = E[k, symbol[x[0]]] / num_states  # Uniform start probability
    
    scaling_factors[0] = np.sum(forward[:, 0])
    forward[:, 0] /= scaling_factors[0]

    for i in range(1, n):
        for k in range(num_states):
            forward[k, i] = sum(forward[j, i-1] * T[j, k] * E[k, symbol[x[i]]] for j in range(num_states))
        
        # Normalize
        scaling_factors[i] = np.sum(forward[:, i])
        if scaling_factors[i] == 0:
            scaling_factors[i] = 1e-10  # Prevent division by zero
        forward[:, i] /= scaling_factors[i]
    
    return forward, scaling_factors

def backward_alg(x, alphabet, states, T, E):
    n = len(x)
    num_states = len(states)
    symbol = {s: i for i, s in enumerate(alphabet)}

    backward = np.zeros((num_states, n))
    backward[:, n-1] = 1  # Base case

    scaling_factors = np.zeros(n)
    scaling_factors[n-1] = 1  # Base case for backward scaling

    for i in range(n-2, -1, -1):
        for k in range(num_states):
            backward[k, i] = sum(backward[j, i+1] * T[k, j] * E[j, symbol[x[i+1]]] for j in range(num_states))
        
        # Normalize
        scaling_factors[i] = np.sum(backward[:, i])
        if scaling_factors[i] == 0:
            scaling_factors[i] = 1e-10  # Prevent division by zero
        backward[:, i] /= scaling_factors[i]

    return backward, scaling_factors


######################################################################

def baum_welch(x, alphabet, states, transition, emission, num_iterations):
    n = len(x)
    num_states = len(states)
    symbol = {s: i for i, s in enumerate(alphabet)}
    T = np.array([[transition[s1][s2] for s2 in states] for s1 in states])
    E = np.array([[emission[s][sym] for sym in alphabet] for s in states])
    
    for _ in tqdm(range(num_iterations)):
        # E-Step
        forward, scaling_factors = forward_alg(x, alphabet, states, T, E)
        backward, backward_scaling_factors = backward_alg(x, alphabet, states, T, E)
        
        # M-Step
        # Use forward and backward probabilities to compute new values
        gamma = np.zeros((num_states, n))
        xi = np.zeros((num_states, num_states, n-1))

        # Gamma
        for i in range(n):
            for k in range(num_states):
                gamma[k, i] = (forward[k, i] * backward[k, i])

        # Xi
        for i in range(n-1):
            for k in range(num_states):
                for j in range(num_states):
                    xi[k, j, i] = (forward[k, i] * T[k, j] * E[j, symbol[x[i+1]]] * backward[j, i+1])
            xi[:, :, i] /= scaling_factors[i+1]  # Normalize xi using next scaling factor

        # Update Transition and Emission matrices
        for k in range(num_states):
            T[k, :] = np.sum(xi[k, :, :], axis=1) / np.sum(gamma[k, :-1])
            T[k, :] /= np.sum(T[k, :])  # Normalize the transition probabilities

            for sym in range(len(alphabet)):
                E[k, sym] = np.sum(gamma[k, np.array([symbol[xi] for xi in x]) == sym]) / np.sum(gamma[k, :])
            E[k, :] /= np.sum(E[k, :])  # Normalize the emission probabilities

    # Convert results to output format
    updated_transition = {states[k]: {states[j]: T[k, j] for j in range(num_states)} for k in range(num_states)}
    updated_emission = {states[k]: {alphabet[sym]: E[k, sym] for sym in range(len(alphabet))} for k in range(num_states)}

    return updated_transition, updated_emission

######################################################################

def posterior_decode(x, alphabet, states, transition, emission):
    n = len(x)
    num_states = len(states)
    T = np.array([[transition[s1][s2] for s2 in states] for s1 in states])
    E = np.array([[emission[s][sym] for sym in alphabet] for s in states])
    
    forward, scaling_factors = forward_alg(x, alphabet, states, T, E)
    backward, backward_scaling_factors = backward_alg(x, alphabet, states, T, E)
    
    posterior = np.zeros((num_states, n))
    for i in range(n):
        for k in range(num_states):
            posterior[k, i] = (forward[k, i] * backward[k, i])

    return posterior

######################################################################

def main(intervals, output_file):
    alphabet = ["x", "y", "z", "n"]
    states = ["S1", "S2"] # S1 is the gene state, S2 is the Non-Gene state
    initial_transition = {
        "S1": {"S1": 0.9, "S2": 0.1},
        "S2": {"S1": 0.1, "S2": 0.9}
    }
    initial_emission = {
        "S1": {"x": 0.3, "y": 0.3, "z": 0.3, "n": 0.1},
        "S2": {"x": 0.1, "y": 0.1, "z": 0.1, "n": 0.7}
    }
    num_interations = 50
    updated_transition, updated_emission = baum_welch(intervals[0], alphabet, states, initial_transition, initial_emission, num_interations)
    print(updated_emission, updated_transition)
    probabilities = posterior_decode(intervals[0], alphabet, states, updated_transition, updated_emission)
    S_probs = probabilities[0]
    indexed_S_probs = list(enumerate(S_probs))
    res = sorted(indexed_S_probs, key=lambda x: x[1], reverse=True)

    with open(output_file, 'w') as file:
        for i in range(50000):
            file.write(f"{res[i][0] + 1}\n")

    print("Success!")

if __name__ == "__main__":
    # Handle CLI
    parser = argparse.ArgumentParser(description="Project 4 Solution")
    parser.add_argument("-i", "--input", required=True, dest="input", help="Input file path")
    parser.add_argument("-o", "--output", required=True, dest="output", help="Output file path (should be a txt)")

    args = parser.parse_args()
    input_file = args.input
    output_file = args.output

    if not os.path.exists(input_file):
        print("Input file does not exist")
        exit(1)

    intervals = parse_reads(input_file)
    main(intervals, "predictions.txt")