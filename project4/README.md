### Roland Yang, UID: 506053914, C122 Winter 2025

# Project 4 C122

This project identifies genomicc intervals overalpping genes based on histone modification data. Given data of 500,000 intervals consisting of 200bp intervals each assigned four characters representing whether a certain amount of marks for X and Y from our emissions pass a certain threshold, this project identifies the top 50,000 intervals that overlap the body of an annotated protein coding gene.

## Usage

```bash
python project4.py -i [input-filepath] -o [output-filepath]
```

- -i, --input: The input data consisting of the intervals.
- -o, --output: The output filepath, should be a .txt file

### Score

Overall: 48138/50000

### Implementation

First we initialize the parameters to reasonable values. We make it such that the transition matrices represent that each state is more likely to stay to itself and that the emission matrix is more likely to emit 'n' much more in one state compared to the other state which is more likely to emit the other alphabet values.

We use the Baum-Welch EM Algorithm at 50 iterations to obtain updated transition and emission matrices. Using these matrices we then run it through another forward-backward algorithm to calculate the probability of each state for each interval in the input data. In both algorithms the forward and backward calculations are scaled using the sum of their values along the way to handle the large input. When calculating xi we normalize it using the scaling factors. When updating the emission and transition matrices, we normalize them using their respective sums. These probabilities are then sorted so we can obtain the 50,000 intervals that have the highest probability of being in the "gene" state.