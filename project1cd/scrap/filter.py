import shutil
from tqdm import tqdm

valid_genomes = []
with open("valid_genomes.txt", "r") as file:
    for line in file:
        valid_genomes.append(line.strip())

for i in tqdm(valid_genomes):
    shutil.copy(f"project1d/project1d_genome_{i}.fasta", f"project1d-filtered/project1d-filtered_genome_{i}.fasta")