from collections import Counter
import ast

def load_counter_from_file(filename):
    with open(filename, "r") as f:
        content = f.read().strip()  # Read the file and remove extra spaces/newlines

    # Convert the string back into a Counter object
    return ast.literal_eval(content[8:-1])  # Remove the "Counter(" and ")" from the string

res = list(sorted(load_counter_from_file("counters.txt").items(), key=lambda item: item[1], reverse=True))
res2 = list(sorted(load_counter_from_file("counters2.txt").items(), key=lambda item: item[1], reverse=True))

output = ""
output2 = ""

valid_genomes = []

for z in range(len(res2)):
    if res2[z][1] > 1000:
        valid_genomes.append(res2[z][0])

for i in range(len(res)):
    # print(res[i][0], res[i][1])
    output += f"{res[i][0]} {res[i][1]}\n"

for j in range(len(res2)):
    output2 += f"{res2[j][0]} {res2[j][1]}\n"

with open("counters_sorted.txt", "w") as file:
    file.write(output)

with open("counters2_sorted.txt", "w") as file:
    file.write(output2)

print(len(valid_genomes))
with open("valid_genomes.txt", "w") as file:
    for genome in valid_genomes:
        file.write(f"{genome}\n")