# Help Credited: https://github.com/eeskin/CM122_starter_code/blob/master/HP1/basic_aligner.py
    
def parse_reads(file_path):
    reads = []
    with open(file_path, 'r') as file:
        current_read = ""
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if current_read:
                    reads.append(current_read)
                current_read = ""
            else:
                current_read += line
        if current_read:
            reads.append(current_read)
    return reads

def parse_reference_genome(file_path):
    with open(file_path, 'r') as gFile:
        first_line = True
        ref_genome = ''
        for line in gFile:
            if first_line:
                first_line = False
                continue
            ref_genome += line.strip()
    return ref_genome