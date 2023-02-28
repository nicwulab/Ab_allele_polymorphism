import ast
import json


def parse_result_file(filename):
    file1 = open(filename, 'r')
    string = file1.readline()
    json_acceptable_string = string.replace("'", "\"")
    d = json.loads(json_acceptable_string)
    return d

def parse_result_with_tuple(filename):
    # Open the file for reading
    with open(filename, 'r') as f:
        # Read the contents of the file as a string
        contents = f.read()

    # Use ast.literal_eval to convert the string into a dictionary
    return ast.literal_eval(contents)

def parse_pyir_output(file):
    # Open the JSON file
    list_of_dict = []

    file1 = open(file, 'r')
    Lines = file1.readlines()

    count = 0
    # Strips the newline character
    for line in Lines:
        count += 1
        aDict = json.loads(line)
        list_of_dict.append(aDict)

    return list_of_dict

def write_list_to_file(l, file_name):
    import csv

    with open(file_name, "w") as f:
        wr = csv.writer(f)
        wr.writerows(l)

def save_fasta_as_csv(fasta_file):
    from Bio import SeqIO
    import csv

    # Open the FASTA file for reading
    with open(fasta_file, 'r') as fasta_file:
        # Parse the FASTA file using Biopython
        records = SeqIO.parse(fasta_file, 'fasta')

        # Open a CSV file for writing
        with open('output.csv', 'w', newline='') as csv_file:
            # Create a CSV writer object
            writer = csv.writer(csv_file)

            # Write the header row
            writer.writerow(['id', 'sequence'])

            # Loop over the sequences and write each one to the CSV file
            for record in records:
                writer.writerow([record.id, str(record.seq)])


