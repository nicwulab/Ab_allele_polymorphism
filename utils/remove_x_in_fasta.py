from Bio import SeqIO

def remove_x_in_fasta(input_file, output_file):
    # Read in the input FASTA file
    input_file = input_file
    records = list(SeqIO.parse(input_file, "fasta"))

    # Edit the sequences
    for record in records:
        # Change the first character of the sequence to "X"
        sequence = record.seq
        new_sequence = sequence.replace("X", "")
        record.seq = new_sequence

    # Write out the new FASTA file
    output_file = output_file
    with open(output_file, "w") as handle:
        SeqIO.write(records, handle, "fasta")
