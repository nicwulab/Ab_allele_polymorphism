from Bio import SeqIO

def get_sequence(pdb_list, dir, df, output_file):
    file_format = '.pdb'
    count = 0

    for pdb in pdb_list:
        pdbfile = dir+pdb+file_format

        pdb_row = df.loc[pdb]
        lchain = pdb_row['Lchain']
        hchain = pdb_row['Hchain']

        with open(pdbfile, 'r') as pdb_file:
            for record in SeqIO.parse(pdbfile, 'pdb-atom'):
                # print('>' + record.id)
                # print(record.seq)
                chain = record.id.split(":")[1]
                if chain == lchain or chain == hchain:
                    with open(output_file, "a") as output_handle:
                        SeqIO.write(record, output_handle, "fasta")
                        count += 1
    print(str(count) + ' pdbs saved to fasta')

