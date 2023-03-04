from Bio import SeqIO
import pandas as pd
import sys

def parse_summary_file(filename):
    df = pd.read_table(filename)
    df = df.drop_duplicates(subset=['pdb'], keep=False)
    df_HL = df[['pdb', 'Hchain', 'Lchain', 'compound']]

    return df_HL

def get_sequence(dir, output_file, summary_file):
    file_format = '.pdb'
    count = 0

    df = parse_summary_file(summary_file)
    # print(df)
    pdb_list = df['pdb'].tolist()
    df.set_index("pdb", inplace=True)

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


def main(dir, output_file, summary_file):
    # get positions of interacting amino acids for each pdb with only one chain pairings
    get_sequence(dir, output_file, summary_file)

if __name__ == "__main__":
    pdb_dir = sys.argv[1]
    output_file = sys.argv[2]
    summary_file = sys.argv[3]
    main(pdb_dir, output_file, summary_file)