import pandas as pd
def parse_summary_file(filename, output):
    df = pd.read_table(filename)
    df = df.drop_duplicates(subset=['pdb'], keep=False)
    df_HL = df[['pdb', 'Hchain', 'Lchain', 'antigen_chain', 'antigen_name', 'antigen_species', 'compound', 'resolution']]
    df_HL.to_csv(output, sep="\t")

    return df_HL

parse_summary_file('../data/20230217_0084705_summary.tsv', 'pdb_with_only_one_chain_pairing.tsv')