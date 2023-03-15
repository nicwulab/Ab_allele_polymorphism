import pandas as pd
import collections


def get_allele_percentage(all_paired_antibodies_df, gene_type, result_df):
    new_all_paired_antibodies_df = all_paired_antibodies_df.loc[:, ['Name', gene_type]]

    new_all_paired_antibodies_df[gene_type+'_gene_one'] = new_all_paired_antibodies_df[gene_type].str.split(',').str[
        0]

    new_all_paired_antibodies_df[[gene_type+'_gene_one', 'allele']] = new_all_paired_antibodies_df[
        gene_type+'_gene_one'].str.split('*', expand=True)


    grouped_df = new_all_paired_antibodies_df.groupby(gene_type+'_gene_one')

    # for each group, do pairwise comparison
    for group_id, group_df in grouped_df:

        gene_id = group_df[gene_type+'_gene_one'].unique()[0]
        value_counts = group_df['allele'].value_counts(normalize=True)

        percentages = value_counts * 100

        for index, row in percentages.iteritems():
            value = index
            percentage = row

            new_row = {'gene_id': gene_id,
                       'allele': value,
                       'percentage': percentage}

            result_df = result_df.append(new_row, ignore_index=True)

    return result_df


def igv_csv_to_df(filename):
    df = pd.read_csv(filename)
    df[['Id', 'allele']] = df['Id'].str.split('*', expand=True)
    df = df.drop('domain_no', axis=1)
    df = df.drop('hmm_species', axis=1)
    df = df.drop('e-value', axis=1)
    df = df.drop('identity_species', axis=1)
    df = df.drop('v_gene', axis=1)
    df = df.drop('v_identity', axis=1)
    df = df.drop('j_gene', axis=1)
    df = df.drop('j_identity', axis=1)
    df = df.drop('seqstart_index', axis=1)
    df = df.drop('score', axis=1)
    df = df.drop('seqend_index', axis=1)
    df = df.drop('chain_type', axis=1)
    return df

def get_igv_percentage(igv_H_df, igv_KL_df, epitope_identification_with_group_id_df, table_one_df):
    result_df = pd.DataFrame(columns=['gene_id', 'location', 'variant', 'percentage'])

    for index, row_data in epitope_identification_with_group_id_df.iterrows():
        gene_id = row_data['gene']
        location = row_data['location']
        epitope_amino_acid_original = row_data['amino_acid_original']

        if gene_id[:3] == 'IGH':
            igv = igv_H_df
        elif gene_id[:3] == 'IGK' or gene_id[:3] == 'IGL':
            igv = igv_KL_df

        gene_id_df = igv[igv['Id'] == gene_id]

        d = collections.defaultdict(int)

        for gene_index, gene_row_data in gene_id_df.iterrows():
            allele = gene_row_data['allele']
            ac_name = gene_row_data[location]

            selected_row = table_one_df.loc[(table_one_df['allele'] == allele) & (table_one_df['gene_id'] == gene_id)]
            if selected_row.empty:
                continue
            d[ac_name] += selected_row['percentage'].item()

        for k, v in d.items():
            new_row = {'gene_id': gene_id,
                       'location': location,
                       'variant': k,
                       'percentage': v}
            result_df = result_df.append(new_row, ignore_index=True)

        epitope_identification_with_group_id_df.loc[index, 'baseline_freq_original'] = d[epitope_amino_acid_original]

    return result_df, epitope_identification_with_group_id_df


def main():
    all_paired_antibodies_df = pd.read_csv('data/all_paired_antibodies_from_GB_v6.csv')
    result_df = pd.DataFrame(columns=['gene_id', 'allele', 'percentage'])

    df = get_allele_percentage(all_paired_antibodies_df, 'Heavy_V_gene', result_df)
    result_df = result_df.append(df, ignore_index=True)
    df = get_allele_percentage(all_paired_antibodies_df, 'Light_V_gene', result_df)
    result_df = result_df.append(df, ignore_index=True)
    result_df = result_df.drop_duplicates()
    result_df = result_df.reset_index(drop=True)
    result_df.to_csv('baseline_variation_table_1.csv')

    igv_H_df = igv_csv_to_df('data/anarci_igv_output.csv_H.csv')
    igv_KL_df = igv_csv_to_df('data/anarci_igv_output.csv_KL.csv')
    epitope_identification_with_group_id_df = pd.read_csv('results/epitope_idenfication/epitope_identification_with_group_id.csv', index_col=0)

    df, epitope_identification_with_group_id_df = get_igv_percentage(igv_H_df, igv_KL_df, epitope_identification_with_group_id_df, result_df)
    df = df.drop_duplicates()
    df = df.reset_index(drop=True)
    df.to_csv('baseline_variation_table_2.csv')
    epitope_identification_with_group_id_df = epitope_identification_with_group_id_df[['pdb_id', 'chain_id', 'location',
                                                                                       'dssp_location', 'gene', 'amino_acid_original',
                                                                                       'list_amino_acid_variants', 'variants', 'baseline_freq_original',
                                                                                       'description', 'antigen_name', 'antigen_chain',
                                                                                       'antigen_species', 'resolution',
                                                                                       'DDG','antigen_loc_list',
                                                                                       'antigen_species_id','groupID']]
    epitope_identification_with_group_id_df['groupID'] = epitope_identification_with_group_id_df['groupID'].astype(str).apply(lambda x: x.replace('.0', ''))
    epitope_identification_with_group_id_df['groupID'] = epitope_identification_with_group_id_df['groupID'].astype(str).apply(lambda x: x.replace('nan', ''))

    epitope_identification_with_group_id_df.to_csv('epitope_identification_with_baseline_variation.csv')


if __name__ == "__main__":
    main()