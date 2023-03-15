import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from Ab_allele.utils.sinaplot import sinaplot


def plot_histogram(df, column, graph_title, xlabel, ylabel):
    df['floatddg'] = df[column].astype('float64')

    q = df['floatddg'].quantile(0.99)
    qdf = df[df['floatddg'] < q]
    q_low = df['floatddg'].quantile(0.01)
    q_hi = df['floatddg'].quantile(0.99)

    df_filtered = df[(df["floatddg"] < q_hi) & (df["floatddg"] > q_low)]
    list_filtered = df_filtered["floatddg"].tolist()
    sns.histplot(list_filtered, stat='count')
    plt.title(graph_title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    plt.savefig('ddgs_histogram_0.99.png')
    plt.clf()

    l = df[column].tolist()
    intl = [float(i) for i in l if 'nan' not in i]

    sns.histplot(intl, stat='count')
    plt.title(graph_title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.savefig('ddgs_histogram.png')
    plt.clf()

    # DDG cap = 20 and floor = -10
    l = df[column].tolist()
    intl = [float(i) for i in l if 'nan' not in i]
    for i in range(len(intl)):
        if intl[i] > 20:
            intl[i] = 20
        if intl[i] < -10:
            intl[i] = -10
    sns.histplot(intl, stat='count')
    plt.title(graph_title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.savefig('ddgs_histogram20.png')
    plt.clf()


def plot_baseline_freq_ddgs_scatter_plot(baseline_freq_ddgs_df, column_name):
    baseline_freq_ddgs_df['baseline_freq_original'] = baseline_freq_ddgs_df['baseline_freq_original'].astype('float64')
    baseline_freq_ddgs_df['DDG'] = baseline_freq_ddgs_df['DDG'].astype('float64')

    sns.scatterplot(data=baseline_freq_ddgs_df, x="baseline_freq_original", y="DDG")
    plt.savefig(column_name+'baseline_freq_ddgs_scatterplot.png')
    plt.clf()

    baseline_freq_ddgs_df.loc[baseline_freq_ddgs_df['DDG'] < -10, 'DDG'] = -10
    baseline_freq_ddgs_df.loc[baseline_freq_ddgs_df['DDG'] > 20, 'DDG'] = 20
    sns.scatterplot(data=baseline_freq_ddgs_df, x="baseline_freq_original", y="DDG")
    plt.savefig(column_name+'baseline_freq_ddgs_scatterplot_20.png')
    plt.clf()

    r = baseline_freq_ddgs_df.corr(method="spearman")
    print(r)

def plot_baseline_freq_species_sinaplot(baseline_freq_species):

    plt.gcf().set_size_inches(25, 25)

    baseline_freq_species['baseline_freq_original'] = baseline_freq_species['baseline_freq_original'].astype('float64')
    ddgs_species = baseline_freq_species.dropna()

    sinaplot(x='baseline_freq_original', y='antigen_species_id', data=ddgs_species)
    plt.savefig('baseline_freq_species_sinaplot.png')
    plt.clf()


def plot_ddgs_species_sinaplot(ddgs_species):

    plt.gcf().set_size_inches(25,25)

    ddgs_species['DDG'] = ddgs_species['DDG'].astype('float64')
    ddgs_species = ddgs_species.dropna()

    list_filtered = ddgs_species['DDG'].tolist()

    sinaplot(x='DDG', y='antigen_species_id', data=ddgs_species)
    plt.savefig('ddgs_species_sinaplot.png')
    plt.clf()

    q_low = ddgs_species['DDG'].quantile(0.01)
    q_hi = ddgs_species['DDG'].quantile(0.99)

    df_filtered = ddgs_species[(ddgs_species["DDG"] < q_hi) & (ddgs_species["DDG"] > q_low)]
    list_filtered = df_filtered['DDG'].tolist()

    sinaplot(x='DDG', y='antigen_species_id', data=df_filtered)
    plt.savefig('ddgs_species_sinaplot0.99.png')
    plt.clf()

    # DDG cap = 20 and floor = -10
    ddgs_species.loc[ddgs_species['DDG'] < -10, 'DDG'] = -10
    ddgs_species.loc[ddgs_species['DDG'] > 20, 'DDG'] = 20

    sinaplot(x='DDG', y='antigen_species_id', data=ddgs_species)
    plt.savefig('ddgs_species_sinaplot20.png')
    plt.clf()

def main():

    epitope_identification_with_baseline_variation_df = pd.read_csv('epitope_identification_with_baseline_variation.csv')
    ddgs = epitope_identification_with_baseline_variation_df[['DDG']]
    baseline_freq_ddgs = epitope_identification_with_baseline_variation_df[['baseline_freq_original', 'DDG']]
    ddgs_species = epitope_identification_with_baseline_variation_df[['DDG', 'antigen_species_id']]
    baseline_freq_species = epitope_identification_with_baseline_variation_df[['baseline_freq_original', 'antigen_species_id']]


    plot_histogram(ddgs, 'DDG', 'Distribution of DDG', 'DDG value', 'Count')
    plot_baseline_freq_ddgs_scatter_plot(baseline_freq_ddgs, '')
    HIV_df = epitope_identification_with_baseline_variation_df[epitope_identification_with_baseline_variation_df['antigen_species_id'] == 'HIV']
    HIV_df = HIV_df[['baseline_freq_original', 'DDG']]
    plot_baseline_freq_ddgs_scatter_plot(HIV_df, 'HIV')

    influenza_df = epitope_identification_with_baseline_variation_df[epitope_identification_with_baseline_variation_df['antigen_species_id'] == 'influenza']
    influenza_df = influenza_df[['baseline_freq_original', 'DDG']]
    plot_baseline_freq_ddgs_scatter_plot(influenza_df, 'influenza')

    plot_ddgs_species_sinaplot(ddgs_species)
    plot_baseline_freq_species_sinaplot(baseline_freq_species)


if __name__ == "__main__":
    main()
